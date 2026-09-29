#!/usr/bin/env python3
"""Seek vs. avoid monochromatic lines: a value-rule experiment on
Hales-Jewett CNFs, run through CaDiCaL's decision callback (IPASIR-UP).

v4 = v2 + the "saved" arm (phase saving) + the "dynseek" variable order.

Every decision is made by our code:
  variable = first unassigned point in a fixed order (incidence, index,
             random), or, for the dynseek order, the unassigned point with
             the largest S(x) = s0(x) + s1(x), ties by incidence order;
  value    = chosen by the arm's rule.

Arms
  cadical   CaDiCaL decides (decide() returns 0); reference only
  saved     value = the variable's last assigned value (initially True)
  constant  value always True
  random    value by seeded coin flip
  seek      value = colour that extends monochromatic runs most
  avoid     value = colour that extends monochromatic runs least

Score used by seek/avoid, for an unassigned point x and colour c:
  s_c(x) = sum over lines L through x whose assigned points all have
           colour c (and at least one is assigned) of 2^(number assigned).
dynseek keeps S(y) = sum over lines L through y of w(L) incrementally, where
w(L) = 2^a if the a assigned points of L all have one colour, else 0.

Requires: pip install python-sat   (needs CaDiCaL 1.9.5 = 'cadical195')

Usage (variable orders: --orders incidence,index,random,dynseek; --resume
skips finished runs; --summary reprints the summary of --out)
  python3 hj_value_rule_experiment_v4.py --self-test
  python3 hj_value_rule_experiment_v4.py --instances 4:5 --orders dynseek --seeds 2
  python3 hj_value_rule_experiment_v4.py --summary --out results.csv
"""

import argparse
import csv
import datetime as dt
import heapq
import itertools
import multiprocessing as mp
import random
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

from pysat.engines import Propagator
from pysat.solvers import Solver

ARMS = ["cadical", "saved", "constant", "random", "seek", "avoid"]

# (k, n): k letters, dimension n, 2 colours. HJ(4;2) is SAT for these n;
# HJ(3;2) is UNSAT for every n >= 4.
DEFAULT_INSTANCES = [(4, 3), (4, 4), (4, 5), (4, 6), (3, 4), (3, 5), (3, 6)]


# --------------------------------------------------------------------------
# 1. Encoding
# --------------------------------------------------------------------------

def build_instance(k, n):
    """Points of [k]^n, numbered like the group's encoder: variable v is the
    point whose base-k digits (most significant first) spell v-1.
    Returns (points, lines, through, inc) with 0-based point indices."""
    points = list(itertools.product(range(k), repeat=n))
    index = {p: i for i, p in enumerate(points)}
    lines = []
    for root in itertools.product(range(-1, k), repeat=n):   # -1 = star
        if -1 not in root:
            continue
        lines.append([index[tuple(a if r == -1 else r for r in root)]
                      for a in range(k)])
    through = [[] for _ in points]
    for li, line in enumerate(lines):
        for p in line:
            through[p].append(li)
    inc = [len(t) for t in through]        # = number of lines through p
    return points, lines, through, inc


def clauses_for(lines):
    """Two clauses per line: not all colour 1, not all colour 0."""
    cls = []
    for line in lines:
        vs = [p + 1 for p in line]
        cls.append([-v for v in vs])      # forbids all True
        cls.append(vs)                    # forbids all False
    return cls


# --------------------------------------------------------------------------
# 2. The propagator: tracks the assignment, makes every decision
# --------------------------------------------------------------------------

class ValueRule(Propagator):
    def __init__(self, arm, order, through, lines, seed, dynamic=False,
                 check=False):
        super().__init__()
        self.arm = arm
        self.order = order                   # variables (1-based), fixed
        self.rank = {v: r for r, v in enumerate(order)}
        self.through = through
        self.nlines = len(lines)
        self.rng = random.Random(seed)
        nvars = len(order)
        self.value = [0] * (nvars + 1)        # 0 unassigned, +1 True, -1 False
        self.saved = [1] * (nvars + 1)        # last value, for the saved arm
        self.cnt = [[0] * self.nlines, [0] * self.nlines]  # [False, True]
        self.trail = []                      # assigned variables, in order
        self.levels = []                     # trail length at each new level
        self.ptr = 0                         # all order[:ptr] are assigned
        self.decisions_ours = 0
        self.mismatch = None
        # Decision audit: our literal must be the first assignment reported
        # on the level it opens; if another literal comes first, CaDiCaL
        # rejected ours and decided by itself.
        self.pending = self.expect = 0
        self.verified = self.fallbacks = 0
        # dynseek: variable = argmax S, kept incrementally (_shift) and
        # selected with a lazy max-heap (_pick_dynamic).
        self.dynamic = dynamic
        self.check = check                   # brute-force audit at each pick
        self.check_fail = 0
        if dynamic:
            self.lines1 = [[p + 1 for p in L] for L in lines]
            self.rk = [0] * (nvars + 1)       # incidence rank, as a list
            for r, v in enumerate(order):
                self.rk[v] = r
            self.W = [0] * self.nlines        # w(L)
            self.S = [0] * (nvars + 1)        # S(y), for every point
            self.dirty = [False] * (nvars + 1)
            self.dl = []                      # points to (re)push on the heap
            self.heap = list(range(nvars))    # keys rank - S*N; all S = 0
            self.heap_limit = 4 * nvars
            self.rebuilds = 0

    # --- bookkeeping ---------------------------------------------------
    def _assign(self, lit):
        v = abs(lit)
        s = 1 if lit > 0 else -1
        if self.value[v] == s:
            return                           # duplicate notification
        if self.value[v] == -s:              # should not happen; be safe
            self._unassign(v)
        self.value[v] = s
        self.saved[v] = s
        c = 1 if s > 0 else 0
        if self.dynamic:
            self._shift(v, c, 1)
        else:
            row = self.cnt[c]
            for li in self.through[v - 1]:
                row[li] += 1
        self.trail.append(v)

    def _unassign(self, v):
        s = self.value[v]
        if s == 0:
            return
        c = 1 if s > 0 else 0
        if self.dynamic:
            self._shift(v, c, -1)
            if not self.dirty[v]:
                self.dirty[v] = True
                self.dl.append(v)
        else:
            row = self.cnt[c]
            for li in self.through[v - 1]:
                row[li] -= 1
        self.value[v] = 0
        r = self.rank[v]
        if r < self.ptr:
            self.ptr = r

    def _shift(self, v, c, step):
        """Add step (+1 assign, -1 unassign) to the colour-c count of every
        line L through v; where w(L) changes, add the change to S(y) for all
        y on L and mark the unassigned ones for a heap push."""
        c0, c1 = self.cnt
        row = c1 if c else c0
        W, S, value = self.W, self.S, self.value
        dirty, dl, lines1 = self.dirty, self.dl, self.lines1
        for li in self.through[v - 1]:
            row[li] += step
            a0 = c0[li]
            a1 = c1[li]
            w = 1 << (a0 + a1) if (a0 == 0) != (a1 == 0) else 0
            d = w - W[li]
            if d:
                W[li] = w
                for y in lines1[li]:
                    S[y] += d
                    if not value[y] and not dirty[y]:
                        dirty[y] = True
                        dl.append(y)

    # --- IPASIR-UP callbacks -------------------------------------------
    def on_assignment(self, lit, fixed=False):
        if self.expect:                      # first report on a new level
            if lit == self.expect:
                self.verified += 1
            else:
                self.fallbacks += 1
            self.expect = 0
        self._assign(lit)

    def on_new_level(self):
        self.levels.append(len(self.trail))
        self.expect, self.pending = self.pending, 0

    def on_backtrack(self, to):
        self.expect = 0                      # decision undone before reported
        while len(self.levels) > to:
            start = self.levels.pop()
            while len(self.trail) > start:
                self._unassign(self.trail.pop())

    def check_model(self, model):
        # All constraints are CNF clauses, so CaDiCaL has already checked
        # them. We only record whether our tracked state agrees with it.
        self.mismatch = sum(1 for lit in model
                            if self.value[abs(lit)] != (1 if lit > 0 else -1))
        return True

    def propagate(self):
        return []

    def provide_reason(self, lit):
        return []

    def add_clause(self):
        return []

    def decide(self):
        self.pending = 0
        if self.arm == "cadical":
            return 0
        v = self._pick_dynamic() if self.dynamic else self._pick_fixed()
        if not v:
            return 0
        self.decisions_ours += 1
        self.pending = lit = self._polarity(v)
        return lit

    # --- variable choice -----------------------------------------------
    def _pick_fixed(self):
        """Next unassigned variable in the fixed order."""
        order, value = self.order, self.value
        while self.ptr < len(order) and value[order[self.ptr]] != 0:
            self.ptr += 1
        return order[self.ptr] if self.ptr < len(order) else 0

    def _pick_dynamic(self):
        """Unassigned variable with the largest S, ties by incidence rank.
        Heap keys are rank - S*N; an entry is live iff its variable is
        unassigned and its key matches the current S."""
        heap, S, value = self.heap, self.S, self.value
        order, rk, N = self.order, self.rk, len(self.order)
        dirty, dl = self.dirty, self.dl
        if len(heap) + len(dl) > self.heap_limit:
            heap[:] = [r - S[v] * N for r, v in enumerate(order)
                       if not value[v]]
            heapq.heapify(heap)
            self.rebuilds += 1
            for y in dl:
                dirty[y] = False
        else:
            for y in dl:
                dirty[y] = False
                if not value[y]:
                    heapq.heappush(heap, rk[y] - S[y] * N)
        dl.clear()
        v = 0
        while heap:
            key = heap[0]
            r = key % N
            if not value[order[r]] and key == r - S[order[r]] * N:
                v = order[r]
                break
            heapq.heappop(heap)
        if self.check:
            self._audit_pick(v)
        return v

    def _polarity(self, v):
        if self.arm == "saved":
            return v if self.saved[v] > 0 else -v
        if self.arm == "constant":
            return v
        if self.arm == "random":
            return v if self.rng.random() < 0.5 else -v
        s0, s1 = self._scores(v)
        if s0 == s1:
            return v if self.rng.random() < 0.5 else -v
        if self.arm == "seek":
            return v if s1 > s0 else -v
        return v if s1 < s0 else -v             # avoid

    def _scores(self, v):
        """s_c(x): weighted count of lines through x whose assigned points
        all have colour c."""
        c0, c1 = self.cnt
        s0 = s1 = 0
        for li in self.through[v - 1]:
            a0, a1 = c0[li], c1[li]
            if a1 == 0 and a0:
                s0 += 1 << a0
            elif a0 == 0 and a1:
                s1 += 1 << a1
        return s0, s1

    # --- brute-force audit (dynseek, --self-test) ------------------------
    def audit(self):
        """Recompute line counts, w(L) and S from self.value alone.
        Returns (S, number of disagreements with the incremental state)."""
        value = self.value
        S = [0] * len(value)
        bad = 0
        for li, L in enumerate(self.lines1):
            a0 = sum(1 for y in L if value[y] < 0)
            a1 = sum(1 for y in L if value[y] > 0)
            w = 1 << (a0 + a1) if (a0 == 0) != (a1 == 0) else 0
            bad += ((a0 != self.cnt[0][li]) + (a1 != self.cnt[1][li])
                    + (w != self.W[li]))
            for y in L:
                S[y] += w
        bad += sum(1 for y in range(1, len(S)) if S[y] != self.S[y])
        return S, bad

    def _audit_pick(self, v):
        """The pick must be the brute-force argmax of S (ties by rank) and
        S must equal s0 + s1 as computed by _scores."""
        S, bad = self.audit()
        free = [y for y in range(1, len(S)) if not self.value[y]]
        best = min(free, key=lambda y: (-S[y], self.rk[y])) if free else 0
        if bad or v != best or (v and S[v] != sum(self._scores(v))):
            self.check_fail += 1


# --------------------------------------------------------------------------
# 3. One run
# --------------------------------------------------------------------------

def variable_order(mode, inc, seed):
    """Fixed decision order over variables 1..N.
    incidence: decreasing inc(x), ties by increasing index (default)
    index:     increasing variable index (lexicographic order of points)
    random:    uniformly random, fixed per seed
    dynseek:   the incidence order, used only to break ties in S"""
    nvars = len(inc)
    if mode in ("incidence", "dynseek"):
        return sorted(range(1, nvars + 1), key=lambda v: (-inc[v - 1], v))
    if mode == "index":
        return list(range(1, nvars + 1))
    if mode == "random":
        order = list(range(1, nvars + 1))
        random.Random(f"order-{seed}").shuffle(order)
        return order
    raise ValueError(f"unknown order {mode!r}")


def run_one(task, check=False):
    k, n, arm, seed, order_mode = task
    points, lines, through, inc = build_instance(k, n)
    nvars = len(points)
    cls = clauses_for(lines)
    random.Random(seed).shuffle(cls)          # clause-order permutation

    order = variable_order(order_mode, inc, seed)
    dynamic = order_mode == "dynseek" and arm != "cadical"

    prop = ValueRule(arm, order, through, lines, seed, dynamic, check)
    t0 = time.perf_counter()
    with Solver(name="cadical195", bootstrap_with=cls) as s:
        s.connect_propagator(prop)
        for v in range(1, nvars + 1):
            s.observe(v)
        sat = s.solve()
        stats = s.accum_stats()
        model = s.get_model() if sat else None
    elapsed = time.perf_counter() - t0

    mono = ""
    if sat:
        col = {abs(l): l > 0 for l in model}
        mono = sum(1 for L in lines if len({col[p + 1] for p in L}) == 1)

    return {
        "k": k, "n": n, "order": order_mode, "arm": arm, "seed": seed,
        "status": "SAT" if sat else "UNSAT",
        "conflicts": stats.get("conflicts", ""),
        "decisions": stats.get("decisions", ""),
        "propagations": stats.get("propagations", ""),
        "restarts": stats.get("restarts", ""),
        "our_decisions": prop.decisions_ours,
        "verified": prop.verified,            # our literal seen on its level
        "fallbacks": prop.fallbacks,          # must be 0
        "seconds": round(elapsed, 3),
        "mono_lines_in_model": mono,          # must be 0 for SAT
        "state_mismatch": prop.mismatch if prop.mismatch is not None else "",
        "check_fail": prop.check_fail if check else "",
        "vars": nvars, "lines": len(lines),
    }


def worker(task, timeout):
    """Runs one task in a child process so it can be killed on timeout."""
    q = mp.Queue()
    p = mp.Process(target=_child, args=(q, task))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        p.join()
        k, n, arm, seed, order_mode = task
        return {"k": k, "n": n, "order": order_mode, "arm": arm, "seed": seed,
                "status": "TIMEOUT", "seconds": timeout}
    if not q.empty():
        return q.get()
    k, n, arm, seed, order_mode = task
    return {"k": k, "n": n, "order": order_mode, "arm": arm, "seed": seed,
            "status": "CRASH"}


def _child(q, task):
    q.put(run_one(task))


# --------------------------------------------------------------------------
# 4. Self-test: incremental S against brute force
# --------------------------------------------------------------------------

def self_test(steps=3000):
    """(a) Random sequences of decisions, propagations, duplicate reports and
    backtracks (restarts included) on HJ(4;2),4 and HJ(3;2),4: after every
    step the incremental counts, w(L) and S must equal a recomputation from
    the assignment, and every pick must be the brute-force argmax.
    (b) Real CaDiCaL runs with the same audit at every decision."""
    rng = random.Random(20260929)
    failures = 0
    print("(a) random assign/backtrack sequences")
    for (k, n) in [(4, 4), (3, 4)]:
        _, lines, through, inc = build_instance(k, n)
        N = len(inc)
        prop = ValueRule("avoid", variable_order("dynseek", inc, 0), through,
                         lines, 1, dynamic=True, check=True)
        bad = picks = backtracks = 0
        for _ in range(steps):
            level = len(prop.levels)
            free = [v for v in range(1, N + 1) if not prop.value[v]]
            r = rng.random()
            if not free or (level and r < 0.05):
                prop.on_backtrack(rng.randrange(level) if level else 0)
                backtracks += 1
            elif (level == 0 and r < 0.97) or r < 0.3:
                lit = prop.decide()              # pick audited inside
                picks += 1
                prop.on_new_level()
                if rng.random() < 0.1:           # conflict before the report
                    prop.on_backtrack(rng.randrange(level + 1))
                    backtracks += 1
                else:
                    prop.on_assignment(lit)
            else:
                v = rng.choice(free)
                lit = v if rng.random() < 0.5 else -v
                prop.on_assignment(lit, fixed=level == 0)
                if rng.random() < 0.05:
                    prop.on_assignment(lit)      # duplicate report
            bad += prop.audit()[1]
        f = bad + prop.check_fail + prop.fallbacks
        failures += f
        print(f"  HJ({k};2),{n}: {steps} steps, {picks} picks, {backtracks} "
              f"backtracks, {prop.rebuilds} heap rebuilds: state mismatches "
              f"{bad}, bad picks {prop.check_fail}, fallbacks "
              f"{prop.fallbacks}  {'OK' if f == 0 else 'FAIL'}")

    print("(b) CaDiCaL runs, audited at every decision")
    for (k, n) in [(4, 4), (4, 5), (3, 4), (3, 5)]:
        for arm in ("avoid", "saved"):
            row = run_one((k, n, arm, 20260927, "dynseek"), check=True)
            f = (row["check_fail"] + row["fallbacks"]
                 + (row["mono_lines_in_model"] or 0)
                 + (row["state_mismatch"] or 0)
                 + (row["status"] != ("SAT" if k == 4 else "UNSAT"))
                 + (row["decisions"] != row["our_decisions"]))
            failures += f
            print(f"  HJ({k};2),{n} {arm:<6} {row['status']:<6} conflicts "
                  f"{row['conflicts']:>6}  decisions {row['decisions']:>6} "
                  f"(ours {row['our_decisions']}, verified {row['verified']}, "
                  f"fallbacks {row['fallbacks']})  bad picks "
                  f"{row['check_fail']}  {'OK' if f == 0 else 'FAIL'}")
    print("SELF-TEST", "PASSED" if failures == 0 else f"FAILED ({failures})")
    return 0 if failures == 0 else 1


# --------------------------------------------------------------------------
# 5. Driver
# --------------------------------------------------------------------------

FIELDS = ["k", "n", "order", "arm", "seed", "status", "conflicts", "decisions",
          "propagations", "restarts", "our_decisions", "verified", "fallbacks",
          "seconds", "mono_lines_in_model", "state_mismatch", "vars", "lines"]


def summarise(rows):
    groups = {}
    for r in rows:
        groups.setdefault((int(r["k"]), int(r["n"]), r.get("order", "incidence"),
                           r["arm"]), []).append(r)
    print(f"\n{'instance':<12}{'order':<11}{'arm':<10}{'done':>6}"
          f"{'median conf':>14}{'IQR':>12}{'median GLR':>12}"
          f"{'median props':>15}{'median sec':>12}")
    for (k, n, om, arm) in sorted(groups, key=lambda g: (-g[0], g[1], g[2],
                                                         ARMS.index(g[3]))):
        g = groups[(k, n, om, arm)]
        ok = [r for r in g if r["status"] in ("SAT", "UNSAT")]
        conf = sorted(int(r["conflicts"]) for r in ok)
        if conf:
            q = statistics.quantiles(conf, n=4) if len(conf) >= 4 else \
                [conf[0], 0, conf[-1]]
            glr = [int(r["conflicts"]) / max(int(r["decisions"]), 1) for r in ok]
            props = statistics.median(int(r["propagations"]) for r in ok)
            secs = statistics.median(float(r["seconds"]) for r in ok)
            print(f"HJ({k};2),{n:<4}{om:<11}{arm:<10}{len(ok):>3}/{len(g):<2}"
                  f"{statistics.median(conf):>14,.0f}{q[2] - q[0]:>12,.0f}"
                  f"{statistics.median(glr):>12.4f}{props:>15,.0f}"
                  f"{secs:>12.1f}")
        else:
            print(f"HJ({k};2),{n:<4}{om:<11}{arm:<10}{0:>3}/{len(g):<2}"
                  f"   (none finished)")


def report_checks(rows):
    """Warnings for failed sanity checks (rows may come from the CSV)."""
    def nonzero(r, col):
        return str(r.get(col, "") or "") not in ("", "0")
    bad = [r for r in rows if nonzero(r, "mono_lines_in_model")]
    if bad:
        print(f"\nWARNING: {len(bad)} SAT models failed the line check")
    drift = [r for r in rows if nonzero(r, "state_mismatch")]
    if drift:
        print(f"NOTE: {len(drift)} runs had tracked-state drift at the model "
              f"(see state_mismatch column)")
    fb = [r for r in rows if nonzero(r, "fallbacks")]
    if fb:
        print(f"WARNING: {len(fb)} runs had decisions CaDiCaL made itself "
              f"(see fallbacks column)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instances", default="",
                    help="comma list of k:n, e.g. 4:5,3:4 (default: built-in)")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--orders", default="incidence",
                    help="comma list of variable orders: "
                         "incidence,index,random,dynseek")
    ap.add_argument("--resume", action="store_true",
                    help="append to --out and skip runs already recorded there")
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--first-seed", type=int, default=20260927)
    ap.add_argument("--timeout", type=int, default=180,
                    help="seconds per run (s03 kills processes near 225 s)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default="value_rule_results.csv")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--self-test", action="store_true",
                    help="check incremental S against brute force and exit")
    ap.add_argument("--summary", action="store_true",
                    help="print the summary of every row in --out and exit")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.summary:
        with open(args.out) as fh:
            rows = list(csv.DictReader(fh))
        summarise(rows)
        report_checks(rows)
        return 0

    instances = ([tuple(int(x) for x in s.split(":"))
                  for s in args.instances.split(",")]
                 if args.instances else DEFAULT_INSTANCES)
    arms = [a for a in args.arms.split(",") if a]
    seeds = [args.first_seed + i for i in range(args.seeds)]
    orders = [o for o in args.orders.split(",") if o]
    tasks = [(k, n, arm, s, o) for (k, n) in instances for o in orders
             for arm in arms for s in seeds]
    done = set()
    fields = FIELDS
    if args.resume and Path(args.out).exists():
        with open(args.out) as fh:
            reader = csv.DictReader(fh)
            fields = reader.fieldnames or FIELDS   # keep an older file's header
            for r in reader:
                done.add((int(r["k"]), int(r["n"]), r["arm"], int(r["seed"]),
                          r.get("order", "incidence")))
        tasks = [t for t in tasks if t not in done]

    print("Value-rule experiment (CaDiCaL 1.9.5 via PySAT IPASIR-UP)")
    for (k, n) in instances:
        _, lines, _, _ = build_instance(k, n)
        print(f"  HJ({k};2),{n}: {k**n:,} variables, {len(lines):,} lines, "
              f"{'SAT' if k == 4 else 'UNSAT'} expected")
    print(f"  arms: {', '.join(arms)}   orders: {', '.join(orders)}")
    if done:
        print(f"  resuming: {len(done)} runs already in {args.out}")
    print(f"  {len(tasks)} runs, {args.workers} in parallel, "
          f"timeout {args.timeout}s each")
    if args.dry_run:
        return 0

    rows = []
    append = args.resume and Path(args.out).exists()
    with open(args.out, "a" if append else "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if not append:
            w.writeheader()
        with mp.pool.ThreadPool(args.workers) as pool:
            for i, row in enumerate(pool.imap_unordered(
                    lambda t: worker(t, args.timeout), tasks), 1):
                w.writerow(row)
                fh.flush()
                rows.append(row)
                print(f"[{i}/{len(tasks)}] HJ({row['k']};2),{row['n']} "
                      f"{row['order']:<9} {row['arm']:<8} seed {row['seed']}  "
                      f"{row['status']:<7} "
                      f"conflicts {row.get('conflicts', '-')}", flush=True)
    summarise(rows)
    report_checks(rows)
    print(f"\nResults: {Path(args.out).resolve()}")
    return 0


if __name__ == "__main__":
    import multiprocessing.pool  # noqa: F401  (for ThreadPool)
    sys.exit(main())
