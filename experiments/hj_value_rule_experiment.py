#!/usr/bin/env python3
"""Seek vs. avoid monochromatic lines: a value-rule experiment on
Hales-Jewett CNFs, run through CaDiCaL's decision callback (IPASIR-UP).

Every decision is made by our code:
  variable = first unassigned point in a fixed order (descending incidence,
             ties by index) -- identical in every arm;
  value    = chosen by the arm's rule.

Arms
  cadical   CaDiCaL decides (decide() returns 0); reference only
  constant  fixed order, value always True
  random    fixed order, value by seeded coin flip
  seek      fixed order, value = colour that extends monochromatic runs most
  avoid     fixed order, value = colour that extends monochromatic runs least

Score used by seek/avoid, for an unassigned point x and colour c:
  s_c(x) = sum over lines L through x whose assigned points all have
           colour c (and at least one is assigned) of 2^(number assigned).

Requires: pip install python-sat   (needs CaDiCaL 1.9.5 = 'cadical195')

Usage
  python3 hj_value_rule_experiment.py --dry-run
  python3 hj_value_rule_experiment.py
  python3 hj_value_rule_experiment.py --instances 4:3,4:4 --seeds 3
"""

import argparse
import csv
import datetime as dt
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

ARMS = ["cadical", "constant", "random", "seek", "avoid"]

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
    def __init__(self, arm, order, through, lines, seed):
        super().__init__()
        self.arm = arm
        self.order = order                   # variables (1-based), fixed
        self.rank = {v: r for r, v in enumerate(order)}
        self.through = through
        self.nlines = len(lines)
        self.rng = random.Random(seed)
        nvars = len(order)
        self.value = [0] * (nvars + 1)        # 0 unassigned, +1 True, -1 False
        self.cnt = [[0] * self.nlines, [0] * self.nlines]  # [False, True]
        self.trail = []                      # assigned variables, in order
        self.levels = []                     # trail length at each new level
        self.ptr = 0                         # all order[:ptr] are assigned
        self.decisions_ours = 0
        self.mismatch = None

    # --- bookkeeping ---------------------------------------------------
    def _assign(self, lit):
        v = abs(lit)
        s = 1 if lit > 0 else -1
        if self.value[v] == s:
            return                           # duplicate notification
        if self.value[v] == -s:              # should not happen; be safe
            self._unassign(v)
        self.value[v] = s
        c = 1 if s > 0 else 0
        row = self.cnt[c]
        for li in self.through[v - 1]:
            row[li] += 1
        self.trail.append(v)

    def _unassign(self, v):
        s = self.value[v]
        if s == 0:
            return
        c = 1 if s > 0 else 0
        row = self.cnt[c]
        for li in self.through[v - 1]:
            row[li] -= 1
        self.value[v] = 0
        r = self.rank[v]
        if r < self.ptr:
            self.ptr = r

    # --- IPASIR-UP callbacks -------------------------------------------
    def on_assignment(self, lit, fixed=False):
        self._assign(lit)

    def on_new_level(self):
        self.levels.append(len(self.trail))

    def on_backtrack(self, to):
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
        if self.arm == "cadical":
            return 0
        # next unassigned variable in the fixed order
        order, value = self.order, self.value
        while self.ptr < len(order) and value[order[self.ptr]] != 0:
            self.ptr += 1
        if self.ptr == len(order):
            return 0
        v = order[self.ptr]
        self.decisions_ours += 1
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


# --------------------------------------------------------------------------
# 3. One run
# --------------------------------------------------------------------------

def run_one(task):
    k, n, arm, seed = task
    points, lines, through, inc = build_instance(k, n)
    nvars = len(points)
    cls = clauses_for(lines)
    random.Random(seed).shuffle(cls)          # clause-order permutation

    # fixed variable order: descending incidence, ties by index
    order = sorted(range(1, nvars + 1), key=lambda v: (-inc[v - 1], v))

    prop = ValueRule(arm, order, through, lines, seed)
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
        "k": k, "n": n, "arm": arm, "seed": seed,
        "status": "SAT" if sat else "UNSAT",
        "conflicts": stats.get("conflicts", ""),
        "decisions": stats.get("decisions", ""),
        "propagations": stats.get("propagations", ""),
        "restarts": stats.get("restarts", ""),
        "our_decisions": prop.decisions_ours,
        "seconds": round(elapsed, 3),
        "mono_lines_in_model": mono,          # must be 0 for SAT
        "state_mismatch": prop.mismatch if prop.mismatch is not None else "",
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
        k, n, arm, seed = task
        return {"k": k, "n": n, "arm": arm, "seed": seed, "status": "TIMEOUT"}
    return q.get() if not q.empty() else {
        "k": task[0], "n": task[1], "arm": task[2], "seed": task[3],
        "status": "CRASH"}


def _child(q, task):
    q.put(run_one(task))


# --------------------------------------------------------------------------
# 4. Driver
# --------------------------------------------------------------------------

FIELDS = ["k", "n", "arm", "seed", "status", "conflicts", "decisions",
          "propagations", "restarts", "our_decisions", "seconds",
          "mono_lines_in_model", "state_mismatch", "vars", "lines"]


def summarise(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["k"], r["n"], r["arm"]), []).append(r)
    print(f"\n{'instance':<12}{'arm':<10}{'done':>6}{'median conf':>14}"
          f"{'IQR':>12}{'median dec/conf':>17}")
    for (k, n, arm) in sorted(groups, key=lambda g: (-g[0], g[1],
                                                      ARMS.index(g[2]))):
        g = groups[(k, n, arm)]
        ok = [r for r in g if r["status"] in ("SAT", "UNSAT")]
        conf = sorted(r["conflicts"] for r in ok)
        if conf:
            q = statistics.quantiles(conf, n=4) if len(conf) >= 4 else \
                [conf[0], 0, conf[-1]]
            dpc = [r["decisions"] / max(r["conflicts"], 1) for r in ok]
            print(f"HJ({k};2),{n:<4}{arm:<10}{len(ok):>3}/{len(g):<2}"
                  f"{statistics.median(conf):>14,.0f}{q[2] - q[0]:>12,.0f}"
                  f"{statistics.median(dpc):>17.2f}")
        else:
            print(f"HJ({k};2),{n:<4}{arm:<10}{0:>3}/{len(g):<2}   (none finished)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instances", default="",
                    help="comma list of k:n, e.g. 4:5,3:4 (default: built-in)")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--first-seed", type=int, default=20260927)
    ap.add_argument("--timeout", type=int, default=180,
                    help="seconds per run (s03 kills processes near 225 s)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default="value_rule_results.csv")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    instances = ([tuple(int(x) for x in s.split(":"))
                  for s in args.instances.split(",")]
                 if args.instances else DEFAULT_INSTANCES)
    arms = [a for a in args.arms.split(",") if a]
    seeds = [args.first_seed + i for i in range(args.seeds)]
    tasks = [(k, n, arm, s) for (k, n) in instances for arm in arms
             for s in seeds]

    print("Value-rule experiment (CaDiCaL 1.9.5 via PySAT IPASIR-UP)")
    for (k, n) in instances:
        _, lines, _, _ = build_instance(k, n)
        print(f"  HJ({k};2),{n}: {k**n:,} variables, {len(lines):,} lines, "
              f"{'SAT' if k == 4 else 'UNSAT'} expected")
    print(f"  arms: {', '.join(arms)}")
    print(f"  {len(tasks)} runs, {args.workers} in parallel, "
          f"timeout {args.timeout}s each")
    if args.dry_run:
        return 0

    rows = []
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        with mp.pool.ThreadPool(args.workers) as pool:
            for i, row in enumerate(pool.imap_unordered(
                    lambda t: worker(t, args.timeout), tasks), 1):
                w.writerow(row)
                fh.flush()
                rows.append(row)
                print(f"[{i}/{len(tasks)}] HJ({row['k']};2),{row['n']} "
                      f"{row['arm']:<8} seed {row['seed']}  {row['status']:<7} "
                      f"conflicts {row.get('conflicts', '-')}", flush=True)
    summarise(rows)

    bad = [r for r in rows if r.get("mono_lines_in_model", "") not in ("", 0)]
    if bad:
        print(f"\nWARNING: {len(bad)} SAT models failed the line check")
    drift = [r for r in rows if r.get("state_mismatch", "") not in ("", 0)]
    if drift:
        print(f"NOTE: {len(drift)} runs had tracked-state drift at the model "
              f"(see state_mismatch column)")
    print(f"\nResults: {Path(args.out).resolve()}")
    return 0


if __name__ == "__main__":
    import multiprocessing.pool  # noqa: F401  (for ThreadPool)
    sys.exit(main())