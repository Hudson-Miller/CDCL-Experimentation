#!/usr/bin/env python3
"""Symmetric HJ(4,2) calibration: which heuristic to run on Mouhib's open n=14 instance.

Instance (Mouhib, arXiv 2607.02226): colourings of [4]^n invariant under coordinate
permutations. Variables = cells = multiplicity vectors m in N^4 with sum n (C(n+3,3) of
them). Each combinatorial line with s stars and fixed part f becomes the quadruple
{f + s*e_a : a = 0..3}; constraint: not all four cells the same colour (C(n+3,4) quads).
n <= 13 is satisfiable (Mouhib's weight colouring); n = 14 is open.

Arms (same instance, same clause shuffle per seed); default = kissat,inc_avoid,dyn_avoid:
  kissat     plain Kissat binary on the same CNF (baseline; skipped if not found)
  inc_avoid  static incidence order + avoid value rule (IPASIR-UP propagator in CaDiCaL)
  dyn_avoid  dynseek variable rule + avoid value rule (IPASIR-UP propagator in CaDiCaL)
  cadical    plain CaDiCaL 1.9.5 (optional, not run by default)

Usage:
  python3 hj_sym_check.py --dry-run
  python3 hj_sym_check.py --n 12,13 --seeds 3 --timeout 14400 --workers 24 --kissat /path/to/kissat --out sym_check.csv
"""
import argparse, csv, itertools, multiprocessing as mp, os, random, re, shutil
import statistics, subprocess, sys, tempfile, time
from math import comb
from pathlib import Path

from pysat.engines import Propagator
from pysat.solvers import Solver

ARMS = ["kissat", "inc_avoid", "dyn_avoid", "cadical"]
DEFAULT_ARMS = ["kissat", "inc_avoid", "dyn_avoid"]
FIELDS = ["n", "arm", "seed", "status", "conflicts", "decisions", "seconds", "valid", "cells", "quads"]


# ---------------------------------------------------------------- instance
def comps(t, parts):
    if parts == 1:
        yield (t,)
        return
    for i in range(t + 1):
        for r in comps(t - i, parts - 1):
            yield (i,) + r


def build(n, k=4):
    cells = list(comps(n, k))
    cid = {c: i for i, c in enumerate(cells)}
    quads = set()
    for s in range(1, n + 1):
        for f in comps(n - s, k):
            quads.add(tuple(sorted(cid[tuple(f[j] + (s if j == a else 0) for j in range(k))]
                                   for a in range(k))))
    quads = sorted(quads)
    through = [[] for _ in cells]
    for li, q in enumerate(quads):
        for p in q:
            through[p].append(li)
    return cells, quads, through


def clauses(quads, seed):
    cls = [[p + 1 for p in q] for q in quads] + [[-(p + 1) for p in q] for q in quads]
    random.Random(seed).shuffle(cls)
    return cls


def valid(model_pos, quads):
    """model_pos: set of cells coloured True (0-based)."""
    return all(0 < sum(p in model_pos for p in q) < len(q) for q in quads)


# ---------------------------------------------------------------- propagator
class Rule(Propagator):
    """Makes every decision. Variable: static incidence order (dyn=False) or the
    unassigned cell maximising S = sum of 2^a over quads whose a assigned cells share
    one colour (dyn=True, ties by incidence order). Value: avoid."""

    def __init__(self, order, through, quads, seed, dyn):
        super().__init__()
        self.order, self.through, self.quads, self.dyn = order, through, quads, dyn
        self.rank = {v: r for r, v in enumerate(order)}
        self.rng = random.Random(seed)
        nv = len(order)
        self.value = [0] * (nv + 1)
        self.c0 = [0] * len(quads)
        self.c1 = [0] * len(quads)
        self.S = [0] * (nv + 1)
        self.trail, self.levels, self.ptr = [], [], 0
        self.our_decisions = 0

    @staticmethod
    def _w(a0, a1):
        if a1 == 0 and a0:
            return 1 << a0
        if a0 == 0 and a1:
            return 1 << a1
        return 0

    def _shift(self, v, colour, step):
        cnt = self.c1 if colour else self.c0
        for li in self.through[v - 1]:
            if self.dyn:
                before = self._w(self.c0[li], self.c1[li])
            cnt[li] += step
            if self.dyn:
                d = self._w(self.c0[li], self.c1[li]) - before
                if d:
                    for p in self.quads[li]:
                        self.S[p + 1] += d

    def _assign(self, lit):
        v, s = abs(lit), (1 if lit > 0 else -1)
        if self.value[v] == s:
            return
        if self.value[v] == -s:
            self._unassign(v)
        self.value[v] = s
        self._shift(v, s > 0, +1)
        self.trail.append(v)

    def _unassign(self, v):
        s = self.value[v]
        if s == 0:
            return
        self._shift(v, s > 0, -1)
        self.value[v] = 0
        if self.rank[v] < self.ptr:
            self.ptr = self.rank[v]

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
        return True

    def propagate(self):
        return []

    def provide_reason(self, lit):
        return []

    def add_clause(self):
        return []

    def decide(self):
        value = self.value
        if self.dyn:
            best, bv = -1, 0
            for v in self.order:                       # first strict max = best rank
                if value[v] == 0 and self.S[v] > best:
                    best, bv = self.S[v], v
            if not bv:
                return 0
            v = bv
        else:
            while self.ptr < len(self.order) and value[self.order[self.ptr]] != 0:
                self.ptr += 1
            if self.ptr == len(self.order):
                return 0
            v = self.order[self.ptr]
        self.our_decisions += 1
        s0 = s1 = 0
        for li in self.through[v - 1]:
            a0, a1 = self.c0[li], self.c1[li]
            if a1 == 0 and a0:
                s0 += 1 << a0
            elif a0 == 0 and a1:
                s1 += 1 << a1
        if s0 == s1:
            return v if self.rng.random() < 0.5 else -v
        return v if s1 < s0 else -v                    # avoid


# ---------------------------------------------------------------- runs
def run_one(task, kissat_bin):
    n, arm, seed = task
    cells, quads, through = build(n)
    cls = clauses(quads, seed)
    row = {"n": n, "arm": arm, "seed": seed, "cells": len(cells), "quads": len(quads)}
    t0 = time.perf_counter()
    if arm == "kissat":
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "f.cnf")
            with open(f, "w") as fh:
                fh.write(f"p cnf {len(cells)} {len(cls)}\n")
                fh.writelines(" ".join(map(str, c)) + " 0\n" for c in cls)
            out = subprocess.run([kissat_bin, f"--seed={seed % 2**31}", f],
                                 capture_output=True, text=True).stdout
        row["seconds"] = round(time.perf_counter() - t0, 3)
        m = re.search(r"^c conflicts:\s+(\d+)", out, re.M)
        row["conflicts"] = m.group(1) if m else ""
        m = re.search(r"^c decisions:\s+(\d+)", out, re.M)
        row["decisions"] = m.group(1) if m else ""
        if "s SATISFIABLE" in out:
            lits = [int(x) for l in out.splitlines() if l.startswith("v ") for x in l[2:].split()]
            row["status"], row["valid"] = "SAT", valid({l - 1 for l in lits if l > 0}, quads)
        else:
            row["status"] = "UNSAT" if "s UNSATISFIABLE" in out else "UNKNOWN"
        return row
    with Solver(name="cadical195", bootstrap_with=cls) as s:
        if arm != "cadical":
            inc = [len(t) for t in through]
            order = sorted(range(1, len(cells) + 1), key=lambda v: (-inc[v - 1], v))
            prop = Rule(order, through, quads, seed, dyn=(arm == "dyn_avoid"))
            s.connect_propagator(prop)
            for v in range(1, len(cells) + 1):
                s.observe(v)
        sat = s.solve()
        st = s.accum_stats()
        model = s.get_model() if sat else None
    row["seconds"] = round(time.perf_counter() - t0, 3)
    row["conflicts"], row["decisions"] = st.get("conflicts", ""), st.get("decisions", "")
    row["status"] = "SAT" if sat else "UNSAT"
    if sat:
        row["valid"] = valid({l - 1 for l in model if l > 0}, quads)
    return row


def _child(q, task, kissat_bin):
    q.put(run_one(task, kissat_bin))


def worker(task, timeout, kissat_bin):
    ctx = mp.get_context("spawn")                  # fork from threads can hang children
    q = ctx.Queue()
    p = ctx.Process(target=_child, args=(q, task, kissat_bin))
    p.start()
    p.join(timeout)
    n, arm, seed = task
    if p.is_alive():
        p.terminate()
        p.join()
        return {"n": n, "arm": arm, "seed": seed, "status": "TIMEOUT", "seconds": timeout}
    return q.get() if not q.empty() else {"n": n, "arm": arm, "seed": seed, "status": "CRASH"}


def summarise(rows):
    g = {}
    for r in rows:
        g.setdefault((int(r["n"]), r["arm"]), []).append(r)
    print(f"\n{'n':>3} {'arm':<10}{'solved':>8}{'median conflicts':>18}{'median sec (solved)':>21}{'all valid':>10}")
    for (n, arm) in sorted(g, key=lambda x: (x[0], ARMS.index(x[1]))):
        rs = g[(n, arm)]
        ok = [r for r in rs if r["status"] == "SAT"]
        conf = [int(r["conflicts"]) for r in ok if str(r.get("conflicts", "")).isdigit()]
        sec = [float(r["seconds"]) for r in ok]
        print(f"{n:>3} {arm:<10}{len(ok):>5}/{len(rs):<2}"
              f"{(f'{statistics.median(conf):,.0f}' if conf else '-'):>18}"
              f"{(f'{statistics.median(sec):.1f}' if sec else '-'):>21}"
              f"{str(all(str(r.get('valid')) == 'True' for r in ok)) if ok else '-':>10}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", default="12,13")
    ap.add_argument("--arms", default=",".join(DEFAULT_ARMS))
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--first-seed", type=int, default=20260927)
    ap.add_argument("--timeout", type=int, default=14400)
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--kissat", default="kissat", help="path to kissat binary")
    ap.add_argument("--out", default="sym_check.csv")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--summary", action="store_true", help="only reprint the summary of --out")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if a.summary:
        summarise(list(csv.DictReader(open(a.out))))
        return
    ns = [int(x) for x in a.n.split(",")]
    arms = [x for x in a.arms.split(",") if x]
    kb = shutil.which(a.kissat) or (a.kissat if os.path.exists(a.kissat) else None)
    if "kissat" in arms and not kb:
        print(f"kissat not found at '{a.kissat}': skipping the kissat arm")
        arms.remove("kissat")
    for n in ns:
        cells, quads, _ = build(n)
        assert len(cells) == comb(n + 3, 3) and len(quads) == comb(n + 3, 4)
        print(f"symmetric HJ(4,2) n={n}: {len(cells)} cells, {len(quads)} quads")
    tasks = [(n, arm, a.first_seed + i) for n in ns for arm in arms for i in range(a.seeds)]
    done = set()
    if a.resume and Path(a.out).exists():
        done = {(int(r["n"]), r["arm"], int(r["seed"])) for r in csv.DictReader(open(a.out))}
        tasks = [t for t in tasks if t not in done]
    print(f"{len(tasks)} runs, {a.workers} in parallel, timeout {a.timeout}s; kissat: {kb}")
    if a.dry_run:
        return
    rows = []
    append = a.resume and Path(a.out).exists()
    with open(a.out, "a" if append else "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if not append:
            w.writeheader()
        with mp.pool.ThreadPool(a.workers) as pool:
            for i, row in enumerate(pool.imap_unordered(lambda t: worker(t, a.timeout, kb), tasks), 1):
                w.writerow(row)
                fh.flush()
                rows.append(row)
                print(f"[{i}/{len(tasks)}] n={row['n']} {row['arm']:<10} seed {row['seed']} "
                      f"{row['status']:<8} conflicts {row.get('conflicts', '-')} "
                      f"sec {row.get('seconds', '-')} valid {row.get('valid', '-')}", flush=True)
    summarise(rows)


if __name__ == "__main__":
    import multiprocessing.pool  # noqa: F401
    sys.exit(main())
