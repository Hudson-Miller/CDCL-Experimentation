#!/usr/bin/env python3
"""Symmetric colourings for Hales-Jewett lower bounds.

A colouring of [k]^n that is invariant under permuting the n coordinates is
determined by the multiplicity vector (m_1,...,m_k) of a point rather than the
point itself.  That collapses k^n variables to C(n+k-1, k-1): at k=4, n=12 it
is 455 variables instead of 16.7 million.

A combinatorial line is a fixed word on n-m coordinates plus a moving set of
size m >= 1.  If the fixed part has multiplicity vector f, the line's k points
have vectors f + m*e_a for a = 1..k, so "no monochromatic line" becomes a
not-all-equal constraint over those k types.

Invariance under permuting the k SYMBOLS as well is impossible for every n:
the all-stars line runs through the k constant points, which all share the
multiset {n,0,...,0}, so it is monochromatic by construction.  The script
reports that rather than silently returning UNSAT.

Sanity checks: k=2 should be SAT to n=2 and UNSAT at n=3 (HJ(2;2)=3);
k=3 should be SAT to n=3 and UNSAT at n=4 (HJ(3;2)=4).  Both pass.

Requires: pip install python-sat
"""

import argparse, itertools, sys, time

try:
    from pysat.solvers import Minisat22
except ImportError:
    sys.exit("need python-sat:  pip install python-sat")


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def lines_as_types(n, k):
    """One entry per orbit of combinatorial lines: the k multiplicity vectors
    its points carry."""
    for m in range(1, n + 1):
        for f in compositions(n - m, k):
            pts = []
            for a in range(k):
                c = list(f)
                c[a] += m
                pts.append(tuple(c))
            yield tuple(pts)


def build(n, k, symbol_sym=False):
    idx, clauses = {}, []
    for pts in lines_as_types(n, k):
        keys = [tuple(sorted(p, reverse=True)) if symbol_sym else p for p in pts]
        for key in keys:
            if key not in idx:
                idx[key] = len(idx) + 1
        lits = sorted({idx[key] for key in keys})
        if len(lits) == 1:              # all k points share a type -> forced mono
            return idx, None
        clauses.append(lits)            # not all colour 0
        clauses.append([-l for l in lits])   # not all colour 1
    return idx, clauses


def solve(n, k, symbol_sym=False, timeout=None):
    idx, clauses = build(n, k, symbol_sym)
    if clauses is None:
        return idx, None, "forced monochromatic"
    t0 = time.time()
    with Minisat22(bootstrap_with=clauses) as s:
        ok = s.solve()
        model = s.get_model() if ok else None
    dt = time.time() - t0
    if ok:
        return idx, {abs(l): l > 0 for l in model}, f"SAT ({dt:.1f}s)"
    return idx, None, f"UNSAT ({dt:.1f}s)"


def verify(n, k, idx, assign, symbol_sym=False):
    """Expand back to all k^n points and check every line explicitly.  This is
    the step that catches an error in the orbit argument above, so run it on
    any n whose result you intend to claim."""
    def colour(x):
        comp = tuple(x.count(a) for a in range(k))
        key = tuple(sorted(comp, reverse=True)) if symbol_sym else comp
        return assign[idx[key]]
    bad = 0
    for tmpl in itertools.product(range(-1, k), repeat=n):
        if -1 not in tmpl:
            continue
        if len({colour(tuple(a if t == -1 else t for t in tmpl))
                for a in range(k)}) == 1:
            bad += 1
    return bad


def export(n, k, path):
    idx, clauses = build(n, k)
    if clauses is None:
        print(f"  n={n}: forced monochromatic, nothing to export")
        return
    with open(path, "w") as fh:
        fh.write(f"c S_n-invariant Hales-Jewett, k={k}, n={n}\n")
        fh.write(f"p cnf {len(idx)} {len(clauses)}\n")
        for c in clauses:
            fh.write(" ".join(map(str, c)) + " 0\n")
    print(f"  {path}: {len(idx)} vars, {len(clauses)} clauses")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", type=int, default=4, help="alphabet size")
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=9)
    ap.add_argument("--verify-upto", type=int, default=9,
                    help="expand and check against the full instance up to this n")
    ap.add_argument("--export", type=str, default="12,13,14,15",
                    help="comma-separated n values to write as DIMACS and exit")
    args = ap.parse_args()
    k = args.k

    if args.export:
        print(f"Exporting S_n-invariant instances, k={k}:")
        for n in [int(x) for x in args.export.split(",")]:
            export(n, k, f"sym_HJ_{k}_2_{n}.cnf")
        return

    print(f"k={k}   (colour depends only on the multiplicity vector)")
    print(f"{'n':>3} {'vars':>6} {'points':>12}  result")
    for n in range(args.min_n, args.max_n + 1):
        idx, assign, status = solve(n, k)
        note = ""
        if assign and n <= args.verify_upto:
            note = f"   verified: {verify(n, k, idx, assign)} mono lines in full instance"
        print(f"{n:>3} {len(idx):>6} {k**n:>12}  {status}{note}", flush=True)


if __name__ == "__main__":
    main()