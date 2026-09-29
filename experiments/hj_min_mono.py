#!/usr/bin/env python3
"""Exact minimum number of monochromatic combinatorial lines over all 2-colourings of [3]^n
(MaxSAT: one soft constraint per line, solved with RC2 + CaDiCaL). Also reports how many
monochromatic lines the digit-sum colouring chi(x) = floor(sum(x)/2) mod 2 has, and prints an
optimal colouring's summary. Usage: python3 hj_min_mono.py N   (writes results/min_mono_3_N.txt)"""
import itertools, os, sys, time
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2

def build(n, k=3):
    pts = list(itertools.product(range(k), repeat=n)); idx = {p: i for i, p in enumerate(pts)}
    lines = [[idx[tuple(a if r == -1 else r for r in root)] for a in range(k)]
             for root in itertools.product(range(-1, k), repeat=n) if -1 in root]
    return pts, lines

def main():
    n = int(sys.argv[1]); pts, lines = build(n); N = len(pts)
    chi = sum(1 for L in lines if len({(sum(pts[q]) // 2) % 2 for q in L}) == 1)
    print(f"[3]^{n}: {N} points, {len(lines)} lines; digit-sum colouring has {chi} monochromatic lines", flush=True)
    w = WCNF(); top = N
    for L in lines:
        top += 1; xs = [p + 1 for p in L]
        w.append(xs + [top]); w.append([-x for x in xs] + [top]); w.append([-top], weight=1)
    t = time.time()
    with RC2(w, solver="cd19", verbose=1) as rc2:   # prints "c cost: X" each time the lower bound rises
        mdl = rc2.compute(); opt = rc2.cost
    col = [1 if mdl[i] > 0 else 0 for i in range(N)]
    mono = [L for L in lines if len({col[q] for q in L}) == 1]
    out = (f"[3]^{n}: MINIMUM monochromatic lines = {opt} (verified {len(mono)}); digit-sum colouring: {chi}; "
           f"{time.time() - t:.1f}s\n#true = {sum(col)} of {N}\ncolouring (points in lexicographic order): "
           f"{''.join(map(str, col))}\n")
    print(out, flush=True)
    os.makedirs("experiments/results", exist_ok=True)
    open(f"experiments/results/min_mono_3_{n}.txt", "w").write(out)

if __name__ == "__main__":
    main()
