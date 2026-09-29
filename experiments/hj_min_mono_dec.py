#!/usr/bin/env python3
"""Decision version of the minimum-monochromatic-lines question for 2-colourings of [3]^n:
    "is there a 2-colouring of [3]^n with at most B monochromatic combinatorial lines?"
UNSAT for B = D(n) - 1, where D(n) is the digit-sum colouring's count, proves that the
colouring chi(x) = floor(sum(x)/2) mod 2 is optimal at n.

Encoding: x_p (point p coloured 1), r_L >= [line L monochromatic] via (OR x_L | r_L) and
(OR -x_L | r_L), and sum r_L <= B (sequential counter).
Symmetry breaking (--sb, sound for the full group S_n x S_3 x colour swap): x_{0...0} = 0 and
x <=_lex g(x) for the adjacent coordinate swaps and the digit swaps (0 1), (1 2), all in the
lexicographic point order. The lex-least member of any orbit satisfies all of them, and the
number of monochromatic lines is invariant under the group, so the answer is unchanged.

Usage:
  python3 hj_min_mono_dec.py --n 6 --bound 134 --sb --kissat ~/kissat/build/kissat
  python3 hj_min_mono_dec.py --n 5 --bound 14 --sb --write-only f.cnf
"""
import argparse, itertools, os, re, subprocess, sys, tempfile, time
from pysat.card import CardEnc, EncType


def build(n, k=3):
    pts = list(itertools.product(range(k), repeat=n))
    idx = {p: i for i, p in enumerate(pts)}
    lines = [[idx[tuple(a if r == -1 else r for r in root)] for a in range(k)]
             for root in itertools.product(range(-1, k), repeat=n) if -1 in root]
    return pts, idx, lines


def lex_leader(x, y, top):
    """Clauses for x <=_lex y (lists of DIMACS literals, same length). Returns (clauses, top)."""
    cls, e = [], None                     # e = 'prefix so far equal' (None = true)
    for xi, yi in zip(x, y):
        if xi == yi:
            continue
        pre = [] if e is None else [-e]
        cls.append(pre + [-xi, yi])       # prefix equal -> x_i <= y_i
        top += 1; ne = top                # ne <- prefix equal and x_i == y_i
        cls.append(pre + [-xi, -yi, ne])
        cls.append(pre + [xi, yi, ne])
        e = ne
    return cls, top


def encode(n, bound, sb):
    pts, idx, lines = build(n)
    N, cls = len(pts), []
    R = list(range(N + 1, N + len(lines) + 1))
    for r, L in zip(R, lines):
        xs = [p + 1 for p in L]
        cls.append(xs + [r]); cls.append([-v for v in xs] + [r])
    top = N + len(lines)
    card = CardEnc.atmost(lits=R, bound=bound, top_id=top, encoding=EncType.seqcounter)
    cls += card.clauses; top = max(top, card.nv)
    if sb:
        cls.append([-1])                                   # point 0...0 has colour 0
        X = list(range(1, N + 1))
        perms = []
        for j in range(n - 1):                              # adjacent coordinate swaps
            perms.append([idx[p[:j] + (p[j + 1], p[j]) + p[j + 2:]] for p in pts])
        for a, b in ((0, 1), (1, 2)):                       # digit swaps
            sw = {a: b, b: a}
            perms.append([idx[tuple(sw.get(d, d) for d in p)] for p in pts])
        for pi in perms:                                    # involutions: (g x)_i = x_{pi(i)}
            c, top = lex_leader(X, [pi[i] + 1 for i in range(N)], top)
            cls += c
    return cls, top, N, len(lines)


def write_cnf(path, cls, nv):
    with open(path, "w") as fh:
        fh.write(f"p cnf {nv} {len(cls)}\n")
        fh.writelines(" ".join(map(str, c)) + " 0\n" for c in cls)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--bound", type=int, required=True)
    ap.add_argument("--sb", action="store_true", help="add symmetry-breaking constraints")
    ap.add_argument("--kissat", default=None, help="kissat binary (default: solve with CaDiCaL via PySAT)")
    ap.add_argument("--write-only", default=None, help="only write the CNF to this path")
    a = ap.parse_args()
    t0 = time.time()
    cls, nv, N, nl = encode(a.n, a.bound, a.sb)
    tag = f"[3]^{a.n}, at most {a.bound} monochromatic lines, symmetry breaking {'on' if a.sb else 'off'}"
    print(f"{tag}: {nv} variables, {len(cls)} clauses ({N} points, {nl} lines); encoded in {time.time()-t0:.1f}s", flush=True)
    if a.write_only:
        write_cnf(a.write_only, cls, nv); return
    t0 = time.time()
    if a.kissat:
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "f.cnf"); write_cnf(f, cls, nv)
            out = subprocess.run([os.path.expanduser(a.kissat), f], capture_output=True, text=True).stdout
        status = "SAT" if "s SATISFIABLE" in out else "UNSAT" if "s UNSATISFIABLE" in out else "UNKNOWN"
        m = re.search(r"^c conflicts:\s+(\d+)", out, re.M)
        conf = m.group(1) if m else "?"
    else:
        from pysat.solvers import Solver
        with Solver(name="cadical195", bootstrap_with=cls) as s:
            status = "SAT" if s.solve() else "UNSAT"; conf = s.accum_stats().get("conflicts", "?")
    res = f"{tag}: {status}; conflicts {conf}; {time.time()-t0:.1f}s"
    print(res, flush=True)
    os.makedirs("experiments/results", exist_ok=True)
    with open(f"experiments/results/min_mono_dec_n{a.n}_b{a.bound}_{'sb' if a.sb else 'nosb'}.txt", "w") as fh:
        fh.write(res + "\n")


if __name__ == "__main__":
    sys.exit(main())
