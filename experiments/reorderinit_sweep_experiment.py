#!/usr/bin/env python3
"""Map Kissat's structural-reorder response curve on Hales-Jewett instances.

Two arms:

  sweep   --no-bump --reorderinit=X  over a grid of X, several clause-order
          permutations per point.  Produces a response curve against the
          ~2% noise floor already measured in random_order_HJ_4_2_7.csv.

  tumble  isolation of the 12x reorder effect.  The existing variance data
          compares `--no-bump` against `--no-bump --no-tumble --reorder=0`,
          which differ in two options at once.  The middle arm
          `--no-bump --reorder=0` (tumble left on) is the missing cell.

Per run we record total conflicts, every reorder event's conflict count and
search mode, decisions, decisions-per-conflict, restarts, rephases and mode
switches.  decisions-per-conflict was 7.40 on the good reorderinit=3000 run
and 11.70 on the poor reorderinit=7000 run, so it is worth carrying through.

Wall clock is not a useful figure of merit here: fitting the two runs in
reorder_explainability_HJ_4_2_7 gives ~1.8s of fixed parse and preprocess
against ~0.12s and ~0.66s of actual search.  Kissat's own process-time is
parsed where available; otherwise compare conflicts.

Usage:
    python3 experiments/reorderinit_sweep_experiment.py --dry-run
    python3 experiments/reorderinit_sweep_experiment.py
    python3 experiments/reorderinit_sweep_experiment.py --arms sweep
    python3 experiments/reorderinit_sweep_experiment.py --cnf cnfs/HJ_3_3_9.cnf
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# --------------------------------------------------------------------------
# Defaults.  Paths match the layout on s03.
# --------------------------------------------------------------------------

DEFAULT_KISSAT = Path.home() / "kissat" / "build" / "kissat"
DEFAULT_REPO = Path.home() / "CDCL-Experimentation"
DEFAULT_CNF = DEFAULT_REPO / "cnfs" / "HJ_4_2_7.cnf"
DEFAULT_RESULTS = DEFAULT_REPO / "experiments" / "results"

# 25 thresholds.  Dense below 5000 because Table 5's minimum (5,187 conflicts
# at reorderinit=100) sits there; Kissat's own default is 1e4.
REORDERINIT_GRID = [
    100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500,
    5000, 5500, 6000, 6500, 7000, 8000, 9000, 10000, 11000, 12000, 13000,
    14000, 15000,
]

# Matches the convention in random_order_HJ_4_2_7.csv so the two files merge.
FIRST_SEED = 20260815

TUMBLE_ARMS = [
    ("no_bump", ["--no-bump"]),
    ("no_bump_reorder_0", ["--no-bump", "--reorder=0"]),                    # the missing cell
    ("no_bump_no_tumble_reorder_0", ["--no-bump", "--no-tumble", "--reorder=0"]),
]

SAT, UNSAT = 10, 20


# --------------------------------------------------------------------------
# CNF handling
# --------------------------------------------------------------------------

def read_dimacs(path: Path):
    """Return (n_vars, n_clauses, clauses).  Token-based, so clauses that wrap
    across lines still parse correctly."""
    n_vars = n_clauses = None
    tokens = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("c"):
                continue
            if line.startswith("p"):
                parts = line.split()
                n_vars, n_clauses = int(parts[2]), int(parts[3])
                continue
            tokens.extend(line.split())

    clauses, current = [], []
    for tok in tokens:
        lit = int(tok)
        if lit == 0:
            clauses.append(current)
            current = []
        else:
            current.append(lit)
    if current:
        clauses.append(current)

    if n_clauses is not None and len(clauses) != n_clauses:
        print(f"  warning: header says {n_clauses} clauses, parsed {len(clauses)}",
              file=sys.stderr)
    return n_vars, len(clauses), clauses


def permutation_hash(clauses) -> str:
    """16 hex chars over the clause sequence.

    NOTE: this will not reproduce the permutation_hash values already in
    random_order_HJ_4_2_7.csv unless it happens to use the same function.
    If you need the two files to join on that column, copy the hashing
    routine out of whichever script produced that CSV and drop it in here.
    The seed column joins correctly either way.
    """
    h = hashlib.blake2b(digest_size=8)
    for clause in clauses:
        h.update((" ".join(map(str, clause)) + "\n").encode())
    return h.hexdigest()


def write_dimacs(path: Path, n_vars: int, clauses) -> None:
    with path.open("w") as fh:
        fh.write(f"p cnf {n_vars} {len(clauses)}\n")
        for clause in clauses:
            fh.write(" ".join(map(str, clause)) + " 0\n")


def build_permutations(cnf: Path, workdir: Path, n_random: int):
    """Materialise each clause ordering once and reuse it across every
    configuration.  Shuffling clause order only -- never variable numbering,
    which would destroy the incidence ranking we are measuring."""
    n_vars, n_clauses, clauses = read_dimacs(cnf)
    print(f"Instance: {cnf} ({n_vars:,} variables, {n_clauses:,} clauses)")

    variants = []

    original = list(clauses)
    path = workdir / "perm_original.cnf"
    write_dimacs(path, n_vars, original)
    variants.append({"ordering": "original", "run": 0, "seed": "",
                     "permutation_hash": permutation_hash(original), "path": path})

    for i in range(n_random):
        seed = FIRST_SEED + i
        shuffled = list(clauses)
        random.Random(seed).shuffle(shuffled)
        path = workdir / f"perm_{seed}.cnf"
        write_dimacs(path, n_vars, shuffled)
        variants.append({"ordering": "random", "run": i + 1, "seed": seed,
                         "permutation_hash": permutation_hash(shuffled), "path": path})

    return n_vars, n_clauses, variants


# --------------------------------------------------------------------------
# Running and parsing
# --------------------------------------------------------------------------

STAT_PATTERNS = {
    "conflicts": re.compile(r"^c conflicts:\s+(\d+)"),
    "decisions": re.compile(r"^c decisions:\s+(\d+)"),
    "propagations": re.compile(r"^c propagations:\s+(\d+)"),
    "restarts": re.compile(r"^c restarts:\s+(\d+)"),
    "rephased": re.compile(r"^c rephased:\s+(\d+)"),
    "switched": re.compile(r"^c switched:\s+(\d+)"),
    "reordered": re.compile(r"^c reordered:\s+(\d+)"),
}
PROCESS_TIME = re.compile(r"^c process-time:.*?([\d.]+)\s*seconds")
EVENT = re.compile(r"^c reorder-weight event=(\d+) conflicts=(\d+) mode=(\w+)")


def parse_output(text: str) -> dict:
    out = {k: "" for k in STAT_PATTERNS}
    out["process_time"] = ""
    events, modes, seen = [], [], set()

    for line in text.splitlines():
        for key, pat in STAT_PATTERNS.items():
            if out[key] == "":
                m = pat.match(line)
                if m:
                    out[key] = int(m.group(1))
        if out["process_time"] == "":
            m = PROCESS_TIME.match(line)
            if m:
                out["process_time"] = float(m.group(1))
        m = EVENT.match(line)
        if m:
            ev = int(m.group(1))
            if ev not in seen:          # --reorderlog=1 gives one row per event
                seen.add(ev)
                events.append(int(m.group(2)))
                modes.append(m.group(3))

    out["event_conflicts"] = ";".join(map(str, events))
    out["event_modes"] = ";".join(modes)
    out["first_reorder_conflict"] = events[0] if events else ""
    out["logged_events"] = len(events)
    return out


def run_once(kissat: Path, args, cnf: Path, timeout: int) -> dict:
    cmd = [str(kissat)] + args + ["--reorderlog=1", str(cnf)]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        elapsed = time.perf_counter() - start
        stdout, rc = proc.stdout, proc.returncode
        status = {SAT: "SAT", UNSAT: "UNSAT"}.get(rc, f"OTHER_{rc}")
    except subprocess.TimeoutExpired:
        elapsed = time.perf_counter() - start
        stdout, rc, status = "", None, "TIMEOUT"

    record = {"status": status, "returncode": rc if rc is not None else "",
              "wall_seconds": round(elapsed, 6)}
    record.update(parse_output(stdout))

    conflicts, decisions = record.get("conflicts"), record.get("decisions")
    if isinstance(conflicts, int) and isinstance(decisions, int) and conflicts:
        record["decisions_per_conflict"] = round(decisions / conflicts, 4)
    else:
        record["decisions_per_conflict"] = ""
    record["command"] = " ".join(cmd)
    return record


# --------------------------------------------------------------------------
# Job construction
# --------------------------------------------------------------------------

def build_jobs(args, variants):
    jobs = []
    if "sweep" in args.arms:
        for init in args.grid:
            for v in variants[: args.replicates + 1]:
                jobs.append({"arm": "sweep", "configuration": f"no_bump_reorderinit_{init}",
                             "reorderinit": init, "solver_args": ["--no-bump", f"--reorderinit={init}"],
                             "variant": v})
    if "tumble" in args.arms:
        for name, solver_args in TUMBLE_ARMS:
            for v in variants[: args.tumble_replicates + 1]:
                jobs.append({"arm": "tumble", "configuration": name, "reorderinit": "",
                             "solver_args": list(solver_args), "variant": v})
    return jobs


FIELDS = ["timestamp_utc", "arm", "configuration", "instance", "solver", "solver_version",
          "reorderinit", "ordering", "run", "seed", "permutation_hash", "status",
          "returncode", "wall_seconds", "process_time", "conflicts", "decisions",
          "decisions_per_conflict", "propagations", "restarts", "rephased", "switched",
          "reordered", "logged_events", "first_reorder_conflict", "event_conflicts",
          "event_modes", "command"]


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

def mann_whitney(a, b):
    """Two-sided rank-sum with tie correction and normal approximation.
    Uses scipy when present.  The approximation is unreliable below ~8 per
    group, so it is reported as an indication, not a p-value to quote."""
    if not a or not b:
        return ""
    try:
        from scipy.stats import mannwhitneyu
        return round(float(mannwhitneyu(a, b, alternative="two-sided").pvalue), 8)
    except Exception:
        pass

    import math
    combined = sorted([(x, 0) for x in a] + [(x, 1) for x in b])
    ranks, i, n = [0.0] * len(combined), 0, len(combined)
    tie_term = 0.0
    while i < n:
        j = i
        while j + 1 < n and combined[j + 1][0] == combined[i][0]:
            j += 1
        avg = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[k] = avg
        t = j - i + 1
        if t > 1:
            tie_term += t ** 3 - t
        i = j + 1

    n1, n2 = len(a), len(b)
    r1 = sum(r for r, (_, g) in zip(ranks, combined) if g == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    var = n1 * n2 * (n + 1) / 12.0 - n1 * n2 * tie_term / (12.0 * n * (n - 1))
    if var <= 0:
        return ""
    z = (abs(u1 - mu) - 0.5) / math.sqrt(var)
    return round(math.erfc(z / math.sqrt(2)), 8)


def summarise(rows, out_path: Path) -> None:
    groups = {}
    for r in rows:
        if r["status"] != "SAT" or not isinstance(r["conflicts"], int):
            continue
        groups.setdefault((r["arm"], r["configuration"], r["reorderinit"]), []).append(r)

    baseline = next((v for k, v in groups.items() if k[2] == 10000), None)
    base_conf = [r["conflicts"] for r in baseline] if baseline else []

    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["arm", "configuration", "reorderinit", "n", "timeouts",
                    "median_conflicts", "iqr_conflicts", "min_conflicts", "max_conflicts",
                    "median_decisions_per_conflict", "median_first_reorder_conflict",
                    "median_reorders", "p_vs_default_10000"])
        for key in sorted(groups, key=lambda k: (k[0], k[2] if k[2] != "" else 0, k[1])):
            arm, config, init = key
            g = groups[key]
            conf = sorted(r["conflicts"] for r in g)
            dpc = [r["decisions_per_conflict"] for r in g if r["decisions_per_conflict"] != ""]
            first = [r["first_reorder_conflict"] for r in g
                     if isinstance(r["first_reorder_conflict"], int)]
            reord = [r["reordered"] for r in g if isinstance(r["reordered"], int)]
            timeouts = sum(1 for r in rows
                           if (r["arm"], r["configuration"], r["reorderinit"]) == key
                           and r["status"] == "TIMEOUT")
            q1, q3 = (statistics.quantiles(conf, n=4)[0], statistics.quantiles(conf, n=4)[2]) \
                if len(conf) >= 4 else (conf[0], conf[-1])
            p = mann_whitney(conf, base_conf) if base_conf and init != 10000 else ""
            w.writerow([arm, config, init, len(conf), timeouts,
                        statistics.median(conf), round(q3 - q1, 1), conf[0], conf[-1],
                        round(statistics.median(dpc), 3) if dpc else "",
                        statistics.median(first) if first else "",
                        statistics.median(reord) if reord else "", p])

    print(f"\nSummary written to {out_path}")
    print("\nReminder: the noise floor from random_order_HJ_4_2_7.csv is about +/-2%")
    print("at the default reorderinit=1e4.  Treat any spread inside that band as noise.")


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kissat", type=Path, default=DEFAULT_KISSAT)
    ap.add_argument("--cnf", type=Path, default=DEFAULT_CNF)
    ap.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--arms", default="sweep,tumble",
                    help="comma-separated: sweep, tumble")
    ap.add_argument("--replicates", type=int, default=4,
                    help="random permutations per sweep point (plus the original ordering)")
    ap.add_argument("--tumble-replicates", type=int, default=21,
                    help="random permutations per tumble arm (plus the original ordering)")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--grid", type=str, default="",
                    help="comma-separated reorderinit values; default is the built-in grid")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    args.arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    args.grid = [int(x) for x in args.grid.split(",")] if args.grid else REORDERINIT_GRID

    if not args.kissat.exists():
        print(f"error: kissat not found at {args.kissat}", file=sys.stderr)
        return 1
    if not args.cnf.exists():
        print(f"error: CNF not found at {args.cnf}", file=sys.stderr)
        return 1

    version = subprocess.run([str(args.kissat), "--version"],
                             capture_output=True, text=True).stdout.strip()
    print("Kissat reorderinit sweep")
    print(f"Kissat:   {args.kissat} ({version})")

    workdir = Path(tempfile.mkdtemp(prefix="reorderinit_sweep_"))
    try:
        n_random = max(args.replicates, args.tumble_replicates)
        _, _, variants = build_permutations(args.cnf, workdir, n_random)
        jobs = build_jobs(args, variants)

        print(f"Arms:     {', '.join(args.arms)}")
        print(f"Grid:     {len(args.grid)} thresholds, {args.replicates + 1} orderings each")
        print(f"Jobs:     {len(jobs)} runs | timeout {args.timeout}s per run")

        if args.dry_run:
            print("\nFirst and last few commands:")
            for j in jobs[:3] + jobs[-3:]:
                print("  " + " ".join([str(args.kissat)] + j["solver_args"]
                                      + ["--reorderlog=1", str(j["variant"]["path"])]))
            return 0

        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        outdir = args.results_dir / f"reorderinit_sweep_{args.cnf.stem}_{stamp}"
        outdir.mkdir(parents=True, exist_ok=True)
        csv_path = outdir / "sweep_runs.csv"

        rows = []
        with csv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDS)
            writer.writeheader()
            for i, job in enumerate(jobs, 1):
                v = job["variant"]
                rec = run_once(args.kissat, job["solver_args"], v["path"], args.timeout)
                row = {
                    "timestamp_utc": dt.datetime.now(dt.timezone.utc)
                                       .strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "arm": job["arm"], "configuration": job["configuration"],
                    "instance": args.cnf.name, "solver": "kissat", "solver_version": version,
                    "reorderinit": job["reorderinit"], "ordering": v["ordering"],
                    "run": v["run"], "seed": v["seed"],
                    "permutation_hash": v["permutation_hash"],
                }
                row.update({k: rec.get(k, "") for k in FIELDS if k not in row})
                writer.writerow(row)
                fh.flush()                       # survive a crash mid-sweep
                rows.append(row)
                print(f"[{i}/{len(jobs)}] {job['configuration']:<34} "
                      f"{v['ordering']}/{v['run']:<3} {rec['status']:<8} "
                      f"conflicts {rec.get('conflicts', '?')!s:>8}  "
                      f"d/c {rec.get('decisions_per_conflict', '?')!s:>7}  "
                      f"reorders@ {rec.get('event_conflicts', '') or '-'}")

        summarise(rows, outdir / "sweep_summary.csv")
        print(f"\nRun data: {csv_path}")
        return 0
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())