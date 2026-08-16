#!/usr/bin/env python3

"""
Random-variable-order experiment for Hales-Jewett CNFs.

Goal
----
Determine whether Kissat's dramatic improvement with --no-bump comes from

    (A) simply using a static decision order,

or

    (B) the PARTICULAR initial variable ordering in the Hales-Jewett CNF.

For each run, we randomly relabel every Boolean variable, producing an
isomorphic SAT instance.  We then solve that instance with:

    kissat --no-bump --no-tumble --reorder=0

The SAT problem is mathematically identical; only the variable labels/order
have changed.

All results are appended to a CSV file.  Runs are reproducible by seed and
the script can safely be restarted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
import re
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


# ============================================================
# Defaults for Colgate cluster
# ============================================================

DEFAULT_PROJECT = Path.home() / "CDCL-Experimentation"
DEFAULT_INSTANCE = DEFAULT_PROJECT / "cnfs" / "HJ_4_2_7.cnf"
DEFAULT_KISSAT = Path.home() / "kissat" / "build" / "kissat"

DEFAULT_RESULTS = (
    DEFAULT_PROJECT
    / "experiments"
    / "results"
    / "random_order_static_HJ_4_2_7.csv"
)

DEFAULT_BASE_SEED = 20260815
DEFAULT_RUNS = 20
DEFAULT_TIMEOUT = 120.0


# ============================================================
# Kissat output
# ============================================================

CONFLICT_RE = re.compile(r"c conflicts:\s+([0-9]+)")


# ============================================================
# DIMACS input/output
# ============================================================

def read_dimacs(path: Path):
    """
    Read a DIMACS CNF.

    Returns
    -------
    nvars : int
        Number of Boolean variables.

    clauses : list[list[int]]
        Clauses represented by signed DIMACS literals.
    """

    nvars = None
    expected_clauses = None

    clauses = []
    current_clause = []

    with path.open("r") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith("c"):
                continue

            if line.startswith("p"):
                parts = line.split()

                if len(parts) != 4 or parts[1] != "cnf":
                    raise ValueError(
                        f"Unexpected DIMACS header: {line}"
                    )

                nvars = int(parts[2])
                expected_clauses = int(parts[3])
                continue

            for token in line.split():
                literal = int(token)

                if literal == 0:
                    clauses.append(current_clause)
                    current_clause = []
                else:
                    current_clause.append(literal)

    if nvars is None:
        raise ValueError("No 'p cnf ...' header found.")

    if current_clause:
        raise ValueError(
            "DIMACS file ended before a clause-terminating 0."
        )

    if expected_clauses != len(clauses):
        raise ValueError(
            f"Header says {expected_clauses} clauses, "
            f"but parsed {len(clauses)}."
        )

    return nvars, clauses


def rename_literal(literal: int, permutation: list[int]) -> int:
    """
    Rename one signed DIMACS literal.

    permutation[v] is the new name for original variable v.
    """

    variable = abs(literal)
    new_variable = permutation[variable]

    if literal < 0:
        return -new_variable

    return new_variable


def write_permuted_dimacs(
    path: Path,
    nvars: int,
    clauses: list[list[int]],
    permutation: list[int],
):
    """
    Write an isomorphic CNF obtained by globally relabeling variables.
    """

    with path.open("w") as f:
        f.write(
            f"c Random variable relabeling experiment\n"
            f"p cnf {nvars} {len(clauses)}\n"
        )

        for clause in clauses:
            renamed = [
                rename_literal(lit, permutation)
                for lit in clause
            ]

            f.write(
                " ".join(str(lit) for lit in renamed)
                + " 0\n"
            )


# ============================================================
# Variable permutations
# ============================================================

def identity_permutation(nvars: int) -> list[int]:
    return list(range(nvars + 1))


def reverse_permutation(nvars: int) -> list[int]:
    """
    1 -> n
    2 -> n-1
    ...
    n -> 1
    """

    permutation = [0]

    permutation.extend(
        range(nvars, 0, -1)
    )

    return permutation


def random_permutation(
    nvars: int,
    seed: int,
) -> list[int]:

    rng = random.Random(seed)

    labels = list(range(1, nvars + 1))
    rng.shuffle(labels)

    return [0] + labels


def permutation_hash(
    permutation: list[int],
) -> str:
    """
    Short reproducibility fingerprint.
    """

    text = ",".join(
        str(x)
        for x in permutation[1:]
    )

    return hashlib.sha256(
        text.encode("utf-8")
    ).hexdigest()[:16]


# ============================================================
# Solver
# ============================================================

def run_kissat(
    kissat: Path,
    cnf: Path,
    timeout: float,
):
    """
    Run Kissat with bumping, tumbling, and reordering disabled.
    """

    command = [
        str(kissat),
        "--no-bump",
        "--no-tumble",
        "--reorder=0",
        str(cnf),
    ]

    start = time.perf_counter()

    try:
        process = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )

        elapsed = time.perf_counter() - start

        output = (
            process.stdout
            + "\n"
            + process.stderr
        )

        conflict_match = CONFLICT_RE.search(output)

        if conflict_match:
            conflicts = int(
                conflict_match.group(1)
            )
        else:
            conflicts = None

        if "s SATISFIABLE" in output:
            status = "SAT"

        elif "s UNSATISFIABLE" in output:
            status = "UNSAT"

        else:
            status = "UNKNOWN"

        return {
            "status": status,
            "conflicts": conflicts,
            "time_seconds": elapsed,
            "returncode": process.returncode,
        }

    except subprocess.TimeoutExpired:

        elapsed = time.perf_counter() - start

        return {
            "status": "TIMEOUT",
            "conflicts": None,
            "time_seconds": elapsed,
            "returncode": None,
        }


# ============================================================
# CSV tracking
# ============================================================

FIELDS = [
    "instance",
    "solver",
    "configuration",
    "ordering",
    "run",
    "seed",
    "permutation_hash",
    "time_seconds",
    "conflicts",
    "status",
    "returncode",
]


def initialize_results_file(path: Path):

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if path.exists():
        return

    with path.open(
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=FIELDS,
        )

        writer.writeheader()


def load_completed_runs(path: Path):
    """
    Allows experiment to resume after logout/crash.

    Returns set of keys:
        (ordering, run, seed)
    """

    completed = set()

    if not path.exists():
        return completed

    with path.open("r", newline="") as f:

        reader = csv.DictReader(f)

        for row in reader:

            completed.add(
                (
                    row["ordering"],
                    int(row["run"]),
                    row["seed"],
                )
            )

    return completed


def append_result(
    path: Path,
    row: dict,
):

    with path.open(
        "a",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=FIELDS,
        )

        writer.writerow(row)


# ============================================================
# Console summaries
# ============================================================

def print_current_summary(results_path: Path):

    if not results_path.exists():
        return

    random_conflicts = []
    random_times = []

    original_conflicts = None
    original_time = None

    with results_path.open("r") as f:

        reader = csv.DictReader(f)

        for row in reader:

            if row["status"] != "SAT":
                continue

            if row["conflicts"]:
                conflicts = int(row["conflicts"])
            else:
                conflicts = None

            runtime = float(
                row["time_seconds"]
            )

            if row["ordering"] == "original":

                original_conflicts = conflicts
                original_time = runtime

            elif row["ordering"] == "random":

                if conflicts is not None:
                    random_conflicts.append(
                        conflicts
                    )

                random_times.append(
                    runtime
                )

    print()
    print("=" * 60)
    print("CURRENT SUMMARY")
    print("=" * 60)

    if original_conflicts is not None:

        print(
            f"Original order: "
            f"{original_conflicts:,} conflicts, "
            f"{original_time:.3f} s"
        )

    if random_conflicts:

        print()
        print(
            f"Random permutations completed: "
            f"{len(random_conflicts)}"
        )

        print(
            f"Minimum conflicts: "
            f"{min(random_conflicts):,}"
        )

        print(
            f"Median conflicts:  "
            f"{statistics.median(random_conflicts):,.1f}"
        )

        print(
            f"Mean conflicts:    "
            f"{statistics.mean(random_conflicts):,.1f}"
        )

        print(
            f"Maximum conflicts: "
            f"{max(random_conflicts):,}"
        )

        print(
            f"Median time:       "
            f"{statistics.median(random_times):.3f} s"
        )

    print("=" * 60)
    print()


# ============================================================
# One experiment
# ============================================================

def run_one_ordering(
    *,
    ordering_name,
    run_number,
    seed,
    permutation,
    nvars,
    clauses,
    kissat,
    instance,
    results,
    timeout,
    completed,
):

    seed_string = (
        ""
        if seed is None
        else str(seed)
    )

    key = (
        ordering_name,
        run_number,
        seed_string,
    )

    if key in completed:

        print(
            f"Skipping already-completed "
            f"{ordering_name} run {run_number}"
        )

        return

    perm_hash = permutation_hash(
        permutation
    )

    print(
        f"{ordering_name.upper():8s} "
        f"run={run_number:3d} "
        f"seed={seed_string or '-':10s} "
        f"hash={perm_hash}"
    )

    with tempfile.NamedTemporaryFile(
        suffix=".cnf",
        delete=False,
    ) as temporary_file:

        temporary_path = Path(
            temporary_file.name
        )

    try:

        write_permuted_dimacs(
            temporary_path,
            nvars,
            clauses,
            permutation,
        )

        result = run_kissat(
            kissat,
            temporary_path,
            timeout,
        )

        conflicts_text = (
            "None"
            if result["conflicts"] is None
            else f'{result["conflicts"]:,}'
        )

        print(
            f"    status={result['status']:7s} "
            f"time={result['time_seconds']:8.3f}s "
            f"conflicts={conflicts_text}"
        )

        row = {
            "instance": instance.name,
            "solver": "kissat",
            "configuration": "no_bump_no_tumble_reorder_0",
            "ordering": ordering_name,
            "run": run_number,
            "seed": seed_string,
            "permutation_hash": perm_hash,
            "time_seconds":
                f'{result["time_seconds"]:.6f}',
            "conflicts":
                ""
                if result["conflicts"] is None
                else result["conflicts"],
            "status": result["status"],
            "returncode":
                ""
                if result["returncode"] is None
                else result["returncode"],
        }

        append_result(
            results,
            row,
        )

        completed.add(key)

    finally:

        temporary_path.unlink(
            missing_ok=True
        )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--instance",
        type=Path,
        default=DEFAULT_INSTANCE,
    )

    parser.add_argument(
        "--kissat",
        type=Path,
        default=DEFAULT_KISSAT,
    )

    parser.add_argument(
        "--runs",
        type=int,
        default=DEFAULT_RUNS,
        help="Number of random permutations.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_BASE_SEED,
        help="Seed for first random permutation.",
    )

    parser.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_TIMEOUT,
        help="Timeout per solver run in seconds.",
    )

    parser.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_RESULTS,
    )

    parser.add_argument(
        "--skip-reverse",
        action="store_true",
        help="Do not run deterministic reversed ordering.",
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # Sanity checks
    # --------------------------------------------------------

    if not args.instance.exists():
        raise FileNotFoundError(
            f"CNF not found: {args.instance}"
        )

    if not args.kissat.exists():
        raise FileNotFoundError(
            f"Kissat not found: {args.kissat}"
        )

    print()
    print("Hales-Jewett static-order experiment")
    print("=" * 60)

    print(
        f"Instance: {args.instance}"
    )

    print(
        f"Kissat:   {args.kissat}"
    )

    print(
        f"Runs:     {args.runs}"
    )

    print(
        f"Timeout:  {args.timeout}s"
    )

    print(
        f"Results:  {args.results}"
    )

    print()

    # --------------------------------------------------------
    # Parse CNF once
    # --------------------------------------------------------

    nvars, clauses = read_dimacs(
        args.instance
    )

    print(
        f"Variables: {nvars:,}"
    )

    print(
        f"Clauses:   {len(clauses):,}"
    )

    print()

    initialize_results_file(
        args.results
    )

    completed = load_completed_runs(
        args.results
    )

    # --------------------------------------------------------
    # Original variable ordering
    # --------------------------------------------------------

    run_one_ordering(
        ordering_name="original",
        run_number=0,
        seed=None,
        permutation=identity_permutation(nvars),
        nvars=nvars,
        clauses=clauses,
        kissat=args.kissat,
        instance=args.instance,
        results=args.results,
        timeout=args.timeout,
        completed=completed,
    )

    # --------------------------------------------------------
    # Reverse variable ordering
    # --------------------------------------------------------

    if not args.skip_reverse:

        run_one_ordering(
            ordering_name="reverse",
            run_number=0,
            seed=None,
            permutation=reverse_permutation(
                nvars
            ),
            nvars=nvars,
            clauses=clauses,
            kissat=args.kissat,
            instance=args.instance,
            results=args.results,
            timeout=args.timeout,
            completed=completed,
        )

    # --------------------------------------------------------
    # Random variable orderings
    # --------------------------------------------------------

    for run_number in range(
        1,
        args.runs + 1,
    ):

        seed = (
            args.seed
            + run_number
            - 1
        )

        permutation = random_permutation(
            nvars,
            seed,
        )

        run_one_ordering(
            ordering_name="random",
            run_number=run_number,
            seed=seed,
            permutation=permutation,
            nvars=nvars,
            clauses=clauses,
            kissat=args.kissat,
            instance=args.instance,
            results=args.results,
            timeout=args.timeout,
            completed=completed,
        )

        print_current_summary(
            args.results
        )

    print()
    print("Experiment complete.")

    print_current_summary(
        args.results
    )


if __name__ == "__main__":
    main()
