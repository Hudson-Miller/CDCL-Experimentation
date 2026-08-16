#!/usr/bin/env python3
"""Run a fixed Kissat ablation across Hales--Jewett CNF instances.

The script discovers ``cnfs/HJ_*.cnf`` relative to the repository root,
runs four fixed Kissat configurations, and appends one CSV row immediately
after each run. Existing ``(instance, configuration)`` pairs are skipped,
including runs whose recorded result is TIMEOUT or UNKNOWN.

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
from typing import Iterable, Optional, Sequence, TextIO


DEFAULT_TIMEOUT_SECONDS = 300.0


@dataclass(frozen=True)
class Configuration:
    name: str
    options: tuple[str, ...]


CONFIGURATIONS = (
    Configuration("baseline", ()),
    Configuration("no_bump", ("--no-bump",)),
    Configuration("no_bump_no_reorder", ("--no-bump", "--reorder=0")),
    Configuration("no_bump_stable_reorder", ("--no-bump", "--reorder=1")),
)


CSV_FIELDS = (
    "instance",
    "hj_family",
    "hj_parameters",
    "hj_k",
    "hj_r",
    "hj_n",
    "variables",
    "clauses",
    "configuration",
    "wall_time_seconds",
    "conflicts",
    "number_of_reorders",
    "status",
    "return_code",
    "timeout",
    "timestamp_utc",
)


CONFLICTS_RE = re.compile(
    r"^\s*c\s+conflicts\s*:\s*([0-9][0-9,]*)\b",
    re.IGNORECASE | re.MULTILINE,
)
REORDERS_RE = re.compile(
    r"^\s*c\s+reorder(?:ed|s)\s*:\s*([0-9][0-9,]*)\b",
    re.IGNORECASE | re.MULTILINE,
)


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repository_root = script_path.parent.parent

    parser = argparse.ArgumentParser(
        description=(
            "Run four fixed Kissat configurations over cnfs/HJ_*.cnf and "
            "append resumable results to a CSV file."
        )
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_TIMEOUT_SECONDS,
        metavar="SECONDS",
        help="per-run timeout in seconds (default: 300)",
    )
    parser.add_argument(
        "--filter",
        action="append",
        default=[],
        metavar="PATTERN",
        help=(
            "only run filenames matching this shell-style glob; a value with "
            "no glob characters is treated as a substring (repeatable)"
        ),
    )
    parser.add_argument(
        "--cnf-dir",
        type=Path,
        default=repository_root / "cnfs",
        metavar="PATH",
        help="CNF directory (default: <repository>/cnfs)",
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=repository_root
        / "experiments"
        / "results"
        / "hj_family_benchmark.csv",
        metavar="PATH",
        help=(
            "append-only result CSV (default: "
            "<repository>/experiments/results/hj_family_benchmark.csv)"
        ),
    )
    parser.add_argument(
        "--kissat",
        metavar="PATH",
        help=(
            "Kissat executable; otherwise use KISSAT_BIN, "
            "~/kissat/build/kissat, or kissat from PATH"
        ),
    )

    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error("--timeout must be greater than zero")
    return args


def resolve_executable(explicit: Optional[str]) -> str:
    """Return an executable Kissat path or raise a useful error."""

    requested = explicit or os.environ.get("KISSAT_BIN")
    candidates: list[str] = []

    if requested:
        candidates.append(requested)
    else:
        candidates.extend(("~/kissat/build/kissat", "kissat"))

    checked: list[str] = []

    for candidate in candidates:
        expanded = os.path.expanduser(candidate)
        checked.append(expanded)

        if os.path.sep in expanded or (
            os.path.altsep is not None and os.path.altsep in expanded
        ):
            path = Path(expanded)
            if path.is_file() and os.access(path, os.X_OK):
                return str(path.resolve())
            continue

        located = shutil.which(expanded)
        if located:
            return located

    rendered = ", ".join(repr(item) for item in checked)
    raise FileNotFoundError(
        "Could not find an executable Kissat binary. Checked: "
        f"{rendered}. Pass --kissat PATH or set KISSAT_BIN."
    )


def natural_sort_key(path: Path) -> list[object]:
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", path.name)
    ]


def matches_filters(filename: str, patterns: Iterable[str]) -> bool:
    patterns = list(patterns)
    if not patterns:
        return True

    for pattern in patterns:
        if any(character in pattern for character in "*?["):
            if fnmatch.fnmatchcase(filename, pattern):
                return True
        elif pattern in filename:
            return True

    return False


def discover_instances(cnf_dir: Path, patterns: Iterable[str]) -> list[Path]:
    if not cnf_dir.is_dir():
        raise FileNotFoundError(f"CNF directory does not exist: {cnf_dir}")

    instances = [
        path.resolve()
        for path in cnf_dir.glob("HJ_*.cnf")
        if path.is_file() and matches_filters(path.name, patterns)
    ]
    return sorted(instances, key=natural_sort_key)


def parse_dimacs_header(path: Path) -> tuple[int, int]:
    """Read the variable and clause counts from a DIMACS p cnf line."""

    with path.open("r", encoding="ascii", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("c"):
                continue

            parts = stripped.split()

            if (
                len(parts) >= 4
                and parts[0] == "p"
                and parts[1].lower() == "cnf"
            ):
                try:
                    return int(parts[2]), int(parts[3])
                except ValueError as error:
                    raise ValueError(
                        f"Invalid DIMACS counts in {path} at line {line_number}"
                    ) from error

    raise ValueError(
        f"No valid 'p cnf <variables> <clauses>' header in {path}"
    )


def parse_hj_filename(path: Path) -> dict[str, object]:
    """Parse common names such as HJ_4_2_7.cnf into HJ parameters."""

    suffix = path.stem[3:] if path.stem.startswith("HJ_") else path.stem

    def named_value(letter: str) -> Optional[int]:
        match = re.search(
            rf"(?:^|[_\-;,]){letter}\s*=?\s*(\d+)(?:$|[_\-;,])",
            suffix,
            re.IGNORECASE,
        )
        return int(match.group(1)) if match else None

    k = named_value("k")
    r = named_value("r")
    n = named_value("n")

    # The repository's established form is HJ_k_r_n.cnf.
    numbers = [int(value) for value in re.findall(r"\d+", suffix)]

    if k is None and len(numbers) >= 1:
        k = numbers[0]
    if r is None and len(numbers) >= 2:
        r = numbers[1]
    if n is None and len(numbers) >= 3:
        n = numbers[2]

    family = f"HJ({k};{r})" if k is not None and r is not None else "HJ"

    return {
        "hj_family": family,
        "hj_parameters": suffix,
        "hj_k": "" if k is None else k,
        "hj_r": "" if r is None else r,
        "hj_n": "" if n is None else n,
    }


def parse_counter(
    pattern: re.Pattern[str],
    output: str,
) -> Optional[int]:
    matches = pattern.findall(output)
    if not matches:
        return None

    return int(matches[-1].replace(",", ""))


def parse_status(
    output: str,
    return_code: int,
    timed_out: bool,
) -> str:
    if timed_out:
        return "TIMEOUT"

    if re.search(
        r"^\s*s\s+UNSATISFIABLE\b",
        output,
        re.IGNORECASE | re.MULTILINE,
    ):
        return "UNSAT"

    if re.search(
        r"^\s*s\s+SATISFIABLE\b",
        output,
        re.IGNORECASE | re.MULTILINE,
    ):
        return "SAT"

    if re.search(
        r"^\s*s\s+UNKNOWN\b",
        output,
        re.IGNORECASE | re.MULTILINE,
    ):
        return "UNKNOWN"

    # Standard SAT competition return codes.
    if return_code == 10:
        return "SAT"
    if return_code == 20:
        return "UNSAT"

    return "UNKNOWN"


def kill_process_group(process: subprocess.Popen[str]) -> None:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    except PermissionError:
        try:
            process.kill()
        except ProcessLookupError:
            return


def run_kissat(
    executable: str,
    configuration: Configuration,
    instance: Path,
    timeout_seconds: float,
) -> dict[str, object]:
    command = [
        executable,
        *configuration.options,
        str(instance),
    ]

    started_at = datetime.now(timezone.utc)
    started = time.monotonic()

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
        start_new_session=True,
    )

    timed_out = False

    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        timed_out = True
        kill_process_group(process)
        stdout, stderr = process.communicate()
    except KeyboardInterrupt:
        kill_process_group(process)
        process.communicate()
        raise

    wall_time = time.monotonic() - started
    output = "\n".join((stdout, stderr))
    return_code = process.returncode
    status = parse_status(output, return_code, timed_out)

    return {
        "wall_time_seconds": f"{wall_time:.6f}",
        "conflicts": parse_counter(CONFLICTS_RE, output),
        "number_of_reorders": parse_counter(REORDERS_RE, output),
        "status": status,
        "return_code": return_code,
        "timeout": "true" if timed_out else "false",
        "timestamp_utc": started_at.isoformat(timespec="seconds"),
    }


def read_completed_runs(
    results_path: Path,
) -> set[tuple[str, str]]:
    if not results_path.exists() or results_path.stat().st_size == 0:
        return set()

    with results_path.open(
        "r",
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []

        if tuple(fields) != CSV_FIELDS:
            raise ValueError(
                f"Refusing to append to {results_path}: its CSV header does "
                "not match this script. Existing data was left unchanged. "
                "Choose a new file with --results PATH."
            )

        completed: set[tuple[str, str]] = set()

        for row in reader:
            instance = (row.get("instance") or "").strip()
            configuration = (row.get("configuration") or "").strip()

            if instance and configuration:
                completed.add((instance, configuration))

        return completed


def open_results_writer(
    results_path: Path,
) -> tuple[TextIO, csv.DictWriter]:
    results_path.parent.mkdir(parents=True, exist_ok=True)

    needs_header = (
        not results_path.exists()
        or results_path.stat().st_size == 0
    )

    handle = results_path.open(
        "a",
        newline="",
        encoding="utf-8",
    )
    writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)

    if needs_header:
        writer.writeheader()
        handle.flush()
        os.fsync(handle.fileno())

    return handle, writer


def persist_row(
    handle: TextIO,
    writer: csv.DictWriter,
    row: dict[str, object],
) -> None:
    writer.writerow(
        {field: row.get(field, "") for field in CSV_FIELDS}
    )
    handle.flush()
    os.fsync(handle.fileno())


def printable_count(value: object) -> str:
    if value is None or value == "":
        return "?"

    return f"{int(value):,}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)

    try:
        executable = resolve_executable(args.kissat)
        cnf_dir = args.cnf_dir.expanduser().resolve()
        results_path = args.results.expanduser().resolve()
        instances = discover_instances(cnf_dir, args.filter)
        completed = read_completed_runs(results_path)
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if not instances:
        filter_note = (
            f" matching {args.filter!r}"
            if args.filter
            else ""
        )
        print(
            f"No HJ_*.cnf files found in "
            f"{cnf_dir}{filter_note}."
        )
        return 0

    metadata: dict[Path, dict[str, object]] = {}

    try:
        for instance in instances:
            variables, clauses = parse_dimacs_header(instance)
            metadata[instance] = {
                **parse_hj_filename(instance),
                "variables": variables,
                "clauses": clauses,
            }
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    pending = [
        (instance, configuration)
        for instance in instances
        for configuration in CONFIGURATIONS
        if (instance.name, configuration.name) not in completed
    ]

    total_pairs = len(instances) * len(CONFIGURATIONS)
    skipped = total_pairs - len(pending)

    print(f"Kissat: {executable}")
    print(
        f"Instances: {len(instances)} | "
        f"configurations: {len(CONFIGURATIONS)} | "
        f"pending: {len(pending)} | skipped: {skipped}"
    )
    print(
        f"Timeout: {args.timeout:g}s per run | "
        f"results: {results_path}"
    )

    if not pending:
        print(
            "All discovered (instance, configuration) "
            "runs are already recorded."
        )
        return 0

    status_counts = {
        "SAT": 0,
        "UNSAT": 0,
        "UNKNOWN": 0,
        "TIMEOUT": 0,
    }
    completed_now = 0

    try:
        handle, writer = open_results_writer(results_path)

        with handle:
            for index, (instance, configuration) in enumerate(
                pending,
                start=1,
            ):
                print(
                    f"[{index}/{len(pending)}] "
                    f"{instance.name} | "
                    f"{configuration.name} ...",
                    end=" ",
                    flush=True,
                )

                result = run_kissat(
                    executable,
                    configuration,
                    instance,
                    args.timeout,
                )

                row = {
                    "instance": instance.name,
                    **metadata[instance],
                    "configuration": configuration.name,
                    **result,
                }

                persist_row(handle, writer, row)
                completed_now += 1
                status_counts[str(result["status"])] += 1

                elapsed = float(
                    str(result["wall_time_seconds"])
                )

                print(
                    f"{result['status']} in {elapsed:.2f}s | "
                    f"conflicts "
                    f"{printable_count(result['conflicts'])} | "
                    f"reorders "
                    f"{printable_count(result['number_of_reorders'])}"
                )

    except KeyboardInterrupt:
        print(
            f"\nInterrupted. {completed_now} new run(s) were saved; "
            "rerun the same command to resume.",
            file=sys.stderr,
        )
        return 130
    except OSError as error:
        print(f"\nerror: {error}", file=sys.stderr)
        return 2

    counts = ", ".join(
        f"{status}={status_counts[status]}"
        for status in (
            "SAT",
            "UNSAT",
            "UNKNOWN",
            "TIMEOUT",
        )
        if status_counts[status]
    )

    print(
        f"Finished: saved {completed_now} new run(s), "
        f"skipped {skipped}; "
        f"{counts or 'no new results'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())