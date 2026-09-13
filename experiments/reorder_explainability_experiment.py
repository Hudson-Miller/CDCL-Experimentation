#!/usr/bin/env python3
"""Collect Kissat structural-reorder data for HJ(4;2),7.

This is a data-collection script only.  It runs the two comparison settings

    --no-bump --reorderinit=3000
    --no-bump --reorderinit=7000

with the instrumented Kissat option ``--reorderlog=N``.  It preserves the
complete solver output and writes the raw reorder records to CSV.  It does
not decode Hales--Jewett variables or analyze the results.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


DEFAULT_PROJECT = Path.home() / "CDCL-Experimentation"
DEFAULT_CNF = DEFAULT_PROJECT / "cnfs" / "HJ_4_2_7.cnf"
DEFAULT_KISSAT = Path.home() / "kissat" / "build" / "kissat"
DEFAULT_TIMEOUT = 600.0
DEFAULT_TOP = 100
DEFAULT_REORDERINITS = (3000, 7000)

CONFLICT_RE = re.compile(
    r"^c\s+conflicts:\s*([0-9][0-9,]*)\b", re.MULTILINE
)
REORDERS_RE = re.compile(
    r"^c\s+reordered:\s*([0-9][0-9,]*)\b", re.MULTILINE
)
REORDER_WEIGHT_RE = re.compile(
    r"reorder-weight\s+"
    r"event=(?P<event>\d+)\s+"
    r"conflicts=(?P<conflicts>\d+)\s+"
    r"mode=(?P<mode>stable|focused)\s+"
    r"rank=(?P<rank>\d+)\s+"
    r"variable=(?P<variable>\d+)\s+"
    r"positive=(?P<positive>\S+)\s+"
    r"negative=(?P<negative>\S+)\s+"
    r"combined=(?P<combined>\S+)"
)

SUMMARY_FIELDS = (
    "timestamp_utc",
    "run_name",
    "instance",
    "instance_path",
    "variables",
    "clauses",
    "reorderinit",
    "reorderlog_top",
    "solver_arguments",
    "timeout_seconds",
    "status",
    "time_seconds",
    "total_conflicts",
    "reported_reorders",
    "logged_events",
    "logged_rows",
    "event_conflicts",
    "returncode",
    "solver_path",
    "solver_version",
    "log_file",
    "error",
)

WEIGHT_FIELDS = (
    "timestamp_utc",
    "run_name",
    "instance",
    "reorderinit",
    "reorderlog_top",
    "status",
    "total_conflicts",
    "time_seconds",
    "event",
    "conflicts",
    "mode",
    "rank",
    "dimacs_variable",
    "positive_weight",
    "negative_weight",
    "combined_weight",
    "solver_version",
    "log_file",
)


@dataclass(frozen=True)
class ReorderWeight:
    event: int
    conflicts: int
    mode: str
    rank: int
    variable: int
    positive: str
    negative: str
    combined: str


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect raw Kissat reorder weights for HJ(4;2),7."
    )
    parser.add_argument("--kissat", type=Path, default=DEFAULT_KISSAT)
    parser.add_argument("--cnf", type=Path, default=DEFAULT_CNF)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_PROJECT / "experiments" / "results",
        help="parent directory for a new timestamped result directory",
    )
    parser.add_argument("--top", type=int, default=DEFAULT_TOP)
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    parser.add_argument(
        "--reorderinit",
        dest="reorderinits",
        action="append",
        type=int,
        help="initial reorder threshold; repeat for multiple runs",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate paths and print commands without creating files",
    )
    args = parser.parse_args(argv)
    if args.top < 1:
        parser.error("--top must be positive")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.reorderinits is None:
        args.reorderinits = list(DEFAULT_REORDERINITS)
    if any(value < 0 for value in args.reorderinits):
        parser.error("--reorderinit values must be nonnegative")
    if len(set(args.reorderinits)) != len(args.reorderinits):
        parser.error("duplicate --reorderinit value")
    return args


def absolute(path: Path) -> Path:
    return path.expanduser().resolve()


def dimacs_header(path: Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("c"):
                continue
            fields = line.split()
            if len(fields) == 4 and fields[:2] == ["p", "cnf"]:
                variables, clauses = int(fields[2]), int(fields[3])
                if variables < 1 or clauses < 1:
                    raise ValueError("DIMACS counts must be positive")
                return variables, clauses
            raise ValueError(
                f"{path}:{line_number}: expected a DIMACS p cnf header"
            )
    raise ValueError(f"{path}: missing DIMACS p cnf header")


def run_short_command(command: list[str]) -> tuple[int, str]:
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=15,
        check=False,
    )
    return completed.returncode, completed.stdout + "\n" + completed.stderr


def inspect_solver(kissat: Path) -> str:
    _, version_output = run_short_command([str(kissat), "--version"])
    version_lines = [line.strip() for line in version_output.splitlines() if line.strip()]
    version = version_lines[0] if version_lines else "unknown"

    _, help_output = run_short_command([str(kissat), "--help"])
    if "--reorderlog" not in help_output:
        raise RuntimeError(
            f"{kissat} does not support --reorderlog. "
            "Patch and rebuild Kissat before running this experiment."
        )
    return version


def terminate_process(process: subprocess.Popen[str], force: bool = False) -> None:
    chosen_signal = signal.SIGKILL if force else signal.SIGTERM
    try:
        if os.name == "posix":
            os.killpg(process.pid, chosen_signal)
        elif force:
            process.kill()
        else:
            process.terminate()
    except ProcessLookupError:
        pass


def run_command(command: list[str], timeout: float) -> tuple[str, str, int, bool, float]:
    started = time.perf_counter()
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        start_new_session=(os.name == "posix"),
    )
    timed_out = False
    try:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            terminate_process(process)
            try:
                stdout, stderr = process.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                terminate_process(process, force=True)
                stdout, stderr = process.communicate()
    except KeyboardInterrupt:
        terminate_process(process)
        try:
            process.communicate(timeout=2)
        except subprocess.TimeoutExpired:
            terminate_process(process, force=True)
            process.communicate()
        raise
    elapsed = time.perf_counter() - started
    return stdout, stderr, int(process.returncode), timed_out, elapsed


def classify_status(output: str, returncode: int, timed_out: bool) -> str:
    if timed_out:
        return "TIMEOUT"
    if re.search(r"^s\s+UNSATISFIABLE\s*$", output, re.MULTILINE):
        return "UNSAT"
    if re.search(r"^s\s+SATISFIABLE\s*$", output, re.MULTILINE):
        return "SAT"
    if returncode == 20:
        return "UNSAT"
    if returncode == 10:
        return "SAT"
    if returncode not in (0, 10, 20):
        return "ERROR"
    return "UNKNOWN"


def last_counter(pattern: re.Pattern[str], output: str) -> Optional[int]:
    matches = pattern.findall(output)
    return int(matches[-1].replace(",", "")) if matches else None


def parse_weights(output: str) -> list[ReorderWeight]:
    records = []
    for match in REORDER_WEIGHT_RE.finditer(output):
        records.append(
            ReorderWeight(
                event=int(match.group("event")),
                conflicts=int(match.group("conflicts")),
                mode=match.group("mode"),
                rank=int(match.group("rank")),
                variable=int(match.group("variable")),
                positive=match.group("positive"),
                negative=match.group("negative"),
                combined=match.group("combined"),
            )
        )
    return records


def utc_timestamp() -> str:
    return (
        dt.datetime.now(dt.timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


def make_output_directory(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    base = root / f"reorder_explainability_HJ_4_2_7_{stamp}"
    for suffix in range(10_000):
        candidate = base if suffix == 0 else root / f"{base.name}_{suffix}"
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            continue
    raise RuntimeError(f"could not create a unique directory under {root}")


def initialize_csv(path: Path, fields: tuple[str, ...]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        csv.DictWriter(handle, fieldnames=fields).writeheader()
        handle.flush()
        os.fsync(handle.fileno())


def append_rows(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def write_log(
    path: Path,
    command: list[str],
    stdout: str,
    stderr: str,
    status: str,
    elapsed: float,
    returncode: int,
) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(f"command: {shlex.join(command)}\n")
        handle.write(f"status: {status}\n")
        handle.write(f"elapsed_seconds: {elapsed:.9f}\n")
        handle.write(f"returncode: {returncode}\n")
        handle.write("\n===== STDOUT =====\n")
        handle.write(stdout)
        if stdout and not stdout.endswith("\n"):
            handle.write("\n")
        handle.write("\n===== STDERR =====\n")
        handle.write(stderr)
        if stderr and not stderr.endswith("\n"):
            handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def error_summary(status: str, stderr: str, records: list[ReorderWeight]) -> str:
    if status in {"SAT", "UNSAT"} and not records:
        return "solver finished but emitted no reorder-weight records"
    if status not in {"ERROR", "UNKNOWN"}:
        return ""
    lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    return (lines[-1] if lines else "solver produced no recognizable result")[:500]


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    kissat = absolute(args.kissat)
    cnf = absolute(args.cnf)
    output_root = absolute(args.output_root)

    if not kissat.is_file() or not os.access(kissat, os.X_OK):
        print(f"error: Kissat is not executable: {kissat}", file=sys.stderr)
        return 2
    if not cnf.is_file():
        print(f"error: CNF does not exist: {cnf}", file=sys.stderr)
        return 2

    try:
        variables, clauses = dimacs_header(cnf)
        if cnf.name == "HJ_4_2_7.cnf" and variables != 4**7:
            raise ValueError(
                f"expected 16,384 variables for HJ_4_2_7.cnf, found {variables:,}"
            )
        version = inspect_solver(kissat)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    commands = [
        [
            str(kissat),
            "--no-bump",
            f"--reorderinit={reorderinit}",
            f"--reorderlog={args.top}",
            str(cnf),
        ]
        for reorderinit in args.reorderinits
    ]

    print("Kissat structural-reorder data collection")
    print(f"Kissat:   {kissat} ({version})")
    print(f"Instance: {cnf} ({variables:,} variables, {clauses:,} clauses)")
    print(f"Top rows: {args.top} per reorder | timeout: {args.timeout:g}s per run")

    if args.dry_run:
        print("\nCommands:")
        for command in commands:
            print(f"  {shlex.join(command)}")
        return 0

    try:
        output_dir = make_output_directory(output_root)
        summary_csv = output_dir / "run_summary.csv"
        weights_csv = output_dir / "reorder_weights.csv"
        initialize_csv(summary_csv, SUMMARY_FIELDS)
        initialize_csv(weights_csv, WEIGHT_FIELDS)
    except (OSError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"Output:   {output_dir}\n")
    failed = False

    for index, (reorderinit, command) in enumerate(
        zip(args.reorderinits, commands), 1
    ):
        run_name = f"reorderinit_{reorderinit}"
        print(f"[{index}/{len(commands)}] {run_name} ... ", end="", flush=True)
        try:
            stdout, stderr, returncode, timed_out, elapsed = run_command(
                command, args.timeout
            )
            combined = stdout + "\n" + stderr
            status = classify_status(combined, returncode, timed_out)
            total_conflicts = last_counter(CONFLICT_RE, combined)
            reported_reorders = last_counter(REORDERS_RE, combined)
            records = parse_weights(combined)
            event_conflicts = sorted({record.conflicts for record in records})
            logged_events = len({record.event for record in records})
            timestamp = utc_timestamp()
            log_path = output_dir / f"{run_name}.log"
            write_log(
                log_path,
                command,
                stdout,
                stderr,
                status,
                elapsed,
                returncode,
            )

            summary_row: dict[str, object] = {
                "timestamp_utc": timestamp,
                "run_name": run_name,
                "instance": cnf.name,
                "instance_path": str(cnf),
                "variables": variables,
                "clauses": clauses,
                "reorderinit": reorderinit,
                "reorderlog_top": args.top,
                "solver_arguments": shlex.join(command[1:-1]),
                "timeout_seconds": f"{args.timeout:.9g}",
                "status": status,
                "time_seconds": f"{elapsed:.9f}",
                "total_conflicts": "" if total_conflicts is None else total_conflicts,
                "reported_reorders": (
                    "" if reported_reorders is None else reported_reorders
                ),
                "logged_events": logged_events,
                "logged_rows": len(records),
                "event_conflicts": ";".join(map(str, event_conflicts)),
                "returncode": returncode,
                "solver_path": str(kissat),
                "solver_version": version,
                "log_file": str(log_path),
                "error": error_summary(status, stderr, records),
            }
            append_rows(summary_csv, SUMMARY_FIELDS, [summary_row])

            weight_rows = [
                {
                    "timestamp_utc": timestamp,
                    "run_name": run_name,
                    "instance": cnf.name,
                    "reorderinit": reorderinit,
                    "reorderlog_top": args.top,
                    "status": status,
                    "total_conflicts": (
                        "" if total_conflicts is None else total_conflicts
                    ),
                    "time_seconds": f"{elapsed:.9f}",
                    "event": record.event,
                    "conflicts": record.conflicts,
                    "mode": record.mode,
                    "rank": record.rank,
                    "dimacs_variable": record.variable,
                    "positive_weight": record.positive,
                    "negative_weight": record.negative,
                    "combined_weight": record.combined,
                    "solver_version": version,
                    "log_file": str(log_path),
                }
                for record in records
            ]
            append_rows(weights_csv, WEIGHT_FIELDS, weight_rows)

            conflict_display = (
                "?" if total_conflicts is None else f"{total_conflicts:,}"
            )
            event_display = ", ".join(f"{value:,}" for value in event_conflicts) or "none"
            print(
                f"{status} in {elapsed:.2f}s | total conflicts {conflict_display} "
                f"| reorder events {event_display} | rows {len(records)}"
            )
            if status not in {"SAT", "UNSAT"} or not records:
                failed = True
        except KeyboardInterrupt:
            print("interrupted")
            print(f"Completed-run data remains in {output_dir}")
            return 130
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
            print(f"HARNESS ERROR: {exc}", file=sys.stderr)
            return 2

    print("\nFinished collecting data. Please send back:")
    print(f"  {summary_csv}")
    print(f"  {weights_csv}")
    print(f"  {output_dir}/*.log")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
