#!/usr/bin/env python3
"""Run the fixed Hales--Jewett bump x reorder ablation with Kissat.

Primary configurations (and only these, unless --stable-diagnostic is used):

    baseline                 kissat
    baseline_no_reorder      kissat --reorder=0
    no_bump                  kissat --no-bump
    no_bump_no_reorder       kissat --no-bump --reorder=0

The script discovers plain DIMACS files named HJ_<k>_<r>_<n>.cnf, appends
one result per completed run to a CSV, flushes and fsyncs every row, and can
resume safely from that CSV.  It uses only the Python standard library.
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
from typing import Iterable, Optional


DEFAULT_TIMEOUT = 300.0
DEFAULT_KISSAT = Path.home() / "kissat" / "build" / "kissat"
DEFAULT_PROJECT = Path.home() / "CDCL-Experimentation"
DEFAULT_CNF_DIR = DEFAULT_PROJECT / "cnfs"
DEFAULT_RESULTS = (
    DEFAULT_PROJECT
    / "experiments"
    / "results"
    / "hj_2x2_final_20260815.csv"
)

HJ_NAME_RE = re.compile(r"^HJ_(\d+)_(\d+)_(\d+)\.cnf$")
CONFLICT_RE = re.compile(r"^c\s+conflicts:\s*([0-9][0-9,]*)\b", re.MULTILINE)
REORDER_RE = re.compile(r"^c\s+reordered:\s*([0-9][0-9,]*)\b", re.MULTILINE)


@dataclass(frozen=True)
class Configuration:
    name: str
    bump: str
    reorder: str
    arguments: tuple[str, ...]


PRIMARY_CONFIGURATIONS = (
    Configuration("baseline", "on", "default", ()),
    Configuration("baseline_no_reorder", "on", "off", ("--reorder=0",)),
    Configuration("no_bump", "off", "default", ("--no-bump",)),
    Configuration(
        "no_bump_no_reorder",
        "off",
        "off",
        ("--no-bump", "--reorder=0"),
    ),
)

STABLE_DIAGNOSTIC = Configuration(
    "no_bump_stable_reorder",
    "off",
    "stable-only",
    ("--no-bump", "--reorder=1"),
)

CSV_FIELDS = (
    "timestamp_utc",
    "instance",
    "instance_path",
    "hj_k",
    "hj_r",
    "hj_n",
    "variables",
    "clauses",
    "configuration",
    "bump",
    "reorder",
    "solver_arguments",
    "timeout_seconds",
    "status",
    "time_seconds",
    "conflicts",
    "reorders",
    "returncode",
    "solver_path",
    "solver_version",
    "cnf_bytes",
    "cnf_mtime_ns",
    "log_file",
    "error",
)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the fixed 2x2 Kissat bump x reorder HJ ablation."
    )
    parser.add_argument("--kissat", type=Path, default=DEFAULT_KISSAT)
    parser.add_argument("--cnf-dir", type=Path, default=DEFAULT_CNF_DIR)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=None,
        help="Raw-output directory (default: <results stem>_logs beside CSV).",
    )
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    parser.add_argument(
        "--filter-regex",
        default=r"^HJ_\d+_\d+_\d+\.cnf$",
        help="Regex applied to discovered filenames.",
    )
    parser.add_argument(
        "--stable-diagnostic",
        action="store_true",
        help="Also run the optional --no-bump --reorder=1 diagnostic.",
    )
    parser.add_argument(
        "--rerun-timeouts",
        action="store_true",
        help="Append new attempts for prior TIMEOUT rows instead of skipping them.",
    )
    parser.add_argument(
        "--rerun-failures",
        action="store_true",
        help="Append new attempts for prior ERROR/UNKNOWN rows.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print discovery and commands without creating files or running Kissat.",
    )
    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    try:
        args.compiled_filter = re.compile(args.filter_regex)
    except re.error as exc:
        parser.error(f"invalid --filter-regex: {exc}")
    return args


def absolute(path: Path) -> Path:
    return path.expanduser().resolve()


def discover_instances(cnf_dir: Path, filename_filter: re.Pattern[str]) -> list[Path]:
    found: list[tuple[tuple[int, int, int], Path]] = []
    for path in cnf_dir.glob("HJ_*.cnf"):
        match = HJ_NAME_RE.fullmatch(path.name)
        if not match or not filename_filter.search(path.name) or not path.is_file():
            continue
        found.append((tuple(int(value) for value in match.groups()), path.resolve()))
    found.sort(key=lambda item: (item[0], item[1].name))
    return [path for _, path in found]


def hj_parameters(path: Path) -> tuple[int, int, int]:
    match = HJ_NAME_RE.fullmatch(path.name)
    if not match:
        raise ValueError(f"not an HJ filename: {path.name}")
    return tuple(int(value) for value in match.groups())  # type: ignore[return-value]


def dimacs_header(path: Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("c"):
                continue
            parts = line.split()
            if len(parts) == 4 and parts[:2] == ["p", "cnf"]:
                variables = int(parts[2])
                clauses = int(parts[3])
                if variables < 0 or clauses < 0:
                    raise ValueError("negative value in DIMACS header")
                return variables, clauses
            raise ValueError(f"expected DIMACS header, found: {line[:160]!r}")
    raise ValueError("missing DIMACS 'p cnf' header")


def solver_version(kissat: Path) -> str:
    try:
        completed = subprocess.run(
            [str(kissat), "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable: {exc}"
    lines = (completed.stdout + "\n" + completed.stderr).splitlines()
    return lines[0].strip() if lines else "unknown"


def initialize_csv(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
        if header != list(CSV_FIELDS):
            raise RuntimeError(
                f"refusing to append to {path}: its header does not match this experiment"
            )
        return
    try:
        with path.open("x", newline="", encoding="utf-8") as handle:
            csv.DictWriter(handle, fieldnames=CSV_FIELDS).writeheader()
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError:
        initialize_csv(path)


def append_csv(path: Path, row: dict[str, object]) -> None:
    with path.open("a", newline="", encoding="utf-8") as handle:
        csv.DictWriter(handle, fieldnames=CSV_FIELDS).writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def load_latest_statuses(path: Path) -> dict[tuple[str, str], str]:
    latest: dict[tuple[str, str], str] = {}
    if not path.exists():
        return latest
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(CSV_FIELDS):
            raise RuntimeError(f"unexpected CSV header in {path}")
        for row in reader:
            latest[(row["instance"], row["configuration"])] = row["status"]
    return latest


def should_skip(status: str, rerun_timeouts: bool, rerun_failures: bool) -> bool:
    if status == "TIMEOUT" and rerun_timeouts:
        return False
    if status in {"ERROR", "UNKNOWN"} and rerun_failures:
        return False
    return True


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
                stdout, stderr = process.communicate(timeout=5.0)
            except subprocess.TimeoutExpired:
                terminate_process(process, force=True)
                stdout, stderr = process.communicate()
    except KeyboardInterrupt:
        terminate_process(process)
        try:
            process.communicate(timeout=2.0)
        except subprocess.TimeoutExpired:
            terminate_process(process, force=True)
            process.communicate()
        raise
    elapsed = time.perf_counter() - started
    return stdout, stderr, int(process.returncode), timed_out, elapsed


def last_stat(pattern: re.Pattern[str], text: str) -> Optional[int]:
    matches = pattern.findall(text)
    if not matches:
        return None
    return int(matches[-1].replace(",", ""))


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


def error_summary(status: str, stderr: str) -> str:
    if status not in {"ERROR", "UNKNOWN"}:
        return ""
    lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    return (lines[-1] if lines else "solver produced no recognizable result")[:500]


def unique_log_path(log_dir: Path, instance: str, configuration: str) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    base = log_dir / f"{Path(instance).stem}__{configuration}.log"
    if not base.exists():
        return base
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    for attempt in range(1, 10_000):
        candidate = base.with_name(f"{base.stem}__{stamp}__{attempt}{base.suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"could not choose a unique log filename under {log_dir}")


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


def configurations(include_stable: bool) -> tuple[Configuration, ...]:
    if include_stable:
        return PRIMARY_CONFIGURATIONS + (STABLE_DIAGNOSTIC,)
    return PRIMARY_CONFIGURATIONS


def display_count(value: Optional[int]) -> str:
    return "?" if value is None else f"{value:,}"


def command_for(kissat: Path, config: Configuration, cnf: Path) -> list[str]:
    return [str(kissat), *config.arguments, str(cnf)]


def print_dry_run(
    kissat: Path, instances: Iterable[Path], configs: Iterable[Configuration]
) -> None:
    for instance in instances:
        for config in configs:
            print(f"{instance.name} | {config.name}")
            print(f"  {shlex.join(command_for(kissat, config, instance))}")


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    kissat = absolute(args.kissat)
    cnf_dir = absolute(args.cnf_dir)
    results = absolute(args.results)
    log_dir = absolute(
        args.log_dir
        if args.log_dir is not None
        else results.parent / f"{results.stem}_logs"
    )
    configs = configurations(args.stable_diagnostic)

    if not cnf_dir.is_dir():
        print(f"error: CNF directory does not exist: {cnf_dir}", file=sys.stderr)
        return 2
    if not kissat.is_file() or not os.access(kissat, os.X_OK):
        print(f"error: Kissat is not an executable file: {kissat}", file=sys.stderr)
        return 2

    instances = discover_instances(cnf_dir, args.compiled_filter)
    if not instances:
        print(
            f"error: no matching cnfs/HJ_*.cnf files found under {cnf_dir}",
            file=sys.stderr,
        )
        return 2

    if args.dry_run:
        print(f"Kissat: {kissat}")
        print(f"Instances: {len(instances)} | configurations: {len(configs)}")
        print_dry_run(kissat, instances, configs)
        return 0

    try:
        initialize_csv(results)
        latest = load_latest_statuses(results)
    except (OSError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    pending: list[tuple[Path, Configuration]] = []
    skipped = 0
    for instance in instances:
        for config in configs:
            prior = latest.get((instance.name, config.name))
            if prior is not None and should_skip(
                prior, args.rerun_timeouts, args.rerun_failures
            ):
                skipped += 1
            else:
                pending.append((instance, config))

    version = solver_version(kissat)
    print(f"Kissat: {kissat} ({version})")
    print(
        f"Instances: {len(instances)} | configurations: {len(configs)} | "
        f"pending: {len(pending)} | skipped: {skipped}"
    )
    print(f"Timeout: {args.timeout:g}s per run | results: {results}")
    if not pending:
        print("Nothing to do: all selected instance/configuration pairs are recorded.")
        return 0

    totals: dict[str, int] = {}
    for index, (instance, config) in enumerate(pending, start=1):
        print(
            f"[{index}/{len(pending)}] {instance.name} | {config.name} ... ",
            end="",
            flush=True,
        )
        try:
            k, r, n = hj_parameters(instance)
            variables, clauses = dimacs_header(instance)
            stat = instance.stat()
            command = command_for(kissat, config, instance)
            stdout, stderr, returncode, timed_out, elapsed = run_command(
                command, args.timeout
            )
            combined = stdout + "\n" + stderr
            status = classify_status(combined, returncode, timed_out)
            conflicts = last_stat(CONFLICT_RE, combined)
            reorders = last_stat(REORDER_RE, combined)
            if reorders is None and status in {"SAT", "UNSAT"}:
                reorders = 0
            log_path = unique_log_path(log_dir, instance.name, config.name)
            write_log(
                log_path,
                command,
                stdout,
                stderr,
                status,
                elapsed,
                returncode,
            )
            row: dict[str, object] = {
                "timestamp_utc": dt.datetime.now(dt.timezone.utc)
                .isoformat(timespec="seconds")
                .replace("+00:00", "Z"),
                "instance": instance.name,
                "instance_path": str(instance),
                "hj_k": k,
                "hj_r": r,
                "hj_n": n,
                "variables": variables,
                "clauses": clauses,
                "configuration": config.name,
                "bump": config.bump,
                "reorder": config.reorder,
                "solver_arguments": shlex.join(config.arguments),
                "timeout_seconds": f"{args.timeout:.9g}",
                "status": status,
                "time_seconds": f"{elapsed:.9f}",
                "conflicts": "" if conflicts is None else conflicts,
                "reorders": "" if reorders is None else reorders,
                "returncode": returncode,
                "solver_path": str(kissat),
                "solver_version": version,
                "cnf_bytes": stat.st_size,
                "cnf_mtime_ns": stat.st_mtime_ns,
                "log_file": str(log_path),
                "error": error_summary(status, stderr),
            }
            append_csv(results, row)
        except KeyboardInterrupt:
            print("interrupted; the current pair was not added to the CSV")
            return 130
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
            print(f"HARNESS ERROR: {exc}", file=sys.stderr)
            return 2

        totals[status] = totals.get(status, 0) + 1
        print(
            f"{status} in {elapsed:.2f}s | conflicts {display_count(conflicts)} "
            f"| reorders {display_count(reorders)}"
        )

    summary = ", ".join(f"{key}={totals[key]}" for key in sorted(totals))
    print(
        f"Finished: saved {len(pending)} new run(s), skipped {skipped}"
        + (f"; {summary}" if summary else "")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
