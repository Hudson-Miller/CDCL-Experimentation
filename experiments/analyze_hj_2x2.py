#!/usr/bin/env python3
"""Print family-wise tables for the final Hales--Jewett 2x2 CSV.

This script uses only the Python standard library.  Its ratios are:

    R_VSIDS = baseline_conflicts / no_bump_conflicts

    reorder@bump = baseline_no_reorder_conflicts / baseline_conflicts

    reorder@no-bump =
        no_bump_no_reorder_conflicts / no_bump_conflicts

For the two reorder ratios, values greater than one mean enabling Kissat's
default reordering reduced conflicts.  Ratios are printed only when both
runs solved the instance with the same SAT/UNSAT status.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


DEFAULT_CSV = (
    Path.home()
    / "CDCL-Experimentation"
    / "experiments"
    / "results"
    / "hj_2x2_final_20260815.csv"
)
HJ_NAME_RE = re.compile(r"^HJ_(\d+)_(\d+)_(\d+)\.cnf$")
PRIMARY = (
    "baseline",
    "baseline_no_reorder",
    "no_bump",
    "no_bump_no_reorder",
)


@dataclass(frozen=True)
class Run:
    instance: str
    k: int
    r: int
    n: int
    configuration: str
    status: str
    conflicts: Optional[int]
    seconds: Optional[float]
    reorders: Optional[int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze a final HJ bump x reorder benchmark CSV."
    )
    parser.add_argument("csv", nargs="?", type=Path, default=DEFAULT_CSV)
    parser.add_argument(
        "--family",
        action="append",
        default=[],
        metavar="K,R",
        help="Restrict output to a family; may be repeated (example: --family 3,3).",
    )
    parser.add_argument(
        "--show-diagnostics",
        action="store_true",
        help="Include non-primary configurations in run tables.",
    )
    return parser.parse_args()


def parse_optional_int(value: str) -> Optional[int]:
    value = value.strip().replace(",", "")
    return None if not value else int(value)


def parse_optional_float(value: str) -> Optional[float]:
    value = value.strip()
    return None if not value else float(value)


def parse_parameters(row: dict[str, str]) -> tuple[int, int, int]:
    try:
        return int(row["hj_k"]), int(row["hj_r"]), int(row["hj_n"])
    except (KeyError, TypeError, ValueError):
        match = HJ_NAME_RE.fullmatch(row.get("instance", ""))
        if not match:
            raise ValueError(f"cannot parse HJ parameters from row: {row!r}")
        return tuple(int(value) for value in match.groups())  # type: ignore[return-value]


def load_latest(path: Path) -> tuple[dict[tuple[int, int, int, str], Run], int]:
    latest: dict[tuple[int, int, int, str], Run] = {}
    duplicates = 0
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"instance", "configuration", "status", "conflicts", "time_seconds"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"CSV is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            k, r, n = parse_parameters(row)
            run = Run(
                instance=row["instance"],
                k=k,
                r=r,
                n=n,
                configuration=row["configuration"],
                status=row["status"],
                conflicts=parse_optional_int(row.get("conflicts", "")),
                seconds=parse_optional_float(row.get("time_seconds", "")),
                reorders=parse_optional_int(row.get("reorders", "")),
            )
            key = (k, r, n, run.configuration)
            if key in latest:
                duplicates += 1
            latest[key] = run
    return latest, duplicates


def parse_family_filters(values: list[str]) -> set[tuple[int, int]]:
    result: set[tuple[int, int]] = set()
    for value in values:
        parts = value.split(",")
        if len(parts) != 2:
            raise ValueError(f"invalid family {value!r}; expected K,R")
        result.add((int(parts[0]), int(parts[1])))
    return result


def ratio(numerator: Optional[Run], denominator: Optional[Run]) -> Optional[float]:
    if numerator is None or denominator is None:
        return None
    if numerator.status not in {"SAT", "UNSAT"}:
        return None
    if numerator.status != denominator.status:
        return None
    if numerator.conflicts is None or denominator.conflicts is None:
        return None
    if denominator.conflicts == 0:
        return math.inf if numerator.conflicts > 0 else None
    return numerator.conflicts / denominator.conflicts


def format_integer(value: Optional[int]) -> str:
    return "-" if value is None else f"{value:,}"


def format_seconds(value: Optional[float]) -> str:
    return "-" if value is None else f"{value:.3f}"


def format_ratio(value: Optional[float]) -> str:
    if value is None:
        return "-"
    if math.isinf(value):
        return "inf"
    if value >= 100:
        return f"{value:.1f}"
    if value >= 10:
        return f"{value:.2f}"
    return f"{value:.3f}"


def print_table(headers: tuple[str, ...], rows: list[tuple[str, ...]]) -> None:
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))
    numeric = set(range(len(headers))) - {1, 2}

    def rendered(row: tuple[str, ...]) -> str:
        cells = []
        for index, cell in enumerate(row):
            cells.append(cell.rjust(widths[index]) if index in numeric else cell.ljust(widths[index]))
        return "  ".join(cells)

    print(rendered(headers))
    print("  ".join("-" * width for width in widths))
    for row in rows:
        print(rendered(row))


def transition_flags(
    ratios: list[tuple[int, Optional[float]]]
) -> tuple[dict[int, str], list[str]]:
    flags: dict[int, str] = {}
    messages: list[str] = []
    prior: Optional[tuple[int, float, int]] = None
    for n, value in ratios:
        if value is None or math.isnan(value) or value == 1.0:
            continue
        sign = 1 if value > 1.0 else -1
        if prior is not None and sign != prior[2]:
            old_n, old_value, _ = prior
            arrow = "CROSS UP" if sign > 0 else "CROSS DOWN"
            flags[n] = arrow
            messages.append(
                f"{arrow}: n={old_n} (R={format_ratio(old_value)}) -> "
                f"n={n} (R={format_ratio(value)})"
            )
        prior = (n, value, sign)
    return flags, messages


def main() -> int:
    args = parse_args()
    csv_path = args.csv.expanduser().resolve()
    if not csv_path.is_file():
        print(f"error: CSV does not exist: {csv_path}", file=sys.stderr)
        return 2
    try:
        selected_families = parse_family_filters(args.family)
        latest, duplicates = load_latest(csv_path)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    grouped: dict[tuple[int, int], dict[int, dict[str, Run]]] = {}
    for (k, r, n, configuration), run in latest.items():
        if selected_families and (k, r) not in selected_families:
            continue
        grouped.setdefault((k, r), {}).setdefault(n, {})[configuration] = run
    if not grouped:
        print("error: no rows matched the selected families", file=sys.stderr)
        return 2

    print(f"Source: {csv_path}")
    print("Latest CSV row wins when an instance/configuration was rerun.")
    if duplicates:
        print(f"Rerun rows superseded: {duplicates}")
    print()
    print("R = baseline / no_bump; R>1 means no-bump used fewer conflicts.")
    print("Reorder effects are OFF / ON; effect>1 means reordering used fewer conflicts.")

    total_instances = 0
    for family in sorted(grouped):
        k, r = family
        by_n = grouped[family]
        total_instances += len(by_n)
        print()
        print(f"=== HJ_{k}_{r}_n ===")
        print()
        print("Runs")
        run_rows: list[tuple[str, ...]] = []
        for n in sorted(by_n):
            configs = by_n[n]
            names = list(PRIMARY)
            if args.show_diagnostics:
                names.extend(sorted(set(configs).difference(PRIMARY)))
            for name in names:
                run = configs.get(name)
                if run is None:
                    run_rows.append((str(n), name, "MISSING", "-", "-", "-"))
                else:
                    run_rows.append(
                        (
                            str(n),
                            name,
                            run.status,
                            format_integer(run.conflicts),
                            format_seconds(run.seconds),
                            format_integer(run.reorders),
                        )
                    )
        print_table(
            ("n", "configuration", "status", "conflicts", "time_s", "reorders"),
            run_rows,
        )

        ratios_by_n: list[tuple[int, Optional[float]]] = []
        effect_values: dict[int, tuple[Optional[float], Optional[float], Optional[float]]] = {}
        for n in sorted(by_n):
            configs = by_n[n]
            r_vsids = ratio(configs.get("baseline"), configs.get("no_bump"))
            reorder_bump = ratio(
                configs.get("baseline_no_reorder"), configs.get("baseline")
            )
            reorder_no_bump = ratio(
                configs.get("no_bump_no_reorder"), configs.get("no_bump")
            )
            ratios_by_n.append((n, r_vsids))
            effect_values[n] = (r_vsids, reorder_bump, reorder_no_bump)
        flags, transition_messages = transition_flags(ratios_by_n)

        print()
        print("Effects (conflict ratios)")
        effect_rows = []
        for n in sorted(effect_values):
            r_vsids, reorder_bump, reorder_no_bump = effect_values[n]
            effect_rows.append(
                (
                    str(n),
                    format_ratio(r_vsids),
                    format_ratio(reorder_bump),
                    format_ratio(reorder_no_bump),
                    flags.get(n, ""),
                )
            )
        print_table(
            ("n", "R_vsids", "reorder@bump", "reorder@no-bump", "transition"),
            effect_rows,
        )
        print()
        if transition_messages:
            for message in transition_messages:
                print(f"TRANSITION: {message}")
        else:
            print("TRANSITION: no observed crossing of R=1 among comparable solved rows.")

    print()
    print(f"Families: {len(grouped)} | instances: {total_instances} | latest runs: {len(latest)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
