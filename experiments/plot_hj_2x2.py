#!/usr/bin/env python3
"""Optional matplotlib plots for the HJ 2x2 benchmark CSV.

This file is deliberately separate: neither the benchmark nor the text
analysis requires matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot HJ 2x2 conflict ratios.")
    parser.add_argument("csv", nargs="?", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def optional_int(value: str) -> Optional[int]:
    value = value.strip().replace(",", "")
    return int(value) if value else None


def solved(row: dict[str, str]) -> bool:
    return row.get("status") in {"SAT", "UNSAT"} and optional_int(
        row.get("conflicts", "")
    ) is not None


def params(row: dict[str, str]) -> tuple[int, int, int]:
    try:
        return int(row["hj_k"]), int(row["hj_r"]), int(row["hj_n"])
    except (KeyError, ValueError):
        match = HJ_NAME_RE.fullmatch(row.get("instance", ""))
        if not match:
            raise ValueError(f"cannot parse HJ parameters from {row!r}")
        return tuple(int(value) for value in match.groups())  # type: ignore[return-value]


def load_latest(path: Path) -> dict[tuple[int, int, int], dict[str, dict[str, str]]]:
    grouped: dict[tuple[int, int, int], dict[str, dict[str, str]]] = {}
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            k, r, n = params(row)
            grouped.setdefault((k, r, n), {})[row["configuration"]] = row
    return grouped


def conflict_ratio(
    rows: dict[str, dict[str, str]], numerator: str, denominator: str
) -> Optional[float]:
    top = rows.get(numerator)
    bottom = rows.get(denominator)
    if top is None or bottom is None or not solved(top) or not solved(bottom):
        return None
    if top["status"] != bottom["status"]:
        return None
    a = optional_int(top["conflicts"])
    b = optional_int(bottom["conflicts"])
    assert a is not None and b is not None
    if b == 0:
        return math.inf if a > 0 else None
    value = a / b
    return value if value > 0 and math.isfinite(value) else None


def main() -> int:
    args = parse_args()
    csv_path = args.csv.expanduser().resolve()
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else csv_path.with_name(f"{csv_path.stem}_ratios.png")
    )
    if not csv_path.is_file():
        print(f"error: CSV does not exist: {csv_path}", file=sys.stderr)
        return 2
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "error: this optional script needs matplotlib; install it with "
            "'python3 -m pip install --user matplotlib'",
            file=sys.stderr,
        )
        return 2

    try:
        data = load_latest(csv_path)
    except (OSError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    definitions = (
        ("VSIDS effect: baseline / no-bump", "baseline", "no_bump"),
        (
            "Reorder effect with bumping: off / on",
            "baseline_no_reorder",
            "baseline",
        ),
        (
            "Reorder effect without bumping: off / on",
            "no_bump_no_reorder",
            "no_bump",
        ),
    )
    families = sorted({(k, r) for k, r, _ in data})
    figure, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=False, constrained_layout=True)
    plotted = 0
    for axis, (title, numerator, denominator) in zip(axes, definitions):
        for k, r in families:
            points = []
            for kk, rr, n in sorted(data):
                if (kk, rr) != (k, r):
                    continue
                value = conflict_ratio(data[(kk, rr, n)], numerator, denominator)
                if value is not None:
                    points.append((n, value))
            if points:
                plotted += 1
                axis.plot(
                    [n for n, _ in points],
                    [value for _, value in points],
                    marker="o",
                    label=f"HJ_{k}_{r}_n",
                )
        axis.axhline(1.0, color="black", linewidth=1, linestyle="--")
        axis.set_yscale("log")
        axis.set_title(title)
        axis.set_xlabel("n")
        axis.set_ylabel("conflict ratio (log scale)")
        axis.grid(True, which="both", alpha=0.25)
        if axis.lines:
            axis.legend()
    if not plotted:
        print("error: no comparable positive conflict ratios to plot", file=sys.stderr)
        return 2
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
