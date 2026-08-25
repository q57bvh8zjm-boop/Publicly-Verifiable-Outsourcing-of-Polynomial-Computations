#!/usr/bin/env python3
"""Validate reported-data grids and repository source hygiene."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORTED = ROOT / "results" / "reported"
HAN = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


def rows(name: str) -> list[dict[str, str]]:
    with (REPORTED / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def assert_grid(name: str, actual: set[tuple[str, ...]], expected: set[tuple[str, ...]]) -> None:
    if actual == expected:
        return
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    raise AssertionError(f"{name}: missing={missing[:5]}, extra={extra[:5]}")


def validate_figure_grids() -> None:
    m_large = {str(value) for value in range(200, 2001, 200)}

    figure2 = rows("figure2_low_degree.csv")
    assert_grid(
        "Figure 2",
        {(row["scheme"], row["d"], row["t"], row["m"]) for row in figure2},
        {(scheme, "2", "1", m) for scheme in ("Gamma1", "Gamma2", "Gamma3") for m in m_large},
    )

    figure3 = rows("figure3_thresholds.csv")
    assert_grid(
        "Figure 3",
        {(row["scheme"], row["d"], row["t"], row["m"]) for row in figure3},
        {(scheme, "2", t, m) for scheme in ("Gamma1", "Gamma2") for t in ("2", "3") for m in m_large},
    )

    figure4 = rows("figure4_high_degree.csv")
    assert_grid(
        "Figure 4",
        {(row["scheme"], row["d"], row["t"], row["m"]) for row in figure4},
        {(scheme, d, "1", str(m)) for scheme in ("Gamma1", "Gamma2", "Gamma3") for d in ("4", "8", "10") for m in range(1, 16)},
    )

    figure5 = rows("figure5_gamma1_vs_pi3.csv")
    assert_grid(
        "Figure 5",
        {(row["d"], row["t"], row["m"]) for row in figure5},
        {(d, t, m) for d in ("2", "3") for t in ("1", "2") for m in m_large},
    )

    for filename in (
        "figure2_low_degree.csv",
        "figure3_thresholds.csv",
        "figure4_high_degree.csv",
        "figure5_gamma1_vs_pi3.csv",
        "table_microbenchmark.csv",
        "table_var.csv",
    ):
        for row in rows(filename):
            for key, value in row.items():
                if value is None or not value.strip():
                    raise AssertionError(f"{filename}: empty value in column {key}")


def validate_source_comments() -> None:
    violations: list[str] = []
    for path in sorted((*ROOT.rglob("*.c"), *ROOT.rglob("*.h"))):
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if HAN.search(line):
                violations.append(f"{path.relative_to(ROOT)}:{line_number}")
    if violations:
        raise AssertionError("Chinese characters found in C source: " + ", ".join(violations[:10]))


def main() -> int:
    try:
        validate_figure_grids()
        validate_source_comments()
    except (AssertionError, KeyError, OSError) as error:
        print(f"Artifact validation failed: {error}", file=sys.stderr)
        return 1
    print("Artifact data grids and source-comment checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
