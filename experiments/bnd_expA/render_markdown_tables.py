#!/usr/bin/env python
"""Render the Experiment A CSV tables as a supplementary Markdown appendix."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def render(path: Path) -> str:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        rows = list(reader)
    if not rows:
        return f"## {path.stem}\n\n_(empty table)_\n"
    widths = len(rows[0])
    normalized = [row + [""] * max(0, widths - len(row)) for row in rows]
    header = normalized[0]
    lines = [
        f"## {path.stem}",
        "",
        "| " + " | ".join(cell(x) for x in header) + " |",
        "|" + "|".join("---" for _ in header) + "|",
    ]
    lines.extend("| " + " | ".join(cell(x) for x in row[:widths]) + " |" for row in normalized[1:])
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir", required=True, type=Path)
    args = parser.parse_args()
    tables = args.report_dir / "tables"
    files = sorted(tables.glob("TABLE_A*.csv"))
    if not files:
        raise SystemExit(f"no TABLE_A*.csv files found under {tables}")
    output = [
        "# Experiment A tables",
        "",
        "Generated directly from the machine-readable CSV files. Numeric source of truth: the CSVs under `tables/`.",
        "",
    ]
    output.extend(render(path) for path in files)
    (args.report_dir / "TABLES_EXP_A.md").write_text("\n".join(output), encoding="utf-8")


if __name__ == "__main__":
    main()
