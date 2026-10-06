"""Replace an incomplete seed catalogue only after full contour/root coverage closes."""

from __future__ import annotations

import csv
import math
import shutil
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    assert len(sys.argv) == 2, "usage: augment_root_catalogue.py candidate_id"
    candidate_id = sys.argv[1]
    folder = HERE / "evaluations" / candidate_id
    coverage = read(folder / "ROOT_COVERAGE.csv")
    assert coverage[0]["status"] == "COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR"
    path = folder / "ROOTS.csv"
    if not (folder / "ROOTS_INITIAL.csv").exists():
        shutil.copy2(path, folder / "ROOTS_INITIAL.csv")
    discovered = read(folder / "DISCOVERED_ROOTS.csv")
    expanded = []
    for row in discovered:
        tau = float(row["local_crossing_ms"])
        frequency = float(row["frequency_at_crossing_hz"])
        expanded.append(
            {
                "candidate_id": candidate_id,
                "rank": 0,
                "seed_id": f"frequency_scan_{row['root_id']}",
                "local_crossing_ms": tau,
                "critical_real": -0.05,
                "critical_imag": 2 * math.pi * frequency,
                "frequency_hz": frequency,
                "small_factor_residual": row["crossing_root_residual"],
                "status": "FREQUENCY_SCAN_ROOT_COVERAGE_CLOSED_BY_CONTOUR",
            }
        )
    expanded.sort(key=lambda row: row["local_crossing_ms"])
    for index, row in enumerate(expanded, start=1):
        row["rank"] = index
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(expanded[0]))
        writer.writeheader()
        writer.writerows(expanded)
    print(candidate_id, "roots", len(expanded), "earliest", expanded[0]["local_crossing_ms"])


if __name__ == "__main__":
    main()
