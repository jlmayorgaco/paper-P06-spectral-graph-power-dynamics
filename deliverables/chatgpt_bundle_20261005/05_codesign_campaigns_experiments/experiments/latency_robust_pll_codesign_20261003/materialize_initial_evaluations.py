"""Expose the already executed first corrector in the generic iteration layout."""

from __future__ import annotations

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    summaries = read("L3_INITIAL_CORRECTOR.csv")
    roots = read("L3_INITIAL_CORRECTOR_ROOTS.csv")
    counts = read("L3_INITIAL_CORRECTOR_COUNTS.csv")
    for summary in summaries:
        candidate_id = summary["candidate_id"]
        dest = HERE / "evaluations" / candidate_id
        dest.mkdir(parents=True, exist_ok=True)
        relevant_roots = [r for r in roots if r["candidate_id"] == candidate_id]
        relevant_counts = [r for r in counts if r["candidate_id"] == candidate_id]
        assert len(relevant_roots) >= 2 and len(relevant_counts) == 3
        write(dest / "ROOTS.csv", relevant_roots)
        write(dest / "ROOT_COUNTS.csv", relevant_counts)
        write(dest / "SPECTRAL.csv", [summary])
        print(candidate_id, len(relevant_roots))


if __name__ == "__main__":
    main()
