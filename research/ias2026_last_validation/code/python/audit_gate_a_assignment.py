"""Independent modal reconciliation for the PowerDynamics Gate A paths."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


REPO = Path(__file__).resolve().parents[4]
CAMPAIGN = REPO / "research" / "ias2026_bulletproof"
RAW = CAMPAIGN / "raw" / "powerdynamics"
rows = []
with (RAW / "gate_a_spectra.csv").open(newline="", encoding="utf-8") as handle:
    rows.extend(csv.DictReader(handle))

groups: dict[tuple[str, str], list[complex]] = {}
for row in rows:
    key = (row["path"], row["tolerance"])
    groups.setdefault(key, []).append(complex(float(row["real"]), float(row["imag"])))
reference_key = next((key for key in groups if key[1] == "1.0e-12"), next(iter(groups)))
reference = np.asarray(groups[reference_key], dtype=complex)
out = []
for key, spectrum in groups.items():
    path, tolerance = key
    current = np.asarray(spectrum, dtype=complex)
    if len(current) != len(reference):
        out.append((path, tolerance, len(current), np.nan, np.nan, "FAIL"))
        continue
    cost = np.abs(reference[:, None] - current[None, :])
    reference_idx, current_idx = linear_sum_assignment(cost)
    matched = cost[reference_idx, current_idx]
    out.append((path, tolerance, len(current), float(np.max(matched)), float(np.mean(matched)),
                "PASS" if np.max(matched) <= 1e-6 else "FAIL"))

with (RAW / "gate_a_assignment.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["reference_path", "reference_tolerance", "path", "tolerance", "eigenvalue_count", "max_hungarian_delta", "mean_hungarian_delta", "status"])
    for path, tolerance, count, maximum, mean, status in out:
        writer.writerow([reference_key[0], reference_key[1], path, tolerance, count, maximum, mean, status])

passed = all(row[-1] == "PASS" for row in out)
(RAW / "gate_a_assignment_report.md").write_text(
    "# PowerDynamics Gate A independent modal assignment\n\n"
    f"status: {'PASS' if passed else 'FAIL'}\n"
    "method: scipy.optimize.linear_sum_assignment on pairwise complex-eigenvalue distance\n"
    f"reference_path: {reference_key[0]}\n"
    f"reference_tolerance: {reference_key[1]}\n"
    "lexicographic_entrywise_comparison: not used for this reconciliation\n"
    "threshold: maximum matched complex-eigenvalue distance <= 1e-6\n"
    "source: raw/powerdynamics/gate_a_spectra.csv\n"
    "result: raw/powerdynamics/gate_a_assignment.csv\n",
    encoding="utf-8",
)
print(f"GATE_A_ASSIGNMENT_{'PASS' if passed else 'FAIL'} paths={len(out)}")
raise SystemExit(0 if passed else 1)
