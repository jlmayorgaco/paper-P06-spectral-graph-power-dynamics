"""Compare the canonical Python and Julia 2x2 terminal transfer maps."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from campaign_root import campaign_root_from_argv

CAMPAIGN = campaign_root_from_argv()
RAW = CAMPAIGN / "raw" / "gfl11"


def load(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    frequencies = np.array([float(r["frequency_hz"]) for r in rows])
    maps = np.empty((len(rows), 2, 2), dtype=complex)
    for k, row in enumerate(rows):
        maps[k] = np.array(
            [
                [float(row["h11_re"]) + 1j * float(row["h11_im"]), float(row["h12_re"]) + 1j * float(row["h12_im"])],
                [float(row["h21_re"]) + 1j * float(row["h21_im"]), float(row["h22_re"]) + 1j * float(row["h22_im"])],
            ]
        )
    sigma_min_resolvent = np.array([float(r["sigma_min"]) for r in rows])
    return frequencies, maps, sigma_min_resolvent


def main() -> int:
    f_python, h_python, s_python = load(RAW / "p1_transfer_canonical_python.csv")
    f_julia, h_julia, s_julia = load(RAW / "p1_transfer_julia.csv")
    if len(f_python) != len(f_julia) or np.max(np.abs(f_python - f_julia)) > 1e-10:
        raise RuntimeError("canonical and Julia frequency grids differ")

    abs_frobenius = np.linalg.norm(h_python - h_julia, axis=(1, 2))
    ref_frobenius = np.linalg.norm(h_python, axis=(1, 2))
    rel_frobenius = abs_frobenius / np.maximum(ref_frobenius, 1e-14)
    singular_python = np.linalg.svd(h_python, compute_uv=False)
    singular_julia = np.linalg.svd(h_julia, compute_uv=False)
    singular_abs = np.abs(singular_python - singular_julia)
    singular_rel = singular_abs / np.maximum(singular_python, 1e-14)
    element_phase = np.angle(h_python) - np.angle(h_julia)
    element_phase = np.angle(np.exp(1j * element_phase))
    element_magnitude = np.maximum(np.abs(h_python), np.abs(h_julia))
    phase_error = np.max(np.where(element_magnitude > 1e-10, np.abs(element_phase), 0.0), axis=(1, 2))
    near_singular = np.minimum(s_python, s_julia) < 0.1
    abs_frobenius_near = abs_frobenius[near_singular]
    backward_error = abs_frobenius / np.maximum(1.0, ref_frobenius)

    def percentile(values: np.ndarray, q: float) -> float:
        return float(np.percentile(values, q)) if len(values) else float("nan")

    summary = {
        "status": "PASS" if percentile(rel_frobenius, 50) < 1e-6 and percentile(rel_frobenius, 99) < 1e-4 else "FAIL",
        "points": int(len(f_python)),
        "frequency_min_hz": float(f_python.min()),
        "frequency_max_hz": float(f_python.max()),
        "dense_0p2_2_hz_points": int(np.count_nonzero((f_python >= 0.2) & (f_python <= 2.0))),
        "frobenius_abs_median": percentile(abs_frobenius, 50),
        "frobenius_abs_p99": percentile(abs_frobenius, 99),
        "frobenius_relative_median": percentile(rel_frobenius, 50),
        "frobenius_relative_p99": percentile(rel_frobenius, 99),
        "largest_singular_relative_median": percentile(singular_rel[:, 0], 50),
        "largest_singular_relative_p99": percentile(singular_rel[:, 0], 99),
        "smallest_singular_relative_median": percentile(singular_rel[:, 1], 50),
        "smallest_singular_relative_p99": percentile(singular_rel[:, 1], 99),
        "phase_error_rad_median": percentile(phase_error, 50),
        "phase_error_rad_p99": percentile(phase_error, 99),
        "near_singular_points": int(np.count_nonzero(near_singular)),
        "near_singular_abs_frobenius_p99": percentile(abs_frobenius_near, 99),
        "backward_error_median": percentile(backward_error, 50),
        "backward_error_p99": percentile(backward_error, 99),
        "target_median_relative_lt_1e-6": bool(percentile(rel_frobenius, 50) < 1e-6),
        "target_p99_relative_lt_1e-4": bool(percentile(rel_frobenius, 99) < 1e-4),
    }
    (RAW / "p1_transfer_comparison.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (RAW / "p1_transfer_comparison.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frequency_hz", "frobenius_abs", "frobenius_relative", "largest_singular_relative", "smallest_singular_relative", "phase_error_rad", "near_singular", "backward_error"])
        for i, frequency in enumerate(f_python):
            writer.writerow([frequency, abs_frobenius[i], rel_frobenius[i], singular_rel[i, 0], singular_rel[i, 1], phase_error[i], int(near_singular[i]), backward_error[i]])
    print(
        f"P1_TRANSFER_{summary['status']} points={summary['points']} "
        f"median_rel={summary['frobenius_relative_median']:.3e} "
        f"p99_rel={summary['frobenius_relative_p99']:.3e}"
    )
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
