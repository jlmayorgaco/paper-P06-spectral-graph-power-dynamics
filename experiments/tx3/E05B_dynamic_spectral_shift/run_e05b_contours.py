from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E05B = ROOT / "experiments" / "tx3" / "E05B_dynamic_spectral_shift"
sys.path[:0] = [str(ROOT), str(E05B)]

from dynamic_spectral_shift import (  # noqa: E402
    connected_resolvent_trace,
    contour_for_reference,
    exact_contour_moments,
    mobius_vertices,
    numerical_contour_moments,
    parse_subset_key,
    subset_key,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"


def load_spectra(population: str) -> tuple[list[str], list[str], np.ndarray]:
    with np.load(ARTIFACT / "spectra" / f"{population}_spectra.npz", allow_pickle=False) as data:
        return (
            data["operating_point_ids"].astype(str).tolist(),
            data["subset_keys"].astype(str).tolist(),
            np.asarray(data["eigenvalues"], dtype=np.complex128),
        )


def main() -> int:
    freeze = json.loads(
        (ARTIFACT / "preregistration" / "E05B_CANDIDATE_CONTOUR_FREEZE.json").read_text(encoding="utf-8")
    )
    numerical = json.loads(
        (ARTIFACT / "preregistration" / "E05B_NUMERICAL_FREEZE.json").read_text(encoding="utf-8")
    )
    nodes = int(numerical["contour_nodes"])
    gap_fraction = float(numerical["contour_gap_fraction"])
    rows: list[dict[str, Any]] = []
    xi_rows: list[dict[str, Any]] = []
    for population in ("development", "holdout"):
        op_ids, subset_keys, array = load_spectra(population)
        subset_index = {key: index for index, key in enumerate(subset_keys)}
        for op_position, op_id in enumerate(op_ids):
            baseline = array[op_position, subset_index["EMPTY"]]
            for candidate in freeze["candidates"]:
                coalition = parse_subset_key(candidate["coalition_key"])
                spectra = {
                    subset: array[op_position, subset_index[subset_key(subset)]]
                    for subset in mobius_vertices(coalition)
                }
                reference = complex(
                    candidate["reference_eigenvalue_real_per_s"],
                    candidate["reference_eigenvalue_imag_per_s"],
                )
                center, radius = contour_for_reference(baseline, reference, gap_fraction)
                exact = exact_contour_moments(spectra, coalition, center, radius)
                numerical_mu0, numerical_mu1 = numerical_contour_moments(
                    spectra, coalition, center, radius, nodes
                )
                error = max(abs(numerical_mu0 - exact["mu0"]), abs(numerical_mu1 - exact["mu1"]))
                uncertainty = float(error + (2 ** len(coalition)) * 1e-10)
                dominant = "real" if abs(exact["mu1"].real) >= abs(exact["mu1"].imag) else "imag"
                dominant_value = exact["mu1"].real if dominant == "real" else exact["mu1"].imag
                rows.append(
                    {
                        "operating_point_id": op_id,
                        "population": population,
                        "candidate_id": candidate["candidate_id"],
                        "coalition_key": candidate["coalition_key"],
                        "order": len(coalition),
                        "contour_center_real_per_s": float(center.real),
                        "contour_center_imag_per_s": float(center.imag),
                        "contour_radius_per_s": float(radius),
                        "center_frequency_hz": float(center.imag / (2.0 * math.pi)),
                        "single_pole_each_vertex": bool(exact["single_pole_each_vertex"]),
                        "vertex_counts_json": json.dumps(exact["vertex_counts"], sort_keys=True),
                        "mu0_real": float(exact["mu0"].real),
                        "mu0_imag": float(exact["mu0"].imag),
                        "mu1_real_per_s": float(exact["mu1"].real),
                        "mu1_imag_per_s": float(exact["mu1"].imag),
                        "mu1_absolute_per_s": float(abs(exact["mu1"])),
                        "dominant_component": dominant,
                        "dominant_component_value_per_s": float(dominant_value),
                        "numerical_mu0_real": float(numerical_mu0.real),
                        "numerical_mu0_imag": float(numerical_mu0.imag),
                        "numerical_mu1_real_per_s": float(numerical_mu1.real),
                        "numerical_mu1_imag_per_s": float(numerical_mu1.imag),
                        "numerical_contour_error": float(error),
                        "mu1_uncertainty": uncertainty,
                        "resolved_mu1": bool(exact["single_pole_each_vertex"] and abs(exact["mu1"]) >= 5.0 * uncertainty),
                    }
                )
                theta = 2.0 * math.pi * np.arange(nodes, dtype=float) / nodes
                for grid_index, angle in enumerate(theta):
                    s = center + radius * np.exp(1j * angle)
                    xi = connected_resolvent_trace(spectra, coalition, s)
                    xi_rows.append(
                        {
                            "operating_point_id": op_id,
                            "population": population,
                            "candidate_id": candidate["candidate_id"],
                            "coalition_key": candidate["coalition_key"],
                            "grid_index": grid_index,
                            "s_real_per_s": float(s.real),
                            "s_imag_per_s": float(s.imag),
                            "xi_real_s": float(xi.real),
                            "xi_imag_s": float(xi.imag),
                        }
                    )
    moments = pd.DataFrame(rows).sort_values(["population", "operating_point_id", "candidate_id"])
    xi = pd.DataFrame(xi_rows).sort_values(["population", "operating_point_id", "candidate_id", "grid_index"])
    moments.to_parquet(ARTIFACT / "contours" / "connected_contour_moments.parquet", index=False)
    xi.to_parquet(ARTIFACT / "contours" / "connected_xi_samples.parquet", index=False)
    payload = {
        "status": "SUCCESS",
        "candidate_count": len(freeze["candidates"]),
        "moment_rows": len(moments),
        "xi_rows": len(xi),
        "maximum_numerical_contour_error": float(moments["numerical_contour_error"].max()),
        "single_pole_valid_rows": int(moments["single_pole_each_vertex"].sum()),
        "resolved_mu1_rows": int(moments["resolved_mu1"].sum()),
    }
    (ARTIFACT / "logs" / "contour_compute.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
