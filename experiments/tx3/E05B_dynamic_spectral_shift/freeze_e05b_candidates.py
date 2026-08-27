from __future__ import annotations

import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E05B = ROOT / "experiments" / "tx3" / "E05B_dynamic_spectral_shift"
sys.path[:0] = [str(ROOT), str(E05B)]

from dynamic_spectral_shift import (  # noqa: E402
    contour_for_reference,
    exact_contour_moments,
    mobius_vertices,
    parse_subset_key,
    positive_oscillatory_poles,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"
OUTPUT = ARTIFACT / "preregistration" / "E05B_CANDIDATE_CONTOUR_FREEZE.json"
MANDATORY = ("A7-A8", "A3-A6", "A2-A7-A8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_spectra(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return (
            data["operating_point_ids"].astype(str).tolist(),
            data["subset_keys"].astype(str).tolist(),
            np.asarray(data["eigenvalues"], dtype=np.complex128),
        )


def main() -> int:
    holdout_result = ARTIFACT / "holdout" / "component_externalities.parquet"
    if holdout_result.exists():
        raise RuntimeError("cannot freeze candidates after E05B holdout results exist")
    cases_path = ARTIFACT / "development" / "component_externalities.parquet"
    spectra_path = ARTIFACT / "spectra" / "development_spectra.npz"
    if not cases_path.is_file() or not spectra_path.is_file():
        raise FileNotFoundError("development decomposition must complete before candidate freeze")
    gates = json.loads((ARTIFACT / "preregistration" / "E05B_GATE_FREEZE.json").read_text(encoding="utf-8"))
    numerical = json.loads((ARTIFACT / "preregistration" / "E05B_NUMERICAL_FREEZE.json").read_text(encoding="utf-8"))
    cases = pd.read_parquet(cases_path)
    eligible: list[dict[str, Any]] = []
    for (coalition, order), group in cases.groupby(["coalition_key", "order"]):
        resolved = group[group["resolved_pole_externality"]]
        same_sign = max((resolved["delta_phi_poles"] > 0).mean(), (resolved["delta_phi_poles"] < 0).mean()) if len(resolved) else 0.0
        if len(resolved) < 12 or same_sign < 0.75:
            continue
        eligible.append(
            {
                "coalition_key": str(coalition),
                "order": int(order),
                "resolved_development_points": len(resolved),
                "same_sign_fraction": float(same_sign),
                "median_delta_phi_poles": float(resolved["delta_phi_poles"].median()),
                "median_absolute_delta_phi_poles": float(resolved["delta_phi_poles"].abs().median()),
                "median_algebraic_fraction_absolute": float(resolved["algebraic_fraction_absolute"].median()),
                "median_pole_fraction_absolute": float(resolved["pole_fraction_absolute"].median()),
            }
        )
    ranked = {
        order: sorted(
            [record for record in eligible if record["order"] == order],
            key=lambda record: (-record["median_absolute_delta_phi_poles"], record["coalition_key"]),
        )
        for order in (2, 3)
    }
    selected_keys: list[str] = []
    for order in (2, 3):
        selected_keys.extend(record["coalition_key"] for record in ranked[order][:2])
    selected_keys.extend(MANDATORY)
    selected_keys = list(dict.fromkeys(selected_keys))
    statistics = {record["coalition_key"]: record for record in eligible}

    op_ids, subset_keys, array = load_spectra(spectra_path)
    op_index = op_ids.index("OP00")
    subset_index = {key: index for index, key in enumerate(subset_keys)}
    baseline = array[op_index, subset_index["EMPTY"]]
    reference_modes = positive_oscillatory_poles(
        baseline, tuple(float(value) for value in numerical["contour_reference_band_hz"])
    )
    gap_fraction = float(numerical["contour_gap_fraction"])
    candidates: list[dict[str, Any]] = []
    for sequence, coalition_key in enumerate(selected_keys, start=1):
        coalition = parse_subset_key(coalition_key)
        spectra = {
            subset: array[op_index, subset_index["EMPTY" if not subset else "-".join(f"A{x}" for x in subset)]]
            for subset in mobius_vertices(coalition)
        }
        modal_choices: list[dict[str, Any]] = []
        for mode in reference_modes:
            center, radius = contour_for_reference(baseline, complex(mode), gap_fraction)
            exact = exact_contour_moments(spectra, coalition, center, radius)
            if not exact["single_pole_each_vertex"] or exact["mu0"] != 0:
                continue
            modal_choices.append(
                {
                    "reference_eigenvalue_real_per_s": float(center.real),
                    "reference_eigenvalue_imag_per_s": float(center.imag),
                    "reference_frequency_hz": float(center.imag / (2.0 * np.pi)),
                    "op00_contour_radius_per_s": float(radius),
                    "op00_mu1_real_per_s": float(exact["mu1"].real),
                    "op00_mu1_imag_per_s": float(exact["mu1"].imag),
                    "op00_mu1_absolute_per_s": float(abs(exact["mu1"])),
                }
            )
        if not modal_choices:
            raise RuntimeError(f"no valid isolated OP00 modal contour for {coalition_key}")
        selected_mode = max(modal_choices, key=lambda record: record["op00_mu1_absolute_per_s"])
        candidates.append(
            {
                "candidate_id": f"D{sequence:02d}",
                "coalition_key": coalition_key,
                "actions": [f"A{value}" for value in coalition],
                "order": len(coalition),
                "mandatory_E05_diagnostic": coalition_key in MANDATORY,
                "selected_by_top_dynamic_rule": coalition_key in {
                    record["coalition_key"] for order in (2, 3) for record in ranked[order][:2]
                },
                "development_statistics": statistics.get(coalition_key),
                "modal_selection_valid_choice_count": len(modal_choices),
                "modal_selection_rule": "largest absolute OP00 connected mu1 among isolated 0.1-30 Hz single-pole contours",
                **selected_mode,
            }
        )
    payload = {
        "freeze_id": "TX3-E05B-CANDIDATE-CONTOUR-1.0",
        "created_utc": datetime.now(UTC).isoformat(),
        "development_population": "all 24 already-seen E03 operating points; treated as development only",
        "holdout_inspected_before_freeze": False,
        "eligible_pair_count": len(ranked[2]),
        "eligible_triple_count": len(ranked[3]),
        "top_dynamic_pairs": [record["coalition_key"] for record in ranked[2][:2]],
        "top_dynamic_triples": [record["coalition_key"] for record in ranked[3][:2]],
        "mandatory_diagnostics": list(MANDATORY),
        "candidates": candidates,
        "development_component_sha256": sha256(cases_path),
        "development_spectra_sha256": sha256(spectra_path),
        "selection_rule": gates["development_selection_rule"],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
