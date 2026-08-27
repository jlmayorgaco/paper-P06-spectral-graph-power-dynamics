from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import andes
import numpy as np
import pandas as pd
from scipy.stats import qmc


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05C = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned"
sys.path[:0] = [str(ROOT), str(E03), str(E05C)]

from campaign_core import action_vector  # noqa: E402
from margin_conditioned import (  # noqa: E402
    STRESS_PARTICIPATION,
    SYNCHRONOUS_PV_BUSES,
    biorthogonal_mac,
    complex_externality,
    eigenvector_hash,
    evaluate_stressed_operator,
    generalized_modes,
    mobius_vertices,
    parse_subset_key,
    positive_oscillatory_indices,
    sha256,
    subset_key,
    track_assignment,
    union_vertices,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05C_margin_conditioned"
PREREG = ARTIFACT / "preregistration"
REPORTS = ARTIFACT / "reports"
E05B = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"
ACTION_FREEZE = ROOT / "artifacts" / "tx3" / "E03_E04" / "preregistration" / "ACTION_SET_FREEZE.json"
SEED = 20260829


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def state_names() -> list[str]:
    system = andes.load(str(ROOT / "systems/ieee39_tx3_gfl3/andes_case.xlsx"), setup=True, no_output=True)
    if not system.PFlow.run() or not system.PFlow.converged:
        raise RuntimeError("state-name reference PF failed")
    system.TDS.init()
    return [str(value) for value in system.dae.x_name]


def main() -> int:
    for name in (
        "preregistration", "baseline_calibration", "holdout", "tables", "figures", "logs",
        "manifests", "reports", "cache",
    ):
        (ARTIFACT / name).mkdir(parents=True, exist_ok=True)
    if (ARTIFACT / "holdout" / "e05c_connected_pole_externalities.parquet").exists():
        raise RuntimeError("cannot preregister after E05C coalition outcomes exist")
    git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    source_paths = [
        ROOT / "experiments/tx3/E05C_margin_conditioned/margin_conditioned.py",
        ROOT / "experiments/tx3/E05C_margin_conditioned/calibrate_e05c_stress.py",
        ROOT / "experiments/tx3/E05C_margin_conditioned/run_e05c_campaign.py",
        ROOT / "tests/unit/test_e05c_margin_conditioned.py",
    ]
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing E05C frozen sources: {missing}")

    dynamic_path = E05B / "tables" / "holdout_dynamic_reproducibility.parquet"
    dynamic = pd.read_parquet(dynamic_path)
    dynamic = dynamic[dynamic["reproducible_dynamic_externality"]].sort_values(
        ["order", "median_absolute_delta_phi_poles"], ascending=[True, False]
    )
    if len(dynamic) != 15 or (dynamic["order"] == 2).sum() != 10 or (dynamic["order"] == 3).sum() != 5:
        raise RuntimeError("frozen E05B reproducible population is not the expected 10 pairs + 5 triples")
    candidate_rows = []
    for sequence, row in enumerate(dynamic.itertuples(index=False), start=1):
        candidate_rows.append(
            {
                "candidate_id": f"MC{sequence:02d}",
                "coalition_key": row.coalition_key,
                "actions": row.coalition_key.split("-"),
                "order": int(row.order),
                "e05b_median_delta_phi_poles": float(row.median_delta_phi_poles),
                "e05b_sign": "positive" if row.median_delta_phi_poles > 0 else "negative",
                "e05b_resolved_holdout_points": int(row.resolved_holdout_points),
                "e05b_same_sign_fraction": float(row.same_sign_fraction),
                "e05b_median_absolute_delta_phi_poles": float(row.median_absolute_delta_phi_poles),
            }
        )
    candidates = pd.DataFrame(candidate_rows)
    candidates.to_parquet(PREREG / "e05c_candidates.parquet", index=False)
    candidate_payload = {
        "freeze_id": "TX3-E05C-CANDIDATES-1.0",
        "created_utc": datetime.now(UTC).isoformat(),
        "source_artifact": str(dynamic_path.relative_to(ROOT)).replace("\\", "/"),
        "source_artifact_sha256": sha256(dynamic_path),
        "selection_rule": "every E05B row with reproducible_dynamic_externality == true; no E05C outcome inspected",
        "candidate_count": len(candidate_rows),
        "pair_count": int((candidates["order"] == 2).sum()),
        "triple_count": int((candidates["order"] == 3).sum()),
        "candidates": candidate_rows,
    }
    write_json(PREREG / "E05C_CANDIDATE_FREEZE.json", candidate_payload)

    samples = qmc.Sobol(d=4, scramble=True, seed=SEED).random_base2(m=3)
    lower = np.asarray([0.94, 0.97, 0.97, -10.0])
    upper = np.asarray([1.02, 1.03, 1.03, 10.0])
    scaled = qmc.scale(samples, lower, upper)
    seeds = [
        {
            "seed_id": f"C{index:02d}",
            "operating_point_id": f"C{index:02d}",
            "split": "e05c_blind_holdout",
            "load_scale": float(row[0]),
            "gfl_dispatch_scale": float(row[1]),
            "reactive_load_scale": float(row[2]),
            "redispatch_mw": float(row[3]),
            "role": "independent_E05C_scrambled_Sobol_seed",
        }
        for index, row in enumerate(scaled)
    ]
    pd.DataFrame(seeds).to_parquet(PREREG / "E05C_BASE_HOLDOUT_MANIFEST.parquet", index=False)
    seed_payload = {
        "freeze_id": "TX3-E05C-BASE-HOLDOUT-1.0",
        "seed": SEED,
        "generator": "independent scrambled Sobol base-2 draw of eight points",
        "ranges": {
            "load_scale": [0.94, 1.02],
            "gfl_dispatch_scale": [0.97, 1.03],
            "reactive_load_scale": [0.97, 1.03],
            "redispatch_mw": [-10.0, 10.0],
        },
        "range_basis": "strict interior of the previously feasible E03/E05B region; no E05C outcomes used",
        "old_points_reused": False,
        "outcomes_inspected_before_freeze": False,
        "operating_points": seeds,
    }
    write_json(PREREG / "E05C_BASE_HOLDOUT_FREEZE.json", seed_payload)

    coalitions = [parse_subset_key(value) for value in candidates["coalition_key"]]
    vertices = union_vertices(coalitions)
    op00 = {
        "seed_id": "OP00",
        "operating_point_id": "OP00",
        "load_scale": 1.0,
        "gfl_dispatch_scale": 1.0,
        "reactive_load_scale": 1.0,
        "redispatch_mw": 0.0,
    }
    empty_point = evaluate_stressed_operator(op00, 1.0, {})
    if empty_point.status != "SUCCESS":
        raise RuntimeError(f"development OP00 empty vertex failed: {empty_point.status}")
    baseline = generalized_modes(empty_point)
    reference_indices = positive_oscillatory_indices(baseline, (0.1, 30.0))
    reference_indices = reference_indices[np.argsort(baseline.eigenvalues[reference_indices].imag)]
    reference_ids = {int(index): f"R{sequence:03d}" for sequence, index in enumerate(reference_indices, start=1)}
    endpoint_modes: dict[tuple[int, ...], Any] = {(): baseline}
    endpoint_tracking: dict[tuple[int, ...], dict[int, tuple[int, float]]] = {
        (): {int(index): (int(index), 1.0) for index in reference_indices}
    }
    for subset in vertices:
        if not subset:
            continue
        previous = baseline
        current_map = {int(index): int(index) for index in reference_indices}
        minimum_score = {int(index): 1.0 for index in reference_indices}
        for beta in (0.25, 0.5, 0.75, 1.0):
            point = evaluate_stressed_operator(
                op00, 1.0, action_vector({f"A{value}": beta for value in subset})
            )
            if point.status != "SUCCESS":
                raise RuntimeError(
                    f"development homotopy {subset_key(subset)} beta={beta} failed: {point.status}"
                )
            target = generalized_modes(point)
            previous_indices = positive_oscillatory_indices(previous, (0.05, 45.0))
            assignment = track_assignment(previous, target, previous_indices)
            for reference_index, current_index in list(current_map.items()):
                match = assignment.get(current_index)
                if match is None:
                    current_map.pop(reference_index)
                    minimum_score.pop(reference_index)
                    continue
                current_map[reference_index] = match[0]
                minimum_score[reference_index] = min(minimum_score[reference_index], match[1])
            previous = target
        endpoint_modes[subset] = previous
        endpoint_tracking[subset] = {
            reference_index: (target_index, minimum_score[reference_index])
            for reference_index, target_index in current_map.items()
        }
    names = state_names()
    mode_rows = []
    for candidate in candidate_rows:
        coalition = parse_subset_key(candidate["coalition_key"])
        required = mobius_vertices(coalition)
        eligible = []
        for reference_index in reference_indices:
            values: dict[tuple[int, ...], complex] = {(): complex(baseline.eigenvalues[reference_index])}
            minimum_bmac = 1.0
            complete = True
            for subset in required:
                if not subset:
                    continue
                match = endpoint_tracking[subset].get(int(reference_index))
                if match is None or match[1] < 0.95:
                    complete = False
                    break
                target_index, score = match
                minimum_bmac = min(minimum_bmac, score)
                values[subset] = complex(endpoint_modes[subset].eigenvalues[target_index])
            if not complete:
                continue
            delta = complex_externality(values, coalition)
            eligible.append((abs(delta.real), int(reference_index), minimum_bmac, delta))
        if not eligible:
            mode_rows.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "coalition_key": candidate["coalition_key"],
                    "order": candidate["order"],
                    "lock_status": "NO_LOCKABLE_MODE",
                }
            )
            continue
        _, reference_index, minimum_bmac, delta = max(eligible, key=lambda item: (item[0], -item[1]))
        participation = np.abs(baseline.left[:, reference_index].conj() * baseline.right[:, reference_index])
        participation /= max(float(np.sum(participation)), 1e-30)
        top = np.argsort(participation)[-6:][::-1]
        value = complex(baseline.eigenvalues[reference_index])
        mode_rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "coalition_key": candidate["coalition_key"],
                "order": candidate["order"],
                "lock_status": "LOCKED",
                "modal_family_id": reference_ids[reference_index],
                "reference_mode_index": reference_index,
                "reference_eigenvalue_real_per_s": float(value.real),
                "reference_eigenvalue_imag_per_s": float(value.imag),
                "reference_frequency_hz": float(value.imag / (2.0 * np.pi)),
                "reference_damping_ratio": float(-value.real / abs(value)),
                "minimum_development_vertex_bmac": float(minimum_bmac),
                "development_delta_lambda_real_per_s": float(delta.real),
                "development_delta_lambda_imag_per_s": float(delta.imag),
                "development_delta_lambda_absolute_per_s": float(abs(delta)),
                "left_right_eigenvectors_sha256": eigenvector_hash(baseline, reference_index),
                "dominant_state_participation_json": json.dumps(
                    [
                        {"state_index": int(index), "state_name": names[index], "normalized_participation": float(participation[index])}
                        for index in top
                    ],
                    sort_keys=True,
                ),
            }
        )
    modes = pd.DataFrame(mode_rows)
    modes.to_parquet(PREREG / "e05c_mode_family_freeze.parquet", index=False)
    mode_payload = {
        "freeze_id": "TX3-E05C-MODE-FAMILIES-1.0",
        "reference_operating_point": "OP00",
        "reference_band_hz": [0.1, 30.0],
        "development_vertex_bmac_minimum": 0.95,
        "selection_rule": "largest absolute real connected pole externality at OP00 among modes tracked with BMAC >= 0.95 at every required endpoint vertex",
        "holdout_inspected_before_freeze": False,
        "locked_count": int((modes["lock_status"] == "LOCKED").sum()),
        "no_lockable_count": int((modes["lock_status"] != "LOCKED").sum()),
        "mode_families": mode_rows,
    }
    write_json(PREREG / "E05C_MODE_FAMILY_FREEZE.json", mode_payload)

    stress = {
        "freeze_id": "TX3-E05C-STRESS-RULE-1.0",
        "primary_coordinate": "uniform constant-power-factor scaling of every seeded PQ load",
        "mapping": "P_L(gamma)=gamma*P_L(seed); Q_L(gamma)=gamma*Q_L(seed)",
        "gamma_start": 1.0,
        "coarse_increment": 0.02,
        "gamma_cap": 1.50,
        "boundary_bisection_tolerance": 0.0025,
        "safe_endpoint_fraction": 0.95,
        "tau_grid": [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0],
        "synchronous_participation_buses": list(SYNCHRONOUS_PV_BUSES),
        "synchronous_participation_factors": [float(value) for value in STRESS_PARTICIPATION],
        "participation_basis": "fixed proportional canonical active-power contribution of retained non-reference synchronous PV generators",
        "slack_rule": "slack absorbs only residual PF mismatch",
        "gfl_reference_rule": "fixed along gamma except explicit A7 vertex",
        "hard_stops": ["PF_FAIL", "INIT_FAIL", "voltage_min<0.85", "voltage_max>1.15", "generator_limit", "LIMITER_ACTIVE", "NONSMOOTH_REGIME", "EIG_FAIL", "EIGENPAIR_RESIDUAL_FAIL", "BASELINE_UNSTABLE"],
        "coalitions_inspected_during_endpoint_calibration": False,
    }
    write_json(PREREG / "E05C_STRESS_RULE_FREEZE.json", stress)
    numerical = {
        "freeze_id": "TX3-E05C-NUMERICS-1.0",
        "git_sha_before_preregistration_commit": git_sha,
        "mode_solver": "dense generalized eig of exact algebraic Schur complement (Ar,M), left and right vectors",
        "development_mode_bmac_minimum": 0.95,
        "stress_and_action_tracking_bmac_minimum": 0.90,
        "stress_refinement": "recursive half-step then quarter-step continuation; no manual switch",
        "action_homotopy_nodes": [0.0, 0.25, 0.5, 0.75, 1.0],
        "action_homotopy_refinement_nodes": [9, 17],
        "eigenpair_normalized_residual_maximum": 1e-10,
        "pseudo_mode_filter_radius": 1e-5,
        "target_tracking_band_hz": [0.05, 45.0],
        "dynamic_pole_factor_band_hz": [0.1, 30.0],
        "dynamic_pole_factor_nodes": 96,
        "composition_closure_tolerance": 1e-10,
        "source_hashes": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in source_paths},
        "new_paraemt_runs_permitted": False,
        "E06_execution_permitted": False,
    }
    write_json(PREREG / "E05C_NUMERICAL_FREEZE.json", numerical)
    gates = {
        "freeze_id": "TX3-E05C-GATES-1.0",
        "C3b_A": {
            "valid_seed_paths_minimum": 6,
            "high_stress_levels": [0.875, 1.0],
            "destabilizing_sign_fraction_minimum": 0.75,
            "median_final_rho_comp_minimum": 0.10,
            "median_high_to_nominal_amplification_minimum": 2.0,
            "smooth_only": True,
        },
        "C3b_B_route_1": {"replicated_hard_or_screen_failures_minimum": 2},
        "C3b_B_route_2": {
            "seed_paths_minimum": 4,
            "eligible_tau_levels": [0.875, 1.0],
            "rho_comp_minimum": 0.25,
            "delta_zeta_maximum": -0.0025,
            "destabilizing_sign_fraction_minimum": 0.75,
        },
        "screen_damping_ratio": 0.05,
        "old_E05_absolute_delta_zeta_threshold": 0.0025,
        "C3a_reopened": False,
        "C3_original_assessed": False,
        "E06_consideration_rule": "true only if C3b-B is SUPPORTED; execution remains forbidden",
    }
    write_json(PREREG / "E05C_GATE_FREEZE.json", gates)
    top_panels = {
        "freeze_id": "TX3-E05C-FIGURE-PANELS-1.0",
        "selection_source": "E05B frozen dynamic ranking only",
        "primary_panel_candidates": candidate_rows[:2] + candidate_rows[10:12],
    }
    write_json(PREREG / "E05C_FIGURE_PANEL_FREEZE.json", top_panels)

    prereg_text = f"""# E05C margin-conditioned connected pole externalities preregistration

**Freeze date:** 2026-08-27. **Git SHA before freeze:** `{git_sha}`.

E05C preserves C1/C2 supported, C3a rejected, E05B D1/D2/D3 supported, and original cycle/SCC C3
unassessed. It tests a new margin-conditioned claim using all 15 E05B-reproducible coalitions, one
development-locked mode per coalition, uniform constant-power-factor load stress, and eight new Sobol
seeds. A1--A8 and the old `|Delta_S zeta| >= 0.0025` threshold are immutable.

Stress endpoints are calibrated with the EMPTY coalition only. No action vertex may be evaluated until
`E05C_BASELINE_STRESS_LIMITS.parquet` is frozen and committed. The nine tau levels, redispatch factors,
tracking/refinement rules, failure statuses, C3b-A amplification gate, and C3b-B two routes are fixed in
the accompanying JSON files. E05C uses fully re-equilibrated ANDES only; no ParaEMT or E06 execution is
permitted. Even a positive C3b-B result stops for independent review.
"""
    (PREREG / "E05C_PREREGISTRATION.md").write_text(prereg_text, encoding="utf-8")
    (REPORTS / "E05C_CANDIDATE_FREEZE.md").write_text(
        f"# E05C candidate freeze\n\nAll {len(candidate_rows)} E05B-reproducible coalitions are frozen: 10 pairs and 5 triples. Source SHA-256: `{sha256(dynamic_path)}`. No candidate was removed for small effect.\n",
        encoding="utf-8",
    )
    (REPORTS / "E05C_MODE_FAMILY_FREEZE.md").write_text(
        f"# E05C mode-family freeze\n\n{mode_payload['locked_count']}/15 candidates have one OP00 development mode locked by the preregistered real connected-pole ranking with endpoint BMAC >= 0.95; {mode_payload['no_lockable_count']} are `NO_LOCKABLE_MODE`.\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "FROZEN", "candidates": len(candidate_rows), "locked_modes": mode_payload["locked_count"], "new_seeds": 8}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
