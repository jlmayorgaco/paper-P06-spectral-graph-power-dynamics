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
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E05D)]

from weak_grid import CASE, line_branch_records  # noqa: E402


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
PREREG = ARTIFACT / "preregistration"
E05B = ROOT / "artifacts" / "tx3" / "E05B_dynamic_spectral_shift"
ACTION_FREEZE = ROOT / "artifacts" / "tx3" / "E03_E04" / "preregistration" / "ACTION_SET_FREEZE.json"
SEED = 20260830
PREREG_NAMES = (
    "E05D_PREREGISTRATION.md", "E05D_CANDIDATE_FREEZE.json", "E05D_HOLDOUT_FREEZE.json",
    "E05D_GRID_STRESS_BRANCH_FREEZE.json", "E05D_STRESS_RULE_FREEZE.json",
    "E05D_CRITICAL_MODE_RULE_FREEZE.json", "E05D_TRACKING_FREEZE.json",
    "E05D_C3C_GATE_FREEZE.json", "E05D_NUMERICAL_FREEZE.json",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    for name in ("preregistration", "baseline_calibration", "holdout", "tables", "figures", "logs", "manifests", "cache"):
        (ARTIFACT / name).mkdir(parents=True, exist_ok=True)
    if (ARTIFACT / "holdout" / "e05d_connected_externalities.parquet").exists():
        raise RuntimeError("cannot preregister after E05D coalition outcomes exist")
    git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    source_paths = [
        E05D / "weak_grid.py", E05D / "preregister_e05d.py", E05D / "calibrate_e05d_stress.py",
        E05D / "freeze_e05d_endpoints_modes.py", E05D / "run_e05d_campaign.py",
        E05D / "finalize_e05d.py", E05D / "validate_e05d.py", E05D / "package_e05d.py",
        ROOT / "tests" / "unit" / "test_e05d_weak_grid.py",
    ]
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing: raise FileNotFoundError(f"missing E05D frozen sources: {missing}")

    dynamic_path = E05B / "tables" / "holdout_dynamic_reproducibility.parquet"
    dynamic = pd.read_parquet(dynamic_path)
    dynamic = dynamic[dynamic["reproducible_dynamic_externality"]].sort_values(
        ["order", "median_absolute_delta_phi_poles"], ascending=[True, False]
    )
    if len(dynamic) != 15 or int(dynamic["order"].eq(2).sum()) != 10 or int(dynamic["order"].eq(3).sum()) != 5:
        raise RuntimeError("E05B reproducible population is not exactly 10 pairs + 5 triples")
    candidates = []
    for sequence, row in enumerate(dynamic.itertuples(index=False), start=1):
        candidates.append({
            "candidate_id": f"WG{sequence:02d}", "coalition_key": row.coalition_key,
            "actions": row.coalition_key.split("-"), "order": int(row.order),
            "e05b_median_delta_phi_poles": float(row.median_delta_phi_poles),
            "e05b_sign": "positive" if row.median_delta_phi_poles > 0 else "negative",
            "e05b_same_sign_fraction": float(row.same_sign_fraction),
            "e05b_resolved_count": int(row.resolved_holdout_points),
            "source_sha256": sha256(dynamic_path),
        })
    pd.DataFrame(candidates).to_parquet(PREREG / "e05d_candidates.parquet", index=False)
    write_json(PREREG / "E05D_CANDIDATE_FREEZE.json", {
        "freeze_id": "TX3-E05D-CANDIDATES-1.0", "created_utc": datetime.now(UTC).isoformat(),
        "source_artifact": str(dynamic_path.relative_to(ROOT)).replace("\\", "/"),
        "source_artifact_sha256": sha256(dynamic_path),
        "selection_rule": "all and only E05B rows with reproducible_dynamic_externality == true",
        "candidate_count": 15, "pair_count": 10, "triple_count": 5,
        "coalition_outcomes_inspected_before_freeze": False, "candidates": candidates,
    })

    samples = qmc.Sobol(d=4, scramble=True, seed=SEED).random_base2(m=3)
    lower = np.asarray([0.94, 0.97, 0.97, -10.0])
    upper = np.asarray([1.02, 1.03, 1.03, 10.0])
    scaled = qmc.scale(samples, lower, upper)
    seeds = [{
        "seed_id": f"D{index:02d}", "operating_point_id": f"D{index:02d}",
        "sobol_index": index, "split": "e05d_blind_holdout",
        "load_scale": float(row[0]), "gfl_dispatch_scale": float(row[1]),
        "reactive_load_scale": float(row[2]), "redispatch_mw": float(row[3]),
        "role": "independent_E05D_scrambled_Sobol_seed",
    } for index, row in enumerate(scaled)]
    pd.DataFrame(seeds).to_parquet(PREREG / "e05d_holdout_manifest.parquet", index=False)
    write_json(PREREG / "E05D_HOLDOUT_FREEZE.json", {
        "freeze_id": "TX3-E05D-HOLDOUT-1.0", "seed": SEED,
        "previous_scientific_Sobol_seeds": [20260827, 20260828, 20260829],
        "generator": "independent scrambled Sobol base-2 draw of eight points",
        "ranges": {"load_scale": [0.94, 1.02], "gfl_dispatch_scale": [0.97, 1.03], "reactive_load_scale": [0.97, 1.03], "redispatch_mw": [-10.0, 10.0]},
        "range_basis": "identical established E05C operating-point dimensions and ranges; not widened",
        "canonical_case": str(CASE.relative_to(ROOT)).replace("\\", "/"),
        "canonical_case_sha256": sha256(CASE), "old_points_reused": False,
        "coalition_outcomes_inspected_before_freeze": False, "operating_points": seeds,
    })

    system = andes.load(str(CASE), setup=True, no_output=True)
    branches = line_branch_records(system)
    write_json(PREREG / "E05D_GRID_STRESS_BRANCH_FREEZE.json", {
        "freeze_id": "TX3-E05D-GRID-BRANCHES-1.0", "canonical_case_sha256": sha256(CASE),
        "topology_only_rule": "include every ANDES Line record with trans==0; exclude every fixed-tap transformer record trans==1; record local GFL interfaces and A8 with more specific exclusion reasons",
        "included_mapping": "R_l(kappa)=kappa*R_l0; X_l(kappa)=kappa*X_l0; line shunts unchanged; R/X preserved",
        "included_count": sum(record["included"] for record in branches),
        "excluded_count": sum(not record["included"] for record in branches),
        "alternative_branch_sets_compared": False, "dynamic_results_inspected": False, "branches": branches,
    })
    action_payload = json.loads(ACTION_FREEZE.read_text(encoding="utf-8"))
    write_json(PREREG / "E05D_STRESS_RULE_FREEZE.json", {
        "freeze_id": "TX3-E05D-STRESS-RULE-1.0", "primary_coordinate": "kappa_grid >= 1",
        "coordinate_count": 1, "mapping": "scale R and X of every frozen included external-grid line by kappa_grid",
        "unscaled": ["GFL internal/filter impedance", "synchronous-machine internal impedance", "controller parameters", "transformer tap ratios", "fixed shunts", "loads", "A8 branch 25-37"],
        "kappa_start": 1.0, "coarse_increment": 0.1,
        "coarse_kappa_grid": [round(1.0 + 0.1 * index, 10) for index in range(21)],
        "kappa_cap": 3.0, "boundary_bisection_tolerance": 0.005,
        "safe_endpoint_fraction_dynamic": 0.90, "safe_endpoint_fraction_physical": 0.95,
        "tau_grid": [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0],
        "baseline_calibration_population": "EMPTY coalition only",
        "first_boundary_priority": ["DAMPING_SCREEN_0.05", "BASELINE_UNSTABLE", "PF_FAIL", "VOLTAGE_LIMIT", "GENERATOR_LIMIT", "LIMITER_ACTIVE", "NONSMOOTH_REGIME", "EIG_FAIL", "EIGENPAIR_RESIDUAL_FAIL", "NUMERICAL_UNRESOLVED"],
        "action_source": str(ACTION_FREEZE.relative_to(ROOT)).replace("\\", "/"),
        "action_source_sha256": sha256(ACTION_FREEZE), "actions": action_payload["actions"],
        "coalition_outcomes_inspected_during_endpoint_calibration": False,
    })
    write_json(PREREG / "E05D_CRITICAL_MODE_RULE_FREEZE.json", {
        "freeze_id": "TX3-E05D-CRITICAL-MODE-RULE-1.0",
        "selection_population": "per-seed EMPTY coalition at tau=1 only",
        "selection_rule": "valid positive-imaginary oscillatory mode with minimum damping ratio in 0.1-30 Hz",
        "band_hz": [0.1, 30.0], "excluded": ["zero modes", "algebraic/infinite modes", "numerical artifacts"],
        "stored_fields": ["lambda", "frequency", "damping", "left eigenvector", "right eigenvector", "condition number", "dominant device/state participation", "eigenpair residual"],
        "freeze_before_coalition_results": True, "manual_replacement_permitted": False,
    })
    write_json(PREREG / "E05D_TRACKING_FREEZE.json", {
        "freeze_id": "TX3-E05D-TRACKING-1.0", "algorithm": "global one-to-one Hungarian assignment on biorthogonal MAC",
        "stress_direction": "tau=1 backward to tau=0", "BMAC_minimum": 0.90,
        "continuity_diagnostics": ["BMAC", "eigenvalue continuity", "frequency continuity"],
        "stress_refinement": "recursive half-step then quarter-step; otherwise MODE_TRACK_FAIL",
        "action_homotopy_nodes": [0.0, 0.25, 0.5, 0.75, 1.0],
        "action_refinement_node_counts": [9, 17], "all_active_actions_scaled_together": True,
        "manual_reassignment_permitted": False,
    })
    write_json(PREREG / "E05D_C3C_GATE_FREEZE.json", {
        "freeze_id": "TX3-E05D-C3C-GATES-1.0", "claim": "C3c_WEAK_GRID_MATERIAL_EXTERNALITY",
        "route_1": {"replicated_failure_seed_minimum": 2, "same_coalition_required": True, "eligible_failures": ["hard_composition_failure", "screen_composition_failure"]},
        "route_2": {"replicated_seed_minimum": 4, "same_coalition_required": True, "eligible_tau_levels": [0.875, 1.0], "rho_comp_minimum": 0.25, "delta_zeta_maximum": -0.0025, "destabilizing_sign_fraction_minimum": 0.75},
        "screen_damping_ratio": 0.05, "historical_absolute_delta_zeta_threshold": 0.0025,
        "SCR_is_diagnostic_not_gate": True, "amplification_alone_is_not_success": True,
        "if_supported": {"c3_mechanism_localization_authorized": True, "e06_consideration_authorized": False},
        "if_not_supported": {"c3_mechanism_localization_authorized": False, "e06_consideration_authorized": False, "materiality_search_closed": True, "E05E_authorized": False},
    })
    write_json(PREREG / "E05D_NUMERICAL_FREEZE.json", {
        "freeze_id": "TX3-E05D-NUMERICS-1.0", "git_sha_before_preregistration_commit": git_sha,
        "mode_solver": "dense generalized eig of exact algebraic Schur complement (Ar,M), left and right vectors",
        "stress_and_action_tracking_bmac_minimum": 0.90,
        "eigenpair_normalized_residual_maximum": 1e-10, "pseudo_mode_filter_radius": 1e-5,
        "tracking_band_hz": [0.05, 45.0], "critical_mode_band_hz": [0.1, 30.0],
        "composition_closure_tolerance": 1e-10,
        "valid_statuses": ["SUCCESS", "PF_FAIL", "INIT_FAIL", "EIG_FAIL", "MODE_TRACK_FAIL", "EIGENPAIR_RESIDUAL_FAIL", "NONSMOOTH_REGIME", "LIMITER_ACTIVE", "VOLTAGE_LIMIT", "GENERATOR_LIMIT", "BASELINE_UNSTABLE", "LOWER_ORDER_UNSTABLE", "NUMERICAL_UNRESOLVED"],
        "source_hashes": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in source_paths},
        "new_paraemt_runs_permitted": False, "E06_execution_permitted": False,
    })
    prereg_text = f"""# E05D weak-grid materiality challenge preregistration

**Freeze date:** 2026-08-27. **Git SHA before freeze:** `{git_sha}`.

E05D preserves C1/C2 supported; C3a, C3b-A, and C3b-B rejected; E05B D1/D2/D3 supported;
and original physical cycle/SCC claim C3 unassessed. It creates only
`C3c_WEAK_GRID_MATERIAL_EXTERNALITY`. All 15 E05B-reproducible coalitions, exact A1--A8 endpoints,
eight new Sobol seeds, one topology-only 34-line R/X-preserving grid coordinate, nine tau values,
the actual empty-system margin-setting mode, global Hungarian/BMAC tracking, exact complex Möbius
composition, signed margin metrics, failure statuses, and both C3c routes are frozen in the companion files.

Endpoint calibration evaluates the EMPTY coalition only. Per-seed critical modes are frozen after endpoint
calibration and before any action coalition result. SCR36/37/38 and min SCR are conventional independent
diagnostics and do not determine success. No ParaEMT or E06 run is permitted.

> E05D is the final new physical stress coordinate allowed for TX3.
> If C3c is rejected, no additional load, inertia, controller, line-set,
> voltage-limit, action-amplitude, or mode-selection stress search will be
> conducted to obtain a material connected stability consequence.
"""
    (PREREG / "E05D_PREREGISTRATION.md").write_text(prereg_text, encoding="utf-8")
    (ARTIFACT / "E05D_GRID_STRENGTH_IMPLEMENTATION.md").write_text("""# E05D conventional grid-strength implementation

The diagnostic follows the conventional NERC definition `SCR = S_sc / S_rated`. On system base,
`S_sc,pu = V_POI^2 / |Z_th|`, hence `SCR_i = V_i^2/(|Z_th,i| S_rated,i,pu)`. The passive positive-sequence
Ybus retains lines, transformers, line charging, and fixed shunts; constant-power loads and the three GFL
current sources are opened, while synchronous source buses 30--35 and slack 39 are ideal voltage sources
(incrementally grounded). `Z_th` is the corresponding diagonal of the inverse reduced Ybus. Each GFL rating
is 1000 MVA on the 100-MVA system base. This is a diagnostic, not a claim gate, and conventional SCR can be
optimistic when nearby IBRs interact. Reference: NERC, *Short-Circuit Modeling and System Strength* (2018),
https://www.nerc.com/globalassets/programs/rapa/ra/short_circuit_whitepaper_final_1_26_18.pdf.

Tests in `tests/unit/test_e05d_weak_grid.py` verify the one-source analytical Thevenin formula, inverse
impedance scaling, and the exact 34/12 topology classification.
""", encoding="utf-8")
    hashes = {name: sha256(PREREG / name) for name in PREREG_NAMES}
    write_json(PREREG / "E05D_PREREGISTRATION_SHA256.json", {
        "freeze_id": "TX3-E05D-PREREGISTRATION-HASHES-1.0", "files": list(PREREG_NAMES), "sha256": hashes,
    })
    print(json.dumps({
        "status": "FROZEN", "candidate_count": 15, "pair_count": 10, "triple_count": 5,
        "new_seed": SEED, "seed_count": 8, "included_branch_count": 34, "excluded_branch_count": 12,
        "preregistration_file_count": 9,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
