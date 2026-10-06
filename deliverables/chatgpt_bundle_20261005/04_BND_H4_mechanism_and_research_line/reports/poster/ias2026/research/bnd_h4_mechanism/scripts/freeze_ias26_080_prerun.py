"""Freeze deterministic IAS26-080 holdout IDs and the audited G2 TDS contract."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
REPO = HERE.parents[6]
RESEARCH = REPO / "reports/poster/ias2026/research"
TAU = 1e-8


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_root.resolve()
    selection_path = run / "inputs/IAS26-080_TDS_HOLDOUT_V1.csv"
    config_path = run / "config/IAS26-080_TDS_OPERATING_V1.json"
    if selection_path.exists() or config_path.exists():
        raise RuntimeError("IAS26-080 preregistration outputs already exist; refusing overwrite")
    scenarios = pd.read_csv(run / "derived/SCENARIO_METRICS.csv")
    mc_cases = pd.read_csv(run / "derived/MC_CASES.csv")
    valid = scenarios[scenarios.scenario_valid.astype(bool)].copy()
    picked: set[str] = set()
    rows: list[dict] = []

    def add_stratum(name: str, preferred: pd.DataFrame, sort_columns: list[str], ascending: list[bool], fallback: pd.DataFrame) -> None:
        preferred = preferred[~preferred.scenario_id.isin(picked)].sort_values(sort_columns + ["scenario_id"], ascending=ascending + [True], kind="mergesort")
        chosen = preferred.head(6).copy()
        if len(chosen) < 6:
            pool = fallback[~fallback.scenario_id.isin(picked | set(chosen.scenario_id))].sort_values(sort_columns + ["scenario_id"], ascending=ascending + [True], kind="mergesort")
            chosen = pd.concat([chosen, pool.head(6 - len(chosen))], ignore_index=True)
        for rank, scenario in enumerate(chosen.itertuples(index=False), 1):
            picked.add(scenario.scenario_id)
            actual = "PRIMARY_STRATUM" if scenario.scenario_id in set(preferred.scenario_id) else "PREDECLARED_FALLBACK"
            case_rows = mc_cases[(mc_cases.scenario_id == scenario.scenario_id) & (mc_cases.treatment == "ORIGINAL") & (mc_cases.portfolio_mask != 15)]
            if len(case_rows) != 15:
                raise RuntimeError(f"incomplete proper-subset lattice for {scenario.scenario_id}")
            case_rows = case_rows.sort_values(["alpha_perp", "portfolio_mask"], ascending=[False, True], kind="mergesort")
            critical = case_rows.iloc[0]
            rows.append({
                "tds_stratum": name,
                "rank_within_stratum": rank,
                "scenario_id": scenario.scenario_id,
                "selection_basis": actual,
                "load_scale": scenario.load_scale,
                "epsilon_l2": scenario.epsilon_l2,
                "saturation_count": scenario.saturation_count,
                "h4_minimal_blocker_persists": scenario.h4_minimal_blocker_persists,
                "other_minimal_blocker_ids": scenario.minimal_unstable_blocker_ids,
                "alpha_h4_original": scenario.y_s_h4_alpha_original,
                "alpha_h4_retuned": scenario.alpha_perp_intervention,
                "delta_alpha": scenario.delta_alpha,
                "most_critical_proper_subset_mask": int(critical.portfolio_mask),
                "most_critical_proper_subset_id": str(critical.portfolio_id),
                "most_critical_proper_subset_alpha": float(critical.alpha_perp),
            })

    a = valid[valid.h4_minimal_blocker_persists == "TRUE"]
    add_stratum("A_ROBUST_BLOCKER", a, ["m_s_signed_blocker_margin"], [False], valid)
    b = valid[(valid.h4_minimal_blocker_persists == "TRUE") & (valid.m_s_signed_blocker_margin > 0)]
    add_stratum("B_NEAR_BOUNDARY", b, ["m_s_signed_blocker_margin"], [True], a)
    c = valid[(valid.h4_minimal_blocker_persists == "FALSE") & valid.other_minimal_blocker]
    add_stratum("C_WITNESS_CHANGE", c, ["scenario_id"], [True], valid[valid.h4_minimal_blocker_persists == "FALSE"])
    d_primary = valid[(valid.y_s_h4_alpha_original > TAU) & (valid.alpha_perp_intervention >= -TAU)]
    d_fallback = valid[(valid.y_s_h4_alpha_original > TAU) & (valid.delta_alpha > TAU)]
    add_stratum("D_DIFFICULT_INTERVENTION", d_primary, ["alpha_perp_intervention"], [False], d_fallback if len(d_fallback) else valid)
    selection = pd.DataFrame(rows)
    if len(selection) != 24 or selection.scenario_id.nunique() != 24:
        raise RuntimeError(f"holdout must contain 24 distinct IDs; got {len(selection)} rows / {selection.scenario_id.nunique()} IDs")
    selection.to_csv(selection_path, index=False)

    sources = [
        "reports/poster/ias2026/research/experiments/tx4_tds_final.py",
        "reports/poster/ias2026/research/experiments/G2_tds.py",
        "reports/poster/ias2026/research/experiments/E25_time_domain_validation.py",
        "reports/poster/ias2026/research/experiments/_f7_common.py",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_case.py",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_network.py",
        "reports/poster/ias2026/research/src/ibr_cycles/models/ieee39_devices.py",
    ]
    source_hashes = {rel: sha256(REPO / rel) for rel in sources}
    config = {
        "schema_version": 1,
        "ticket": "IAS26-080",
        "config_id": "IAS26-080_TDS_OPERATING_V1",
        "frozen_after_MC_before_TDS": True,
        "parent_run_id": run.name,
        "model_scope": "frozen IEEE-39 TX4/P4, no-governor, exact GFL11; nonlinear phasor-domain DAE, not EMT",
        "primary_metric": "measured trajectory matrix-pencil alpha_TDS and frequency; eigenvalue used only as comparator and observable-mode locator",
        "decision_threshold_s-1": TAU,
        "canonical_cases": [
            {"case_id": "CANONICAL_PROPER_30_33_35", "members": [30, 33, 35], "g": 0.03625},
            {"case_id": "CANONICAL_H4_ORIGINAL", "members": [30, 33, 35, 37], "g": 0.03625},
            {"case_id": "CANONICAL_H4_RETUNED", "members": [30, 33, 35, 37], "g": 0.25},
        ],
        "canonical_operating_point": "frozen nominal IEEE-39 base network, source scheduled generation and unperturbed loads",
        "canonical_trajectory_count": 18,
        "holdout_trajectory_count": 72,
        "total_trajectory_count": 90,
        "canonical_amplitude_multipliers": [0.5, 1.0, 2.0],
        "historical_disturbances": {
            "D1": {"kind": "rotor-speed kick", "machine_bus": 31, "A0_pu": 0.0001, "duration_s": 0.0},
            "D2": {"kind": "active-load rectangular pulse", "load_bus": 20, "A0_fraction_of_P": 0.02, "duration_s": 0.2},
        },
        "holdout_disturbance": "D2 at A0 (2% active-load pulse, 0.2 s, bus 20), exactly as the frozen historical TX4/G2 common pulse; one trajectory per case per scenario",
        "holdout_strata": {
            "A_ROBUST_BLOCKER": "six valid scenarios with largest positive m_s among H4 persistent blockers",
            "B_NEAR_BOUNDARY": "six remaining valid H4 blockers with smallest positive m_s",
            "C_WITNESS_CHANGE": "six remaining valid scenarios where a different inclusion-minimal unstable portfolio exists in V4; lexicographic scenario-ID tie/order",
            "D_DIFFICULT_INTERVENTION": "six remaining valid original-H4-unstable scenarios where retuned H4 remains at/above -tau, ranked by largest retuned alpha; fallback to material deterioration, then all valid IDs if the class is empty",
            "overlap_policy": "disjoint strata, no repeated IDs; if a preferred class cannot fill six, use only the documented deterministic fallback population, never a new or replacement draw",
            "proper_subset_case": "for each selected ID, evaluate the most critical strict H4 subset: largest original alpha_perp among the 15 proper subsets; portfolio-mask ascending breaks exact ties",
        },
        "integration_and_analysis_contract": {
            "integrator": "SciPy BDF on reduced differential state, algebraic network solved by Newton at every RHS evaluation",
            "rtol": 1e-7,
            "atol": 1e-9,
            "fixed_relative_central_difference_jacobian_step": 1e-6,
            "algebraic_residual_target_inf": 1e-10,
            "disturbance_switch_handling": "separate BDF segments at pulse end; never step across the switch",
            "horizon_s": "clip(3/abs(alpha_eig), 30, 120)",
            "divergence_limit_pu": 0.05,
            "observable": "for each case, difference of the machine-speed pair with maximum separation in the critical right eigenvector, matching historical G2 rule; additionally save the same canonical plot observable omega_sg34-omega_sg39",
            "small_signal_window": "retain prefix while cumulative machine-speed spread <=1e-3 pu; matrix pencil order 14; minimum 200 samples",
            "frequency_diagnostic": "matrix-pencil exponent and 0.05-20 Hz Hann-window FFT on actual sampled trace",
            "outcome": "DIVERGED if speed spread exceeds 0.05 pu or integration/algebraic solve fails; otherwise measured sign from second-half envelope and matrix-pencil; failed fit retained as missing, never replaced",
            "amplitude_linearity": "all amplitudes retained; eligibility is diagnosed from observed spread and fit window, no amplitude is removed post hoc",
            "fit_r2": "report OLS fit R2 for the measured complex-exponential signal on the declared small-signal fit window; diagnostic only, no unregistered threshold",
        },
        "selection_sha256": sha256(selection_path),
        "mc_case_table_sha256": sha256(run / "derived/MC_CASES.csv"),
        "scenario_metrics_sha256": sha256(run / "derived/SCENARIO_METRICS.csv"),
        "historical_and_model_source_sha256": source_hashes,
        "historical_audit": {
            "entry_wrapper": "experiments/tx4_tds_final.py",
            "integrator_contract": "experiments/G2_tds.py",
            "historical_same_model_limitation": "a short common pulse did not visibly excite the RHP mode in the historical blocker trace; this is an explicit falsification risk and will not be overridden by eigenvalue curves",
        },
        "no_scenario_replacement": True,
        "no_retuning_or_controller_change": True,
    }
    config_path.write_text(json.dumps(config, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"selection_path": str(selection_path), "selection_sha256": config["selection_sha256"], "N_selected": 24, "strata": selection.groupby("tds_stratum").scenario_id.nunique().to_dict(), "config_path": str(config_path), "total_trajectories": 90}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
