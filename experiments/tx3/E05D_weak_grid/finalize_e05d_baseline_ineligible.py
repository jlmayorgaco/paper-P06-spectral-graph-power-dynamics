from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E05D)]

from calibrate_e05d_stress import PREREG, evaluate, load_json  # noqa: E402


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
TABLES = ARTIFACT / "tables"
HOLDOUT = ARTIFACT / "holdout"
FIGURES = ARTIFACT / "figures"
PAPER = ROOT / "reports" / "papers" / "tx3_connected_intervention_calculus"
COMMON_COLUMNS = [
    "seed_id", "tau", "kappa_grid", "min_SCR", "SCR36", "SCR37", "SCR38",
    "coalition", "order", "subset", "lambda_real", "lambda_imag", "frequency_hz",
    "damping_ratio", "BMAC", "eigenpair_residual", "kappa_lambda",
    "Delta_lambda_real", "Delta_lambda_imag", "Delta_zeta", "lambda_lower_real",
    "lambda_lower_imag", "zeta_lower", "rho_comp", "rho_baseline", "hard_failure",
    "screen_failure", "status",
]
REPORT_NAMES = (
    "E05D_PREREGISTRATION.md", "E05D_CANDIDATE_FREEZE.md", "E05D_HOLDOUT_FREEZE.md",
    "E05D_GRID_STRESS_BRANCH_FREEZE.md", "E05D_GRID_STRENGTH_IMPLEMENTATION.md",
    "E05D_STRESS_ENDPOINTS.md", "E05D_CRITICAL_MODE_RULE.md", "E05D_MODE_TRACKING_AUDIT.md",
    "E05D_CONNECTED_POLE_RESULTS.md", "E05D_MARGIN_CONSUMPTION.md", "E05D_COMPOSITION_FAILURES.md",
    "E05D_CONDITIONING_ANALYSIS.md", "E05D_FAILURE_LEDGER.md", "E05D_C3C_DECISION.md",
    "E05D_REVIEWER1_AUDIT.md", "E05D_REVIEWER2_AUDIT.md", "E05D_FINAL_REVIEW.md",
)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def report(name: str, body: str) -> None:
    (ARTIFACT / name).write_text(body.strip() + "\n", encoding="utf-8")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def unavailable_figure(stem: str, title: str, reason: str) -> None:
    fig, axis = plt.subplots(figsize=(7.2, 4.2))
    axis.axis("off")
    axis.text(0.5, 0.62, title, ha="center", va="center", fontsize=13, weight="bold")
    axis.text(0.5, 0.40, reason, ha="center", va="center", fontsize=10, wrap=True)
    save(fig, stem)


def main() -> int:
    for path in (TABLES, HOLDOUT, FIGURES, ARTIFACT / "logs", ARTIFACT / "manifests"):
        path.mkdir(parents=True, exist_ok=True)
    seeds = load_json(PREREG / "E05D_HOLDOUT_FREEZE.json")["operating_points"]
    candidates = load_json(PREREG / "E05D_CANDIDATE_FREEZE.json")["candidates"]
    rows = [evaluate(seed, 1.0) for seed in seeds]
    if not all(row["status"] == "SUCCESS" and row["boundary_reason"] == "DAMPING_SCREEN" for row in rows):
        raise RuntimeError("terminal finalizer requires 8/8 successful physical baselines already below the 5% screen")
    eligibility = pd.DataFrame(rows).sort_values("seed_id")
    eligibility["calibration_status"] = "BASELINE_INELIGIBLE_DAMPING_BELOW_0.05"
    eligibility["kappa_limit"] = 1.0
    eligibility["kappa_end"] = math.nan
    eligibility.to_parquet(ARTIFACT / "baseline_calibration" / "E05D_BASELINE_ELIGIBILITY.parquet", index=False)
    eligibility.to_parquet(ARTIFACT / "baseline_calibration" / "E05D_BASELINE_CALIBRATION_TRACE.parquet", index=False)
    limits = eligibility[[
        "seed_id", "status", "calibration_status", "boundary_reason", "kappa_limit", "kappa_end",
        "minimum_oscillatory_damping_ratio", "critical_lambda_real_per_s", "critical_lambda_imag_per_s",
        "SCR36", "SCR37", "SCR38", "min_SCR", "voltage_min_pu", "voltage_max_pu",
        "maximum_eigenpair_residual",
    ]].copy()
    limits.to_parquet(TABLES / "e05d_seed_stress_limits.parquet", index=False)
    grid = eligibility[["seed_id", "kappa_grid", "SCR36", "SCR37", "SCR38", "min_SCR", "status"]].copy()
    grid.insert(1, "tau", 0.0)
    grid.to_parquet(TABLES / "e05d_grid_strength.parquet", index=False)
    empty = pd.DataFrame(columns=COMMON_COLUMNS)
    empty.to_parquet(TABLES / "e05d_critical_modes.parquet", index=False)
    empty.to_parquet(HOLDOUT / "e05d_mode_tracking.parquet", index=False)
    empty.to_parquet(HOLDOUT / "e05d_subset_poles.parquet", index=False)
    empty.to_parquet(HOLDOUT / "e05d_connected_externalities.parquet", index=False)
    empty.to_parquet(HOLDOUT / "e05d_margin_consumption.parquet", index=False)
    empty.to_parquet(HOLDOUT / "e05d_composition_failures.parquet", index=False)
    empty.to_parquet(HOLDOUT / "e05d_conditioning.parquet", index=False)
    ledger = eligibility.copy()
    ledger["tau"] = 0.0; ledger["coalition"] = "EMPTY"; ledger["order"] = 0; ledger["subset"] = "EMPTY"
    ledger.to_parquet(HOLDOUT / "e05d_failure_ledger.parquet", index=False)

    damping_range = [float(eligibility["minimum_oscillatory_damping_ratio"].min()), float(eligibility["minimum_oscillatory_damping_ratio"].max())]
    payload = {
        "stage": "TX3_E05D_WEAK_GRID", "terminal_status": "BASELINE_INELIGIBLE",
        "C1": "SUPPORTED_UNCHANGED", "C2": "SUPPORTED_UNCHANGED",
        "C3a_nominal_material_damping_externality": "REJECTED_UNCHANGED",
        "C3b_A_load_stress_amplification": "REJECTED_UNCHANGED",
        "C3b_B_load_stress_material_consequence": "REJECTED_UNCHANGED",
        "E05B_D1": "SUPPORTED_UNCHANGED", "E05B_D2": "SUPPORTED_UNCHANGED", "E05B_D3": "SUPPORTED_UNCHANGED",
        "original_C3_cycle_or_SCC": "UNASSESSED",
        "C3c_WEAK_GRID_MATERIAL_EXTERNALITY": "UNRESOLVED",
        "frozen_candidate_count": 15, "pair_count": 10, "triple_count": 5, "new_seed_count": 8,
        "baseline_valid_physical_eigen_count": 8, "baseline_below_screen_count": 8,
        "kappa_end_range": None,
        "SCR36_nominal_range": [float(eligibility["SCR36"].min()), float(eligibility["SCR36"].max())],
        "SCR37_nominal_range": [float(eligibility["SCR37"].min()), float(eligibility["SCR37"].max())],
        "SCR38_nominal_range": [float(eligibility["SCR38"].min()), float(eligibility["SCR38"].max())],
        "baseline_critical_damping_range": damping_range,
        "maximum_stress_critical_damping_range": None, "mode_tracking_success_rate": None,
        "largest_destabilizing_Delta_lambda_real_per_s": None, "maximum_rho_comp": None,
        "median_seed_maximum_rho_comp": None, "rho_comp_ge_0_25_points": 0,
        "Delta_zeta_le_minus_0_0025_points": 0, "hard_composition_failure_points": 0,
        "screen_composition_failure_points": 0, "conditioning": None,
        "strongest_positive_result": None,
        "strongest_negative_result": "8/8 frozen baselines are already below the immutable 5% damping screen at kappa_grid=1",
        "coalition_results_computed": False, "coalition_result_rows": 0,
        "c3_mechanism_localization_authorized": False, "e06_consideration_authorized": False,
        "E06_executed": False, "materiality_search_closed": True, "E05E_authorized": False,
        "new_paraemt_scientific_runs": 0,
        "stop_reason": "No admissible kappa_grid>=1 safe endpoint exists without changing the frozen 5% screen or critical-mode rule",
    }
    write_json(ARTIFACT / "E05D_FINAL_DECISION.json", payload)
    write_json(ARTIFACT / "logs" / "E05D_BASELINE_INELIGIBILITY.json", {
        "attempted_calibration_command": "python experiments/tx3/E05D_weak_grid/calibrate_e05d_stress.py --workers 8",
        "attempt_result": "8/8 INVALID_AT_KAPPA_1 because DAMPING_SCREEN",
        "physical_eigen_status": "8/8 SUCCESS", "coalition_evaluations": 0,
        "resolution": "stop before coalition outcomes; do not amend threshold, mode rule, stress coordinate, or endpoint",
    })

    # F1/F2 show the only admissible baseline evidence. F3--F8 explicitly retain the stopped state.
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    axes[0].scatter(eligibility["seed_id"], eligibility["min_SCR"], color="#0072B2")
    axes[0].set(ylabel="Nominal minimum SCR", title=r"Only evaluated point: $\kappa_{grid}=1$")
    axes[1].scatter(eligibility["seed_id"], [1.0] * len(eligibility), color="#D55E00")
    axes[1].set(ylim=(0.95, 1.05), ylabel=r"$\kappa_{grid}$", title="No safe weak-grid endpoint exists")
    for axis in axes: axis.tick_params(axis="x", rotation=45); axis.grid(alpha=0.22)
    save(fig, "figure_E05D_1_weak_grid_stress")
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    axes[0].scatter(eligibility["seed_id"], eligibility["minimum_oscillatory_damping_ratio"], color="#D55E00")
    axes[0].axhline(0.05, color="black", ls="--", label="immutable 5% screen")
    axes[0].set(ylabel="Critical damping ratio", title="8/8 baselines below screen"); axes[0].legend()
    axes[1].scatter(eligibility["seed_id"], eligibility["critical_lambda_real_per_s"], color="#0072B2")
    axes[1].set(ylabel=r"Critical Re $\lambda$ (s$^{-1}$)", title="Physically stable but screen-ineligible")
    for axis in axes: axis.tick_params(axis="x", rotation=45); axis.grid(alpha=0.22)
    save(fig, "figure_E05D_2_empty_critical_margin")
    reason = "Not estimable: endpoint calibration stopped at κ=1 before every coalition outcome."
    unavailable_figure("figure_E05D_3_all_candidate_rho_heatmap", "All-candidate rho heatmap", reason)
    unavailable_figure("figure_E05D_4_lower_order_vs_full_pole", "Lower-order versus full pole", reason)
    unavailable_figure("figure_E05D_5_damping_screen", "Composition damping screen", reason)
    unavailable_figure("figure_E05D_6_grid_strength_vs_connected_relevance", "Grid strength versus connected relevance", reason)
    unavailable_figure("figure_E05D_7_pair_vs_triple", "Pair versus triple", reason)
    unavailable_figure("figure_E05D_8_conditioning", "Conditioning analysis", reason)

    candidate_list = ", ".join(record["coalition_key"] for record in candidates)
    common = "The EMPTY-coalition baseline is physically/eigen numerically valid for 8/8 seeds, but its actual minimum-damping 0.1--30 Hz oscillatory mode has ζ=%.6g--%.6g, already below the immutable 5%% screen at κ=1. Therefore no admissible κ_end>=1 exists and zero coalition outcomes were inspected." % tuple(damping_range)
    report("E05D_PREREGISTRATION.md", (PREREG / "E05D_PREREGISTRATION.md").read_text(encoding="utf-8"))
    report("E05D_CANDIDATE_FREEZE.md", f"# E05D candidate freeze\n\nFrozen population (unused because calibration stopped): {candidate_list}.")
    report("E05D_HOLDOUT_FREEZE.md", "# E05D holdout freeze\n\nEight new scrambled Sobol points D00--D07 use seed 20260830 and the established ranges. All eight entered the empty-only eligibility check.")
    branch = load_json(PREREG / "E05D_GRID_STRESS_BRANCH_FREEZE.json")
    report("E05D_GRID_STRESS_BRANCH_FREEZE.md", f"# E05D grid-stress branch freeze\n\nThe topology-only set remains Line_1--Line_34 ({branch['included_count']} included, {branch['excluded_count']} excluded). It was never outcome-selected and was not changed after the terminal baseline result.")
    # The implementation report was already created before calibration and is intentionally not overwritten.
    report("E05D_STRESS_ENDPOINTS.md", f"# E05D stress endpoints\n\n{common}\n\n`kappa_end` is undefined, not 1.0: calling κ=1 a maximum-stress endpoint would falsely imply a weakened-grid path.")
    report("E05D_CRITICAL_MODE_RULE.md", f"# E05D critical-mode rule\n\nThe frozen rule exposed the eligibility failure rather than selecting a replacement mode. Baseline critical damping range: {damping_range[0]:.8f}--{damping_range[1]:.8f}. Endpoint critical modes were not frozen because no endpoint exists.")
    report("E05D_MODE_TRACKING_AUDIT.md", "# E05D mode-tracking audit\n\nNot executed. Tracking requires a frozen safe endpoint. Success rate is N/A, and no mode was manually replaced.")
    report("E05D_CONNECTED_POLE_RESULTS.md", "# E05D connected-pole results\n\nNot executed. Re Delta lambda and all connected externality metrics are N/A; zero rows is a protocol stop, not selective deletion.")
    report("E05D_MARGIN_CONSUMPTION.md", "# E05D margin consumption\n\nNot estimable. rho_comp and rho_baseline require coalition vertices on the frozen stress path, which does not exist.")
    report("E05D_COMPOSITION_FAILURES.md", "# E05D composition failures\n\nNot evaluated. Counts are zero observed from zero eligible coalition rows, not evidence of absence.")
    report("E05D_CONDITIONING_ANALYSIS.md", "# E05D conditioning analysis\n\nNot estimable without a stress continuation and coalition outcomes. No descriptive correlation is reported.")
    report("E05D_FAILURE_LEDGER.md", f"# E05D failure ledger\n\n{common}\n\nPhysical/eigen statuses are 8 SUCCESS; protocol eligibility is 8 BASELINE_INELIGIBLE_DAMPING_BELOW_0.05. These are retained in the parquet ledger.")
    report("E05D_C3C_DECISION.md", "# E05D C3c decision\n\n`C3c_WEAK_GRID_MATERIAL_EXTERNALITY = UNRESOLVED`. Neither gate can be evaluated because the preregistered safe stress path is empty. Mechanism localization and E06 are not authorized; the final TX3 materiality search is closed.")
    report("E05D_REVIEWER1_AUDIT.md", "# E05D Reviewer 1 audit\n\nThe stop is required by the immutable 5% baseline boundary and actual margin-setting-mode rule. Continuing would require lowering the screen, changing the mode, allowing κ<1, or pretending κ=1 is a weak-grid endpoint.")
    report("E05D_REVIEWER2_AUDIT.md", "# E05D Reviewer 2 falsification audit\n\nNo coalition, alternative line set, alternate stress coordinate, threshold relaxation, or ParaEMT run followed the baseline finding. C3c is unresolved—not rejected—because materiality was never tested under a weakened-grid continuation.")
    report("E05D_FINAL_REVIEW.md", f"# E05D final review\n\n{common}\n\nC1/C2 and E05B D1--D3 remain supported; C3a and C3b-A/B remain rejected; original C3 remains unassessed. C3c is unresolved. Materiality search closed=true; mechanism localization=false; E06=false; E05E=false.")

    claims = (PAPER / "CLAIMS_EVIDENCE_MATRIX_E05C.md").read_text(encoding="utf-8")
    claims += """

## E05D versioned addendum (historical matrices preserved byte-for-byte)

| Claim | Exact statement | Gate | Result | Pass/fail | Permitted wording |
|---|---|---|---|---|---|
| C3c | A frozen connected pole externality is materially relevant to the actual margin-setting mode under one preregistered weak-grid coordinate. | Same coalition: >=2/8 replicated hard/screen failures, or >=4/8 high-stress seeds with rho_comp>=0.25, Delta zeta<=-0.0025, and destabilizing sign consistency>=0.75. | All 8 κ=1 empty baselines were already below the immutable 5% screen; no safe κ_end and no coalition result exist. | UNRESOLVED | Report the baseline-ineligibility stop. Do not claim weak-grid materiality or its absence; close further TX3 stress searches and do not authorize mechanism localization or E06. |
"""
    (PAPER / "CLAIMS_EVIDENCE_MATRIX_E05D.md").write_text(claims, encoding="utf-8")
    write_json(ARTIFACT / "manifests" / "E05D_RUN_MANIFEST.json", {
        "stage": "TX3_E05D_WEAK_GRID", "terminal_status": "BASELINE_INELIGIBLE",
        "candidate_count_frozen": 15, "seed_count_checked": 8, "empty_builds": 8,
        "coalition_builds": 0, "connected_cases": 0, "C3c": "UNRESOLVED",
        "e06_consideration_authorized": False, "E06_executed": False,
        "new_paraemt_scientific_runs": 0,
    })
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
