from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
PREREG = ARTIFACT / "preregistration"
TABLES = ARTIFACT / "tables"
FIGURES = ARTIFACT / "figures"
PAPER = ROOT / "reports" / "papers" / "tx3_connected_intervention_calculus"
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 8.5, "legend.fontsize": 7.2})


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def report(name: str, content: str) -> None:
    (ARTIFACT / name).write_text(content.strip() + "\n", encoding="utf-8")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def safe_ratio(final: float, nominal: float) -> float:
    return float(final / nominal) if nominal != 0.0 and math.isfinite(nominal) else math.nan


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    (ARTIFACT / "manifests").mkdir(parents=True, exist_ok=True)
    candidates = pd.read_parquet(PREREG / "e05d_candidates.parquet")
    limits = pd.read_parquet(ARTIFACT / "baseline_calibration" / "E05D_BASELINE_STRESS_LIMITS.parquet")
    grid = pd.read_parquet(TABLES / "e05d_grid_strength.parquet")
    critical = pd.read_parquet(TABLES / "e05d_critical_modes.parquet")
    empty_tracks = pd.read_parquet(PREREG / "e05d_empty_mode_tracking_freeze.parquet")
    tracks = pd.read_parquet(ARTIFACT / "holdout" / "e05d_mode_tracking.parquet")
    results = pd.read_parquet(ARTIFACT / "holdout" / "e05d_connected_externalities.parquet")
    builds = pd.read_parquet(ARTIFACT / "holdout" / "e05d_failure_ledger.parquet")
    gates = load_json(PREREG / "E05D_C3C_GATE_FREEZE.json")
    compute = load_json(ARTIFACT / "logs" / "e05d_campaign_compute.json")
    valid = results[results["status"] == "SUCCESS"].copy()

    decision_rows: list[dict[str, Any]] = []
    amplification_rows: list[dict[str, Any]] = []
    for candidate in candidates.itertuples(index=False):
        frame = results[results["candidate_id"] == candidate.candidate_id]
        smooth = frame[frame["status"] == "SUCCESS"]
        high = smooth[smooth["tau"].isin(gates["route_2"]["eligible_tau_levels"])]
        per_seed = smooth.groupby("seed_id")
        hard_seed_count = int(per_seed["hard_failure"].any().sum())
        screen_seed_count = int(per_seed["screen_failure"].any().sum())
        combined_seed_count = int(per_seed.apply(lambda group: bool(group["hard_failure"].any() or group["screen_failure"].any()), include_groups=False).sum())
        route1 = combined_seed_count >= int(gates["route_1"]["replicated_failure_seed_minimum"])
        material = high[
            (high["rho_comp"] >= float(gates["route_2"]["rho_comp_minimum"]))
            & (high["Delta_zeta"] <= float(gates["route_2"]["delta_zeta_maximum"]))
        ]
        route2_seed_count = int(material["seed_id"].nunique())
        sign_fraction = float((high["Delta_lambda_real"] > 0.0).mean()) if len(high) else 0.0
        route2 = bool(
            route2_seed_count >= int(gates["route_2"]["replicated_seed_minimum"])
            and sign_fraction >= float(gates["route_2"]["destabilizing_sign_fraction_minimum"])
        )
        decision_rows.append({
            "candidate_id": candidate.candidate_id, "coalition": candidate.coalition_key,
            "coalition_key": candidate.coalition_key, "order": int(candidate.order),
            "valid_seed_paths": int(smooth["seed_id"].nunique()),
            "hard_failure_seed_paths": hard_seed_count, "screen_failure_seed_paths": screen_seed_count,
            "combined_failure_seed_paths": combined_seed_count, "route2_seed_paths": route2_seed_count,
            "destabilizing_high_stress_sign_fraction": sign_fraction,
            "maximum_rho_comp": float(smooth["rho_comp"].max()) if len(smooth) else math.nan,
            "minimum_Delta_zeta": float(smooth["Delta_zeta"].min()) if len(smooth) else math.nan,
            "C3c_route_1": "SUPPORTED" if route1 else "REJECTED",
            "C3c_route_2": "SUPPORTED" if route2 else "REJECTED",
            "C3c": "SUPPORTED" if route1 or route2 else "REJECTED",
        })
        nominal = smooth[smooth["tau"] == 0.0].median(numeric_only=True)
        maximum = smooth[smooth["tau"] == 1.0].median(numeric_only=True)
        amplification_rows.append({
            "candidate_id": candidate.candidate_id, "coalition": candidate.coalition_key, "order": int(candidate.order),
            "nominal_Delta_lambda_real": float(nominal.get("Delta_lambda_real", math.nan)),
            "maximum_stress_Delta_lambda_real": float(maximum.get("Delta_lambda_real", math.nan)),
            "Delta_lambda_real_amplification_ratio": safe_ratio(float(maximum.get("Delta_lambda_real", math.nan)), float(nominal.get("Delta_lambda_real", math.nan))),
            "nominal_Delta_zeta": float(nominal.get("Delta_zeta", math.nan)),
            "maximum_stress_Delta_zeta": float(maximum.get("Delta_zeta", math.nan)),
            "Delta_zeta_amplification_ratio": safe_ratio(float(maximum.get("Delta_zeta", math.nan)), float(nominal.get("Delta_zeta", math.nan))),
            "nominal_rho_comp": float(nominal.get("rho_comp", math.nan)),
            "maximum_stress_rho_comp": float(maximum.get("rho_comp", math.nan)),
            "rho_comp_amplification_ratio": safe_ratio(float(maximum.get("rho_comp", math.nan)), float(nominal.get("rho_comp", math.nan))),
            "nominal_rho_baseline": float(nominal.get("rho_baseline", math.nan)),
            "maximum_stress_rho_baseline": float(maximum.get("rho_baseline", math.nan)),
            "nominal_min_SCR": float(nominal.get("min_SCR", math.nan)),
            "maximum_stress_min_SCR": float(maximum.get("min_SCR", math.nan)),
            "nominal_critical_damping": float(nominal.get("critical_damping_ratio", math.nan)),
            "maximum_stress_critical_damping": float(maximum.get("critical_damping_ratio", math.nan)),
            "nominal_critical_kappa_lambda": float(nominal.get("critical_kappa_lambda", math.nan)),
            "maximum_stress_critical_kappa_lambda": float(maximum.get("critical_kappa_lambda", math.nan)),
        })
    decisions = pd.DataFrame(decision_rows)
    amplification = pd.DataFrame(amplification_rows)
    decisions.to_parquet(TABLES / "e05d_candidate_decisions.parquet", index=False)
    amplification.to_parquet(TABLES / "e05d_secondary_amplification.parquet", index=False)
    c3c = "SUPPORTED" if decisions["C3c"].eq("SUPPORTED").any() else "REJECTED"
    mechanism_authorized = c3c == "SUPPORTED"

    corr_condition = spearmanr(valid["rho_comp"], np.log10(valid["critical_kappa_lambda"]), nan_policy="omit")
    remaining_margin = -valid["lambda_lower_real"]
    corr_margin = spearmanr(valid["rho_comp"], remaining_margin, nan_policy="omit")
    conditioning_growth = float(corr_condition.statistic)
    margin_shrink_growth = -float(corr_margin.statistic)
    if conditioning_growth >= 0.2 and margin_shrink_growth >= 0.2:
        mechanism = "both shrinking decay margin and increasing eigenvalue conditioning"
    elif conditioning_growth >= margin_shrink_growth and conditioning_growth >= 0.2:
        mechanism = "primarily increasing eigenvalue conditioning"
    elif margin_shrink_growth >= 0.2:
        mechanism = "primarily shrinking decay margin"
    else:
        mechanism = "neither factor shows a strong monotone association"
    conditioning = {
        "spearman_rho_comp_vs_log10_critical_kappa_lambda": float(corr_condition.statistic),
        "pvalue_rho_comp_vs_log10_critical_kappa_lambda": float(corr_condition.pvalue),
        "spearman_rho_comp_vs_remaining_lower_order_decay_margin": float(corr_margin.statistic),
        "pvalue_rho_comp_vs_remaining_lower_order_decay_margin": float(corr_margin.pvalue),
        "descriptive_interpretation": mechanism,
        "claim_gate_uses_correlation": False,
    }
    write_json(TABLES / "e05d_conditioning_summary.json", conditioning)

    strongest = valid.loc[valid["rho_comp"].idxmax()].to_dict()
    stabilizing = valid.loc[valid["rho_comp"].idxmin()].to_dict()
    per_seed_maximum = valid.groupby("seed_id")["rho_comp"].max()
    tracking_success_rate = float(tracks["status"].eq("SUCCESS").mean())
    final_payload = {
        "stage": "TX3_E05D_WEAK_GRID",
        "C1": "SUPPORTED_UNCHANGED", "C2": "SUPPORTED_UNCHANGED",
        "C3a_nominal_material_damping_externality": "REJECTED_UNCHANGED",
        "C3b_A_load_stress_amplification": "REJECTED_UNCHANGED",
        "C3b_B_load_stress_material_consequence": "REJECTED_UNCHANGED",
        "E05B_D1": "SUPPORTED_UNCHANGED", "E05B_D2": "SUPPORTED_UNCHANGED", "E05B_D3": "SUPPORTED_UNCHANGED",
        "original_C3_cycle_or_SCC": "UNASSESSED", "C3c_WEAK_GRID_MATERIAL_EXTERNALITY": c3c,
        "frozen_candidate_count": len(candidates), "pair_count": int(candidates["order"].eq(2).sum()),
        "triple_count": int(candidates["order"].eq(3).sum()), "new_seed_count": 8,
        "kappa_end_range": [float(limits["kappa_end"].min()), float(limits["kappa_end"].max())],
        "SCR36_range": [float(grid["SCR36"].min()), float(grid["SCR36"].max())],
        "SCR37_range": [float(grid["SCR37"].min()), float(grid["SCR37"].max())],
        "SCR38_range": [float(grid["SCR38"].min()), float(grid["SCR38"].max())],
        "baseline_critical_damping_range": [float(empty_tracks[empty_tracks["tau"] == 0]["damping_ratio"].min()), float(empty_tracks[empty_tracks["tau"] == 0]["damping_ratio"].max())],
        "maximum_stress_critical_damping_range": [float(critical["damping_ratio"].min()), float(critical["damping_ratio"].max())],
        "mode_tracking_success_rate": tracking_success_rate,
        "largest_destabilizing_Delta_lambda_real_per_s": float(valid["Delta_lambda_real"].max()),
        "maximum_rho_comp": float(valid["rho_comp"].max()),
        "median_seed_maximum_rho_comp": float(per_seed_maximum.median()),
        "rho_comp_ge_0_25_points": int((valid["rho_comp"] >= 0.25).sum()),
        "Delta_zeta_le_minus_0_0025_points": int((valid["Delta_zeta"] <= -0.0025).sum()),
        "hard_composition_failure_points": int(valid["hard_failure"].sum()),
        "screen_composition_failure_points": int(valid["screen_failure"].sum()),
        "conditioning": conditioning,
        "strongest_positive_result": {key: strongest[key] for key in (
            "seed_id", "tau", "kappa_grid", "coalition", "rho_comp", "Delta_lambda_real", "Delta_zeta",
            "lambda_lower_real", "lambda_real", "min_SCR",
        )},
        "strongest_stabilizing_result": {key: stabilizing[key] for key in (
            "seed_id", "tau", "kappa_grid", "coalition", "rho_comp", "Delta_lambda_real", "Delta_zeta", "min_SCR",
        )},
        "c3_mechanism_localization_authorized": mechanism_authorized,
        "e06_consideration_authorized": False, "E06_executed": False,
        "materiality_search_closed": not mechanism_authorized, "E05E_authorized": False,
        "new_paraemt_scientific_runs": 0,
    }
    write_json(ARTIFACT / "E05D_FINAL_DECISION.json", final_payload)

    # F1 — weak-grid stress and conventional SCR.
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for limit in limits.itertuples(index=False):
        local = grid[grid["seed_id"] == limit.seed_id].sort_values("tau")
        axes[0].plot(local["tau"], local["kappa_grid"], "o-", ms=2.5, label=f"{limit.seed_id}: {limit.endpoint_reason}")
        axes[1].plot(local["tau"], local["min_SCR"], "o-", ms=2.5, label=limit.seed_id)
    axes[0].set(xlabel=r"Normalized stress $\tau$", ylabel=r"$\kappa_{grid}$", title="Frozen weak-grid paths")
    axes[1].set(xlabel=r"Normalized stress $\tau$", ylabel="Minimum conventional SCR", title="Independent grid-strength diagnostic")
    for axis in axes: axis.grid(alpha=0.22); axis.legend(ncol=2)
    save(fig, "figure_E05D_1_weak_grid_stress")

    # F2 — empty critical margin.
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for seed, frame in empty_tracks.groupby("seed_id"):
        frame = frame.sort_values("tau")
        axes[0].plot(frame["tau"], frame["damping_ratio"], "o-", ms=2.5, label=seed)
        axes[1].plot(frame["tau"], frame["lambda_real"], "o-", ms=2.5, label=seed)
    axes[0].axhline(0.05, color="black", ls="--", label="5% screen")
    axes[0].set(xlabel=r"$\tau$", ylabel="Critical damping ratio", title="Actual empty-coalition margin-setting mode")
    axes[1].axhline(0.0, color="black", ls="--")
    axes[1].set(xlabel=r"$\tau$", ylabel=r"Critical Re $\lambda$ (s$^{-1}$)", title="Critical decay-rate continuation")
    for axis in axes: axis.grid(alpha=0.22); axis.legend(ncol=2)
    save(fig, "figure_E05D_2_empty_critical_margin")

    # F3 — all candidate rho heatmap, invalid outcomes retained as gray missing cells.
    heat = valid.pivot_table(index="coalition", columns="tau", values="rho_comp", aggfunc="median").reindex(candidates["coalition_key"])
    fig, axis = plt.subplots(figsize=(9.6, 6.0))
    bound = max(abs(float(np.nanmin(heat.to_numpy()))), abs(float(np.nanmax(heat.to_numpy()))), 1e-12)
    cmap = plt.get_cmap("coolwarm").copy(); cmap.set_bad("#888888")
    image = axis.imshow(heat, aspect="auto", cmap=cmap, vmin=-bound, vmax=bound)
    axis.set_yticks(range(len(heat.index)), heat.index); axis.set_xticks(range(len(heat.columns)), [f"{value:g}" for value in heat.columns])
    axis.set(xlabel=r"$\tau$", title=r"All candidates: median $\rho_{comp}$ (invalid = gray)")
    fig.colorbar(image, ax=axis, label=r"median $\rho_{comp}$")
    save(fig, "figure_E05D_3_all_candidate_rho_heatmap")

    # F4/F5 use the objectively closest composition case: global maximum rho_comp.
    focus = valid[valid["candidate_id"] == strongest["candidate_id"]]
    summary = focus.groupby("tau").median(numeric_only=True)
    fig, axis = plt.subplots(figsize=(6.5, 4.6))
    axis.plot(summary["lambda_lower_real"], summary["lambda_lower_imag"], "o-", label=r"$\lambda_{<q}$")
    axis.plot(summary["lambda_real"], summary["lambda_imag"], "s--", label=r"$\lambda_{full}$")
    axis.axvline(0, color="black", lw=0.8); axis.grid(alpha=0.22); axis.legend()
    axis.set(xlabel=r"Re $\lambda$", ylabel=r"Im $\lambda$", title=f"Closest composition approach: {strongest['coalition']}")
    save(fig, "figure_E05D_4_lower_order_vs_full_pole")
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    axis.plot(summary.index, summary["zeta_lower"], "o-", label=r"$\zeta_{<q}$")
    axis.plot(summary.index, summary["damping_ratio"], "s--", label=r"$\zeta_{full}$")
    axis.axhline(0.05, color="black", ls=":", label="5% screen")
    axis.set(xlabel=r"$\tau$", ylabel="Damping ratio", title=f"Damping-screen audit: {strongest['coalition']}")
    axis.grid(alpha=0.22); axis.legend()
    save(fig, "figure_E05D_5_damping_screen")

    fig, axis = plt.subplots(figsize=(6.5, 4.4))
    for order, frame in results.groupby("order"):
        good = frame[frame["status"] == "SUCCESS"]
        axis.scatter(good["min_SCR"], good["rho_comp"], s=13, alpha=0.48, label=f"order {order} SUCCESS")
    failed = results[results["status"] != "SUCCESS"]
    if len(failed): axis.scatter(failed["min_SCR"], np.zeros(len(failed)), marker="x", color="black", label="invalid retained")
    axis.axhline(0.25, color="black", ls="--", label="material gate")
    axis.set(xlabel="Minimum conventional SCR", ylabel=r"$\rho_{comp}$", title="Grid strength versus connected relevance")
    axis.grid(alpha=0.22); axis.legend()
    save(fig, "figure_E05D_6_grid_strength_vs_connected_relevance")

    maxima = valid.groupby(["candidate_id", "coalition", "order", "seed_id"])["rho_comp"].max().reset_index()
    fig, axis = plt.subplots(figsize=(6.2, 4.2))
    axis.boxplot([maxima[maxima["order"] == 2]["rho_comp"], maxima[maxima["order"] == 3]["rho_comp"]], tick_labels=["Pairs", "Triples"], showfliers=True)
    axis.axhline(0.25, color="black", ls="--", label="material gate")
    axis.set(ylabel=r"Maximum frozen-grid $\rho_{comp}$ per seed", title="Pair versus triple interactions")
    axis.grid(axis="y", alpha=0.22); axis.legend()
    save(fig, "figure_E05D_7_pair_vs_triple")

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
    axes[0].scatter(np.log10(valid["critical_kappa_lambda"]), valid["rho_comp"], s=12, alpha=0.45)
    axes[0].set(xlabel=r"$\log_{10}\kappa_\lambda$", ylabel=r"$\rho_{comp}$", title=f"Spearman r={corr_condition.statistic:.3f}")
    axes[1].scatter(remaining_margin, valid["rho_comp"], s=12, alpha=0.45, color="#D55E00")
    axes[1].set(xlabel="Remaining lower-order decay margin", ylabel=r"$\rho_{comp}$", title=f"Spearman r={corr_margin.statistic:.3f}")
    for axis in axes: axis.grid(alpha=0.22)
    save(fig, "figure_E05D_8_conditioning")

    branch_freeze = load_json(PREREG / "E05D_GRID_STRESS_BRANCH_FREEZE.json")
    holdout_freeze = load_json(PREREG / "E05D_HOLDOUT_FREEZE.json")
    candidate_lines = "\n".join(f"- {row.candidate_id}: {row.coalition_key} (order {row.order})" for row in candidates.itertuples(index=False))
    decision_lines = "\n".join(
        f"| {row.coalition} | {row.valid_seed_paths}/8 | {row.hard_failure_seed_paths} | {row.screen_failure_seed_paths} | {row.route2_seed_paths} | {row.maximum_rho_comp:.6g} | {row.C3c} |"
        for row in decisions.itertuples(index=False)
    )
    status_lines = "\n".join(f"- {key}: {value}" for key, value in builds["status"].value_counts(dropna=False).to_dict().items())
    report("E05D_PREREGISTRATION.md", (PREREG / "E05D_PREREGISTRATION.md").read_text(encoding="utf-8"))
    report("E05D_CANDIDATE_FREEZE.md", f"# E05D candidate freeze\n\nAll 15 E05B-reproducible coalitions were frozen without E05D outcomes:\n\n{candidate_lines}")
    report("E05D_HOLDOUT_FREEZE.md", f"# E05D holdout freeze\n\nEight new scrambled Sobol points use seed `{holdout_freeze['seed']}`; no historical operating point is reused. Canonical case SHA-256: `{holdout_freeze['canonical_case_sha256']}`.")
    included = [record for record in branch_freeze["branches"] if record["included"]]
    excluded = [record for record in branch_freeze["branches"] if not record["included"]]
    report("E05D_GRID_STRESS_BRANCH_FREEZE.md", f"# E05D grid-stress branch freeze\n\nThe topology-only rule includes {len(included)} AC lines (`trans=0`) and excludes {len(excluded)} fixed-tap/interface transformer records. Exact IDs: {', '.join(row['line_id'] for row in included)}. The full 46-row decision table is in the JSON freeze.")
    report("E05D_GRID_STRENGTH_IMPLEMENTATION.md", """
# E05D conventional grid-strength implementation

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
""")
    report("E05D_STRESS_ENDPOINTS.md", f"# E05D stress endpoints\n\nEMPTY-only continuation froze kappa_end={limits['kappa_end'].min():.6g}--{limits['kappa_end'].max():.6g}. Boundary reasons: {limits['endpoint_reason'].value_counts().to_dict()}. No coalition was evaluated during calibration.")
    report("E05D_CRITICAL_MODE_RULE.md", f"# E05D critical-mode rule\n\nAt each empty tau=1 endpoint, the positive-imaginary 0.1--30 Hz mode with minimum damping was selected. Its left/right eigenvectors, condition number, residual and dominant state participation were frozen before coalitions. Endpoint damping spans {critical['damping_ratio'].min():.6g}--{critical['damping_ratio'].max():.6g}.")
    report("E05D_MODE_TRACKING_AUDIT.md", f"# E05D mode-tracking audit\n\nBackward empty tracking and action homotopy used global one-to-one Hungarian BMAC with threshold 0.90 and half/quarter or 9/17-node refinement only. Action-subset success rate: {tracking_success_rate:.6%}; minimum successful BMAC: {tracks.loc[tracks.status.eq('SUCCESS'),'BMAC'].min():.6g}. No manual reassignment occurred.")
    report("E05D_CONNECTED_POLE_RESULTS.md", f"# E05D connected-pole results\n\nLargest destabilizing real connected shift: {valid['Delta_lambda_real'].max():.6g} s^-1. Maximum signed rho_comp: {valid['rho_comp'].max():.6g}. All pair/triple complex Möbius identities use the same frozen seed-specific margin-setting mode.")
    report("E05D_MARGIN_CONSUMPTION.md", f"# E05D margin consumption\n\nMaximum rho_comp={final_payload['maximum_rho_comp']:.6g}; median of the eight seed-wise maxima={final_payload['median_seed_maximum_rho_comp']:.6g}; rho_comp>=0.25 points={final_payload['rho_comp_ge_0_25_points']}; destabilizing material damping points={final_payload['Delta_zeta_le_minus_0_0025_points']}.")
    report("E05D_COMPOSITION_FAILURES.md", f"# E05D composition failures\n\nHard failures={final_payload['hard_composition_failure_points']}; 5%-screen failures={final_payload['screen_composition_failure_points']}. Maximum complex closure error is {max(valid['composition_closure_real'].abs().max(), valid['composition_closure_imag'].abs().max()):.3e}.")
    report("E05D_CONDITIONING_ANALYSIS.md", f"# E05D conditioning analysis\n\nSpearman rho_comp vs log10(kappa_lambda)={corr_condition.statistic:.4f} (p={corr_condition.pvalue:.4g}); rho_comp vs remaining lower-order margin={corr_margin.statistic:.4f} (p={corr_margin.pvalue:.4g}). Descriptively: {mechanism}. This correlation does not enter C3c.")
    report("E05D_FAILURE_LEDGER.md", f"# E05D failure ledger\n\nAll physical/eigen evaluations remain in `e05d_failure_ledger.parquet`. Counts:\n\n{status_lines}\n\nScientific rejection, if present, is not converted into missingness.")
    report("E05D_C3C_DECISION.md", f"""
# E05D C3c decision

| Coalition | Valid seeds | Hard seeds | Screen seeds | Route-2 seeds | Maximum rho | C3c |
|---|---:|---:|---:|---:|---:|---|
{decision_lines}

Overall `C3c_WEAK_GRID_MATERIAL_EXTERNALITY = {c3c}`. Mechanism localization authorized:
`{str(mechanism_authorized).lower()}`. E06 remains unauthorized and was not run. If rejected, the TX3
materiality search is closed and there is no E05E.
""")
    report("E05D_REVIEWER1_AUDIT.md", f"# E05D Reviewer 1 audit\n\nThe candidate population, action amplitudes, Sobol seed, topology-only 34-line coordinate, endpoints, critical-mode rule, tracking thresholds, and C3c gates were frozen before coalition results. C3c={c3c}. The conventional SCR diagnostic verifies weakening but is not used as a success criterion.")
    report("E05D_REVIEWER2_AUDIT.md", f"# E05D Reviewer 2 falsification audit\n\nNo alternate line set, stress coordinate, voltage limit, action amplitude, or replacement mode was searched. All statuses are retained. Historical C1/C2 and D1--D3 stay supported; C3a and C3b-A/B stay rejected; original C3 remains unassessed. E06 and ParaEMT were not executed.")
    report("E05D_FINAL_REVIEW.md", f"# E05D final review\n\nC3c={c3c}. The grid-strength coordinate changes min SCR from {grid[grid.tau.eq(0)].min_SCR.min():.6g}--{grid[grid.tau.eq(0)].min_SCR.max():.6g} nominally to {grid[grid.tau.eq(1)].min_SCR.min():.6g}--{grid[grid.tau.eq(1)].min_SCR.max():.6g} at maximum stress. Hard/screen failures={final_payload['hard_composition_failure_points']}/{final_payload['screen_composition_failure_points']}; materiality search closed={str(not mechanism_authorized).lower()}.")

    source_claims = PAPER / "CLAIMS_EVIDENCE_MATRIX_E05C.md"
    claims_text = source_claims.read_text(encoding="utf-8")
    claims_text += f"""

## E05D versioned addendum (historical matrices preserved byte-for-byte)

| Claim | Exact statement | Gate | Result | Pass/fail | Permitted wording |
|---|---|---|---|---|---|
| C3c | The frozen connected pole externality is materially relevant to the actual margin-setting mode under one preregistered weak-grid coordinate. | Same coalition: >=2/8 replicated hard/screen failures, or >=4/8 high-stress seeds with rho_comp>=0.25, Delta zeta<=-0.0025, and destabilizing sign consistency>=0.75. | Maximum rho_comp {final_payload['maximum_rho_comp']:.6g}; hard/screen failures {final_payload['hard_composition_failure_points']}/{final_payload['screen_composition_failure_points']}; material damping points {final_payload['Delta_zeta_le_minus_0_0025_points']}. | {'PASS' if c3c == 'SUPPORTED' else 'FAIL'} | {'A material weak-grid connected interaction was replicated; mechanism localization is the only authorized next question.' if c3c == 'SUPPORTED' else 'No material connected stability consequence was established; close further TX3 stress searches and do not authorize mechanism localization or E06.'} |
"""
    (PAPER / "CLAIMS_EVIDENCE_MATRIX_E05D.md").write_text(claims_text, encoding="utf-8")
    write_json(ARTIFACT / "manifests" / "E05D_RUN_MANIFEST.json", {
        "stage": "TX3_E05D_WEAK_GRID", "candidate_count": 15, "seed_count": 8, "tau_count": 9,
        "physical_builds": compute["physical_build_count"], "subset_tracks": len(tracks),
        "connected_cases": len(results), "C3c": c3c,
        "c3_mechanism_localization_authorized": mechanism_authorized,
        "e06_consideration_authorized": False, "E06_executed": False, "new_paraemt_scientific_runs": 0,
    })
    print(json.dumps(final_payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
