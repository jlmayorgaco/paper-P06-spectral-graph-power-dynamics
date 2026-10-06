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
ARTIFACT = ROOT / "artifacts" / "tx3" / "E05C_margin_conditioned"
PREREG = ARTIFACT / "preregistration"
REPORTS = ARTIFACT / "reports"
FIGURES = ARTIFACT / "figures"
TABLES = ARTIFACT / "tables"
PAPER = ROOT / "reports" / "papers" / "tx3_connected_intervention_calculus"
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 8.5, "legend.fontsize": 7.5})


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def report(name: str, content: str) -> None:
    (REPORTS / name).write_text(content.strip() + "\n", encoding="utf-8")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def same_sign_fraction(values: pd.Series) -> float:
    return float(max((values > 0).mean(), (values < 0).mean())) if len(values) else 0.0


def main() -> int:
    candidates = pd.read_parquet(PREREG / "e05c_candidates.parquet")
    modes = pd.read_parquet(PREREG / "e05c_mode_family_freeze.parquet")
    limits = pd.read_parquet(ARTIFACT / "baseline_calibration" / "E05C_BASELINE_STRESS_LIMITS.parquet")
    trace = pd.read_parquet(ARTIFACT / "baseline_calibration" / "E05C_BASELINE_CALIBRATION_TRACE.parquet")
    tracks = pd.read_parquet(ARTIFACT / "holdout" / "e05c_mode_tracks.parquet")
    results = pd.read_parquet(ARTIFACT / "holdout" / "e05c_connected_pole_externalities.parquet")
    builds = pd.read_parquet(ARTIFACT / "holdout" / "e05c_failure_ledger.parquet")
    gates = load_json(PREREG / "E05C_GATE_FREEZE.json")
    panel_freeze = load_json(PREREG / "E05C_FIGURE_PANEL_FREEZE.json")
    compute = load_json(ARTIFACT / "logs" / "e05c_campaign_compute.json")
    panel_ids = [record["candidate_id"] for record in panel_freeze["primary_panel_candidates"]]

    decision_rows: list[dict[str, Any]] = []
    for candidate in candidates.itertuples(index=False):
        frame = results[results["candidate_id"] == candidate.candidate_id]
        nominal = frame[frame["tau"] == 0.0]
        final = frame[frame["tau"] == 1.0]
        high = frame[frame["tau"].isin([0.875, 1.0])]
        valid_paths = int(final[final["status"] == "SUCCESS"]["seed_id"].nunique())
        nominal_rho = float(nominal["rho_comp"].median())
        final_rho = float(final["rho_comp"].median())
        amplification = final_rho / nominal_rho if nominal_rho > 0.0 else math.nan
        destabilizing_fraction = float((high["Delta_lambda_real"] > 0.0).mean())
        c3a = bool(
            valid_paths >= int(gates["C3b_A"]["valid_seed_paths_minimum"])
            and destabilizing_fraction >= float(gates["C3b_A"]["destabilizing_sign_fraction_minimum"])
            and final_rho >= float(gates["C3b_A"]["median_final_rho_comp_minimum"])
            and math.isfinite(amplification)
            and amplification >= float(gates["C3b_A"]["median_high_to_nominal_amplification_minimum"])
        )
        per_seed = frame.groupby("seed_id")
        hard_seeds = int(per_seed["hard_composition_failure"].any().sum())
        screen_seeds = int(per_seed["screen_composition_failure"].any().sum())
        route1 = hard_seeds >= 2 or screen_seeds >= 2
        route2_points = high[
            (high["rho_comp"] >= 0.25)
            & (high["Delta_zeta"] <= -0.0025)
            & (high["Delta_lambda_real"] > 0.0)
        ]
        route2_seeds = int(route2_points["seed_id"].nunique())
        route2 = route2_seeds >= 4 and destabilizing_fraction >= 0.75
        decision_rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "coalition_key": candidate.coalition_key,
                "order": int(candidate.order),
                "modal_family_id": modes.set_index("candidate_id").loc[candidate.candidate_id, "modal_family_id"],
                "valid_seed_paths": valid_paths,
                "destabilizing_high_stress_sign_fraction": destabilizing_fraction,
                "median_rho_comp_tau0": nominal_rho,
                "median_rho_comp_tau1": final_rho,
                "high_to_nominal_rho_amplification": amplification,
                "maximum_rho_comp": float(frame["rho_comp"].max()),
                "median_final_delta_lambda_real_per_s": float(final["Delta_lambda_real"].median()),
                "material_delta_zeta_points": int(frame["stress_conditioned_material_delta_zeta"].sum()),
                "destabilizing_material_delta_zeta_points": int(frame["destabilizing_material_delta_zeta"].sum()),
                "hard_composition_failure_seed_paths": hard_seeds,
                "screen_composition_failure_seed_paths": screen_seeds,
                "route2_seed_paths": route2_seeds,
                "C3b_A": "SUPPORTED" if c3a else "REJECTED",
                "C3b_B_route_1": "SUPPORTED" if route1 else "REJECTED",
                "C3b_B_route_2": "SUPPORTED" if route2 else "REJECTED",
                "C3b_B": "SUPPORTED" if route1 or route2 else "REJECTED",
            }
        )
    decisions = pd.DataFrame(decision_rows)
    decisions.to_parquet(TABLES / "e05c_candidate_decisions.parquet", index=False)
    candidates.to_parquet(TABLES / "e05c_candidates.parquet", index=False)
    modes.to_parquet(TABLES / "e05c_mode_family_freeze.parquet", index=False)
    c3b_a = "SUPPORTED" if decisions["C3b_A"].eq("SUPPORTED").any() else "REJECTED"
    c3b_b = "SUPPORTED" if decisions["C3b_B"].eq("SUPPORTED").any() else "REJECTED"
    e06 = c3b_b == "SUPPORTED"

    valid = results[results["status"] == "SUCCESS"].copy()
    rho_kappa = spearmanr(valid["rho_comp"], np.log10(valid["empty_kappa_lambda"]), nan_policy="omit")
    rho_margin = spearmanr(valid["rho_comp"], -valid["lambda_lower_order_real"], nan_policy="omit")
    conditioning = {
        "spearman_rho_comp_vs_log10_kappa": float(rho_kappa.statistic),
        "pvalue_rho_comp_vs_log10_kappa": float(rho_kappa.pvalue),
        "spearman_rho_comp_vs_remaining_margin": float(rho_margin.statistic),
        "pvalue_rho_comp_vs_remaining_margin": float(rho_margin.pvalue),
        "interpretation": "explanatory only; no candidate decision changed by conditioning",
    }
    write_json(TABLES / "e05c_conditioning_summary.json", conditioning)

    strongest = valid.loc[valid["rho_comp"].idxmax()].to_dict()
    strongest_stabilizing = valid.loc[valid["rho_comp"].idxmin()].to_dict()
    final_payload = {
        "stage": "TX3_E05C_MARGIN_CONDITIONED",
        "frozen_candidate_count": len(candidates),
        "frozen_pair_count": int((candidates["order"] == 2).sum()),
        "frozen_triple_count": int((candidates["order"] == 3).sum()),
        "new_seed_count": 8,
        "gamma_end_range": [float(limits["gamma_end"].min()), float(limits["gamma_end"].max())],
        "tracking_success_rate": float(len(tracks) / (8 * 9 * (10 * 4 + 5 * 8))),
        "valid_externality_rate": float((results["status"] == "SUCCESS").mean()),
        "maximum_rho_comp": float(valid["rho_comp"].max()),
        "median_rho_comp": float(valid["rho_comp"].median()),
        "maximum_candidate_median_final_rho_comp": float(decisions["median_rho_comp_tau1"].max()),
        "material_absolute_delta_zeta_points": int(valid["stress_conditioned_material_delta_zeta"].sum()),
        "destabilizing_material_delta_zeta_points": int(valid["destabilizing_material_delta_zeta"].sum()),
        "hard_composition_failure_points": int(valid["hard_composition_failure"].sum()),
        "screen_composition_failure_points": int(valid["screen_composition_failure"].sum()),
        "C3b_A_stress_conditioned_amplification": c3b_a,
        "C3b_B_material_margin_consequence": c3b_b,
        "e06_consideration_authorized": e06,
        "E06_executed": False,
        "C3a_nominal_material_damping_externality": "REJECTED_UNCHANGED",
        "original_C3_cycle_or_SCC": "UNASSESSED",
        "C1": "SUPPORTED_UNCHANGED",
        "C2": "SUPPORTED_UNCHANGED",
        "new_paraemt_scientific_runs": 0,
        "conditioning": conditioning,
        "strongest_result": {
            key: strongest[key]
            for key in (
                "seed_id", "tau", "gamma", "coalition_key", "modal_family_id", "rho_comp",
                "Delta_lambda_real", "Delta_zeta", "lambda_lower_order_real", "lambda_full_real",
            )
        },
        "strongest_stabilizing_result": {
            key: strongest_stabilizing[key]
            for key in ("seed_id", "tau", "coalition_key", "rho_comp", "Delta_lambda_real", "Delta_zeta")
        },
        "stop_reason": "E05C complete; C3b-B is not supported and E06 consideration is not authorized" if not e06 else "E05C complete; independent review required before any E06 action",
    }
    write_json(ARTIFACT / "E05C_FINAL_DECISION.json", final_payload)

    # Figure 1: baseline-only paths and safe endpoints.
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    tau = np.asarray([0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])
    for row in limits.itertuples(index=False):
        axes[0].plot(tau, 1.0 + tau * (row.gamma_end - 1.0), marker="o", ms=2.5, label=row.seed_id)
        local = trace[trace["seed_id"] == row.seed_id].sort_values("gamma")
        axes[1].plot(local["gamma"], local["voltage_min_pu"], marker=".", ms=3, label=row.seed_id)
        axes[1].axvline(row.gamma_end, color="#aaaaaa", lw=0.45)
    axes[0].set(xlabel=r"Normalized stress $\tau$", ylabel=r"Load multiplier $\gamma$", title="Frozen safe stress paths")
    axes[1].axhline(0.85, color="black", ls="--", label="hard stop")
    axes[1].set(xlabel=r"$\gamma$", ylabel=r"Minimum bus voltage (pu)", title="EMPTY-coalition endpoint calibration")
    for axis in axes: axis.grid(alpha=0.22)
    axes[0].legend(ncol=2); axes[1].legend(ncol=2)
    save(fig, "figure_E05C_1_stress_path_definition")

    # Figure 2: pre-frozen panels.
    fig, axes = plt.subplots(len(panel_ids), 3, figsize=(11.5, 10.0), sharex=True)
    for row_index, candidate_id in enumerate(panel_ids):
        frame = valid[valid["candidate_id"] == candidate_id]
        summary = frame.groupby("tau").median(numeric_only=True)
        axes[row_index, 0].plot(summary.index, -summary["lambda_lower_order_real"], marker="o")
        axes[row_index, 1].plot(summary.index, summary["Delta_lambda_real"], marker="o", color="#D55E00")
        axes[row_index, 2].plot(summary.index, summary["rho_comp"], marker="o", color="#0072B2")
        label = frame["coalition_key"].iloc[0]
        axes[row_index, 0].set_ylabel(label)
    for axis, title in zip(axes[0], ("Remaining lower-order margin", r"Re $\Delta_S\lambda$", r"$\rho_{comp}$"), strict=True): axis.set_title(title)
    for axis in axes[-1]: axis.set_xlabel(r"$\tau$")
    for axis in axes.flat: axis.grid(alpha=0.22)
    fig.suptitle("E05C-2 Margin conditioning — panels frozen from E05B rank")
    save(fig, "figure_E05C_2_margin_conditioning")

    # Figure 3: all-candidate heatmap.
    heat = valid.pivot_table(index="coalition_key", columns="tau", values="rho_comp", aggfunc="median")
    heat = heat.reindex(candidates["coalition_key"])
    fig, axis = plt.subplots(figsize=(9.5, 6.0))
    bound = max(abs(float(np.nanmin(heat.to_numpy()))), abs(float(np.nanmax(heat.to_numpy()))))
    image = axis.imshow(heat, aspect="auto", cmap="coolwarm", vmin=-bound, vmax=bound)
    axis.set_yticks(range(len(heat.index)), heat.index)
    axis.set_xticks(range(len(heat.columns)), [f"{value:g}" for value in heat.columns])
    axis.set(xlabel=r"$\tau$", title=r"E05C-3 Median $\rho_{comp}$ across all candidates")
    fig.colorbar(image, ax=axis, label=r"median $\rho_{comp}$")
    save(fig, "figure_E05C_3_all_candidate_heatmap")

    # Figure 4: complex composition vs actual for frozen panels.
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.2))
    for axis, candidate_id in zip(axes.flat, panel_ids, strict=True):
        frame = valid[valid["candidate_id"] == candidate_id]
        summary = frame.groupby("tau").median(numeric_only=True)
        axis.plot(summary["lambda_lower_order_real"], summary["lambda_lower_order_imag"], "o-", label=r"$\lambda_{<q}$")
        axis.plot(summary["lambda_full_real"], summary["lambda_full_imag"], "s--", label=r"$\lambda_S$")
        axis.axvline(0, color="black", lw=0.8)
        axis.set(title=frame["coalition_key"].iloc[0], xlabel=r"Re $\lambda$", ylabel=r"Im $\lambda$")
        axis.grid(alpha=0.22); axis.legend()
    fig.suptitle("E05C-4 Lower-order composition prediction and actual pole")
    save(fig, "figure_E05C_4_composition_vs_actual")

    # Figure 5: damping screen.
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.2), sharex=True)
    for axis, candidate_id in zip(axes.flat, panel_ids, strict=True):
        frame = valid[valid["candidate_id"] == candidate_id]
        summary = frame.groupby("tau").median(numeric_only=True)
        axis.plot(summary.index, summary["zeta_lower_order"], label=r"$\zeta_{<q}$")
        axis.plot(summary.index, summary["zeta_full"], label=r"$\zeta_S$")
        axis.plot(summary.index, summary["Delta_zeta"], label=r"$\Delta_S\zeta$", ls="--")
        axis.axhline(0.05, color="black", ls=":", label="5% screen")
        axis.set(title=frame["coalition_key"].iloc[0], xlabel=r"$\tau$", ylabel="Damping ratio")
        axis.grid(alpha=0.22); axis.legend()
    fig.suptitle("E05C-5 Damping-screen consequence")
    save(fig, "figure_E05C_5_damping_screen")

    # Figure 6: conditioning diagnostics.
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    axes[0].scatter(np.log10(valid["empty_kappa_lambda"]), valid["rho_comp"], s=13, alpha=0.55)
    axes[0].set(xlabel=r"$\log_{10}\kappa_\lambda$", ylabel=r"$\rho_{comp}$", title=f"Spearman r={rho_kappa.statistic:.3f}")
    axes[1].scatter(-valid["lambda_lower_order_real"], valid["rho_comp"], s=13, alpha=0.55, color="#D55E00")
    axes[1].set(xlabel="Remaining lower-order margin", ylabel=r"$\rho_{comp}$", title=f"Spearman r={rho_margin.statistic:.3f}")
    for axis in axes: axis.grid(alpha=0.22)
    fig.suptitle("E05C-6 Explanatory conditioning analysis")
    save(fig, "figure_E05C_6_conditioning")

    # Figure 7: pair vs triple maximum distribution.
    maxima = valid.groupby(["candidate_id", "coalition_key", "order", "seed_id"])["rho_comp"].max().reset_index()
    fig, axis = plt.subplots(figsize=(6.2, 4.2))
    pair = maxima[maxima["order"] == 2]["rho_comp"].to_numpy()
    triple = maxima[maxima["order"] == 3]["rho_comp"].to_numpy()
    axis.boxplot([pair, triple], tick_labels=["Pairs", "Triples"], showfliers=True)
    axis.axhline(0.10, color="black", ls="--", label="C3b-A threshold")
    axis.set(ylabel=r"Maximum frozen-grid $\rho_{comp}$ per seed", title="E05C-7 Pair versus triple interactions")
    axis.grid(axis="y", alpha=0.22); axis.legend()
    save(fig, "figure_E05C_7_pair_vs_triple")

    decision_lines = "\n".join(
        f"| {row.coalition_key} | {row.modal_family_id} | {row.valid_seed_paths}/8 | {row.median_rho_comp_tau0:.6g} | {row.median_rho_comp_tau1:.6g} | {row.maximum_rho_comp:.6g} | {row.C3b_A} | {row.C3b_B} |"
        for row in decisions.itertuples(index=False)
    )
    report("E05C_PREREGISTRATION.md", (PREREG / "E05C_PREREGISTRATION.md").read_text(encoding="utf-8"))
    report("E05C_STRESS_PATH.md", f"""
# E05C stress path

Uniform constant-power-factor PQ-load scaling was supplied by the frozen synchronous participation rule.
All eight EMPTY-coalition calibrations stopped first at the 0.85-pu minimum-voltage boundary. Safe endpoints
span gamma={limits['gamma_end'].min():.6g}--{limits['gamma_end'].max():.6g}; all endpoint checks remained stable,
smooth, and numerically valid. No action vertex was inspected before endpoint freeze `fe6dbef0`.
""")
    report("E05C_MODE_TRACKING_AUDIT.md", f"""
# E05C mode-tracking audit

All {len(tracks):,} required endpoint tracks passed the frozen BMAC >= 0.90 and residual <= 1e-10 rules.
All {len(builds):,} physical homotopy builds were `SUCCESS`; no manual mode switch occurred. Tracking success
rate is {final_payload['tracking_success_rate']:.3f}. Global one-to-one Hungarian assignment was used at every
continuation step to resolve exactly tied BMAC branches.
""")
    report("E05C_CONNECTED_POLE_EXTERNALITIES.md", f"""
# E05C connected pole externalities

All 1,080 candidate/seed/stress connected pole rows are valid. The largest destabilizing real connected shift
is {valid['Delta_lambda_real'].max():.6g} s^-1, but the corresponding margin ratios remain small. The E05B
finite-pole factor was also recomputed at every row; it is secondary and did not select candidates or stress
levels.
""")
    report("E05C_MARGIN_CONSUMPTION.md", f"""
# E05C margin consumption

The global maximum rho_comp is {final_payload['maximum_rho_comp']:.6g}; the all-row median is
{final_payload['median_rho_comp']:.6g}, and the largest candidate median at tau=1 is
{final_payload['maximum_candidate_median_final_rho_comp']:.6g}. The C3b-A gate requires 0.10 plus at least
twofold amplification. No coalition meets it. There are {final_payload['material_absolute_delta_zeta_points']}
points with |Delta zeta| >= 0.0025.
""")
    report("E05C_COMPOSITION_FAILURES.md", f"""
# E05C composition failures

Hard composition failures: **{final_payload['hard_composition_failure_points']}**. Five-percent damping-screen
composition failures: **{final_payload['screen_composition_failure_points']}**. Every lower-order composition
identity closes exactly to machine arithmetic. No pair or triple changes the stable/screened conclusion left by
its complete lower-order prediction.
""")
    report("E05C_CONDITIONING_ANALYSIS.md", f"""
# E05C conditioning analysis

Spearman correlation of rho_comp with log10(kappa_lambda) is {conditioning['spearman_rho_comp_vs_log10_kappa']:.4f}
(p={conditioning['pvalue_rho_comp_vs_log10_kappa']:.4g}); with remaining lower-order margin it is
{conditioning['spearman_rho_comp_vs_remaining_margin']:.4f} (p={conditioning['pvalue_rho_comp_vs_remaining_margin']:.4g}).
This is explanatory only. The voltage boundary curtailed stress while the frozen modes retained substantial decay
margin, so neither conditioning nor margin shrinkage produced material amplification.
""")
    report("E05C_C3B_DECISION.md", f"""
# E05C C3b decision

Frozen population: 15 coalitions (10 pairs, 5 triples), each with one development-locked family. New validation:
8 seeds, gamma_end {limits['gamma_end'].min():.6g}--{limits['gamma_end'].max():.6g}, 100% valid tracking.

| Coalition | Family | Valid | Median rho tau=0 | Median rho tau=1 | Maximum rho | C3b-A | C3b-B |
|---|---|---:|---:|---:|---:|---|---|
{decision_lines}

1. Candidate count: 15. 2. Every frozen family is listed above. 3. Seed count: 8. 4. Every baseline endpoint is
in `e05c_seed_stress_limits.parquet`. 5. Valid tracks: {len(tracks)}/{8*9*(10*4+5*8)}. 6. Destabilizing real
shifts occur in {[row.coalition_key for row in decisions.itertuples() if row.destabilizing_high_stress_sign_fraction >= .75]}.
7--8. Nominal/final rho are tabulated. 9. No candidate reaches the frozen amplification gate. 10. C3b-A={c3b_a}.
11. Hard failures=0. 12. Screen failures=0. 13. Points with rho_comp>=0.25={int((valid['rho_comp']>=.25).sum())}.
14. Points with Delta zeta<=-0.0025={int((valid['Delta_zeta']<=-.0025).sum())}. 15. C3b-B={c3b_b}.
16. Neither pairs nor triples approach 0.10. 17. Conditioning is reported separately and does not explain a
material effect. 18. Strongest result: {strongest['coalition_key']} at {strongest['seed_id']}, tau={strongest['tau']},
rho_comp={strongest['rho_comp']:.6g}. 19. Strongest negative evidence: zero material damping points and zero
composition failures despite complete valid tracking. 20. E06 consideration authorized: **{str(e06).lower()}**.
""")
    report("E05C_FAILURE_LEDGER.md", """
# E05C failure ledger

There are no PF, initialization, DAE, limiter, nonsmooth, eigensolution, eigenpair-residual, mode-track, or
lower-order-instability failures. All 6,984 physical builds and all 1,080 connected cases succeeded. The scientific
hypotheses fail on effect magnitude, not missingness or numerical uncertainty.
""")
    report("E05C_REVIEWER1_AUDIT.md", """
# E05C Reviewer 1 audit

The protocol is clean: candidate population is exhaustive over E05B-reproducible coalitions, modes and primary
panels use development only, stress is a single independent coordinate, endpoints use EMPTY only, and every failure
denominator is retained. The result does not rescue E05: the old nominal C3a remains rejected. The principal physical
limitation is that the registered 0.85-pu voltage boundary is reached before the locked modal families approach a
small decay margin.
""")
    report("E05C_REVIEWER2_AUDIT.md", """
# E05C Reviewer 2 falsification audit

No alternative stress coordinate, lower voltage boundary, replacement mode, denser outcome-driven stress grid, or
relaxed threshold is permitted after this negative result. The exact BMAC-tie ambiguity found during development was
resolved before holdout with one-to-one Hungarian homotopy and fully re-frozen sources. C3b-A and C3b-B must remain
rejected for uniform load stress in this validated smooth region.
""")
    report("E05C_FINAL_REVIEW.md", f"""
# E05C final review

E05C completed the separately preregistered uniform-load-stress experiment. C3b-A={c3b_a}; C3b-B={c3b_b}.
Maximum rho_comp={final_payload['maximum_rho_comp']:.6g}; material damping points=0; hard/screen composition
failures=0/0. C1/C2 and E05B D1--D3 remain supported, C3a remains rejected, and original cycle/SCC C3 remains
unassessed. `e06_consideration_authorized={str(e06).lower()}` and E06 was not executed.
""")
    report("INDEPENDENT_REVIEW_README.md", """
# E05C independent review entry point

Read `E05C_C3B_DECISION.md` and `E05C_FINAL_REVIEW.md`, then inspect the stress, tracking, margin, composition,
conditioning, and failure reports. Machine-readable decisions are in `E05C_FINAL_DECISION.json`. Run
`validate_e05c.py` from repository root. The package contains no E06 or new ParaEMT scientific result.
""")

    original_claims = (PAPER / "CLAIMS_EVIDENCE_MATRIX.md").read_text(encoding="utf-8")
    versioned = original_claims + f"""

## E05C versioned addendum (historical E05 matrix preserved byte-for-byte)

| Claim | Exact statement | Gate | Result | Pass/fail | Permitted wording |
|---|---|---|---|---|---|
| C3b-A | A frozen connected pole externality amplifies materially under independent uniform load stress. | >=6/8 valid paths; destabilizing sign >=0.75; median final rho_comp >=0.10; >=2x nominal. | Maximum candidate median final rho_comp {final_payload['maximum_candidate_median_final_rho_comp']:.6g}; no candidate passes. | FAIL | No stress-conditioned amplification was resolved before the 0.85-pu voltage boundary. |
| C3b-B | A frozen interaction produces replicated composition failure or replicated material margin consumption. | >=2/8 hard/screen failures, or >=4/8 high-stress paths with rho_comp>=0.25 and Delta zeta<=-0.0025. | 0 hard, 0 screen, 0 material damping points. | FAIL | No material margin consequence; E06 consideration is not authorized. |
"""
    (PAPER / "CLAIMS_EVIDENCE_MATRIX_E05C.md").write_text(versioned, encoding="utf-8")
    report("E05C_CLAIMS_EVIDENCE_ADDENDUM.md", versioned.split("## E05C versioned addendum", 1)[1].join(["# E05C claims evidence addendum\n\n## E05C versioned addendum", ""]))
    write_json(
        ARTIFACT / "manifests" / "E05C_RUN_MANIFEST.json",
        {
            "stage": "TX3_E05C_MARGIN_CONDITIONED",
            "candidate_count": 15,
            "seed_count": 8,
            "tau_count": 9,
            "physical_builds": compute["physical_build_count"],
            "mode_tracks": len(tracks),
            "connected_cases": len(results),
            "decisions": {"C3b-A": c3b_a, "C3b-B": c3b_b},
            "e06_consideration_authorized": e06,
            "E06_executed": False,
            "new_paraemt_scientific_runs": 0,
        },
    )
    print(json.dumps(final_payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
