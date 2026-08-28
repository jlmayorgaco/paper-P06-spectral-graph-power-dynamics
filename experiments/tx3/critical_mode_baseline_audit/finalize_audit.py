from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import andes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT / "artifacts" / "tx3" / "critical_mode_baseline_audit"
FIGURES = AUDIT / "figures"
TABLES = AUDIT / "tables"
PAPER = ROOT / "reports" / "papers" / "tx3_connected_intervention_calculus"
E05C = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned"
AUDIT_SOURCE = ROOT / "experiments" / "tx3" / "critical_mode_baseline_audit"
sys.path[:0] = [str(ROOT), str(E05C), str(AUDIT_SOURCE)]

from margin_conditioned import damping_ratio, evaluate_stressed_operator, generalized_modes  # noqa: E402
from run_audit import CASE, mode_anatomy, state_names  # noqa: E402


plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 8.5, "legend.fontsize": 7.2})


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    modes = pd.read_parquet(TABLES / "critical_modes.parquet")
    anatomy = pd.read_parquet(TABLES / "critical_mode_state_anatomy.parquet")
    subsets = pd.read_parquet(TABLES / "critical_mode_subset_poles.parquet")
    results = pd.read_parquet(TABLES / "critical_mode_externalities.parquet")
    builds = pd.read_parquet(TABLES / "audit_failure_ledger.parquet")
    compute = load_json(AUDIT / "logs" / "audit_compute.json")
    valid = results[results["status"] == "SUCCESS"].copy()

    # Reconstruct the development modes used for the externality-dominant E05C locks.
    op00 = {"seed_id": "OP00", "load_scale": 1.0, "gfl_dispatch_scale": 1.0, "reactive_load_scale": 1.0, "redispatch_mw": 0.0}
    point = evaluate_stressed_operator(op00, 1.0, {})
    if point.status != "SUCCESS": raise RuntimeError(f"OP00 comparison baseline failed: {point.status}")
    development = generalized_modes(point); names = state_names()
    locks = pd.read_parquet(ROOT / "artifacts/tx3/E05C_margin_conditioned/preregistration/e05c_mode_family_freeze.parquet")
    comparison_rows = []
    for lock in locks.itertuples(index=False):
        index = int(lock.reference_mode_index); value = complex(development.eigenvalues[index])
        _, summary = mode_anatomy(development, index, names)
        comparison_rows.append({
            "coalition": lock.coalition_key, "externality_mode_index": index,
            "externality_mode_family_id": lock.modal_family_id,
            "externality_mode_frequency_hz": float(value.imag / (2.0 * math.pi)),
            "externality_mode_damping_ratio": damping_ratio(value),
            "externality_mode_classification": summary["classification"],
            "externality_mode_synchronous_participation": summary["participation_scores"]["synchronous_electromechanical"],
            "externality_mode_pll_participation": summary["participation_scores"]["PLL_dominated"],
            "externality_mode_current_control_participation": summary["participation_scores"]["current_control_dominated"],
            "critical_mode_median_frequency_hz": float(modes["frequency_hz"].median()),
            "critical_mode_median_damping_ratio": float(modes["damping_ratio"].median()),
            "critical_mode_modal_family": modes["classification"].mode().iat[0],
            "critical_mode_median_synchronous_participation": float(modes["synchronous_participation"].median()),
            "critical_mode_median_pll_participation": float(modes["pll_participation"].median()),
            "critical_mode_median_current_control_participation": float(modes["current_control_participation"].median()),
        })
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_parquet(TABLES / "externality_vs_critical_mode_comparison.parquet", index=False)

    candidate_rows = []
    for coalition, frame in valid.groupby("coalition", sort=False):
        candidate_rows.append({
            "coalition": coalition, "order": int(frame["order"].iat[0]), "valid_seeds": frame["seed_id"].nunique(),
            "median_Delta_lambda_real": float(frame["Delta_lambda_real"].median()),
            "minimum_Delta_lambda_real": float(frame["Delta_lambda_real"].min()),
            "maximum_Delta_lambda_real": float(frame["Delta_lambda_real"].max()),
            "median_Delta_zeta": float(frame["Delta_zeta"].median()),
            "minimum_Delta_zeta": float(frame["Delta_zeta"].min()),
            "maximum_Delta_zeta": float(frame["Delta_zeta"].max()),
            "median_rho_comp": float(frame["rho_comp"].median()),
            "minimum_rho_comp": float(frame["rho_comp"].min()),
            "maximum_rho_comp": float(frame["rho_comp"].max()),
            "destabilizing_real_fraction": float((frame["Delta_lambda_real"] > 0.0).mean()),
            "hard_failure_count_descriptive": int(frame["hard_failure"].sum()),
            "screen_failure_count_descriptive": int(frame["screen_failure"].sum()),
        })
    summary = pd.DataFrame(candidate_rows)
    order = load_json(ROOT / "artifacts/tx3/E05D_weak_grid/preregistration/E05D_CANDIDATE_FREEZE.json")["candidates"]
    ordering = {record["coalition_key"]: index for index, record in enumerate(order)}
    summary["sort"] = summary["coalition"].map(ordering); summary = summary.sort_values("sort").drop(columns="sort")
    summary.to_parquet(TABLES / "critical_mode_coalition_summary.parquet", index=False)

    colors = {"synchronous_participation": "#0072B2", "pll_participation": "#D55E00", "current_control_participation": "#009E73", "other_participation": "#999999"}
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.5), sharex=True)
    bottoms = np.zeros(len(modes))
    for column, label in (("synchronous_participation", "Synchronous"), ("pll_participation", "PLL"), ("current_control_participation", "Converter controls"), ("other_participation", "Other")):
        axes[0].bar(modes["seed_id"], modes[column], bottom=bottoms, color=colors[column], label=label)
        bottoms += modes[column].to_numpy()
    bottoms = np.zeros(len(modes))
    for column, label, color in (("synchronous_endogeneity", "Synchronous", "#0072B2"), ("pll_endogeneity", "PLL", "#D55E00"), ("current_control_endogeneity", "Converter controls", "#009E73"), ("other_endogeneity", "Other", "#999999")):
        axes[1].bar(modes["seed_id"], modes[column], bottom=bottoms, color=color, label=label)
        bottoms += modes[column].to_numpy()
    axes[0].set(ylabel="Normalized |left×right|", title="Critical-mode participation anatomy")
    axes[1].set(ylabel="Normalized descriptor endogeneity", xlabel="Frozen E05D seed", title="Critical-mode endogeneity anatomy")
    for axis in axes: axis.set_ylim(0, 1); axis.grid(axis="y", alpha=0.2); axis.legend(ncol=4, loc="upper center")
    save(fig, "figure_1_critical_mode_participation_anatomy")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    axes[0].scatter(modes["seed_id"], modes["frequency_hz"], color="#0072B2")
    axes[0].set(ylabel="Frequency (Hz)", title="Margin-setting mode frequency")
    axes[1].scatter(modes["seed_id"], modes["damping_ratio"], color="#D55E00")
    axes[1].axhline(0.05, color="black", ls="--", label="5% screen")
    axes[1].set(ylabel="Damping ratio", title="Margin-setting mode damping"); axes[1].legend()
    for axis in axes: axis.tick_params(axis="x", rotation=45); axis.grid(alpha=0.2)
    save(fig, "figure_2_critical_mode_frequency_damping")

    fig, axis = plt.subplots(figsize=(10.2, 4.8))
    positions = np.arange(len(summary)); data = [valid[valid["coalition"] == coalition]["rho_comp"].to_numpy() for coalition in summary["coalition"]]
    box = axis.boxplot(data, positions=positions, widths=0.65, patch_artist=True, showfliers=True)
    for patch, candidate_order in zip(box["boxes"], summary["order"], strict=True): patch.set_facecolor("#56B4E9" if candidate_order == 2 else "#E69F00")
    for position, values in zip(positions, data, strict=True): axis.scatter(np.full(len(values), position), values, s=10, color="#333333", alpha=0.5)
    axis.axhline(0, color="black", lw=0.8); axis.axhline(0.25, color="black", ls="--", label="historical materiality reference (not an audit gate)")
    axis.set_xticks(positions, summary["coalition"], rotation=55, ha="right")
    axis.set(ylabel=r"Descriptive $\rho_{comp}$ at $\kappa_{grid}=1$", title="All 15 frozen coalitions on the actual margin-setting mode")
    axis.grid(axis="y", alpha=0.2); axis.legend()
    save(fig, "figure_3_critical_mode_rho_all_coalitions")

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4))
    axes[0].scatter(comparison["coalition"], comparison["externality_mode_damping_ratio"], color="#D55E00", label="Externality-dominant development mode")
    axes[0].axhspan(float(modes["damping_ratio"].min()), float(modes["damping_ratio"].max()), color="#0072B2", alpha=0.25, label="Actual critical-mode range")
    axes[0].set(ylabel="Damping ratio", title="Damping separation"); axes[0].tick_params(axis="x", rotation=60); axes[0].legend()
    axes[1].scatter(comparison["externality_mode_frequency_hz"], comparison["externality_mode_current_control_participation"], color="#009E73", label="Externality-dominant")
    axes[1].scatter([modes["frequency_hz"].median()], [modes["current_control_participation"].median()], marker="*", s=180, color="#0072B2", label="Actual critical median")
    axes[1].set(xlabel="Frequency (Hz)", ylabel="Converter-control participation", title="Physical-family separation")
    for axis in axes: axis.grid(alpha=0.2)
    axes[1].legend()
    save(fig, "figure_4_externality_vs_margin_mode")

    strongest = valid.loc[valid["rho_comp"].idxmax()].to_dict()
    weakest = valid.loc[valid["rho_comp"].idxmin()].to_dict()
    classification_counts = modes["classification"].value_counts().to_dict()
    decision = {
        "scope": "DESCRIPTIVE_POST_MORTEM_ONLY", "claim_gate_created": False,
        "statuses_unchanged": {"C3a": "REJECTED", "C3b-A": "REJECTED", "C3b-B": "REJECTED", "C3c": "UNRESOLVED"},
        "critical_mode_classification_counts": classification_counts,
        "critical_damping_range": [float(modes["damping_ratio"].min()), float(modes["damping_ratio"].max())],
        "critical_frequency_range_hz": [float(modes["frequency_hz"].min()), float(modes["frequency_hz"].max())],
        "critical_synchronous_participation_range": [float(modes["synchronous_participation"].min()), float(modes["synchronous_participation"].max())],
        "critical_converter_control_participation_range": [float(modes["current_control_participation"].min()), float(modes["current_control_participation"].max())],
        "externality_mode_damping_range": [float(comparison["externality_mode_damping_ratio"].min()), float(comparison["externality_mode_damping_ratio"].max())],
        "tracking_success_rate": float(subsets["status"].eq("SUCCESS").mean()),
        "connected_case_success_rate": float(results["status"].eq("SUCCESS").mean()),
        "maximum_rho_comp": float(valid["rho_comp"].max()), "minimum_rho_comp": float(valid["rho_comp"].min()),
        "median_absolute_rho_comp": float(valid["rho_comp"].abs().median()),
        "maximum_absolute_Delta_zeta": float(valid["Delta_zeta"].abs().max()),
        "hard_failure_count_descriptive": int(valid["hard_failure"].sum()),
        "screen_failure_count_descriptive": int(valid["screen_failure"].sum()),
        "strongest_destabilizing": {key: strongest[key] for key in ("seed_id", "coalition", "rho_comp", "Delta_lambda_real", "Delta_zeta", "lambda_lower_real", "lambda_full_real")},
        "strongest_stabilizing": {key: weakest[key] for key in ("seed_id", "coalition", "rho_comp", "Delta_lambda_real", "Delta_zeta")},
        "interpretation": "the actual margin-setting mode is synchronous-electromechanical, whereas E05B-population externality-dominant development modes are converter-control/PLL families with much larger damping",
        "mechanism_surgery_authorized": False, "E06_authorized": False,
    }
    write_json(AUDIT / "CRITICAL_MODE_BASELINE_AUDIT_DECISION.json", decision)

    table_lines = "\n".join(
        f"| {row.coalition} | {row.order} | {row.valid_seeds}/8 | {row.median_Delta_lambda_real:.4e} | {row.median_Delta_zeta:.4e} | {row.median_rho_comp:.4e} | {row.minimum_rho_comp:.4e} | {row.maximum_rho_comp:.4e} |"
        for row in summary.itertuples(index=False)
    )
    report = f"""# TX3 Critical-Mode Baseline Audit

## Scope

This is a descriptive post-mortem at `kappa_grid=1`. It uses the eight frozen E05D seeds, all 15 frozen
E05B dynamic coalitions, and the original A1--A8 endpoints. It creates no claim gate and cannot change
`C3a=REJECTED`, `C3b-A=REJECTED`, `C3b-B=REJECTED`, or `C3c=UNRESOLVED`.

## Actual margin-setting mode

All eight seeds select the same physical family under the frozen rule: {classification_counts}. Frequencies
span {modes['frequency_hz'].min():.5f}--{modes['frequency_hz'].max():.5f} Hz and damping spans
{modes['damping_ratio'].min():.6f}--{modes['damping_ratio'].max():.6f}. Synchronous-machine plus governor/
excitation participation is {modes['synchronous_participation'].min():.4f}--{modes['synchronous_participation'].max():.4f};
PLL participation is {modes['pll_participation'].min():.4f}--{modes['pll_participation'].max():.4f}; converter-control
participation is {modes['current_control_participation'].min():.4f}--{modes['current_control_participation'].max():.4f}.
The descriptor-weighted endogeneity audit gives the same family-level interpretation. Classification uses a
fixed 0.60 dominance threshold on normalized absolute left/right participation; all state-level values and the
left/right vectors are released.

## Connected effects on that mode

Every one of the {len(subsets)} unique seed/subset tracks and {len(results)} candidate/seed connected cases is
retained. Tracking success is {compute['mode_track_success_rate']:.3%}; the minimum successful BMAC is
{subsets.loc[subsets.status.eq('SUCCESS'), 'BMAC'].min():.6f}. The global rho_comp range is
{valid['rho_comp'].min():.6g}--{valid['rho_comp'].max():.6g}; median absolute rho_comp is
{valid['rho_comp'].abs().median():.6g}. The largest absolute damping externality is
{valid['Delta_zeta'].abs().max():.6g}. Hard/screen composition failures are
{int(valid['hard_failure'].sum())}/{int(valid['screen_failure'].sum())}. These values are descriptive and do
not retroactively activate any historical gate.

| Coalition | Order | Valid | Median Re Delta lambda | Median Delta zeta | Median rho_comp | Min rho | Max rho |
|---|---:|---:|---:|---:|---:|---:|---:|
{table_lines}

## Why the externality and margin stories differ

The E05C development locks, which selected the largest real pole externality from the E05B population, have
damping {comparison['externality_mode_damping_ratio'].min():.4f}--{comparison['externality_mode_damping_ratio'].max():.4f}
and are dominated by PLL/current-command/electrical-control states. The true margin-setting family has roughly
{modes['damping_ratio'].median():.4f} damping and is synchronous-electromechanical. Thus the poles carrying the
demonstrated controller externalities are not the poles setting system damping margin in this benchmark.

The audit supports an explanation, not a new hypothesis test: finite controller interactions can be real and
reproducible without materially consuming the decay margin of a different, synchronous-dominated critical mode.
No cycle/SCC causality, mechanism surgery, E06, or new ParaEMT claim is authorized.
"""
    report_path = AUDIT / "TX3_CRITICAL_MODE_BASELINE_AUDIT.md"
    report_path.write_text(report, encoding="utf-8")
    (PAPER / "TX3_CRITICAL_MODE_BASELINE_AUDIT.md").write_text(report, encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
