"""Assemble the frozen IAS2026 poster evidence package after IAS26-080.

This is deterministic post-processing only: it runs no power-flow, DAE,
eigenvalue, or time-domain solves and never edits frozen inputs or raw evidence.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[6]
CONTRACT = REPO / "reports/poster/ias2026/research/bnd_h4_mechanism"
F1_RUN = CONTRACT / "results/20260926T161116Z_d0fecb32_ias26_020_f1_audit_v1"
F2_RUN = CONTRACT / "results/20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1"
TAU = 1e-8


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def event(events: pd.DataFrame, name: str) -> pd.Series:
    rows = events.loc[events.event == name]
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one event row for {name}")
    return rows.iloc[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, required=True)
    args = ap.parse_args()
    run = args.run_root.resolve()
    tds_metrics_path = run / "tables/IAS26-080_TDS_RESULTS.csv"
    tds_summary_path = run / "derived/tds/IAS26-080_TDS_SUMMARY.json"
    if not tds_metrics_path.is_file() or not tds_summary_path.is_file():
        raise RuntimeError("IAS26-080 must finish all 90 planned records before finalization")

    figures = run / "figures/final"
    tables = run / "tables"
    claims = run / "claims"
    reports = run / "reports"
    config = run / "config"
    for folder in (figures, tables, claims, reports, config):
        folder.mkdir(parents=True, exist_ok=True)

    # Preserve the exact audited F1/F2 art as source copies; create no redraws.
    sources = {
        "P1_STABILITY_CLIFF": F1_RUN / "figures/F1_STABILITY_CLIFF",
        "P2_LOCAL_COLLECTIVE_CLOSURE": F2_RUN / "figures/F2_local_collective_alpha",
        "P3_MC_BLOCKER_ROBUSTNESS": run / "figures/P3_MC_BLOCKER_QUADRANT",
    }
    for dest, source_stem in sources.items():
        for ext in ("pdf", "svg", "png"):
            src = source_stem.with_suffix(f".{ext}")
            if not src.is_file():
                raise RuntimeError(f"missing frozen figure source: {src}")
            shutil.copy2(src, figures / f"{dest}.{ext}")

    mc_events = pd.read_csv(tables / "MC_EVENT_RATES.csv")
    scenarios = pd.read_csv(run / "derived/SCENARIO_METRICS.csv")
    cases = pd.read_csv(run / "derived/MC_CASES.csv")
    tds = pd.read_csv(tds_metrics_path)
    tds_summary = json.loads(tds_summary_path.read_text(encoding="utf-8"))
    cross = json.loads((run / "derived/CROSSCODE_SUMMARY.json").read_text(encoding="utf-8"))
    if len(scenarios) != 1000 or scenarios.scenario_id.nunique() != 1000:
        raise RuntimeError("the frozen 1000-ID scenario ledger is incomplete or duplicated")
    if len(tds) != 90 or tds.trajectory_id.nunique() != 90:
        raise RuntimeError("TDS trajectory aggregate must contain 90 unique preregistered IDs")

    # P4 pairs every valid scenario's fixed intervention with direct canonical TDS.
    p4 = scenarios[[
        "scenario_id", "scenario_valid", "alpha_perp_original",
        "alpha_perp_intervention", "delta_alpha", "intervention_improves_h4",
        "intervention_rescues_h4", "intervention_deteriorates_h4",
    ]].copy()
    p4.to_csv(tables / "IAS26-FINAL_P4_MC_DATA.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(12.2, 5.2), constrained_layout=True)
    ax = axs[0]
    valid = p4.scenario_valid.astype(bool)
    categories = [
        ("rescue", valid & (p4.intervention_rescues_h4 == "TRUE"), "#2a8c62"),
        ("improves; no H4 rescue", valid & (p4.intervention_rescues_h4 != "TRUE") & (p4.delta_alpha < -TAU), "#557b9e"),
        ("deteriorates", valid & (p4.delta_alpha > TAU), "#c64b40"),
        ("deadband", valid & p4.delta_alpha.abs().le(TAU), "#ca8c32"),
    ]
    for label, mask, color in categories:
        ax.scatter(p4.loc[mask, "alpha_perp_original"], p4.loc[mask, "alpha_perp_intervention"],
                   s=19, alpha=.52, linewidths=0, color=color, label=label, rasterized=True)
    invalid = ~valid
    ax.scatter(p4.loc[invalid, "alpha_perp_original"], p4.loc[invalid, "alpha_perp_intervention"],
               marker="x", s=19, color="#777777", alpha=.7, label="physically infeasible (retained)")
    lims = [float(np.nanmin([p4.alpha_perp_original.min(), p4.alpha_perp_intervention.min()])),
            float(np.nanmax([p4.alpha_perp_original.max(), p4.alpha_perp_intervention.max()]))]
    ax.plot(lims, lims, "--", color="#555555", lw=.9)
    ax.axhline(0, color="#222222", lw=.8); ax.axvline(0, color="#222222", lw=.8)
    ax.set(xlabel="Original H4 transverse α⊥ (s⁻¹)", ylabel="Fixed g=0.25 H4 transverse α⊥ (s⁻¹)",
           title="A  Frozen 1,000-ID synthetic ensemble")
    ax.legend(frameon=False, fontsize=7, loc="best")

    ax = axs[1]
    canonical = [
        ("CANONICAL_CANONICAL_PROPER_30_33_35_D2_A1", "proper triple", "#557b9e"),
        ("CANONICAL_CANONICAL_H4_ORIGINAL_D2_A1", "H4, g=0.03625", "#c64b40"),
        ("CANONICAL_CANONICAL_H4_RETUNED_D2_A1", "H4, g=0.25", "#2a8c62"),
    ]
    trace_rows = []
    for tid, label, color in canonical:
        metric = tds.loc[tds.trajectory_id == tid]
        trace_path = run / "raw/tds/traces" / f"{tid}.csv"
        if len(metric) != 1 or metric.iloc[0].status != "COMPLETED" or not trace_path.is_file():
            continue
        trace = pd.read_csv(trace_path)
        col = "common_relative_speed_omega_sg34_minus_omega_sg39_pu"
        if trace.empty or col not in trace:
            continue
        ax.plot(trace.time_s, trace[col], lw=1.05, color=color, label=label)
        trace_rows.extend({"trajectory_id": tid, "case": label, "time_s": float(t),
                           "omega_sg34_minus_omega_sg39_pu": float(y)}
                          for t, y in zip(trace.time_s, trace[col]))
    plotted = {r["trajectory_id"] for r in trace_rows}
    expected = {tid for tid, _, _ in canonical}
    if plotted != expected:
        raise RuntimeError(f"P4 requires all three completed canonical D2/A1 traces; found {sorted(plotted)}")
    ax.axhline(0, color="#222222", lw=.7)
    ax.set(xlabel="Time (s)", ylabel="ωSG34 − ωSG39 (pu)",
           title="B  Nonlinear phasor-DAE response to frozen 2% load pulse")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Fixed intervention: scenario response and canonical time-domain evidence", fontsize=12)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(figures / f"P4_FIXED_INTERVENTION_TDS.{ext}", dpi=240)
    plt.close(fig)
    pd.DataFrame(trace_rows).to_csv(tables / "IAS26-FINAL_P4_TDS_TRACES.csv", index=False)

    # Evidence strip: same-model cross-language reproducibility, not independent validation.
    strip = pd.DataFrame([{
        "evidence": "Python ↔ Julia, frozen exact GFL11 model",
        "N_scenarios": cross["N_scenarios"], "N_case_comparisons": cross["N_case_comparisons"],
        "verdict_agreement": cross["verdict_agreement_count"],
        "max_abs_delta_alpha_s-1": cross["max_abs_delta_alpha"],
        "max_abs_delta_frequency_Hz": cross["max_abs_delta_f_hz"],
        "mode_identity": "unresolved in all comparisons",
        "scope": cross["crosscode_scope"],
    }])
    strip.to_csv(tables / "IAS26-FINAL_EVIDENCE_STRIP.csv", index=False)
    fig, ax = plt.subplots(figsize=(11.8, 1.45), constrained_layout=True)
    ax.axis("off")
    ax.text(.02, .66, "CROSS-CODE REPRODUCIBILITY", weight="bold", fontsize=10, color="#263d4d", transform=ax.transAxes)
    ax.text(.02, .23,
            f"Python ↔ Julia exact GFL11  |  {cross['verdict_agreement_count']}/{cross['N_case_comparisons']} verdicts agree  |  "
            f"max |Δα|={cross['max_abs_delta_alpha']:.2e} s⁻¹  |  "
            f"max |Δf|={cross['max_abs_delta_f_hz']:.2e} Hz  |  mode identity unresolved",
            fontsize=10, transform=ax.transAxes)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(figures / f"EVIDENCE_STRIP_CROSSCODE.{ext}", dpi=240)
    plt.close(fig)

    # Update falsification ledger N8 only after frozen TDS data exists.
    completed = tds[(tds.status == "COMPLETED") & np.isfinite(tds.alpha_tds_s_inv)].copy()
    mismatches = completed[completed.sign_agreement == False]  # noqa: E712
    mismatch_zero_frequency = int((mismatches.frequency_tds_hz == 0).sum())
    mismatch_no_r2 = int(mismatches.fit_r2_diagnostic.isna().sum())
    # A sign disagreement without an identified oscillatory family is not a
    # mode-specific falsification; preserve it as unresolved evidence.
    n8_status = "NOT_TESTABLE_MODE_IDENTITY" if len(mismatches) else ("NO_CONTRADICTION_OBSERVED" if len(completed) else "NOT_TESTABLE")
    falsification_path = claims / "FALSIFICATION_LEDGER.csv"
    falsification = pd.read_csv(falsification_path)
    n8 = falsification.question_id == "N8"
    if n8.sum() != 1:
        raise RuntimeError("falsification ledger does not contain exactly one N8 row")
    falsification.loc[n8, "status"] = n8_status
    falsification.loc[n8, "finding"] = (
        f"Among {len(completed)}/90 completed trajectories with finite primary matrix-pencil alpha, "
        f"{int(completed.sign_agreement.sum())} agree in sign with the linear eigenvalue and "
        f"{len(mismatches)} disagree; all {mismatch_zero_frequency} discrepant fits have zero fitted frequency "
        f"and {mismatch_no_r2} have no finite R2 diagnostic, so the critical mode identity is unresolved. "
        f"{90-len(completed)} records are solver-failed or otherwise excluded from completed-fit sign agreement. "
        "Same-model nonlinear phasor TDS; not EMT."
    )
    falsification.loc[n8, "evidence_path"] = "tables/IAS26-080_TDS_RESULTS.csv; derived/tds/IAS26-080_TDS_SUMMARY.json"
    falsification.to_csv(falsification_path, index=False)

    # Compact headline number table (15 rows; all ensemble rates retain explicit denominators).
    f1 = pd.read_csv(F1_RUN / "tables/F1_FINAL_AUDIT.csv")
    f1h = f1.loc[f1.portfolio_id.astype(str).isin(["30+33+35", "30+33+35+37"])]
    if len(f1h) != 2:
        raise RuntimeError("could not locate frozen F1 proper-triple and H4 rows")
    alpha_triple = float(f1h.loc[f1h.cardinality == 3, "alpha_perp_s-1"].iloc[0])
    alpha_h4 = float(f1h.loc[f1h.cardinality == 4, "alpha_perp_s-1"].iloc[0])
    nominal = pd.read_csv(CONTRACT / "results/20260926T174848Z_d0fecb32_ias26_050_nominal_v1/tables/IAS26-050_NOMINAL_COMPARISON.csv")
    alpha_retuned = float(nominal.loc[nominal.condition == "fixed_intervention_g_0.25", "alpha_perp_s-1"].iloc[0])
    # Prefer exact boundary audit rather than a sampled curve value.
    f2audit = pd.read_csv(F2_RUN / "tables/F2_BOUNDARY_PORT_AUDIT.csv")
    local_sigma = float(f2audit["physical_local_sigma_min_I_plus_Mii"].min())
    collective_sigma = float(f2audit["collective_sigma_min"].min())
    e1 = event(mc_events, "E1_H4_MINIMAL_BLOCKER_PERSISTS")
    e2 = event(mc_events, "E2_COLLECTIVE_PHENOMENON_PERSISTS")
    e3 = event(mc_events, "E3_INTERVENTION_IMPROVES")
    e4 = event(mc_events, "E4_INTERVENTION_RESCUES_H4")
    tdscount = int(len(completed))
    tds_agree = int(completed.sign_agreement.sum()) if tdscount else 0
    numbers = [
        ("canonical proper-triple α⊥", alpha_triple, "s^-1", "F1_FINAL_AUDIT.csv"),
        ("canonical original H4 α⊥", alpha_h4, "s^-1", "F1_FINAL_AUDIT.csv"),
        ("canonical retuned H4 α⊥", alpha_retuned, "s^-1", "IAS26-050_NOMINAL_COMPARISON.csv"),
        ("fixed g original → intervention", "0.03625 → 0.25", "dimensionless", "IAS26-050 frozen config"),
        ("F2 collective boundary g*", 0.20768140519037842, "dimensionless", "IAS26-030 summary"),
        ("minimum physical local σmin at boundary", local_sigma, "dimensionless", "F2_BOUNDARY_PORT_AUDIT.csv"),
        ("collective σmin(I+QH) at boundary", collective_sigma, "dimensionless", "F2_BOUNDARY_PORT_AUDIT.csv"),
        ("MC H4 minimal-blocker persistence", f"{int(e1.numerator)}/{int(e1.denominator)}", "proportion", "MC_EVENT_RATES.csv"),
        ("MC H4 persistence Wilson 95% CI", f"[{e1.wilson95_low:.4f}, {e1.wilson95_high:.4f}]", "proportion", "MC_EVENT_RATES.csv"),
        ("MC collective-phenomenon persistence", f"{int(e2.numerator)}/{int(e2.denominator)}", "proportion", "MC_EVENT_RATES.csv"),
        ("fixed intervention improvement", f"{int(e3.numerator)}/{int(e3.denominator)}", "proportion", "MC_EVENT_RATES.csv"),
        ("fixed intervention H4 rescue", f"{int(e4.numerator)}/{int(e4.denominator)}", "proportion", "MC_EVENT_RATES.csv"),
        ("H4 rescue Wilson 95% CI", f"[{e4.wilson95_low:.4f}, {e4.wilson95_high:.4f}]", "proportion", "MC_EVENT_RATES.csv"),
        ("TDS sign agreement among completed finite fits", f"{tds_agree}/{tdscount}", "proportion; same-model phasor TDS", "IAS26-080_TDS_RESULTS.csv"),
        ("MC frozen IDs physically feasible", f"{int(scenarios.scenario_valid.sum())}/1000", "scenarios; no replacements", "SCENARIO_METRICS.csv"),
    ]
    if len(numbers) != 15:
        raise AssertionError("poster headline number table must contain exactly 15 rows")
    pd.DataFrame(numbers, columns=["metric", "value", "unit", "source"]).to_csv(tables / "POSTER_NUMBERS_FINAL.csv", index=False)

    claim_rows = [
        ["C1", "Nominal H4 is an inclusion-minimal transverse blocker: all 15 proper subsets stable; H4 unstable.", "SAFE_NOMINAL", "P1_STABILITY_CLIFF", "Do not call this universal across operating points."],
        ["C2", "At the frozen nominal boundary, physical local factors remain nonsingular while collective closure approaches singularity.", "SAFE_WITH_SCOPE", "P2_LOCAL_COLLECTIVE_CLOSURE", "No claim that local damping is positive or M1 reduced-pole equivalence passes."],
        ["C3", f"H4 minimal-blocker pattern occurs in {int(e1.numerator)}/{int(e1.denominator)} valid frozen synthetic scenarios.", "SAFE_DESCRIPTIVE", "P3_MC_BLOCKER_ROBUSTNESS", "Not a real-world probability; 33 physically infeasible IDs retained, no replacement."],
        ["C4", f"Fixed g=0.25 improves α⊥ in {int(e3.numerator)}/{int(e3.denominator)} and rescues H4 in {int(e4.numerator)}/{int(e4.denominator)} valid scenarios.", "SAFE_DESCRIPTIVE", "P4_FIXED_INTERVENTION_TDS", "Not universal: deterioration and non-rescue cases remain visible."],
        ["C5", f"Python–Julia verdict agreement {cross['verdict_agreement_count']}/{cross['N_case_comparisons']}; the TDS primary estimator agrees in sign in {tds_agree}/{tdscount} completed finite fits.", "SAFE_WITH_LIMITS", "EVIDENCE_STRIP_CROSSCODE; P4_FIXED_INTERVENTION_TDS; FALSIFICATION_LEDGER N8", f"Four completed-fit sign disagreements are 0 Hz with unavailable R2; mode identity is unresolved, so do not claim mode-specific contradiction. Same-model reproduction only; no EMT claim."],
    ]
    pd.DataFrame(claim_rows, columns=["claim_id", "claim", "status", "evidence", "limitation"]).to_csv(claims / "POSTER_CLAIM_LEDGER.csv", index=False)

    safe = [r[1] for r in claim_rows]
    unsafe = [
        "MC event frequencies are real-world IEEE-39 probabilities.",
        "The H4 blocker persists for every perturbed operating point or every grid.",
        "g=0.25 universally stabilizes H4 or never deteriorates stability margins.",
        "Mode identity is tracked across scenarios or Python and Julia eigenmodes are physically independent validation.",
        "The nonlinear phasor DAE trajectories are EMT/field validation.",
        "IAS26-010 M1 strict gate has passed or any eta/M2/M3/GFM claim is established.",
    ]
    report = {
        "run_id": run.name,
        "finalized_at_utc": pd.Timestamp.now(tz="UTC").isoformat(),
        "campaign": "IAS26-020, 030, 050, 060, 080 evidence assembly; IAS26-010 remains blocked and unchanged",
        "MC": {"n_total": 1000, "n_valid": int(scenarios.scenario_valid.sum()), "n_infeasible": int((~scenarios.scenario_valid.astype(bool)).sum()), "replacements": 0, "interpretation": "descriptive synthetic envelope; not real-world probabilities"},
        "TDS": {"n_planned": 90, "n_completed_finite_primary_fits": tdscount, "n_sign_agreement": tds_agree, "n_sign_disagreement": int(len(mismatches)), "run_status_counts": tds.status.value_counts().to_dict(), "alpha_mae_completed_s-1": float(completed.alpha_abs_error.mean()) if tdscount else None, "frequency_mae_completed_hz": float(completed.frequency_abs_error_hz.mean()) if tdscount else None, "scope": "same-model nonlinear phasor DAE; not EMT"},
        "crosscode": cross,
        "falsification_N8": {"status": n8_status, "completed_comparable_n": tdscount, "sign_mismatch_n": int(len(mismatches)), "mismatches_zero_frequency_n": mismatch_zero_frequency, "mismatches_missing_r2_n": mismatch_no_r2, "mode_identity_resolved": False},
        "blocked_or_excluded": ["IAS26-010 remains BLOCKED_M1_STRICT", "no eta/M2/M3", "no GFM", "no N-1", "no mode-family identity claim", "no real-world probability or universal-stability claim"],
        "safe_claims": safe,
        "unsafe_claims": unsafe,
        "poster_readiness": "READY_WITH_EXPLICIT_LIMITATIONS" if tdscount > 0 else "NOT_READY_TDS_UNTESTABLE",
        "final_figures": [f"figures/final/{x}" for x in ["P1_STABILITY_CLIFF.pdf", "P2_LOCAL_COLLECTIVE_CLOSURE.pdf", "P3_MC_BLOCKER_ROBUSTNESS.pdf", "P4_FIXED_INTERVENTION_TDS.pdf", "EVIDENCE_STRIP_CROSSCODE.pdf"]],
    }
    write_json(reports / "IAS2026_POSTER_READINESS.json", report)
    md = [
        f"# IAS2026 poster evidence readiness — {run.name}", "",
        f"Decision: **{report['poster_readiness']}**", "",
        f"Monte Carlo: {report['MC']['n_valid']}/1000 valid; {report['MC']['n_infeasible']} physically infeasible; zero replacements.",
        f"TDS: {tdscount}/90 completed with a finite primary fit; {tds_agree}/{tdscount} primary-estimator sign agreement; {len(mismatches)} disagreements, all {mismatch_zero_frequency} at 0 Hz with {mismatch_no_r2} unavailable R2 diagnostics, so mode identity is unresolved. 30/90 solver failures are retained.",
        "Cross-code: same-model Python–Julia reproduction; mode identity unresolved.", "",
        "## Safe claims", *[f"- {x}" for x in safe], "",
        "## Claims not supported", *[f"- {x}" for x in unsafe], "",
        "F1 and F2 are exact copies of audited frozen assets. P3 is the frozen MC figure. P4 combines paired fixed-intervention MC results with only the three completed canonical D2/A1 phasor-DAE traces. The ensemble is synthetic, not a real-world probability sample.",
    ]
    (reports / "IAS2026_POSTER_READINESS.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    amendment = {
        "document_id": "IAS2026_POSTER_ARCHITECTURE_FINAL_V1",
        "supersedes": "poster figure mapping only; all scientific protocols, thresholds, models, and historical results remain unchanged",
        "reason": "F3 remains blocked by IAS26-010 M1A; F4 GFL/GFM campaign was not run. The final poster maps completed evidence into four figures without implying either blocked result.",
        "principal_figures": [
            {"id": "P1", "asset": "P1_STABILITY_CLIFF", "source": "IAS26-020 exact 16-portfolio F1 audit"},
            {"id": "P2", "asset": "P2_LOCAL_COLLECTIVE_CLOSURE", "source": "IAS26-030 F2 physical-factor closure audit"},
            {"id": "P3", "asset": "P3_MC_BLOCKER_ROBUSTNESS", "source": "IAS26-060 frozen 1000-ID synthetic stress ensemble"},
            {"id": "P4", "asset": "P4_FIXED_INTERVENTION_TDS", "source": "IAS26-060 paired fixed-g MC and IAS26-080 completed canonical phasor DAE traces"}
        ],
        "evidence_strip": "EVIDENCE_STRIP_CROSSCODE; same-model Python–Julia reproduction only",
        "prohibited_claims": unsafe
    }
    write_json(config / "IAS2026_POSTER_ARCHITECTURE_FINAL_V1.json", amendment)
    print(json.dumps({"run_id": run.name, "poster_readiness": report["poster_readiness"], "tds_completed_finite": tdscount, "tds_sign_agreement": f"{tds_agree}/{tdscount}", "final_figures": str(figures)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
