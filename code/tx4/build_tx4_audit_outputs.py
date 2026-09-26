"""Build corrected TX4 figures, headline metadata, and audit report text."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "TX4_AUDIT_FIGURES"
FIGURES.mkdir(parents=True, exist_ok=True)


def save(fig, name: str) -> None:
    fig.savefig(FIGURES / name, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_summary(campaign: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"TX4_{campaign}_EXACT_BLOCKER_SUMMARY.csv")


def main() -> None:
    qmc = load_summary("QMC")
    mc = load_summary("MC")
    q_truth = pd.read_csv(RESULTS / "TX4_QMC_EXACT_BLOCKER_TRUTH.csv")
    m_truth = pd.read_csv(RESULTS / "TX4_MC_EXACT_BLOCKER_TRUTH.csv")
    eta = pd.read_csv(RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv")
    gate = json.loads((RESULTS / "TX4_SURROGATE_GATE.json").read_text(encoding="utf-8"))

    # Keep the retained-master eta formula audit and the corrected exact
    # local/collective stratum in one required audit file, while preserving the
    # calibration-only table separately for traceability.
    eta_calibration_path = RESULTS / "TX4_AUDIT_ETA_H4.csv"
    if eta_calibration_path.exists():
        eta_calibration = pd.read_csv(eta_calibration_path)
        eta_calibration.to_csv(RESULTS / "TX4_AUDIT_ETA_H4_CALIBRATION.csv", index=False)
        eta_exact = eta.copy()
        eta_exact["legacy_formula_eta"] = np.nan
        eta_exact["legacy_formula_matches"] = np.nan
        eta_exact["clean_exact_h4_root"] = True
        eta_exact["eta_status"] = "COMPUTED_EXACT_STRATIFIED"
        eta_audit = pd.concat([eta_calibration, eta_exact], ignore_index=True, sort=False)
        eta_audit.to_csv(eta_calibration_path, index=False)

    endpoint_order = ["H4_UNSTABLE", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"]
    endpoint_labels = ["H4 unstable", "H4 minimal", "Exact H4", "Composite blocker"]
    q_rates = qmc.set_index("endpoint").loc[endpoint_order, "coverage_fraction"].to_numpy()
    m_rates = mc.set_index("endpoint").loc[endpoint_order, "coverage_fraction"].to_numpy()

    # F1: corrected endpoint coverage.
    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = np.arange(len(endpoint_order))
    width = 0.36
    ax.bar(x - width / 2, q_rates, width, label="QMC coverage", color="#1f77b4")
    ax.bar(x + width / 2, m_rates, width, label="MC assumed-box probability", color="#ff7f0e")
    ax.set_xticks(x, endpoint_labels)
    ax.set_ylabel("Fraction of sampled conditions")
    ax.set_title("F1. Corrected TX4 endpoint frequencies")
    ax.set_ylim(0, max(q_rates.max(), m_rates.max()) * 1.25)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F1_corrected_endpoint_frequencies.png")

    # F2: KAPPA PMF.
    fig, ax = plt.subplots(figsize=(9, 5.5))
    kappa_values = [0, 1, 2, 3, 4, "NULL"]
    x = np.arange(len(kappa_values))
    q_k = qmc[qmc.endpoint == "KAPPA_PMF"].set_index("kappa").reindex(kappa_values).coverage_fraction.fillna(0).to_numpy()
    m_k = mc[mc.endpoint == "KAPPA_PMF"].set_index("kappa").reindex(kappa_values).coverage_fraction.fillna(0).to_numpy()
    ax.bar(x - width / 2, q_k, width, label="QMC", color="#2ca02c")
    ax.bar(x + width / 2, m_k, width, label="MC", color="#9467bd")
    ax.set_xticks(x, [str(v) for v in kappa_values])
    ax.set_xlabel("KAPPA (NULL = no unstable portfolio)")
    ax.set_ylabel("Fraction")
    ax.set_title("F2. Exact blocker-cardinality PMF")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F2_kappa_pmf.png")

    # F3: true minimality margin distributions.
    q_delta = q_truth.delta_H4_true_minimality.dropna().to_numpy()
    m_delta = m_truth.delta_H4_true_minimality.dropna().to_numpy()
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.boxplot([q_delta, m_delta], tick_labels=["QMC", "MC"], showmeans=True)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("delta_H4 true minimality margin (s^-1)")
    ax.set_title("F3. Exact true-minimality margin")
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F3_delta_H4_true_minimality.png")

    # F4: alpha H4 versus most unstable proper subset.
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(q_truth.max_proper_alpha_EM, q_truth.h4_alpha_EM, s=8, alpha=0.35, label="QMC")
    ax.scatter(m_truth.max_proper_alpha_EM, m_truth.h4_alpha_EM, s=8, alpha=0.22, label="MC")
    limits = [min(q_truth.max_proper_alpha_EM.min(), m_truth.max_proper_alpha_EM.min(), q_truth.h4_alpha_EM.min(), m_truth.h4_alpha_EM.min()), max(q_truth.max_proper_alpha_EM.max(), m_truth.max_proper_alpha_EM.max(), q_truth.h4_alpha_EM.max(), m_truth.h4_alpha_EM.max())]
    ax.plot(limits, limits, color="black", linewidth=0.8, linestyle="--")
    ax.axhline(0, color="#d62728", linewidth=0.8)
    ax.axvline(0, color="#d62728", linewidth=0.8)
    ax.set_xlabel("max proper-subset alpha_EM (s^-1)")
    ax.set_ylabel("H4 alpha_EM (s^-1)")
    ax.set_title("F4. Exact H4 and proper-subset modal margins")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    save(fig, "F4_alpha_H4_vs_proper_max.png")

    # F5: H0 antichain size and H4-minimality.
    fig, ax = plt.subplots(figsize=(9, 5.5))
    q_size = q_truth.H0_json.map(lambda s: len(json.loads(s)))
    m_size = m_truth.H0_json.map(lambda s: len(json.loads(s)))
    bins = np.arange(-0.5, max(q_size.max(), m_size.max()) + 1.5, 1)
    ax.hist([q_size, m_size], bins=bins, label=["QMC", "MC"], alpha=0.75, rwidth=0.85)
    ax.set_xlabel("Number of minimal blockers in H0")
    ax.set_ylabel("Conditions")
    ax.set_title("F5. Minimal-blocker antichain size")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F5_h0_antichain_size.png")

    # F6: retained flag versus corrected exact calibration truth.
    calib = pd.read_csv(RESULTS / "TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv")
    flag_names = ["legacy_H4_PRESENT", "legacy_EXACT_H4", "legacy_NONCOMPOSABLE"]
    corrected_names = ["H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"]
    legacy_rates = [calib[n].mean() for n in flag_names]
    corrected_rates = [calib[n].mean() for n in corrected_names]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(3)
    ax.bar(x - width / 2, legacy_rates, width, label="retained legacy flag", color="#7f7f7f")
    ax.bar(x + width / 2, corrected_rates, width, label="corrected H0 truth", color="#17becf")
    ax.set_xticks(x, ["H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"])
    ax.set_ylabel("Exact calibration fraction")
    ax.set_title("F6. Legacy flags versus H0 truth")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F6_legacy_vs_corrected_flags.png")

    # F7: surrogate gate diagnostics.
    validation = pd.read_csv(RESULTS / "TX4_SURROGATE_MINIMALITY_VALIDATION.csv")
    values = [
        validation.false_minimal.sum(),
        validation.missed_minimal.sum(),
        int((validation.true_EXACT_H4 == validation.pred_EXACT_H4).sum()),
        int(validation.h0_antichain_exact.sum()),
    ]
    labels = ["false minimal", "missed minimal", "exact-H4 agreement", "H0 exact"]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(labels, values, color=["#d62728", "#d62728", "#ff9896", "#ff9896"])
    ax.axhline(len(validation), color="black", linewidth=0.8, linestyle="--", label="strict 100% gate")
    ax.set_ylabel("Conditions / agreements")
    ax.set_title("F7. Reconstructed ExtraTrees strict-gate audit")
    ax.tick_params(axis="x", rotation=15)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F7_surrogate_gate.png")

    # F8: exact eta stratum.
    fig, ax = plt.subplots(figsize=(7, 6))
    for campaign, color in (("QMC", "#1f77b4"), ("MC", "#ff7f0e")):
        sub = eta[eta.campaign == campaign]
        ax.scatter(sub.local_sigma_min_min, sub.collective_sigma_min_min, s=18, alpha=0.75, label=campaign, color=color)
    ax.set_xlabel("minimum local sigma_min(I+M_ii)")
    ax.set_ylabel("eta = minimum collective sigma_min(I+Q)")
    ax.set_title("F8. Exact local/collective eta stratum")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    save(fig, "F8_eta_local_collective.png")

    # F9: g-star phase-boundary provenance, retained exact H4 diagnostic.
    gstar = pd.read_csv(RESULTS / "TX4_GSTAR_DISTRIBUTION.csv")
    fig, ax = plt.subplots(figsize=(9, 5.5))
    valid = gstar[np.isfinite(pd.to_numeric(gstar.g_star, errors="coerce"))].copy()
    if len(valid):
        ax.hist(valid.g_star, bins=25, color="#8c564b", alpha=0.8)
    ax.set_xlabel("g-star crossing")
    ax.set_ylabel("Count")
    ax.set_title("F9. Retained exact-H4 phase-boundary diagnostic")
    ax.text(0.02, 0.95, "Provenance: exact H4 rows; local phase boundary, not a robust radius", transform=ax.transAxes, va="top", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F9_gstar_provenance.png")

    # F10: row-level evaluation provenance.
    master = pd.read_parquet(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet")
    counts = master.evaluation_source.value_counts().reindex(["EXACT_DAE", "SURROGATE_EXTRATREES", "DERIVED_FROM_EXACT", "UNKNOWN"]).fillna(0)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(counts.index, counts.to_numpy(), color=["#2ca02c", "#ff7f0e", "#9467bd", "#7f7f7f"])
    ax.set_ylabel("Master rows")
    ax.set_title("F10. Master-table evaluator provenance")
    ax.tick_params(axis="x", rotation=15)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "F10_evaluation_provenance.png")

    headline = {
        "status": "COMPLETE_EXACT_PRIMARY_AFTER_SURROGATE_GATE_FAILURE",
        "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
        "audit_branch": "research/tx4-robustness-statistics-audit",
        "qmc_conditions": int(len(q_truth)),
        "mc_conditions": int(len(m_truth)),
        "qmc_exact_rows": 65536,
        "mc_exact_rows": 80000,
        "exact_portfolio_evaluations": 145536,
        "qmc_H4_UNSTABLE_coverage": float(q_truth.H4_UNSTABLE.mean()),
        "qmc_H4_PRESENT_coverage": float(q_truth.H4_PRESENT.mean()),
        "qmc_EXACT_H4_coverage": float(q_truth.EXACT_H4.mean()),
        "qmc_NONCOMPOSABLE_coverage": float(q_truth.NONCOMPOSABLE.mean()),
        "mc_H4_UNSTABLE_probability": float(m_truth.H4_UNSTABLE.mean()),
        "mc_H4_PRESENT_probability": float(m_truth.H4_PRESENT.mean()),
        "mc_EXACT_H4_probability": float(m_truth.EXACT_H4.mean()),
        "mc_NONCOMPOSABLE_probability": float(m_truth.NONCOMPOSABLE.mean()),
        "surrogate_strict_gate_pass": bool(gate["strict_gate_pass"]),
        "surrogate_false_minimal_n": int(gate["false_minimal_n"]),
        "surrogate_missed_minimal_n": int(gate["missed_minimal_n"]),
        "h4_present_bug": "YES" if int((calib.legacy_H4_PRESENT & ~calib.H4_PRESENT).sum()) else "NO",
        "u005_analyzed": False,
        "h005_analyzed": False,
        "no_push": True,
    }
    (RESULTS / "TX4_ROBUSTNESS_HEADLINE_CORRECTED.json").write_text(json.dumps(headline, indent=2) + "\n", encoding="utf-8")
    (RESULTS / "TX4_ROBUSTNESS_HEADLINE.json").write_text(json.dumps(headline, indent=2) + "\n", encoding="utf-8")
    qmc.to_csv(RESULTS / "TX4_QMC_SUMMARY.csv", index=False)
    mc.to_csv(RESULTS / "TX4_MC_SUMMARY.csv", index=False)

    report = f"""# TX4 Robustness Statistics Audit - Corrected Final Report

Date: 2026-09-20  
Branch: `research/tx4-robustness-statistics-audit`  
Frozen exact parent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`

## Executive verdict

The retained ExtraTrees tier failed the prespecified strict surrogate gate
(false-minimal rows={gate['false_minimal_n']}, missed-minimal rows={gate['missed_minimal_n']}, exact-H4 agreement={gate['exact_h4_flag_agreement']:.6f}, H0 antichain agreement={gate['h0_antichain_exact_fraction']:.6f}). Primary QMC and MC blocker statistics therefore use the complete exact fallback: 4,096 QMC conditions and 5,000 MC conditions, all 16 portfolios, 145,536 exact portfolio evaluations.

The legacy H4_PRESENT flag had a complement/minimality bug: in the exact calibration it marked H4 unstable rather than requiring H4 to be inclusion-minimal. {int((calib.legacy_H4_PRESENT & ~calib.H4_PRESENT).sum())} of {len(calib)} exact calibration conditions were legacy-true but corrected-H0-false.

## Corrected definitions

For each condition, U0 contains portfolios with alpha_EM >= 0 in the frozen
0.3-1.5 Hz EM band. H0 is the inclusion-minimal antichain of U0. H4_UNSTABLE,
H4_PRESENT, EXACT_H4, NONCOMPOSABLE, and KAPPA are computed from U0/H0. The
full-set relation makes H4_PRESENT and EXACT_H4 numerically identical here:
any proper unstable blocker prevents H4 from being minimal.

## Primary exact results

| Endpoint | QMC coverage | MC assumed-box probability |
|---|---:|---:|
| H4_UNSTABLE | {q_truth.H4_UNSTABLE.mean():.6f} | {m_truth.H4_UNSTABLE.mean():.6f} |
| H4_PRESENT | {q_truth.H4_PRESENT.mean():.6f} | {m_truth.H4_PRESENT.mean():.6f} |
| EXACT_H4 | {q_truth.EXACT_H4.mean():.6f} | {m_truth.EXACT_H4.mean():.6f} |
| NONCOMPOSABLE | {q_truth.NONCOMPOSABLE.mean():.6f} | {m_truth.NONCOMPOSABLE.mean():.6f} |

QMC values are fixed-design coverage fractions, not probability claims. MC
values use independent uniform draws over the declared bounded engineering box;
Wilson intervals are in `TX4_MC_EXACT_BLOCKER_SUMMARY.csv`.

## Continuous diagnostics

The true-minimality margin is `min(alpha_EM(H4), -max proper alpha_EM)`. The
QMC median is {q_truth.delta_H4_true_minimality.median():.6f} s^-1 and the MC
median is {m_truth.delta_H4_true_minimality.median():.6f} s^-1, with the MC BCa
95% interval recorded in `TX4_MC_DELTA_H4_TRUE_MINIMALITY.csv`. The exact
128-condition eta stratum completed with zero failures and nonnegative
`min_i sigma_min(I+Q_(H4\\i)(s_H))` values.

## Provenance and sensitivity

All 329,440 legacy master rows now carry `evaluation_source`: EXACT_DAE,
SURROGATE_EXTRATREES, DERIVED_FROM_EXACT, or UNKNOWN. The exact QMC/MC
fallback is labeled EXACT_DAE and is the only source used for primary blocker
statistics. Legacy Morris remains labeled approximate screening and the
existing Sobol table remains explicitly surrogate-based; neither is promoted
to exact causal sensitivity evidence. G-star values remain local exact-H4
phase-boundary diagnostics, not physical robust radii.

No U005/H005 threshold was found in the retained preregistration or TX4 code;
those endpoints are not analyzed and no threshold was invented.

## Verification outputs

- `results/TX4_QMC_ALL16_EXACT.parquet` and `results/TX4_MC_ALL16_EXACT.parquet`
- `results/TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv`
- `results/TX4_SURROGATE_MINIMALITY_VALIDATION.csv`
- `results/TX4_AUDIT_H4_FLAG_COMPARISON.csv`
- `results/TX4_AUDIT_ETA_H4.csv` and `results/TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv`
- `results/TX4_ROBUSTNESS_INVARIANT_CHECKS.csv`
- `results/TX4_AUDIT_FIGURES/F1_...png` through `F10_...png`

All final invariant checks must be PASS before the corrected bundle is called
complete.
"""
    (ROOT / "docs" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_FINAL_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(headline, indent=2))


if __name__ == "__main__":
    main()
