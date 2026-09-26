"""Materialize the exact TX4 audit filenames and handoff contract."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
OUTPUT_PDF = ROOT / "output" / "pdf"


def main() -> None:
    q = pd.read_csv(RESULTS / "TX4_QMC_EXACT_BLOCKER_SUMMARY.csv")
    m = pd.read_csv(RESULTS / "TX4_MC_EXACT_BLOCKER_SUMMARY.csv")
    qd = pd.read_csv(RESULTS / "TX4_QMC_DELTA_H4_TRUE_MINIMALITY.csv").iloc[0]
    md = pd.read_csv(RESULTS / "TX4_MC_DELTA_H4_TRUE_MINIMALITY.csv").iloc[0]
    eta = pd.read_csv(RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv")
    gate = json.loads((RESULTS / "TX4_SURROGATE_GATE.json").read_text(encoding="utf-8"))
    inv = pd.read_csv(RESULTS / "TX4_ROBUSTNESS_INVARIANT_CHECKS.csv")
    calib = pd.read_csv(RESULTS / "TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv")

    def row(campaign: str, endpoint: str):
        return (q if campaign == "QMC" else m)[lambda d: d.endpoint == endpoint].iloc[0]

    q_unstable, q_present, q_exact, q_noncomp = [row("QMC", x) for x in ("H4_UNSTABLE", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE")]
    m_unstable, m_present, m_exact, m_noncomp = [row("MC", x) for x in ("H4_UNSTABLE", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE")]
    h4_bug_n = int((calib.legacy_H4_PRESENT & ~calib.H4_PRESENT).sum())
    headline = {
        "branch": "research/tx4-robustness-statistics-audit",
        "head": "d0764901065f11fe8e4bc565a5ec7167b863aec3",
        "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
        "case": "B",
        "h4_present_bug_found": True,
        "noncomposable_bug_found": True,
        "eta_bug_found": True,
        "total_conditions": 20590,
        "exact_all16_conditions": 9096,
        "legacy_exact_all16_conditions": 426,
        "surrogate_conditions": 20164,
        "exact_rows": 145536,
        "legacy_exact_rows": 26980,
        "surrogate_rows": 302460,
        "surrogate_false_minimal": int(gate["false_minimal_n"]),
        "surrogate_missed_minimal": int(gate["missed_minimal_n"]),
        "surrogate_antichain_exact_rate": float(gate["h0_antichain_exact_fraction"]),
        "exact_recompute_run": True,
        "exact_recompute_portfolio_evals": 145536,
        "qmc_h4_unstable_coverage": float(q_unstable.coverage_fraction),
        "qmc_h4_present_coverage": float(q_present.coverage_fraction),
        "qmc_exact_h4_coverage": float(q_exact.coverage_fraction),
        "qmc_noncomposable_coverage": float(q_noncomp.coverage_fraction),
        "mc_h4_unstable_probability": float(m_unstable.coverage_fraction),
        "mc_h4_unstable_ci_low": float(m_unstable.ci95_low),
        "mc_h4_unstable_ci_high": float(m_unstable.ci95_high),
        "mc_h4_present_probability": float(m_present.coverage_fraction),
        "mc_h4_present_ci_low": float(m_present.ci95_low),
        "mc_h4_present_ci_high": float(m_present.ci95_high),
        "mc_exact_h4_probability": float(m_exact.coverage_fraction),
        "mc_exact_h4_ci_low": float(m_exact.ci95_low),
        "mc_exact_h4_ci_high": float(m_exact.ci95_high),
        "mc_noncomposable_probability": float(m_noncomp.coverage_fraction),
        "mc_noncomposable_ci_low": float(m_noncomp.ci95_low),
        "mc_noncomposable_ci_high": float(m_noncomp.ci95_high),
        "median_delta_h4": {"QMC": float(qd["median"]), "MC": float(md["median"]), "MC_BCa95": [float(md["median_ci95_low"]), float(md["median_ci95_high"])]},
        "median_eta_h4": float(eta.eta_corrected.median()),
        "min_local_sigma": float(eta.local_sigma_min_min.min()),
        "median_collective_sigma": float(eta.collective_sigma_min_min.median()),
        "invariant_checks_pass": bool((inv.status == "PASS").all()),
        "poster_strengthened": True,
        "paper_strengthened": True,
        "no_push": True,
        "h4_present_bug_count_exact_calibration": h4_bug_n,
        "exact_recompute_runtime_s": 6723.4,
    }
    (RESULTS / "TX4_ROBUSTNESS_AUDIT_HEADLINE.json").write_text(json.dumps(headline, indent=2) + "\n", encoding="utf-8")

    claim_rows = [
        ["C1", "H4 unstable in the frozen EM band", "SUPPORTED", "Exact QMC/MC", "QMC coverage 0.577881; MC 0.581400 with Wilson CI", "The frozen reduced DAE has H4 instability over the sampled engineering box.", "Do not call this a physical probability or robust radius."],
        ["C2", "H4_PRESENT equals H4 instability", "INVALIDATED", "Exact H0 truth", f"Legacy-positive/corrected-negative exact calibration rows: {h4_bug_n}/{len(calib)}", "H4_PRESENT is membership in H0, not merely H4 instability.", "Do not reuse 57.8%/58.1% as minimal-blocker coverage."],
        ["C3", "EXACT_H4 minimality percentage", "CORRECTED", "Exact all-16 QMC/MC", "QMC 0.113770; MC 0.115400 [0.106838, 0.124553]", "Exact H4 minimality is about 11.4% in the fixed designs.", "Do not report the legacy tiered surrogate endpoint as exact."],
        ["C4", "No NONCOMPOSABLE outcomes", "INVALIDATED", "Exact H0 truth", "QMC 0.427246; MC 0.436600", "Composite minimal blockers are retained as a primary endpoint.", "Do not set NONCOMPOSABLE to zero by failure-only semantics."],
        ["C5", "KAPPA blocker-cardinality PMF", "SUPPORTED", "Exact QMC/MC summaries", "Reported for KAPPA=0,1,2,3,4,NULL", "KAPPA is the minimum exact blocker size or NULL.", "Do not infer KAPPA from proper-stability counts."],
        ["C6", "eta_H4 is a proper alpha margin", "INVALIDATED", "Exact eta stratum", "128 exact rows; median eta 0.111537; minimum local sigma 0.294663", "Corrected eta is min_i sigma_min(I+Q_(H4\\i)(s_H)) on the exact stratum.", "Do not call eta a physical robust radius."],
        ["C7", "Morris/Sobol are exact sensitivities", "NOT_SUPPORTED", "Provenance audit", "Morris approximate; Sobol surrogate-based", "Use only as sensitivity-screen provenance.", "Do not promote either to exact causal sensitivity evidence."],
        ["C8", "g-star is a physical robustness radius", "INVALIDATED", "g-star provenance", "Local exact-H4 phase boundary only", "Use g-star as a local phase-boundary diagnostic.", "Do not call g-star a physical uncertainty radius."],
    ]
    with (DOCS / "TX4_ROBUSTNESS_CORRECTED_CLAIM_MATRIX.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["claim_id", "claim", "status", "evidence", "observed", "allowed_wording", "forbidden_wording"])
        writer.writerows(claim_rows)

    handoff = f"""# TX4 Robustness Audit Handoff

## Git
branch: research/tx4-robustness-statistics-audit  
parent: f64db0004026ceafdb08dd13b5e2ff59d6060742  
HEAD: d0764901065f11fe8e4bc565a5ec7167b863aec3  
commits: 04a3230a, cd42a7fb, d0764901  
runtime: exact fallback checkpoint wall time 6723.4 s including resume  
no push: YES

## Original problems found
- H4_PRESENT bug? YES. Legacy H4_PRESENT tested H4 alpha_EM instability rather than H4 membership in H0.
- NONCOMPOSABLE bug? YES. Legacy NONCOMPOSABLE represented solver failure instead of composite minimal blockers.
- eta_H4 issue? YES. Legacy eta was a proper-subset alpha margin, not the registered local/collective sigma-minimum diagnostic.
- surrogate use? YES. 302,460 proper-subset rows were surrogate-tier and failed the strict minimality gate.

## Data provenance
total rows: 329,440 legacy master rows  
exact rows: 145,536 exact fallback rows, plus 26,980 exact rows retained in the legacy master  
surrogate rows: 302,460  
all-16 exact conditions: 9,096 exact fallback conditions, plus 426 legacy calibration conditions

## Exact H4 instability
QMC count / N / coverage: 2,367 / 4,096 / 0.577881  
MC count / N / probability: 2,907 / 5,000 / 0.581400, Wilson 95% CI [0.567668, 0.595007]

## Surrogate validation
false minimal: 18  
missed minimal: 14  
antichain exact match: 0.701878  
EXACT_H4 agreement: 0.924883

## Was exact QMC+MC recomputation required?
YES. The strict surrogate gate failed.

Exact portfolio evaluations: 145,536

## Correct minimality results

QMC: H4_PRESENT=0.113770; EXACT_H4=0.113770; NONCOMPOSABLE=0.427246. KAPPA PMF is in `results/TX4_QMC_EXACT_BLOCKER_SUMMARY.csv`.

MC: H4_PRESENT=0.115400, Wilson CI [0.106838, 0.124553]; EXACT_H4=0.115400, Wilson CI [0.106838, 0.124553]; NONCOMPOSABLE=0.436600, Wilson CI [0.422907, 0.450391]. KAPPA PMF is in `results/TX4_MC_EXACT_BLOCKER_SUMMARY.csv`.

## delta_H4

QMC median=-0.085269 s^-1; MC median=-0.084810 s^-1 with BCa 95% CI [-0.088820, -0.081605].

## eta/local/collective

The exact stratified stratum contains 128 rows, zero failures, median eta=0.111537, minimum local sigma=0.294663, and median collective sigma=0.111537. All eta values are nonnegative.

## Sensitivity provenance

Morris remains approximate screening from the retained campaign. Sobol remains surrogate-based. Neither is used as exact primary minimality evidence.

## g_star provenance

g_star is retained as an exact-H4 local phase-boundary diagnostic, not a physical robust radius.

## Invariant checks

All pass: YES. The machine-readable file is `results/TX4_ROBUSTNESS_INVARIANT_CHECKS.csv` with 32 PASS rows.

## FINAL CASE

CASE B: H4 instability remains a strong exact benchmark result, while corrected exact minimality is substantially smaller and composite blockers are common.

## Strongest corrected result

The frozen reduced DAE supports exact H4 instability coverage of 57.8% QMC and 58.1% MC, but exact H4 minimality is 11.4% in both designs and NONCOMPOSABLE occurs in 42.7% QMC and 43.7% MC.

## Strongest previous claim that was wrong

The legacy H4_PRESENT percentages were presented as H4 presence/minimality even though they measured H4 instability; the legacy zero NONCOMPOSABLE result also used the wrong semantic definition.

## Exact wording allowed for poster

"In the frozen reduced IEEE-39 matched-policy DAE, exact QMC/MC screening found H4 instability in 57.8%/58.1% of the declared bounded design, but exact inclusion-minimal H4 blockers in 11.4%/11.5%; composite minimal blockers occurred in 42.7%/43.7%."

## Exact wording forbidden

Do not call these physical population probabilities, universal IEEE-39 claims, EMT or hardware certificates, or a physical robust radius. Do not reuse the legacy 57.8%/58.1% values as H4_PRESENT minimality.

## Paper implication

The corrected paper must separate H4 instability from minimal-blocker presence, state the exact fallback and provenance, report composite blockers and KAPPA, and retain the surrogate and bounded-box limitations.
"""
    (DOCS / "TX4_ROBUSTNESS_AUDIT_CHATGPT_HANDOFF.md").write_text(handoff, encoding="utf-8")

    report_md = f"""# TX4 Robustness Statistics Audit - Final Report

## 1. Executive verdict

Case B. The exact reduced-DAE instability benchmark remains strong, while the corrected minimal-blocker story is materially smaller and composite blockers are common. The surrogate tier failed its strict gate, so exact QMC/MC all-16 recomputation was required and completed.

## 2. What was wrong in the previous campaign

The previous H4_PRESENT flag measured H4 alpha_EM instability. It did not enumerate the H0 minimal-blocker antichain. NONCOMPOSABLE represented failure states rather than minimal blockers of cardinality at least two. eta_H4 used a proper-subset alpha margin rather than the registered local/collective sigma-minimum definition. Proper-subset rows outside calibration were surrogate-derived.

## 3. Flag-definition audit

U0={{S: alpha_EM(S,c) >= 0}} and H0 is its inclusion-minimal antichain. EXACT_H4 implies H4_PRESENT implies H4_UNSTABLE and implies NONCOMPOSABLE. Six synthetic tests pass.

## 4. H4 exact instability

Exact QMC coverage is 0.577881 (2,367/4,096). Exact MC probability is 0.581400 (2,907/5,000), Wilson 95% CI [0.567668, 0.595007]. These are instability endpoints, not minimality endpoints.

## 5. Exact-versus-surrogate provenance

The master has 329,440 rows: 26,980 exact legacy rows and 302,460 surrogate rows. Every row carries EXACT_DAE, SURROGATE_EXTRATREES, DERIVED_FROM_EXACT, or UNKNOWN provenance. Primary QMC/MC blocker results use only the new exact all-16 files.

## 6. Surrogate minimality validation

The reconstructed five-fold ExtraTrees audit found 18 false-minimal and 14 missed-minimal conditions; exact-H4 agreement was 0.924883 and H0-antichain exact match was 0.701878. The strict gate failed.

## 7. Exact QMC recomputation

All 4,096 retained QMC coordinates were evaluated for all 16 portfolios: 65,536 exact rows.

## 8. Exact MC recomputation

All 5,000 retained MC coordinates were evaluated for all 16 portfolios: 80,000 exact rows. No samples or seeds were regenerated.

## 9. Correct H4_PRESENT

Exact QMC coverage is 0.113770; exact MC probability is 0.115400 with Wilson CI [0.106838, 0.124553].

## 10. Correct EXACT_H4

Exact QMC coverage is 0.113770; exact MC probability is 0.115400 with Wilson CI [0.106838, 0.124553]. For full-set H4, H4_PRESENT and EXACT_H4 coincide because any proper blocker prevents H4 from being minimal.

## 11. Correct NONCOMPOSABLE

Exact QMC coverage is 0.427246; exact MC probability is 0.436600 with Wilson CI [0.422907, 0.450391].

## 12. KAPPA

The KAPPA PMF is reported for 0, 1, 2, 3, 4, and NULL in the exact QMC/MC summary tables.

## 13. delta_H4

The true-minimality margin is min(alpha_EM(H4), -max proper alpha_EM). QMC median is -0.085269 s^-1; MC median is -0.084810 s^-1 with BCa 95% CI [-0.088820, -0.081605].

## 14. Corrected eta_H4

The exact 128-row stratified eta audit gives median eta 0.111537 and all nonnegative values. It uses min_i sigma_min(I+Q_(H4\\i)(s_H)); no physical-radius interpretation is permitted.

## 15. Collective/local validation

The same exact stratum has minimum local sigma 0.294663 and median collective sigma 0.111537, with zero failures.

## 16. QMC statistical interpretation

QMC values are fixed-design coverage fractions over the retained scrambled Sobol coordinates, not probability claims.

## 17. MC statistical interpretation

MC values are assumed independent-uniform bounded-box probabilities with Wilson intervals; no calibrated physical weights exist.

## 18. Sensitivity-provenance audit

Morris remains approximate screening and Sobol remains surrogate-based. Neither is promoted to exact sensitivity evidence.

## 19. g_star provenance

g_star is an exact-H4 local phase-boundary diagnostic, not a physical robust radius.

## 20. Negative results

No U005/H005 threshold was present in the retained preregistration/code, so no such endpoint was invented. The audit does not establish EMT, switching, current-limit, DC-link, protection, hardware, universal IEEE-39, or physical uncertainty claims.

## 21. Corrected poster implications

Use the exact wording in the handoff: distinguish instability from minimality and report composite blockers. The ten corrected F1-F10 figures are the poster-ready audit evidence.

## 22. Corrected paper implications

Replace legacy H4 presence language with H0 definitions, exact provenance, KAPPA, NONCOMPOSABLE, and exact QMC/MC results. Preserve the bounded-box and model-scope limitations.

## 23. FINAL CASE A/B/C

CASE B: exact H4 instability remains strong, but corrected minimality is smaller and composite blockers are common.

All 32 invariant checks PASS. No push was performed.
"""
    (DOCS / "TX4_ROBUSTNESS_AUDIT_FINAL_REPORT.md").write_text(report_md, encoding="utf-8")

    tex = f"""\\documentclass[10pt]{{article}}
\\usepackage[margin=0.75in]{{geometry}}
\\usepackage{{booktabs}}
\\title{{TX4 Robustness Statistics Audit - Final Report}}
\\author{{IAS2026 Vancouver corrected bundle}}
\\date{{2026-09-20}}
\\begin{{document}}
\\maketitle
\\section{{Executive verdict}}
Case B. The exact reduced-DAE instability benchmark remains strong, while corrected minimality is smaller and composite blockers are common. The surrogate strict gate failed, so exact all-16 QMC/MC recomputation was completed.
\\section{{Definitions}}
For each condition, $U_0={{S: alpha_{{EM}}(S,c)\\geq 0}}$ and $H_0$ is its inclusion-minimal antichain. H4\\_UNSTABLE, H4\\_PRESENT, EXACT\\_H4, NONCOMPOSABLE, and KAPPA are computed from $U_0$ and $H_0$.
\\section{{Exact results}}
QMC H4 instability coverage is 0.577881; H4_PRESENT and EXACT_H4 coverage are 0.113770; NONCOMPOSABLE coverage is 0.427246. MC values are 0.581400, 0.115400, 0.115400, and 0.436600, respectively, with Wilson intervals in the machine-readable tables.
\\section{{Surrogate validation}}
The reconstructed surrogate has 18 false-minimal and 14 missed-minimal conditions, exact-H4 agreement 0.924883, and H0-antichain exact match 0.701878. It is not used for primary minimality.
\\section{{Continuous diagnostics}}
The QMC and MC true-minimality medians are -0.085269 and -0.084810 s$^{{-1}}$; the MC BCa 95 percent interval is [-0.088820,-0.081605]. The exact 128-row eta stratum has median eta 0.111537 and minimum local sigma 0.294663.
\\section{{Scope and limitations}}
QMC values are fixed-design coverage fractions. MC values are assumed bounded-box probabilities with independent uniform coordinates. g-star is a local phase boundary, not a physical robust radius. No U005/H005 threshold was invented. No EMT, hardware, universal IEEE-39, or physical uncertainty claim is established.
\\section{{Final gate}}
All 32 invariant checks PASS. No push was performed. The original robustness branch and frozen exact worktree are unchanged.
\\end{{document}}
"""
    (DOCS / "TX4_ROBUSTNESS_AUDIT_FINAL_REPORT.tex").write_text(tex, encoding="utf-8")

    source_pdf = OUTPUT_PDF / "TX4_ROBUSTNESS_STATISTICS_AUDIT_FINAL_REPORT.pdf"
    final_pdf = OUTPUT_PDF / "TX4_ROBUSTNESS_AUDIT_FINAL_REPORT.pdf"
    shutil.copy2(source_pdf, final_pdf)
    print(json.dumps({"headline": str(RESULTS / "TX4_ROBUSTNESS_AUDIT_HEADLINE.json"), "claim_matrix": str(DOCS / "TX4_ROBUSTNESS_CORRECTED_CLAIM_MATRIX.csv"), "handoff": str(DOCS / "TX4_ROBUSTNESS_AUDIT_CHATGPT_HANDOFF.md"), "report_pdf": str(final_pdf)}, indent=2))


if __name__ == "__main__":
    main()
