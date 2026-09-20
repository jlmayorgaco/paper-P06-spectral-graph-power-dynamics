"""Build the TX4 robustness report set and PDFs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
RESULTS = ROOT / "results"
H4 = "30+33+35+37"


def load():
    headline = json.loads((RESULTS / "TX4_ROBUSTNESS_HEADLINE.json").read_text(encoding="utf-8"))
    qmc = pd.read_csv(RESULTS / "TX4_QMC_SUMMARY.csv")
    mc = pd.read_csv(RESULTS / "TX4_MC_SUMMARY.csv")
    method = pd.read_csv(RESULTS / "TX4_METHOD_COMPARISON.csv")
    sobol = pd.read_csv(RESULTS / "TX4_SOBOL.csv")
    logistic = pd.read_csv(RESULTS / "TX4_LOGISTIC_MODEL.csv")
    return headline, qmc, mc, method, sobol, logistic


def fmt(x, digits=4):
    try:
        return f"{float(x):.{digits}g}"
    except (TypeError, ValueError):
        return str(x)


def summary_text():
    headline, qmc, mc, method, sobol, logistic = load()
    qh = qmc[qmc.endpoint == "H4_PRESENT"].iloc[0]
    qe = qmc[qmc.endpoint == "EXACT_H4"].iloc[0]
    mh = mc[mc.endpoint == "H4_PRESENT"].iloc[0]
    me = mc[mc.endpoint == "EXACT_H4"].iloc[0]
    return headline, qh, qe, mh, me, method, sobol, logistic


def sections():
    headline, qh, qe, mh, me, method, sobol, logistic = summary_text()
    return [
        ("1. Executive disposition", f"The TX4 campaign is complete as an explicitly tiered robustness analysis. It covers {headline['conditions']:,} conditions and {headline['portfolio_rows']:,} portfolio rows. Every H4 row is an exact reduced-DAE evaluation. The complete proper-subset census outside the calibration skeleton is a calibrated ExtraTrees surrogate and is never presented as exact DAE evidence."),
        ("2. Frozen-parent discipline", "The branch starts at frozen exact P4/GFL11 parent f64db0004026ceafdb08dd13b5e2ff59d6060742. The earlier exact worktree was kept separate and no push was performed. Frozen nominal evidence remains read-only input."),
        ("3. Archaeology result", "The prior closure marked the common uncertainty envelope, robust V4 census, and g-by-epsilon atlas NOT_TESTED. The retained uncertainty manifest had no calibrated weights and forbade portfolio-specific refitting and calling a closure-space radius physical."),
        ("4. Model contract", "The campaign uses the frozen matched-reactive-power reduced semi-explicit IEEE-39 DAE, central-difference Jacobians, and index-one Schur reduction. It does not add EMT switching, current limits, DC-link dynamics, protection, or hardware behavior."),
        ("5. Parameter box", "The seven coordinates are g in [0.020,0.250], k and t in [0.75,1.75], h in [0.50,1.50], epsilon in [0.90,1.10], damping in [0,0.20], and inertia in [0.80,1.20]. These are bounded engineering assumptions, not a calibrated physical distribution."),
        ("6. Sampling design", "The declared design contains 201-point one-dimensional sweeps, six 41x41 two-dimensional maps, scrambled Sobol QMC N=4096, uniform MC N=5000, and the preregistered sensitivity seeds. All IDs, seeds, and parameter values are in the master table."),
        ("7. Endpoint definitions", "Raw alpha_all is the maximum real part of the complete reduced spectrum. alpha_EM is the maximum real part in the 0.3-1.5 Hz transverse band. Because numerical neutral modes contaminate some proper-subset alpha_all values, blocker flags are normalized from alpha_EM; raw alpha_all is preserved."),
        ("8. Evaluator tiers", "The master table contains 6,816 EXACT_16_DAE rows, 20,164 EXACT_H4_DAE rows for non-calibration conditions, and 302,460 SURROGATE_CALIBRATED proper-subset rows. The calibration source is 426 exact all-portfolio conditions."),
        ("9. One-dimensional sweeps", "The 201-point sweeps expose sign changes and gauge behavior across every primary coordinate. The full table is in TX4_1D_SWEEPS.csv. These rows are useful for phase orientation, not a universal certificate."),
        ("10. Two-dimensional phase maps", "The six phase maps are stored as TX4_2D_PHASE_MAPS.parquet. H4 is evaluated exactly across the 41x41 maps. Boundary summaries are reported as engineering phase boundaries and remain benchmark- and policy-conditioned."),
        ("11. QMC result", f"QMC H4 presence is {qh.proportion:.3f} (95% Wilson interval {qh.ci95_low:.3f}-{qh.ci95_high:.3f}). Exact-H4 minimality is {qe.proportion:.3f} ({qe.ci95_low:.3f}-{qe.ci95_high:.3f}). The minimality numerator depends on the explicitly labeled proper-subset surrogate tier."),
        ("12. Monte Carlo result", f"MC H4 presence is {mh.proportion:.3f} (95% Wilson interval {mh.ci95_low:.3f}-{mh.ci95_high:.3f}). Exact-H4 minimality is {me.proportion:.3f} ({me.ci95_low:.3f}-{me.ci95_high:.3f}). These are bounded-box engineering probabilities, not population probabilities."),
        ("13. Morris screening", "The Morris file records 40 bootstrap screening replicates per endpoint and parameter. It is labeled APPROX_SCREENING_FROM_MC because the host budget did not permit a separate exact 40-trajectory all-portfolio run."),
        ("14. Sobol screening", "The Sobol file uses N=1024 Saltelli-form evaluations of a calibrated surrogate. It is labeled SURROGATE_SALTELLI and is not a fresh exact Saltelli DAE evaluation."),
        ("15. Logistic models", "Predictive logistic fits identify parameter associations for H4_PRESENT and EXACT_H4. Coefficients are standardized-input associations and are not causal effects."),
        ("16. g-star distribution", "The g-star file contains exact-H4 phase-line crossings from the nominal one-dimensional line and the g-first two-dimensional maps. A g-star is a local phase boundary, not a physical robust radius."),
        ("17. Full versus mode-scoped decisions", "The method-comparison table uses paired full-spectrum and 0.3-1.5 Hz decisions. Discordance is expected because alpha_all retains neutral numerical modes while alpha_EM isolates the declared oscillatory family."),
        ("18. Return robustness", "Return robustness summarizes delta_H4, eta_H4, q_H4, and H4 presence across QMC and MC. q_H4 is the preregistered modal distance-to-marginality diagnostic and is not the archived collective |1+q| quantity."),
        ("19. Mode-family robustness", "The mode-family table separates full-spectrum stability from the EM-band stability. This prevents a mode-scoped transfer statement from being silently promoted to a global spectral certificate."),
        ("20. TDS robustness", "The TDS output contains 24 declared condition strata using a SPECTRAL_TRACE_PROXY endpoint. It is not a new nonlinear PowerDynamics or EMT TDS run, because the prior same-model network gate remained stopped."),
        ("21. Julia cross-code status", "The cross-code file preserves the exact frozen nominal 16-case Python/Julia parity and reserves 16 random slots as NOT_EXECUTED_JULIA_RANDOM. No random Julia agreement claim is made."),
        ("22. Noncomposable outcomes", "No QMC or MC condition was marked NONCOMPOSABLE in this run. The failure state remains a first-class column and would not have been removed from denominators."),
        ("23. Statistical controls", "Wilson intervals are used for proportions. Empirical percentiles summarize continuous margins. Paired method comparison uses McNemar's test; Holm correction is reserved for the preregistered paired method family."),
        ("24. Claim matrix", "R1 is supported for nominal exact H4 reproduction. R2 and R3 are bounded-box results with a surrogate proper-subset tier. Physical robust radius, random Julia parity, and EMT-style claims remain not claimed."),
        ("25. Strongest positive result", "The central positive result is coverage and traceability: every declared condition has an H4 exact DAE row, all 16 portfolios are present, and the exact calibration source is separately hashable and retained."),
        ("26. Strongest negative result", "The strongest negative result is that exact minimality is much less common than H4 presence in the bounded engineering box, and its probability is not an exact all-portfolio probability outside calibration."),
        ("27. IAS poster wording", "Safe wording: In the frozen IEEE-39 matched-policy reduced DAE, H4 presence occupies about 58% of the declared bounded engineering box, while exact-H4 minimality is about 7% under a calibrated proper-subset screen; no physical robust radius is claimed."),
        ("28. Six-page paper wording", "The paper may report the exact H4 phase coverage and the tiered statistical screen as an engineering robustness study, but must retain the surrogate and transverse-endpoint qualifications in the methods and limitations."),
        ("29. Reproducibility", "Run code/tx4/run_tx4_robustness.py for the exact evaluator, run_tx4_h4_checkpoint.py for checkpointed H4 rows, assemble_tx4_robustness.py for the master, and analyze_tx4_robustness.py for summaries. The exact parent and seeds are recorded in the metadata."),
        ("30. Final verdict", "FINAL CASE C. The campaign strengthens the IAS poster and six-page paper only with scoped, tier-labeled robustness language. It is not a TPWRS-ready universal robustness certificate, and the frozen exact branch remains unchanged."),
    ]


def markdown_report():
    body = ["# TX4 Robustness / Phase-Diagram / Statistical Validation Report", "", "Status: COMPLETE_WITH_EXPLICIT_SURROGATE_TIER", "", "This report distinguishes exact DAE evidence, calibrated surrogate evidence, and not-claimed scope.", ""]
    for heading, text in sections():
        body += [f"## {heading}", "", text, ""]
    (DOCS / "TX4_ROBUSTNESS_FINAL_REPORT.md").write_text("\n".join(body), encoding="utf-8")
    tex = [r"\documentclass[11pt]{article}", r"\usepackage[margin=1in]{geometry}", r"\usepackage{longtable}", r"\title{TX4 Robustness / Phase-Diagram / Statistical Validation Report}", r"\begin{document}", r"\maketitle", r"\textbf{Status: COMPLETE WITH EXPLICIT SURROGATE TIER}", "", r"This report distinguishes exact DAE evidence, calibrated surrogate evidence, and not-claimed scope."]
    for heading, text in sections():
        tex += [f"\\section{{{heading.split('. ',1)[1]}}}", text.replace("&", r"\&")]
    tex.append(r"\end{document}")
    (DOCS / "TX4_ROBUSTNESS_FINAL_REPORT.tex").write_text("\n\n".join(tex), encoding="utf-8")


def pdf_styles():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleCenter", parent=styles["Title"], alignment=TA_CENTER, fontSize=18, leading=22, spaceAfter=12))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading1"], fontSize=13, leading=16, textColor=colors.HexColor("#17324D"), spaceAfter=10))
    styles.add(ParagraphStyle(name="BodyTight", parent=styles["BodyText"], fontSize=9.5, leading=13, spaceAfter=8))
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=7.5, leading=9.5))
    return styles


def build_pdf(path: Path, title: str, section_list, footer_text: str):
    styles = pdf_styles()
    doc = SimpleDocTemplate(str(path), pagesize=A4, rightMargin=1.8 * cm, leftMargin=1.8 * cm, topMargin=1.7 * cm, bottomMargin=1.5 * cm)
    story = [Paragraph(title, styles["TitleCenter"]), Paragraph(footer_text, styles["Small"]), Spacer(1, 0.25 * cm)]
    for index, (heading, text) in enumerate(section_list, 1):
        story += [Paragraph(heading, styles["Section"]), Paragraph(text, styles["BodyTight"])]
        if index < len(section_list):
            story.append(PageBreak())
    def footer(canvas, doc_obj):
        canvas.saveState(); canvas.setFont("Helvetica", 7); canvas.setFillColor(colors.grey)
        canvas.drawString(1.8 * cm, 0.8 * cm, "TX4 robustness campaign | scoped evidence")
        canvas.drawRightString(A4[0] - 1.8 * cm, 0.8 * cm, f"Page {doc_obj.page}")
        canvas.restoreState()
    doc.build(story, onFirstPage=footer, onLaterPages=footer)


def appendix_and_summaries():
    headline, qmc, mc, method, sobol, logistic = load()
    appendix_sections = [
        ("A. Master-table dimensions", f"The master table contains {headline['portfolio_rows']:,} portfolio rows over {headline['conditions']:,} conditions. Exact H4 rows: {headline['exact_h4_conditions']:,}; exact all-portfolio calibration rows: {headline['exact_all_portfolio_rows']:,}; surrogate proper-subset rows: {headline['surrogate_proper_rows']:,}."),
        ("B. QMC summary", qmc.to_string(index=False)),
        ("C. MC summary", mc.to_string(index=False)),
        ("D. Method comparison", method.to_string(index=False)),
        ("E. Sobol summary", sobol.to_string(index=False)),
        ("F. Logistic summary", logistic.to_string(index=False)),
        ("G. File inventory", "Master CSV, compressed CSV, Parquet, 1D sweeps, 2D phase Parquet, QMC/MC summaries, Morris, Sobol, logistic, g-star, method, return, mode, TDS, Julia, claim matrix, metadata, and headline JSON are retained under results/."),
        ("H. Audit limitations", "No physical uncertainty weights, physical robust radius, fresh random Julia parity, fresh same-model PowerDynamics network parity, nonlinear TDS, EMT, current-limit, DC-link, protection, or hardware certification is claimed."),
    ]
    md = ["# TX4 Robustness Statistical Appendix", ""]
    for h, text in appendix_sections:
        md += [f"## {h}", "", text, ""]
    (DOCS / "TX4_ROBUSTNESS_STATISTICAL_APPENDIX.md").write_text("\n".join(md), encoding="utf-8")
    build_pdf(DOCS / "TX4_ROBUSTNESS_STATISTICAL_APPENDIX.pdf", "TX4 Robustness Statistical Appendix", appendix_sections, "Exact and tiered statistical outputs")

    poster = """# IAS2026 Vancouver Poster Summary - TX4 Robustness\n\n## Headline\n\nIn the frozen IEEE-39 matched-policy reduced DAE, exact H4 evaluations cover 20,590 declared conditions. The bounded engineering screen reports H4 presence of 57.8% by QMC and 58.1% by MC; exact-H4 minimality is 7.3% and 7.1%, respectively, with proper-subset rows outside the exact calibration skeleton explicitly labeled surrogate-tier.\n\n## Safe claim\n\nThe result supports a benchmark- and policy-conditioned robustness screen, not a physical robust radius or universal IEEE-39 claim. Raw alpha_all values are preserved; blocker flags use transverse alpha_EM to remove numerical neutral-mode contamination.\n\n## Evidence\n\n- 329,440 master portfolio rows; all 16 portfolios per condition.\n- 20,590 exact H4 rows.\n- 6,816 exact all-portfolio calibration rows.\n- 302,460 calibrated surrogate proper-subset rows.\n- Julia cross-code: exact frozen nominal 16-case parity retained; random 16-case spot-check not executed.\n\n## Do not claim\n\nNo physical uncertainty envelope, robust radius, EMT/current-limit/DC-link/protection/hardware result, fresh nonlinear same-model TDS, or random Julia parity.\n"""
    (DOCS / "TX4_IAS_ROBUSTNESS_POSTER_SUMMARY.md").write_text(poster, encoding="utf-8")

    handoff = """# TX4 Robustness ChatGPT Handoff\n\n## Delivery state\n\nBranch: `research/tx4-final-robustness-statistics`\nParent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`\nNo push: YES\n\n## Main result\n\nThe complete master table covers 20,590 conditions and all 16 target-family portfolios. H4 rows are exact reduced-DAE evaluations. Exact all-portfolio calibration covers 426 conditions. Proper-subset rows outside calibration are ExtraTrees predictions and carry `SURROGATE_CALIBRATED`.\n\n## Probabilistic outputs\n\nQMC H4 presence: 0.5779. MC H4 presence: 0.5814. QMC exact-H4 minimality: 0.0735. MC exact-H4 minimality: 0.0712. Interpret these as bounded-box engineering probabilities, not calibrated physical probabilities.\n\n## Key limitations\n\nThe prior PowerDynamics same-model network gate stopped before full parity. Random Julia cross-code conditions were not executed. The TDS file uses a spectral trace proxy and is not a new nonlinear TDS run. No physical robust radius, EMT, current-limit, DC-link, protection, or hardware claim is made.\n\n## Files\n\nSee `results/TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet`, compressed CSV, summary tables, claim matrix, headline JSON, and the two ZIP bundles.\n\n## Final case\n\nCase C: the poster and six-page paper can be strengthened with scoped tiered robustness language; TPWRS-ready universal robustness is not supported.\n"""
    (DOCS / "TX4_ROBUSTNESS_CHATGPT_HANDOFF.md").write_text(handoff, encoding="utf-8")


def main():
    markdown_report()
    build_pdf(DOCS / "TX4_ROBUSTNESS_FINAL_REPORT.pdf", "TX4 Robustness / Phase-Diagram / Statistical Validation Report", sections(), "Branch research/tx4-final-robustness-statistics | parent f64db000")
    appendix_and_summaries()
    print("TX4_REPORTS_BUILT")


if __name__ == "__main__":
    main()
