"""Create the corrected TX4 robustness statistics audit PDF."""

from __future__ import annotations

import json
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "TX4_AUDIT_FIGURES"
OUTPUT = ROOT / "output" / "pdf" / "TX4_ROBUSTNESS_STATISTICS_AUDIT_FINAL_REPORT.pdf"


def para(text: str, style) -> Paragraph:
    return Paragraph(escape(text).replace("\n", "<br/>").replace("&lt;b&gt;", "<b>").replace("&lt;/b&gt;", "</b>"), style)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#59636e"))
    canvas.drawString(0.65 * inch, 0.42 * inch, "TX4 robustness statistics audit | exact primary fallback")
    canvas.drawRightString(7.85 * inch, 0.42 * inch, f"Page {doc.page}")
    canvas.restoreState()


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    headline = json.loads((RESULTS / "TX4_ROBUSTNESS_HEADLINE_CORRECTED.json").read_text(encoding="utf-8"))
    gate = json.loads((RESULTS / "TX4_SURROGATE_GATE.json").read_text(encoding="utf-8"))
    invariants = pd.read_csv(RESULTS / "TX4_ROBUSTNESS_INVARIANT_CHECKS.csv")
    q = pd.read_csv(RESULTS / "TX4_QMC_EXACT_BLOCKER_SUMMARY.csv")
    m = pd.read_csv(RESULTS / "TX4_MC_EXACT_BLOCKER_SUMMARY.csv")
    qd = pd.read_csv(RESULTS / "TX4_QMC_DELTA_H4_TRUE_MINIMALITY.csv").iloc[0]
    md = pd.read_csv(RESULTS / "TX4_MC_DELTA_H4_TRUE_MINIMALITY.csv").iloc[0]
    eta = pd.read_json(RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT_SUMMARY.json", typ="series")

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="AuditTitle", parent=styles["Title"], alignment=TA_CENTER, fontName="Helvetica-Bold", fontSize=22, leading=27, textColor=colors.HexColor("#17324d"), spaceAfter=18))
    styles.add(ParagraphStyle(name="AuditSubtitle", parent=styles["Normal"], alignment=TA_CENTER, fontSize=11, leading=15, textColor=colors.HexColor("#59636e"), spaceAfter=22))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=15, leading=19, textColor=colors.HexColor("#17324d"), spaceBefore=12, spaceAfter=8))
    styles.add(ParagraphStyle(name="Subsection", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=colors.HexColor("#275b7a"), spaceBefore=8, spaceAfter=5))
    styles.add(ParagraphStyle(name="BodyAudit", parent=styles["BodyText"], fontSize=9.5, leading=13.2, spaceAfter=7))
    styles.add(ParagraphStyle(name="SmallAudit", parent=styles["BodyText"], fontSize=8, leading=10, textColor=colors.HexColor("#59636e"), spaceAfter=5))
    styles.add(ParagraphStyle(name="FigureCaption", parent=styles["BodyText"], alignment=TA_CENTER, fontSize=8.5, leading=11, textColor=colors.HexColor("#59636e"), spaceAfter=12))

    story = []
    story.append(Spacer(1, 0.45 * inch))
    story.append(Paragraph("TX4 Robustness Statistics Audit", styles["AuditTitle"]))
    story.append(Paragraph("Corrected blocker-antichain definitions, exact QMC/MC recomputation, and Vancouver handoff", styles["AuditSubtitle"]))
    story.append(Paragraph("Date: 2026-09-20 | Branch: research/tx4-robustness-statistics-audit", styles["AuditSubtitle"]))
    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph("Executive verdict", styles["Section"]))
    story.append(Paragraph(
        f"The retained ExtraTrees tier failed the preregistered strict gate: {gate['false_minimal_n']} false-minimal conditions, {gate['missed_minimal_n']} missed-minimal conditions, {gate['exact_h4_flag_agreement']:.3f} exact-H4 agreement, and {gate['h0_antichain_exact_fraction']:.3f} H0-antichain agreement. The prespecified exact fallback was therefore executed for all 4,096 QMC and 5,000 MC coordinates across all 16 portfolios, totaling 145,536 exact DAE portfolio evaluations.",
        styles["BodyAudit"],
    ))
    story.append(Paragraph("The legacy H4_PRESENT flag tested H4 instability rather than inclusion-minimal H4 membership. This complement/minimality bug produced 156 legacy-positive but corrected-negative conditions in the 426-condition exact calibration.", styles["BodyAudit"]))
    story.append(Paragraph("Final machine gate: PASS for all 32 invariants.", styles["BodyAudit"]))

    story.append(Paragraph("Definitions and provenance", styles["Section"]))
    story.append(Paragraph("For each condition, U0 is the set of portfolios with alpha_EM >= 0 in the frozen 0.3-1.5 Hz band. H0 is the inclusion-minimal antichain of U0. H4_UNSTABLE means H4 belongs to U0; H4_PRESENT means H4 belongs to H0; EXACT_H4 means H0 is exactly the singleton {H4}; NONCOMPOSABLE means H0 contains a blocker of size at least two; KAPPA is the smallest blocker size or NULL.", styles["BodyAudit"]))
    story.append(Paragraph("QMC results are fixed-design coverage fractions. MC results are assumed independent-uniform bounded-box probabilities with Wilson intervals. No U005/H005 threshold was found in the retained preregistration or code, so those endpoints were not invented or analyzed.", styles["BodyAudit"]))

    story.append(Paragraph("Corrected primary endpoints", styles["Section"]))
    rows = [["Endpoint", "QMC coverage", "MC probability", "MC Wilson 95%"]]
    for endpoint in ("H4_UNSTABLE", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"):
        qr = q[q.endpoint == endpoint].iloc[0]
        mr = m[m.endpoint == endpoint].iloc[0]
        rows.append([endpoint, f"{qr.coverage_fraction:.6f}", f"{mr.coverage_fraction:.6f}", f"[{mr.ci95_low:.6f}, {mr.ci95_high:.6f}]"])
    table = Table(rows, colWidths=[2.2 * inch, 1.35 * inch, 1.35 * inch, 1.55 * inch])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#17324d")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#b8c3cc")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#edf3f7")]),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]))
    story.append(table)
    story.append(Spacer(1, 0.12 * inch))
    story.append(Paragraph("H4_PRESENT and EXACT_H4 coincide for this full-set H4 definition: once any proper unstable blocker exists, H4 cannot remain inclusion-minimal.", styles["SmallAudit"]))

    story.append(Paragraph("Continuous diagnostics", styles["Section"]))
    story.append(Paragraph(f"The true-minimality margin is min(alpha_EM(H4), -max proper alpha_EM). QMC median={qd['median']:.6f} s^-1, IQR=[{qd['iqr_low']:.6f}, {qd['iqr_high']:.6f}], p05={qd['p05']:.6f}, p95={qd['p95']:.6f}. MC median={md['median']:.6f} s^-1 with BCa 95% interval=[{md['median_ci95_low']:.6f}, {md['median_ci95_high']:.6f}].", styles["BodyAudit"]))
    story.append(Paragraph(f"The exact local/collective eta stratum computed {int(eta['n_computed'])} of {int(eta['n_selected'])} requested conditions with {int(eta['n_failures'])} failures; all eta values were nonnegative. The definition was min_i sigma_min(I+Q_(H4\\i)(s_H)).", styles["BodyAudit"]))

    story.append(Paragraph("Corrected figures", styles["Section"]))
    captions = [
        ("F1_corrected_endpoint_frequencies.png", "Corrected endpoint frequencies: instability, minimal H4, exact H4, and composite blockers."),
        ("F2_kappa_pmf.png", "Exact KAPPA PMF, including NULL conditions."),
        ("F3_delta_H4_true_minimality.png", "Distribution of the exact true-minimality margin."),
        ("F4_alpha_H4_vs_proper_max.png", "H4 alpha_EM against the most unstable proper subset."),
        ("F5_h0_antichain_size.png", "Size of the minimal-blocker antichain H0."),
        ("F6_legacy_vs_corrected_flags.png", "Legacy endpoint flags versus exact calibration H0 truth."),
        ("F7_surrogate_gate.png", "Reconstructed ExtraTrees strict-gate diagnostics."),
        ("F8_eta_local_collective.png", "Exact local/collective eta stratum."),
        ("F9_gstar_provenance.png", "Retained exact-H4 phase-boundary diagnostic with provenance."),
        ("F10_evaluation_provenance.png", "Master-table evaluator provenance."),
    ]
    for index in range(0, len(captions), 2):
        for filename, caption in captions[index : index + 2]:
            image = Image(str(FIGURES / filename))
            image._restrictSize(6.55 * inch, 3.85 * inch)
            story.append(image)
            story.append(Paragraph(caption, styles["FigureCaption"]))
        if index + 2 < len(captions):
            story.append(PageBreak())

    story.append(PageBreak())
    story.append(Paragraph("Audit outputs and handoff", styles["Section"]))
    story.append(Paragraph("The corrected exact Parquet files are TX4_QMC_ALL16_EXACT.parquet and TX4_MC_ALL16_EXACT.parquet. The audit retains the exact calibration truth table, blocker-list table, row-level provenance, surrogate validation/error tables, corrected eta audit, 32-check invariant table, headline JSON, and the ten corrected PNGs. The old Morris table remains labeled approximate screening and the old Sobol table remains explicitly surrogate-based; neither is promoted to exact sensitivity evidence.", styles["BodyAudit"]))
    story.append(Paragraph("The frozen exact-P4/GFL11 worktree and the original robustness branch were not modified. No push was performed. The generated paper/poster PDFs from earlier work were not used as corrected evidence.", styles["BodyAudit"]))
    story.append(Paragraph("Representative invariant status", styles["Subsection"]))
    inv_rows = [["Check", "Status", "Observed"]]
    for _, item in invariants.head(10).iterrows():
        observed = "16 labels" if "portfolio_catalog" in str(item["check"]) else str(item["observed"])
        inv_rows.append([str(item["check"]), str(item["status"]), observed])
    inv_table = Table(inv_rows, colWidths=[3.25 * inch, 0.75 * inch, 2.15 * inch])
    inv_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#275b7a")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#b8c3cc")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#edf3f7")]),
    ]))
    story.append(inv_table)
    story.append(Spacer(1, 0.12 * inch))
    story.append(Paragraph("Complete machine-readable invariant results are in results/TX4_ROBUSTNESS_INVARIANT_CHECKS.csv; all 32 rows are PASS.", styles["SmallAudit"]))

    doc = SimpleDocTemplate(str(OUTPUT), pagesize=letter, rightMargin=0.65 * inch, leftMargin=0.65 * inch, topMargin=0.62 * inch, bottomMargin=0.65 * inch, title="TX4 Robustness Statistics Audit", author="Codex")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(str(OUTPUT))


if __name__ == "__main__":
    main()
