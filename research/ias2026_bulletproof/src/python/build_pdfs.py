"""Build the three closure PDFs with reportlab."""
from __future__ import annotations

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

HERE = Path(__file__).resolve()
CAMPAIGN = HERE.parents[2]
OUT = CAMPAIGN / "paper"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = colors.HexColor("#17365D")
GREEN = colors.HexColor("#1F7A5C")
RED = colors.HexColor("#B42318")
GRAY = colors.HexColor("#5B6573")

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="CampaignTitle", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=21, leading=25, textColor=BLUE, alignment=TA_CENTER, spaceAfter=12))
styles.add(ParagraphStyle(name="Subtitle", parent=styles["Normal"], fontName="Helvetica", fontSize=10, leading=14, textColor=GRAY, alignment=TA_CENTER, spaceAfter=16))
styles.add(ParagraphStyle(name="H1C", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=14, leading=18, textColor=BLUE, spaceBefore=10, spaceAfter=6))
styles.add(ParagraphStyle(name="H2C", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=GREEN, spaceBefore=8, spaceAfter=4))
styles.add(ParagraphStyle(name="BodyC", parent=styles["BodyText"], fontName="Helvetica", fontSize=9.2, leading=13, spaceAfter=5))
styles.add(ParagraphStyle(name="SmallC", parent=styles["BodyText"], fontName="Helvetica", fontSize=7.6, leading=10, textColor=GRAY, spaceAfter=3))
styles.add(ParagraphStyle(name="BoxC", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=10, leading=14, textColor=BLUE, backColor=colors.HexColor("#EAF2F8"), borderColor=colors.HexColor("#B8CDE3"), borderWidth=0.5, borderPadding=8, spaceAfter=8))


def p(text: str, style: str = "BodyC") -> Paragraph:
    return Paragraph(text.replace("&", "&amp;"), styles[style])


def bullet(text: str) -> Paragraph:
    return Paragraph("&#8226; " + text.replace("&", "&amp;"), styles["BodyC"])


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D5DCE5"))
    canvas.line(doc.leftMargin, 0.48 * inch, A4[0] - doc.rightMargin, 0.48 * inch)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(GRAY)
    canvas.drawString(doc.leftMargin, 0.30 * inch, "IAS2026 Bulletproof Closure | 2026-09-18")
    canvas.drawRightString(A4[0] - doc.rightMargin, 0.30 * inch, f"Page {doc.page}")
    canvas.restoreState()


def title(story, title: str, subtitle: str):
    story.extend([p(title, "CampaignTitle"), p(subtitle, "Subtitle")])


def summary_table(rows, widths=(1.75 * inch, 4.8 * inch)):
    data = [[p(str(a), "SmallC"), p(str(b), "SmallC")] for a, b in rows]
    table = Table(data, colWidths=list(widths), hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#F2F5F8")),
        ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#B8CDE3")),
        ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D5DCE5")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6), ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return table


def build_executive():
    path = OUT / "FINAL_EXECUTIVE_REPORT.pdf"
    doc = SimpleDocTemplate(str(path), pagesize=A4, rightMargin=0.62 * inch, leftMargin=0.62 * inch, topMargin=0.58 * inch, bottomMargin=0.62 * inch, title="IAS2026 Bulletproof Closure Executive Report")
    s = []
    title(s, "IAS2026 Bulletproof Closure", "Final executive report | branch research/ias2026-bulletproof-closure-v1")
    s.append(p("<b>Disposition:</b> the frozen IEEE-39 nominal closure survives, the stronger historical V9 transfer expectation is refuted, and robustness, second-model, same-mechanism PowerDynamics, and new blind-holdout claims are explicitly untested or stopped by scope.", "BoxC"))
    s.append(p("This report distinguishes document-specified expectations from measured results. No expected value was retuned after observation, and the user's pre-existing dirty checkout was preserved in a separate linked worktree.", "BodyC"))
    s.append(p("1. Results that survived", "H1C"))
    for text in [
        "Gate 0 reproduces the frozen P4 H4 nominal alpha: 0.1270064666836382 s^-1 versus 0.1270065 s^-1.",
        "The frozen P4 structure reports 15 stable proper subsets.",
        "Four randomized theorem/property-test modules pass, including the port bridge and transverse construction.",
        "The frozen boundary and nonlinear audit reproduce retrospectively: g*=0.2076814044, f=0.7064247879 Hz, and 32/32 declared TDS agreements.",
        "The documented governor changes P4 H4 alpha to -0.07454098428813066 s^-1.",
        "PowerDynamics 5.0.0 reaches its official IEEE-39 tutorial equilibrium: 39 buses and 46 branches.",
    ]:
        s.append(bullet(text))
    s.append(Spacer(1, 4))
    fig = CAMPAIGN / "figures" / "png" / "F1_p4_policy_dependence.png"
    if fig.exists():
        s.append(Image(str(fig), width=5.8 * inch, height=3.3 * inch))
        s.append(p("Figure 1. Frozen P4 H4 policy dependence. The governor result is policy-conditioned.", "SmallC"))
    s.append(p("2. Negative and stopped results", "H1C"))
    for text in [
        "V9: 395/512 correct, 116 false-safe, and 1 false-unstable. The plan's 511/512 and zero-error expectation is REFUTED.",
        "The exact V9 blocker antichain is not supported.",
        "PowerDynamics does not reproduce the frozen GFL mechanism because the available tutorial is a different synchronous-machine/governor model.",
        "No second converter model, new blind topology/model holdout, physical robust radius, common uncertainty envelope, new Julia TDS, EMT, current-limit, DC-link, or hardware certification is claimed.",
    ]:
        s.append(bullet(text))
    s.append(p("3. Safe claim", "H1C"))
    s.append(p("In the frozen IEEE-39 benchmark and stated policy, the spectral/collective closure audit reproduces a P4 minimal incompatibility and the declared frozen nonlinear verdicts. The historical V9 transfer expectation is refuted, the documented governor changes the P4 H4 sign, and robustness, second-model, same-mechanism PowerDynamics, and new blind-holdout claims remain untested or stopped by scope.", "BoxC"))
    s.append(p("4. Reproducibility", "H1C"))
    s.append(summary_table([
        ("Gate 0", "reports/GATE0_REPORT.md"),
        ("Theory tests", "logs/theory_tests.log; 4 modules passed"),
        ("PowerDynamics", "raw/powerdynamics/pd39_equilibrium_gate.md"),
        ("Claim ledger", "docs/CLAIM_LEDGER.md"),
        ("Dashboard", "reports/GATE_DASHBOARD.md and .csv"),
        ("Base", "IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE at b9f274e..."),
    ]))
    doc.build(s, onFirstPage=header_footer, onLaterPages=header_footer)


def build_paper():
    path = OUT / "FINAL_PAPER.pdf"
    doc = SimpleDocTemplate(str(path), pagesize=A4, rightMargin=0.65 * inch, leftMargin=0.65 * inch, topMargin=0.62 * inch, bottomMargin=0.62 * inch, title="IAS2026 Closure Paper Draft")
    s = []
    title(s, "Spectral Graph Power Dynamics: A Scoped Closure Audit", "IAS2026 final paper draft | negative evidence retained")
    s.append(p("Abstract", "H1C"))
    s.append(p("We audit a frozen IEEE-39 SG-to-GFL benchmark using spectral and collective closure identities, transverse structure checks, and archived nonlinear evidence. The nominal P4 H4 value is reproduced at 0.1270064667 s^-1, the frozen P4 structure has 15 stable proper subsets, and four randomized algebraic test modules pass. The archived nonlinear table reports 32/32 declared TDS agreements. However, the historical V9 transfer result is 395/512 rather than the stronger expected 511/512, with 116 false-safe and 1 false-unstable cases. An isolated PowerDynamics 5.0.0 run validates only its official IEEE-39 tutorial equilibrium, not parity with the frozen GFL model. The conclusion is therefore benchmark- and policy-conditioned, with robustness and second-model generalization left untested.", "BodyC"))
    s.append(p("1. Scope and protocol", "H1C"))
    s.append(p("The execution layer was created from the IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE tag in a linked worktree. Historical source tables were copied read-only into raw/gate0. Expectations were treated as hypotheses: a mismatch is reported as REFUTED, and an unavailable dependent experiment is reported as NOT_TESTED or STOPPED_BY_GATE.", "BodyC"))
    s.append(summary_table([
        ("Primary benchmark", "Frozen IEEE-39 SG-to-GFL closure model"),
        ("Nominal gate", "P4 H4 alpha absolute error <= 1e-6 s^-1"),
        ("Theory tests", "local/collective factorization, bridge, contextual return, transverse construction"),
        ("Independent software", "PowerDynamics 5.0.0 official IEEE-39 tutorial equilibrium"),
        ("Labels", "IEEE39_VALIDATED, THEOREM_VALIDATED, POWERDYNAMICS_VALIDATED, RETROSPECTIVE, REFUTED, NOT_TESTED, STOPPED_BY_GATE"),
    ]))
    s.append(p("2. Frozen IEEE-39 closure", "H1C"))
    s.append(p("Gate 0 reproduces the nominal and governed P4 values from the immutable final-closure tables. The nominal H4 alpha is 0.1270064666836382 s^-1; the documented governor value is -0.07454098428813066 s^-1. The sign change is evidence of policy dependence and a local repair effect, not a universal certificate. The F7A boundary witness is g*=0.2076814044 at f=0.7064247879 Hz, with boundary Q sign -1 and contextual-return sign +1.", "BodyC"))
    s.append(p("3. Algebra and nonlinear evidence", "H1C"))
    s.append(p("The property suite uses randomized block and transverse constructions. It tests the port identities, determinant factorization, the bridge without inversion of the disturbance matrix, the contextual Schur identity under Q_ii=0, and the transverse rotation/drift residuals. All four modules pass. The frozen G2 table reports 32 declared TDS verdict agreements and two P4 H4 unstable rows.", "BodyC"))
    s.append(p("4. Transfer and independent software", "H1C"))
    s.append(p("The archived V9 reveal contains 512 portfolios: 395 correct, 116 false-safe, and 1 false-unstable. This refutes the stronger near-perfect transfer expectation and exact antichain claim. PowerDynamics reaches its official IEEE-39 tutorial equilibrium, but that tutorial uses a different model class. Same-model parity and same-mechanism reproduction are stopped by gate.", "BodyC"))
    fig2 = CAMPAIGN / "figures" / "png" / "F2_v9_transfer_refutation.png"
    if fig2.exists():
        s.append(Image(str(fig2), width=5.8 * inch, height=3.3 * inch))
        s.append(p("Figure 2. Historical V9 transfer expectation versus revealed result.", "SmallC"))
    s.append(p("5. Limitations and conclusion", "H1C"))
    s.append(p("No physical common uncertainty envelope, robust LFT radius, second converter model, new blind topology/model holdout, new Julia TDS, EMT, current-limit, DC-link, protection, switching-path, or hardware certification is claimed. The safe conclusion is a benchmark-specific, policy-conditioned minimal incompatibility with exact frozen nominal closure and negative evidence against stronger transfer/generalization language.", "BodyC"))
    doc.build(s, onFirstPage=header_footer, onLaterPages=header_footer)


def build_poster():
    path = OUT / "FINAL_POSTER_DRAFT.pdf"
    page = landscape(A4)
    doc = SimpleDocTemplate(str(path), pagesize=page, rightMargin=0.42 * inch, leftMargin=0.42 * inch, topMargin=0.38 * inch, bottomMargin=0.48 * inch, title="IAS2026 Poster Draft")
    s = []
    title(s, "Spectral Graph Power Dynamics", "IAS2026 closure poster draft | scoped claim with negative transfer evidence")
    columns = [
        [p("QUESTION", "H2C"), p("Does the frozen IEEE-39 benchmark admit a compact spectral/collective closure explanation of the P4 blocker, and does it transfer to the archived V9 portfolio set?", "BodyC"), p("METHOD", "H2C"), bullet("Frozen-artifact Gate 0 reproduction"), bullet("Randomized theorem/property tests"), bullet("Retrospective F7A and G2 nonlinear audit"), bullet("Independent PowerDynamics tutorial equilibrium gate"), p("LABELS", "H2C"), p("IEEE39_VALIDATED<br/>THEOREM_VALIDATED<br/>NONLINEAR_TDS_VALIDATED<br/>POWERDYNAMICS_VALIDATED", "BodyC")],
        [p("SURVIVED", "H2C"), p("P4 H4 nominal alpha", "BodyC"), p("<b>0.1270064667 s^-1</b>", "BoxC"), p("P4 proper stable subsets: <b>15</b>", "BodyC"), p("Boundary witness: g*=<b>0.2076814</b>, f=<b>0.7064248 Hz</b>", "BodyC"), p("Frozen TDS agreement: <b>32/32</b>", "BodyC"), p("Governor policy: 0.1270065 -> <b>-0.074541 s^-1</b>", "BodyC"), p("THEOREM CHECK", "H2C"), p("Four randomized test modules pass, with the contextual-return assumption Q_ii=0 explicit.", "BodyC")],
        [p("TRANSFER RESULT", "H2C"), p("Historical V9 reveal", "BodyC"), p("<b>395/512 correct</b>", "BoxC"), p("116 false-safe; 1 false-unstable", "BodyC"), p("The stronger 511/512, zero-error expectation is <b>REFUTED</b>.", "BodyC"), p("DO NOT OVERCLAIM", "H2C"), bullet("No universal IEEE-39 claim"), bullet("No physical robust radius"), bullet("No second converter model"), bullet("No new blind holdout"), bullet("No same-mechanism PowerDynamics parity"), bullet("No EMT/current-limit/DC-link certification")],
    ]
    table = Table([[columns[0], columns[1], columns[2]]], colWidths=[3.35 * inch, 3.35 * inch, 3.35 * inch], hAlign="CENTER")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor("#B8CDE3")),
        ("INNERGRID", (0, 0), (-1, -1), 0.6, colors.HexColor("#B8CDE3")),
        ("BACKGROUND", (0, 0), (-1, -1), colors.white),
        ("LEFTPADDING", (0, 0), (-1, -1), 10), ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("TOPPADDING", (0, 0), (-1, -1), 9), ("BOTTOMPADDING", (0, 0), (-1, -1), 9),
    ]))
    s.append(table)
    s.append(Spacer(1, 8))
    s.append(p("SAFE POSTER CLAIM: In the frozen IEEE-39 benchmark and stated policy, the spectral/collective closure audit reproduces a P4 minimal incompatibility and the declared frozen nonlinear verdicts. The historical V9 transfer expectation is refuted, the documented governor changes the P4 H4 sign, and robustness, second-model, same-mechanism PowerDynamics, and new blind-holdout claims remain untested or stopped by scope.", "BoxC"))
    doc.build(s, onFirstPage=header_footer, onLaterPages=header_footer)


if __name__ == "__main__":
    build_executive()
    build_paper()
    build_poster()
    print("built FINAL_EXECUTIVE_REPORT.pdf, FINAL_PAPER.pdf, FINAL_POSTER_DRAFT.pdf")
