"""Build the TX4 final modal-scope audit PDF from committed artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
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

HERE = Path(__file__).resolve()
ROOT = HERE.parents[5]
PDF = ROOT / "docs" / "TX4_FINAL_MODAL_SCOPE_REPORT.pdf"
FIG = ROOT / "figures" / "final_audit"
RES = ROOT / "results"


def para(text: str, style: ParagraphStyle) -> Paragraph:
    return Paragraph(text.replace("&", "&amp;"), style)


def image(stem: str, width: float = 6.8 * inch) -> Image:
    path = FIG / f"{stem}.png"
    probe = Image(str(path))
    height = width * probe.imageHeight / probe.imageWidth
    img = Image(str(path), width=width, height=height)
    return img


def table(data, widths=None, font=7.4):
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#203040")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#b8c2cc")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), font),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef3f6")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return t


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#b8c2cc"))
    canvas.line(0.65 * inch, 0.52 * inch, 7.85 * inch, 0.52 * inch)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#53606d"))
    canvas.drawString(0.65 * inch, 0.34 * inch, "TX4 final modal-scope correction audit | no push")
    canvas.drawRightString(7.85 * inch, 0.34 * inch, f"page {doc.page}")
    canvas.restoreState()


def build():
    styles = getSampleStyleSheet()
    title = ParagraphStyle("Title", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=22, leading=27, alignment=TA_CENTER, textColor=colors.HexColor("#203040"), spaceAfter=14)
    subtitle = ParagraphStyle("Subtitle", parent=styles["Normal"], fontSize=10.5, leading=15, alignment=TA_CENTER, textColor=colors.HexColor("#53606d"), spaceAfter=18)
    h1 = ParagraphStyle("H1", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=15, leading=19, textColor=colors.HexColor("#203040"), spaceBefore=5, spaceAfter=7)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=11.5, leading=14, textColor=colors.HexColor("#277da1"), spaceBefore=4, spaceAfter=5)
    body = ParagraphStyle("Body", parent=styles["BodyText"], fontName="Helvetica", fontSize=9, leading=13, textColor=colors.HexColor("#20262b"), spaceAfter=7)
    small = ParagraphStyle("Small", parent=body, fontSize=7.8, leading=10.5)
    callout = ParagraphStyle("Callout", parent=body, fontName="Helvetica-Bold", fontSize=11, leading=15, backColor=colors.HexColor("#eef3f6"), borderColor=colors.HexColor("#277da1"), borderWidth=0.7, borderPadding=8, spaceBefore=5, spaceAfter=10)

    story = []
    story += [para("TX4 final modal-scope correction report", title), para("Finite audit of the frozen IEEE-39 phasor-domain blind-prediction campaign", subtitle)]
    story += [para("FINAL CASE: A", callout)]
    story += [para("Executive verdict", h1)]
    story += [para("This correction audit fixes the closure sign convention, replaces invalid normalized local-factor evidence with physical I+M_ii factors, and separates complete-spectrum V9 truth from the targeted 0.3-1.5 Hz task. The V4 collective mechanism remains supported at P4. The targeted EM task is nearly exact, but the global V9 task remains limited and the archived FC10/default table is not same-policy evidence.", body)]
    story += [para("Provenance", h2), para("Branch: research/tx4-final-modal-scope-audit. Starting commit: 2cbb860eb4ba053edee99c1a040c8024f3f6828e. Parent freeze: 69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6. No push. No radius, repair, planner, weak-node, weak-link, IEEE-68, EMT, or new-controller campaign was run.", body)]
    story += [para("Key headline", h2), para("Q_H approaches -1, contextual return approaches +1, and physical local factors remain regular while the collective factor closes. The IAS poster is ready with corrected scope labels. The six-page paper and TPWRS submission are not ready.", body)]
    story += [PageBreak()]

    story += [para("1. Sign correction and physical local factors", h1)]
    story += [para("The exact identity is det(I+Q_H)=det(I+Q_RR) det(I-R_i|R). Therefore the collective closure condition is -1 in sigma(Q_H), while the contextual return condition is +1 in sigma(R_i|R). A previous handoff incorrectly gave both objects the same sign.", body)]
    sign = pd.read_csv(RES / "TX4_Q_VS_RETURN_SIGN_AUDIT.csv")
    sign_rows = [["bus", "abs(Q+1)", "abs(return-1)", "Schur residual"]]
    for _, r in sign.iterrows():
        sign_rows.append([str(int(r.device_i)), f"{r['abs(q_eigenvalue_plus_one)']:.3e}", f"{r['abs(return_eigenvalue_minus_one)']:.3e}", f"{r.schur_residual:.3e}"])
    story += [table(sign_rows, [0.65*inch, 1.45*inch, 1.55*inch, 1.45*inch]), Spacer(1, 8)]
    story += [PageBreak(), para("Sign trajectories", h2), image("F_Q_MINUS1_RETURN_PLUS1", 6.7*inch)]
    story += [para("The physical audit uses M=D K before normalization. The minimum physical local singular value over the frozen gain sweep is 0.2973; the collective singular value reaches 1.26e-8. The minimum proper-subset collective factor at the boundary frequency is 0.2945.", body)]
    story += [PageBreak()]

    story += [para("2. V4 finite blocker and collective interpretation", h1)]
    story += [para("The V4 census is unchanged: 15 proper portfolios are stable and H4={30,33,35,37} is unstable at P4. H4 has alpha=+0.1270065 s^-1 at 0.6222797 Hz. Thus H4 is the sole true blocker in this finite V4 census, not a universal blocker for arbitrary systems.", body)]
    story += [PageBreak(), para("Physical local versus collective singular values", h2), image("F_TRUE_LOCAL_VS_COLLECTIVE", 6.4*inch)]
    story += [PageBreak(), para("V4 Hasse lattice", h2), image("F_V4_TRUE_LATTICE", 6.4*inch)]
    story += [PageBreak()]

    story += [para("3. V9 modal-scope truth", h1)]
    story += [para("All 512 V9 portfolios were reclassified under the same frozen P4 policy as the predictor. The complete transverse spectrum and the 0.3-1.5 Hz band are reported separately.", body)]
    v9_rows = [
        ["task", "n", "stable", "unstable", "correct", "false-safe", "false-unstable"],
        ["global spectrum", "512", "332", "180", "402/512", "110", "0"],
        ["0.3-1.5 Hz band", "512", "442", "70", "511/512", "1", "0"],
    ]
    story += [table(v9_rows, [1.45*inch, .45*inch, .65*inch, .7*inch, .8*inch, .8*inch, .9*inch]), Spacer(1, 8)]
    story += [para("The archived 395/512 and 116 false-safe values came from an FC10/default policy path and are retained only as historical frozen data. They are not same-policy P4 evidence. The corrected same-policy global comparison has 110 false-safe cases and 0 false-unstable cases.", body)]
    taxonomy = pd.read_csv(RES / "TX4_V9_FALSE_SAFE_TAXONOMY.csv")
    story += [para("False-safe taxonomy", h2), table([["global false-safe class", "count"], ["APERIODIC_REAL", str(len(taxonomy))], ["SLOW_OUTSIDE_BAND", "0"], ["TARGET_EM_BAND", "0"], ["FAST_OSCILLATORY", "0"]], [2.5*inch, 1*inch])]
    story += [para("The one targeted EM false-safe is a separate set: 30+31+32+33+34+35+36+37 has alpha_EM=+0.04835 s^-1 at 0.7715 Hz, while its global critical mode is aperiodic. There are no eligible slow false-safe cases, so optional slow-case TDS was skipped.", body)]
    story += [PageBreak()]

    story += [para("4. P7 scope and nonlinear traces", h1)]
    story += [para("P7 is SUPPLEMENT ONLY. The first member, 30+31+32+34+36+38, has global alpha=+422.8355 s^-1 at zero frequency but alpha_EM=-0.5313 s^-1 at 1.3774 Hz. The comparison member is stable at alpha=-0.1511 s^-1 and 1.3203 Hz. This is not a clean test of the validated 0.62 Hz mechanism.", body)]
    story += [PageBreak(), para("Equal-horizon TDS figure", h2), image("F_TDS_FINAL_EQUAL_HORIZON", 6.4*inch)]
    story += [para("The equal-horizon TDS panel reuses the three existing phasor-domain traces over 0-30 s. It is not EMT validation and no new slow-case TDS was run.", body)]
    story += [PageBreak()]

    story += [para("5. Reviewer attack and claim hierarchy", h1)]
    story += [para("Reviewer 1: novelty 5/10, correctness 8/10, model adequacy 5/10, evidence 8/10, usefulness 7/10, IAS impact 7/10, TPWRS readiness 3/10. The algebra is classical; the useful contribution is the finite blind mechanism packaging with explicit scope limits.", body)]
    story += [para("Reviewer 2: novelty 5/10, correctness 8/10, model adequacy 4/10, evidence 7/10, usefulness 6/10, IAS impact 6/10, TPWRS readiness 2/10. The mode is interpretable for the frozen model, but credibility is limited by one custom phasor-domain GFL, no EMT, and no second-model check.", body)]
    story += [para("Allowed headline", h2), para("In the frozen TX4 IEEE-39 phasor-domain model at P4, the V4 H4 boundary is a collective closure loss: physical local factors remain regular while Q_H reaches -1 and its contextual return reaches +1. The targeted 0.3-1.5 Hz V9 task is nearly exact, but the complete-spectrum V9 task contains 110 same-policy aperiodic false-safe cases.", callout)]
    story += [para("Forbidden: universal local regularity, universal V9 safety, probabilistic generalization, a new theorem, EMT validation, six-page-paper readiness, or TPWRS readiness.", body)]
    story += [PageBreak()]

    story += [para("6. Final decision", h1)]
    story += [para("IAS poster redesign: YES. Use the corrected signs, physical local-factor figure, V4 lattice, and explicit global-versus-EM scope. Remove the old P7 main panel. Keep P7 in the supplement with its aperiodic/global qualification.", body)]
    story += [para("Six-page paper ready: NO. TPWRS ready: NO. Recommendation: NARROW. Proceed with the corrected IAS poster and external ChatGPT review of the ZIP. Defer broader journal claims until a second documented dynamic model, more operating points, controller-limit validation, and independent nonlinear validation are available.", body)]
    story += [para("Final answers", h2)]
    answers = [
        ["question", "answer"],
        ["true V4 blocker", "H0={H4} in the finite V4 census"],
        ["Q / return signs", "Q_H -> -1; return -> +1"],
        ["V9 global", "402/512 correct; 110 false-safe, all aperiodic"],
        ["V9 EM band", "511/512 correct; 1 false-safe, 0 false-unstable"],
        ["P7", "SUPPLEMENT ONLY"],
        ["TDS", "equal-horizon replot only; no slow cases eligible"],
        ["poster / paper / TPWRS", "YES / NO / NO"],
    ]
    story += [table(answers, [1.75*inch, 4.5*inch], font=7.8)]
    doc = SimpleDocTemplate(str(PDF), pagesize=letter, rightMargin=.65*inch, leftMargin=.65*inch, topMargin=.62*inch, bottomMargin=.7*inch, title="TX4 final modal-scope correction report", author="Codex")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(PDF)


if __name__ == "__main__":
    build()
