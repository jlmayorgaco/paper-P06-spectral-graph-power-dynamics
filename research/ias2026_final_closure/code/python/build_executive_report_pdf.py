"""Create the executive audit PDF for the final closure root."""

from __future__ import annotations

import csv
import html

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
OUTPUT = ROOT / "reports" / "FINAL_EXECUTIVE_REPORT.pdf"


def paragraph(text: str, style: ParagraphStyle) -> Paragraph:
    return Paragraph(html.escape(text).replace("`", ""), style)


def read_sections() -> list[tuple[str, list[str]]]:
    lines = (ROOT / "reports" / "FINAL_EXECUTIVE_REPORT.md").read_text(encoding="utf-8").splitlines()
    sections: list[tuple[str, list[str]]] = []
    current = "Executive decision"
    body: list[str] = []
    for line in lines:
        if line.startswith("# "):
            continue
        if line.startswith("## "):
            if body:
                sections.append((current, body))
            current = line[3:].strip()
            body = []
        elif line.startswith("|") or line.startswith("---") or not line.strip():
            continue
        else:
            body.append(line[2:].strip() if line.startswith("- ") else line.strip())
    if body:
        sections.append((current, body))
    return sections


def dashboard_rows() -> list[list[str]]:
    with (ROOT / "reports" / "GATE_DASHBOARD.csv").open(newline="", encoding="utf-8") as handle:
        return [[row["gate"], row["status"], row["question"], row["observed"]] for row in csv.DictReader(handle)]


def decorate(canvas, document) -> None:
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#52606D"))
    canvas.drawString(document.leftMargin, 0.35 * inch, "IAS 2026 Vancouver final scientific closure - audit report")
    canvas.drawRightString(letter[0] - document.rightMargin, 0.35 * inch, f"Page {document.page}")
    canvas.restoreState()


styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="TitleCenter", parent=styles["Title"], alignment=TA_CENTER, fontSize=18, leading=22, textColor=colors.HexColor("#102A43"), spaceAfter=8))
styles.add(ParagraphStyle(name="Subtitle", parent=styles["Normal"], alignment=TA_CENTER, fontSize=9, leading=12, textColor=colors.HexColor("#52606D"), spaceAfter=16))
styles.add(ParagraphStyle(name="Section", parent=styles["Heading2"], fontSize=12, leading=15, textColor=colors.HexColor("#0B7285"), spaceBefore=10, spaceAfter=5))
styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], fontSize=9, leading=12, spaceAfter=5))
styles.add(ParagraphStyle(name="BulletSmall", parent=styles["BodyText"], fontSize=8.5, leading=11, leftIndent=12, firstLineIndent=-8, spaceAfter=3))
styles.add(ParagraphStyle(name="TableSmall", parent=styles["BodyText"], fontSize=6.5, leading=8))

story = [
    paragraph("IAS 2026 Vancouver final scientific closure", styles["TitleCenter"]),
    paragraph("Executive audit report | 2026-09-19 | not a submission paper or poster", styles["Subtitle"]),
]
for heading, lines in read_sections():
    story.append(paragraph(heading, styles["Section"]))
    for line in lines:
        story.append(paragraph(line, styles["BulletSmall"] if line.startswith("-") else styles["BodySmall"]))

story.append(paragraph("Gate adjudication", styles["Section"]))
table_data = [[paragraph("Gate", styles["TableSmall"]), paragraph("Status", styles["TableSmall"]), paragraph("Question", styles["TableSmall"]), paragraph("Observed", styles["TableSmall"])]]
for gate, status, question, observed in dashboard_rows():
    table_data.append([paragraph(gate, styles["TableSmall"]), paragraph(status, styles["TableSmall"]), paragraph(question, styles["TableSmall"]), paragraph(observed, styles["TableSmall"])])
table = Table(table_data, colWidths=[0.42 * inch, 1.38 * inch, 1.7 * inch, 3.25 * inch], repeatRows=1)
table.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#102A43")),
    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
    ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#BCCCDC")),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F5F7FA")]),
    ("LEFTPADDING", (0, 0), (-1, -1), 4),
    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ("TOPPADDING", (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
]))
story.append(table)
story.append(Spacer(1, 10))
story.append(paragraph("The final ZIP deliberately excludes FINAL_PAPER.pdf and FINAL_POSTER_DRAFT.pdf. The true same-model gate is stopped, so downstream mechanism, robustness, mixed-holdout, and scaling claims are not promoted.", styles["BodySmall"]))

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
doc = SimpleDocTemplate(str(OUTPUT), pagesize=letter, rightMargin=0.55 * inch, leftMargin=0.55 * inch, topMargin=0.55 * inch, bottomMargin=0.55 * inch, title="IAS 2026 Vancouver final scientific closure")
doc.build(story, onFirstPage=decorate, onLaterPages=decorate)
print(OUTPUT)
