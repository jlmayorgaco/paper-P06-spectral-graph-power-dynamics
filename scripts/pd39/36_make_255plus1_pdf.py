from __future__ import annotations

import html
from pathlib import Path

from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer

ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "docs" / "PD39_255PLUS1_FINAL_REPORT.md"
target = ROOT / "docs" / "PD39_255PLUS1_FINAL_REPORT.pdf"

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="TitlePD39", parent=styles["Title"], alignment=TA_CENTER, fontSize=19, leading=23, spaceAfter=16))
styles.add(ParagraphStyle(name="H1PD39", parent=styles["Heading1"], fontSize=14, leading=17, spaceBefore=12, spaceAfter=7))
styles.add(ParagraphStyle(name="H2PD39", parent=styles["Heading2"], fontSize=11, leading=14, spaceBefore=8, spaceAfter=5))
styles.add(ParagraphStyle(name="BodyPD39", parent=styles["BodyText"], fontSize=9, leading=12, spaceAfter=6))
styles.add(ParagraphStyle(name="BulletPD39", parent=styles["BodyText"], fontSize=8.7, leading=11, leftIndent=14, firstLineIndent=-8, spaceAfter=3))
styles.add(ParagraphStyle(name="MonoPD39", parent=styles["Code"], fontSize=7.4, leading=9, leftIndent=12, spaceAfter=2))


def clean(text: str) -> str:
    # The PDF skill requires ASCII hyphens in generated PDF content.
    text = text.replace("–", "-").replace("—", "-").replace("−", "-").replace("→", "->").replace("↔", "<->")
    return html.escape(text, quote=False)


story = []
for raw in source.read_text(encoding="utf-8").splitlines():
    line = raw.rstrip()
    if not line:
        story.append(Spacer(1, 3))
    elif line.startswith("# "):
        story.append(Paragraph(clean(line[2:]), styles["TitlePD39"]))
    elif line.startswith("## "):
        story.append(Paragraph(clean(line[3:]), styles["H1PD39"]))
    elif line.startswith("### "):
        story.append(Paragraph(clean(line[4:]), styles["H2PD39"]))
    elif line.startswith("- "):
        story.append(Paragraph("- " + clean(line[2:]), styles["BulletPD39"]))
    elif line.startswith("`") and line.endswith("`"):
        story.append(Paragraph(clean(line.strip("`")), styles["MonoPD39"]))
    else:
        story.append(Paragraph(clean(line), styles["BodyPD39"]))


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.drawString(0.65 * inch, 0.42 * inch, "PD39 255+1 mechanism validation - final report")
    canvas.drawRightString(7.85 * inch, 0.42 * inch, f"Page {doc.page}")
    canvas.restoreState()


doc = SimpleDocTemplate(
    str(target), pagesize=letter, rightMargin=0.65 * inch, leftMargin=0.65 * inch,
    topMargin=0.62 * inch, bottomMargin=0.62 * inch,
    title="PD39 255+1 - Final mechanism and generalization validation",
    author="PD39 research campaign",
)
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(target)
