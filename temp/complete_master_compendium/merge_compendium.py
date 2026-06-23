from pathlib import Path

from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[1]

parts = [
    ROOT / "complete_master_compendium" / "build" / "cover_reports.pdf",
    ROOT / "internal_unified_theory" / "main.pdf",
    ROOT / "paper_ieee_transactions" / "main.pdf",
]

out = ROOT / "complete_master_compendium" / "complete_master_compendium.pdf"

writer = PdfWriter()
page_counts = []
for part in parts:
    reader = PdfReader(str(part))
    page_counts.append((part.relative_to(ROOT).as_posix(), len(reader.pages)))
    for page in reader.pages:
        writer.add_page(page)

with out.open("wb") as f:
    writer.write(f)

print(f"Wrote {out}")
for name, pages in page_counts:
    print(f"{name}: {pages} pages")
print(f"Total pages: {sum(pages for _, pages in page_counts)}")
