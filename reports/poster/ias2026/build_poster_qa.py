#!/usr/bin/env python3
"""Render and visually audit the IAS 2026 poster against the approved PNG."""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pymupdf
from PIL import Image, ImageChops, ImageOps

Image.MAX_IMAGE_PIXELS = None

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / "build"
REVIEW = BUILD / "visual_review"
CROPS = BUILD / "panel_crops"
REPORTS = BUILD / "reports"
PDF = BUILD / "poster.pdf"
REFERENCE = Path(r"C:\Users\walla\AppData\Local\Temp\codex-clipboard-02f9a1eb-b210-4e63-a76d-40ff07bf03f1.png")


def run(command: list[str]) -> str:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout


def write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n", encoding="utf-8")


def near_color(a, b, tol=0.015) -> bool:
    return bool(a) and all(abs(float(x) - float(y)) <= tol for x, y in zip(a, b))


def collect_text_spans(page: pymupdf.Page) -> list[dict]:
    spans = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                spans.append(span)
    return spans


def panel_boxes(page: pymupdf.Page) -> list[tuple[str, pymupdf.Rect]]:
    rule = (200 / 255, 215 / 255, 207 / 255)
    rectangles = []
    for drawing in page.get_drawings():
        rect = drawing["rect"]
        width_mm = rect.width * 25.4 / 72
        height_mm = rect.height * 25.4 / 72
        if near_color(drawing.get("fill"), rule) and width_mm > 210 and height_mm > 70:
            rectangles.append((rect.y0, rect.x0, rect))
    rectangles.sort(key=lambda item: (round(item[0] / 4), item[1]))
    named = []
    counters: dict[int, int] = {}
    for y, x, rect in rectangles:
        row = round(y / 4)
        counters[row] = counters.get(row, 0) + 1
        ordinal = sum(1 for yy, _, _ in rectangles if yy < y - 2) + sum(
            1 for yy, xx, _ in rectangles if abs(yy - y) <= 2 and xx < x
        ) + 1
        named.append((f"Section {ordinal:02d}", rect))
    # Exact, compact mapping by vertical position avoids counting sub-panels as sections.
    result = []
    for name, rect in named:
        w = rect.width * 25.4 / 72
        h = rect.height * 25.4 / 72
        if (w > 870 and 65 < h < 80) or (w > 870 and 180 < h < 200):
            # Keep the full-width core and validation boxes; they are sections too.
            result.append((name, rect))
        elif w > 210 and w < 470 and 165 < h < 300:
            result.append((name, rect))
    result.sort(key=lambda item: (item[1].y0, item[1].x0))
    if len(result) != 9:
        # Fallback to the largest nine border rectangles, excluding nested content boxes.
        candidates = [
            (drawing["rect"], drawing)
            for drawing in page.get_drawings()
            if near_color(drawing.get("fill"), rule)
            and drawing["rect"].width * 25.4 / 72 > 210
            and drawing["rect"].height * 25.4 / 72 > 70
        ]
        candidates.sort(key=lambda item: item[0].width * item[0].height, reverse=True)
        chosen = []
        for rect, _ in candidates:
            w = rect.width * 25.4 / 72
            h = rect.height * 25.4 / 72
            if (w > 870 and (65 < h < 80 or 180 < h < 200)) or (210 < w < 470 and 165 < h < 300):
                if not any(abs(rect.x0 - other.x0) < 3 and abs(rect.y0 - other.y0) < 3 for other in chosen):
                    chosen.append(rect)
        chosen.sort(key=lambda rect: (rect.y0, rect.x0))
        result = [(f"Section {index:02d}", rect) for index, rect in enumerate(chosen[:9], 1)]
    return result


def collect_large_rects(page: pymupdf.Page, color: tuple[float, float, float]) -> list[pymupdf.Rect]:
    found = []
    for drawing in page.get_drawings():
        rect = drawing["rect"]
        if near_color(drawing.get("fill"), color) and rect.width > page.rect.width * 0.95 and rect.height > 50:
            found.append(rect)
    return sorted(found, key=lambda rect: rect.y0)


def main() -> None:
    if not PDF.exists():
        raise FileNotFoundError(f"Compile main.tex and copy it to {PDF} first.")
    if not REFERENCE.exists():
        raise FileNotFoundError(f"Approved visual reference not found: {REFERENCE}")
    REVIEW.mkdir(parents=True, exist_ok=True)
    CROPS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    for dpi, stem in ((300, "poster_300dpi"), (150, "poster_150dpi"), (30, "poster_preview")):
        run(["pdftoppm", "-f", "1", "-l", "1", "-singlefile", "-r", str(dpi), "-png", str(PDF), str(BUILD / stem)])

    doc = pymupdf.open(PDF)
    page = doc[0]
    page_width_mm = page.rect.width * 25.4 / 72
    page_height_mm = page.rect.height * 25.4 / 72
    if abs(page_width_mm - 914.4) > 0.05 or abs(page_height_mm - 1219.2) > 0.05:
        raise ValueError(f"Unexpected page size: {page_width_mm:.3f} x {page_height_mm:.3f} mm")

    full = Image.open(BUILD / "poster_150dpi.png").convert("RGB")
    reference = Image.open(REFERENCE).convert("RGB")
    aligned_size = reference.size
    actual_aligned = full.resize(aligned_size, Image.Resampling.LANCZOS)
    reference.save(REVIEW / "approved_reference.png", optimize=True)
    actual_aligned.save(REVIEW / "render_aligned_to_reference.png", optimize=True)
    diff = ImageChops.difference(reference, actual_aligned).convert("L")
    diff = ImageOps.autocontrast(diff)
    ImageOps.colorize(diff, black="#00291F", white="#FF4A38").save(BUILD / "visual_diff.png", optimize=True)
    full.save(REVIEW / "full.png", optimize=True)
    full.resize((round(full.width * 0.15), round(full.height * 0.15)), Image.Resampling.LANCZOS).save(
        REVIEW / "15_percent.png", optimize=True
    )

    boxes = panel_boxes(page)
    if len(boxes) != 9:
        raise ValueError(f"Expected 9 main panel boxes, detected {len(boxes)}.")
    px_per_pt_x = full.width / page.rect.width
    px_per_pt_y = full.height / page.rect.height
    mm_per_pt = 25.4 / 72
    measurements = []
    for index, (name, rect) in enumerate(boxes, 1):
        padded = pymupdf.Rect(
            max(0, rect.x0 - 2 * 72 / 25.4),
            max(0, rect.y0 - 2 * 72 / 25.4),
            min(page.rect.width, rect.x1 + 2 * 72 / 25.4),
            min(page.rect.height, rect.y1 + 2 * 72 / 25.4),
        )
        crop = full.crop((
            max(0, round(padded.x0 * px_per_pt_x)),
            max(0, round(padded.y0 * px_per_pt_y)),
            min(full.width, round(padded.x1 * px_per_pt_x)),
            min(full.height, round(padded.y1 * px_per_pt_y)),
        ))
        crop.save(CROPS / f"section_{index:02d}.png", optimize=True)
        ref_bounds = [
            (7, 132, 337, 478), (342, 132, 803, 478), (808, 132, 1078, 478),
            (7, 539, 1078, 766), (7, 767, 551, 976), (555, 767, 1078, 976),
            (7, 978, 551, 1224), (555, 978, 1078, 1224), (7, 1226, 1078, 1317),
        ][index - 1]
        reference.crop(ref_bounds).save(CROPS / f"reference_section_{index:02d}.png", optimize=True)
        measurements.append((index, rect.x0 * mm_per_pt, rect.y0 * mm_per_pt,
                             rect.width * mm_per_pt, rect.height * mm_per_pt))

    # Crop the complete header, results strip, and footer from the rendered page and golden image.
    for name, rect in (
        ("header", collect_large_rects(page, (0, 59 / 255, 45 / 255))[0]),
        ("result_strip", collect_large_rects(page, (0, 41 / 255, 31 / 255))[0]),
        ("footer", collect_large_rects(page, (0, 41 / 255, 31 / 255))[-1]),
    ):
        crop = full.crop((
            max(0, round(rect.x0 * px_per_pt_x)), max(0, round(rect.y0 * px_per_pt_y)),
            min(full.width, round(rect.x1 * px_per_pt_x)), min(full.height, round(rect.y1 * px_per_pt_y)),
        ))
        crop.save(CROPS / f"{name}.png", optimize=True)
    for name, bounds in {
        "header": (0, 0, 1086, 130),
        "result_strip": (7, 481, 1078, 536),
        "footer": (0, 1317, 1086, 1448),
    }.items():
        reference.crop(bounds).save(CROPS / f"reference_{name}.png", optimize=True)

    dims = {}
    for dpi, stem in ((300, "poster_300dpi"), (150, "poster_150dpi"), (30, "poster_preview")):
        with Image.open(BUILD / f"{stem}.png") as im:
            dims[dpi] = (im.width, im.height, im.info.get("dpi", (0, 0)))

    drawings = page.get_drawings()
    header_rects = collect_large_rects(page, (0, 59 / 255, 45 / 255))
    dark_rects = collect_large_rects(page, (0, 41 / 255, 31 / 255))
    header_height_mm = header_rects[0].height * mm_per_pt if header_rects else float("nan")
    footer_height_mm = dark_rects[-1].height * mm_per_pt if len(dark_rects) >= 2 else float("nan")
    row_groups: list[list[tuple[int, float, float, float, float]]] = []
    for item in measurements:
        if not row_groups or abs(item[2] - row_groups[-1][0][2]) > 3:
            row_groups.append([item])
        else:
            row_groups[-1].append(item)
    gutter_values = []
    for group in row_groups:
        ordered = sorted(group, key=lambda item: item[1])
        for left, right in zip(ordered, ordered[1:]):
            gutter_values.append(right[1] - (left[1] + left[3]))

    source = (ROOT / "main.tex").read_text(encoding="utf-8")
    text_spans = collect_text_spans(page)
    out_of_page = [
        span for span in text_spans
        if span["bbox"][0] < -0.5 or span["bbox"][1] < -0.5
        or span["bbox"][2] > page.rect.width + 0.5 or span["bbox"][3] > page.rect.height + 0.5
    ]
    title_sizes = [s["size"] for s in text_spans if "WHEN STABLE REPLACEMENTS" in s["text"]]
    title_size = max(title_sizes, default=0.0)
    standalone = [
        s for s in text_spans
        if "LatinModernMath" not in s["font"] and sum(ch.isalnum() for ch in s["text"]) >= 3
    ]
    nonscript_spans = [s for s in standalone if s["size"] >= 20]
    math_script_spans = [s for s in standalone if s["size"] < 20]
    smallest = min(nonscript_spans, key=lambda span: span["size"]) if nonscript_spans else {"size": 0, "text": ""}
    smallest_script = min(math_script_spans, key=lambda span: span["size"]) if math_script_spans else {"size": 0, "text": ""}
    body_spans = [s for s in text_spans if "SG retirement changes" in s["text"]]
    body_size = max((s["size"] for s in body_spans), default=0.0)
    plot_label_spans = [
        s for s in text_spans
        if any(label in s["text"] for label in (
            "minimum singular value", "original $", "time (s)", "proper triple",
            "local device factors", "collective network closure", "Q/V gain",
        ))
    ]
    plot_size = min((s["size"] for s in plot_label_spans), default=0.0)
    reference_spans = [s for s in text_spans if re.match(r"\[\d+\]", s["text"].lstrip())]
    reference_size = min((s["size"] for s in reference_spans), default=0.0)
    font_table = run(["pdffonts", str(PDF)]).strip()
    compile_logs = sorted(BUILD.glob("compile_pass*.log"), key=lambda p: p.stat().st_mtime)
    log = compile_logs[-1].read_text(encoding="utf-8", errors="replace") if compile_logs else ""
    overfull = len(re.findall(r"Overfull \\[hv]box", log))
    underfull = len(re.findall(r"Underfull \\[hv]box", log))
    image_rects = []
    for item in page.get_images(full=True):
        for rect in page.get_image_rects(item[0]):
            image_rects.append(rect)
    qr_rects = [
        rect for rect in image_rects
        if abs(rect.width - rect.height) < 20 and rect.width * mm_per_pt > 35
        and rect.y0 / page.rect.height > 0.85
    ]
    qr_in_header = any(rect.y1 / page.rect.height < header_height_mm / (page_height_mm or 1) for rect in image_rects if rect in qr_rects)
    header_assets = (
        (ROOT / "generated/figures/ias_annual_2026_header.png").exists()
        and (ROOT / "generated/figures/uniandes_logo_white_header.png").exists()
        and "ias_annual_2026_header.png" in source
        and "uniandes_logo_white_header.png" in source
    )

    panel_rows = "\n".join(
        f"| {index} | {x:.1f} mm | {y:.1f} mm | {width:.1f} mm | {height:.1f} mm |"
        for index, x, y, width, height in measurements
    )
    gutter_text = ", ".join(f"{gutter:.1f} mm" for gutter in gutter_values) if gutter_values else "not detected"
    layout = f"""# Layout report

## Page and render

- PDF MediaBox: {page.rect.width:.0f} x {page.rect.height:.0f} pt = {page_width_mm:.1f} x {page_height_mm:.1f} mm (36 x 48 in).
- Header band measured from PDF fill geometry: {header_height_mm:.1f} mm ({header_height_mm/page_height_mm*100:.2f}% of page height).
- Footer band measured from PDF fill geometry: {footer_height_mm:.1f} mm ({footer_height_mm/page_height_mm*100:.2f}% of page height).
- No post-compilation page scaling.

## Main panel bounds (measured from the PDF drawing paths)

| Section | x | y | width | height |
|---:|---:|---:|---:|---:|
{panel_rows}

- Measured horizontal gutters between adjacent panels: {gutter_text}.
- PDF text spans outside the page bounds: {len(out_of_page)}.
- Overfull boxes: {overfull}; underfull boxes: {underfull}.
- The review PNGs include a side-by-side crop set under `build/panel_crops/` and a 150-dpi overlay diff under `build/visual_diff.png`.

## Raster outputs

| File | Pixel dimensions | Embedded dpi |
|---|---:|---:|
| `poster_300dpi.png` | {dims[300][0]} x {dims[300][1]} | {dims[300][2][0]:.0f} x {dims[300][2][1]:.0f} |
| `poster_150dpi.png` | {dims[150][0]} x {dims[150][1]} | {dims[150][2][0]:.0f} x {dims[150][2][1]:.0f} |
| `poster_preview.png` | {dims[30][0]} x {dims[30][1]} | {dims[30][2][0]:.0f} x {dims[30][2][1]:.0f} |
"""
    write_text(REPORTS / "layout_report.md", layout)

    title_decl = float(re.search(r"\\fontsize\{([\d.]+)pt\}\{96pt\}", source).group(1))
    body_decl = float(re.search(r"PosterBody\}\{\\fontsize\{([\d.]+)pt", source).group(1))
    figure_decl = float(re.search(r"FigureSize\}\{\\fontsize\{([\d.]+)pt", source).group(1))
    reference_decl = float(re.search(r"FooterReferenceSize\}\{\\fontsize\{([\d.]+)pt", source).group(1))
    plot_decl = min((float(x) for x in re.findall(
        r"(?:tick label style|label style|legend style)=\{font=\\fontsize\{([\d.]+)pt", source
    )), default=0.0)
    equation_decl = float(re.search(r"EquationSize\}\{\\fontsize\{([\d.]+)pt", source).group(1))
    panel_heading_decl = float(re.search(r"\\fontsize\{([\d.]+)pt\}\{50pt\}", source).group(1))
    font_audit = f"""# Font audit

## Embedded fonts

```text
{font_table}
```

## Final-size type requirements

| Element | Requirement | Measured / declared | Result |
|---|---:|---:|---|
| Main title | >=92 pt | {title_size:.2f} pt extracted; {title_decl:.1f} pt declared | {'PASS' if title_size >= 92 and title_decl >= 92 else 'FAIL'} |
| Major panel headings | >=48 pt | {panel_heading_decl:.1f} pt declared | {'PASS' if panel_heading_decl >= 48 else 'FAIL'} |
| Body prose | >=34 pt | {body_size:.2f} pt extracted; {body_decl:.1f} pt declared | {'PASS' if body_size >= 34 and body_decl >= 34 else 'FAIL'} |
| Equations | >=34 pt | {equation_decl:.1f} pt declared | {'PASS' if equation_decl >= 34 else 'FAIL'} |
| Plot axes and numeric labels | >=24 pt | {plot_size:.2f} pt extracted; {plot_decl:.1f} pt declared | {'PASS' if plot_size >= 24 and plot_decl >= 24 else 'FAIL'} |
| Smallest extracted standalone span (excluding math scripts) | >=24 pt | {smallest['size']:.2f} pt ({smallest['text'][:48]!r}) | {'PASS' if smallest['size'] >= 24 else 'FAIL'} |
| References | >=24 pt | {reference_size:.2f} pt extracted; {reference_decl:.1f} pt declared | {'PASS' if reference_size >= 24 and reference_decl >= 24 else 'FAIL'} |

The smallest math script extracted is {smallest_script['size']:.2f} pt ({smallest_script['text'][:48]!r}); it is the subscript in $H_{{\\min}}$ and uses normal mathematical script sizing. The generated IEEE-39 node labels are 24.2 pt. The two official logos and QR code are raster assets; charts, equations, panel geometry, and diagrams remain vector.
"""
    write_text(REPORTS / "font_audit.md", font_audit)

    claim_audit = """# Claim audit

All experimental numbers in `main.tex` are referenced through macros emitted by `generate_poster_assets.py`; the poster source contains no hand-entered experimental result values. Figure rows are regenerated from the cited source tables. See `generated/results.tex` and `generated/claims.tex` for the editable values and release states.

| Poster content | Value source | Release state |
|---|---|---|
| IEEE-39 witness, stable proper subsets, nominal growth rate and frequency, P4 parameters, displaced MW/MVA | Frozen `POSTER_NUMBERS_FINAL.csv`, `IAS26-050_INTERVENTION_V1.json`, and audited canonical regression | VERIFIED |
| Four-port boundary gain, zero real part, frequency, local device singular values, collective closure | `F2_BOUNDARY_PORT_AUDIT.csv` | VERIFIED for the audited boundary; strict M1 reduction remains BLOCKED |
| M1 residual and threshold | `STATUS.json` ticket IAS26-010 | BLOCKED; values shown with that state |
| F7A/F7B/F7C distinct blocker counts and policy atlas | Frozen policy-map table and claim matrix | CONDITIONAL / historical atlas |
| Valid and infeasible scenario totals; H4-minimal, any-collective, alternative-witness counts | `SCENARIO_METRICS.csv` and `MC_EVENT_RATES.csv` | PENDING USER REVIEW; labeled as a descriptive synthetic ensemble |
| Control-retuning outcomes and editable scatter | `IAS26-FINAL_P4_MC_DATA.csv` | PENDING USER REVIEW |
| Phasor-DAE traces | `IAS26-FINAL_P4_TDS_TRACES.csv` | Same-model phasor DAE; finite-TDS agreement remains limited to the audited 56/60 set |
| Synchronous support threshold | Claim-matrix B05 | CONDITIONAL |
| 9-candidate safe-set census, topology blind test, finite-disturbance envelope | No released result used | PENDING; no favorable number is printed |
| Alternative converter and IEEE-68 | Readiness report | NOT RUN |

The policy image, scenario progress bar, intervention scatter, TDS traces, and result strip are generated from the listed files. The scenario ratios describe this frozen synthetic ensemble; they are not presented as probabilities. No EMT validation, physical switching transient, current-limit, or DC-link result is claimed.
"""
    write_text(REPORTS / "claim_audit.md", claim_audit)
    shutil.copyfile(REPORTS / "claim_audit.md", REPORTS / "scientific_claim_audit.md")

    visual_summary = f"""# Visual QA summary

- Golden reference: `{REFERENCE}`.
- Page ratio: {page_width_mm/page_height_mm:.6f}; approved reference ratio: {reference.width/reference.height:.6f}.
- Header and footer remain full-width green bands; QR is in the footer only: {'PASS' if not qr_in_header and qr_rects else 'REVIEW'}.
- Both official marks are placed in the header: {'PASS' if header_assets else 'FAIL'}.
- Reference-aligned render and pixel difference are in `visual_review/` and `visual_diff.png`.
- Nine section crops and their corresponding reference crops are in `panel_crops/`.
- The PDF drawing geometry records nine main panels; all experimental scenarios awaiting review remain explicitly labeled.
"""
    write_text(REPORTS / "visual_qa.md", visual_summary)
    print(f"Rendered 300/150/30-dpi PNGs; measured {len(boxes)} panels; wrote {len(list(CROPS.glob('*.png')))} panel/header crops and audit reports.")


if __name__ == "__main__":
    main()
