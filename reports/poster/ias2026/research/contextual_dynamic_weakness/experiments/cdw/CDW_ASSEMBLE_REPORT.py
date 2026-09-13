"""Phase 25: assemble docs/report_parts/*.md into one ordered report, by section number.

Splits each part file into top-level '## N. ...' blocks (a leading '# Title' line, if
present, is kept as the document title once), sorts all blocks by N across every file,
and concatenates. This lets each part file be written/extended independently while the
final document is strictly numbered 1..N regardless of which file a section lives in.
"""

from __future__ import annotations

import re
from pathlib import Path

PARTS = Path(__file__).resolve().parents[2] / "docs" / "report_parts"
OUT_MD = Path(__file__).resolve().parents[2] / "docs" / "CDW_FULL_CAMPAIGN_REPORT.md"

ORDER = [
    "00_exec_summary.md", "01_front.md", "02_E1.md", "04_E2_E6_E8.md", "03_E34.md",
    "10_E05_E11_E23.md", "09_E07.md", "06_E9.md", "05_E12.md", "07_E13_E14.md",
    "08_E16_E17.md", "26_negative_pubcandidates_next.md", "28_limitations_repro.md",
]

HEAD_RE = re.compile(r"^## (\d+)\. ", re.M)


def blocks_of(text: str):
    """Yield (number, block_text) for every top-level '## N. ...' section."""

    idx = [m.start() for m in HEAD_RE.finditer(text)] + [len(text)]
    nums = [int(m.group(1)) for m in HEAD_RE.finditer(text)]
    for n, (a, b) in zip(nums, zip(idx, idx[1:], strict=False), strict=True):
        yield n, text[a:b].rstrip() + "\n"


def main():
    title = None
    all_blocks = []
    for name in ORDER:
        text = (PARTS / name).read_text(encoding="utf-8")
        m = re.match(r"^# (.+)\n", text)
        if m and title is None:
            title = m.group(1)
        for n, block in blocks_of(text):
            all_blocks.append((n, name, block))
    all_blocks.sort(key=lambda x: x[0])
    seen = [n for n, _, _ in all_blocks]
    gaps = [n for n in range(1, max(seen) + 1) if n not in seen]
    dups = sorted({n for n in seen if seen.count(n) > 1})
    if gaps or dups:
        raise SystemExit(f"assembly numbering error: missing={gaps} duplicated={dups}")
    out = [f"# {title}\n"]
    out += [b for _, _, b in all_blocks]
    OUT_MD.write_text("\n".join(out), encoding="utf-8")
    print(f"wrote {OUT_MD} ({len(all_blocks)} sections, 1..{max(seen)})")


if __name__ == "__main__":
    main()
