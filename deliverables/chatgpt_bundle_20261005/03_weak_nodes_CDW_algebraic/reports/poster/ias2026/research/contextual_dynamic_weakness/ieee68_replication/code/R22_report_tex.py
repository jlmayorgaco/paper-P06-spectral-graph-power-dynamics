# ruff: noqa: E501
"""R22: docs/CDW68_FINAL_REPLICATION_REPORT.md -> .tex with the CDW converter (experiments/cdw/md2tex.py, used unchanged);
only the title block is replaced. Build: lualatex twice in docs/."""

from __future__ import annotations

import _r68 as R  # noqa: I001

import sys

sys.path.insert(0, str(R.CDW / "experiments" / "cdw"))
import md2tex  # noqa: E402

src = R.DOCS / "CDW68_FINAL_REPLICATION_REPORT.md"
tex = md2tex.convert(src.read_text(encoding="utf-8"))
tex = tex.replace(r"\date{2026-09-12}", r"\date{2026-09-13}").replace(
    r"\author{CDW campaign (branch research/contextual-dynamic-weakness)}", r"\author{CDW68 replication (branch research/cdw-ieee68-replication)}")
extra = "\n".join([r"\newunicodechar{⊊}{\ensuremath{\subsetneq}}", r"\newunicodechar{∅}{\ensuremath{\emptyset}}", r"\newunicodechar{⁻}{\textsuperscript{\textminus}}",
                   r"\newunicodechar{¹}{\textsuperscript{1}}", r"\newunicodechar{²}{\textsuperscript{2}}", r"\newunicodechar{⁶}{\textsuperscript{6}}"])
tex = tex.replace(r"\renewcommand{\arraystretch}{1.15}", extra + "\n" + r"\renewcommand{\arraystretch}{1.15}", 1)
(R.DOCS / "CDW68_FINAL_REPLICATION_REPORT.tex").write_text(tex, encoding="utf-8")
print("wrote", R.DOCS / "CDW68_FINAL_REPLICATION_REPORT.tex")
