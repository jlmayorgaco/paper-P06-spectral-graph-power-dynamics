# ruff: noqa: E501
"""Minimal Markdown -> LaTeX converter for the CDW report (restricted syntax, lualatex + fontspec).

Supported: '# Title' (first line), '## ' section, '### ' subsection, '#### ' paragraph,
paragraphs with **bold**, *italic*, `code`, $math$ (passed through), '- ' lists (2-space nesting),
'1. ' enumerations, pipe tables, '![caption](path)' figures, '```' code blocks, '> ' quotes,
'---' page break.
Usage: python md2tex.py in.md out.tex
"""

from __future__ import annotations

import re
import sys

SPECIAL = {"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}


def esc(s: str) -> str:
    out = []
    for ch in s:
        if ch == "\\":
            out.append(r"\textbackslash{}")
        else:
            out.append(SPECIAL.get(ch, ch))
    return "".join(out)


def strip_num(heading: str) -> str:
    """Drop a leading manual 'N.' / 'N.N' numbering — LaTeX numbers sections itself."""

    return re.sub(r"^\d+(\.\d+)*\.\s+", "", heading)


def inline(s: str) -> str:
    # split out math and code first
    parts = re.split(r"(\$[^$]+\$|`[^`]+`)", s)
    res = []
    for p in parts:
        if p.startswith("$") and p.endswith("$") and len(p) > 1:
            res.append(p)
        elif p.startswith("`") and p.endswith("`"):
            code = esc(p[1:-1])
            for ch in ("/", "\\_", "="):
                code = code.replace(ch, ch + r"\allowbreak{}")
            res.append(r"\texttt{" + code + "}")
        else:
            t = esc(p)
            t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
            t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", t)
            res.append(t)
    return "".join(res)


PREAMBLE = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{fontspec}
\setmainfont{Cambria}
\setmonofont{Consolas}[Scale=0.85]
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{float}
\usepackage{longtable,booktabs,array}
\usepackage{enumitem}
\setlist{nosep,leftmargin=1.4em}
\usepackage[hidelinks]{hyperref}
\usepackage{xcolor}
\usepackage{fancyvrb}
\usepackage{newunicodechar}
\newunicodechar{∪}{\ensuremath{\cup}}
\newunicodechar{∩}{\ensuremath{\cap}}
\newunicodechar{⊗}{\ensuremath{\otimes}}
\newunicodechar{⊥}{\ensuremath{\perp}}
\newunicodechar{⊆}{\ensuremath{\subseteq}}
\newunicodechar{⊂}{\ensuremath{\subset}}
\renewcommand{\arraystretch}{1.15}
\setlength{\parskip}{3pt}
\setlength{\parindent}{0pt}
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}
"""


def split_row(r: str) -> list[str]:
    """Split a markdown table row on unescaped '|'; '\\|' becomes a literal pipe."""

    r = r.strip()
    if r.startswith("|"):
        r = r[1:]
    if r.endswith("|") and not r.endswith(r"\|"):
        r = r[:-1]
    cells = re.split(r"(?<!\\)\|", r)
    return [c.strip().replace(r"\|", "|") for c in cells]


def table(rows):
    cells = [split_row(r) for r in rows]
    head, body = cells[0], [c for c in cells[2:]]
    n = len(head)
    lens = [max(len(r[j]) if j < len(r) else 0 for r in [head] + body) for j in range(n)]
    tot = sum(max(x, 3) for x in lens)
    widths = [max(0.06, 0.94 * max(x, 3) / tot) for x in lens]
    spec = "".join(f"L{{{w:.3f}\\linewidth}}" for w in widths)
    out = [r"{\small", r"\begin{longtable}{" + spec + "}", r"\toprule",
           " & ".join(r"\textbf{" + inline(h) + "}" for h in head) + r" \\", r"\midrule", r"\endhead"]
    for r in body:
        r = r + [""] * (n - len(r))
        out.append(" & ".join(inline(c) for c in r[:n]) + r" \\")
    out += [r"\bottomrule", r"\end{longtable}}"]
    return "\n".join(out)


def convert(md: str) -> str:
    lines = md.splitlines()
    out = [PREAMBLE]
    title = lines[0].lstrip("# ").strip() if lines and lines[0].startswith("# ") else "Report"
    out += [r"\title{" + inline(title) + "}", r"\date{2026-09-12}", r"\author{CDW campaign (branch research/contextual-dynamic-weakness)}",
            r"\begin{document}", r"\maketitle", r"\tableofcontents", r"\newpage"]
    i = 1 if lines and lines[0].startswith("# ") else 0
    list_stack = []
    para = []

    def flush_para():
        if para:
            out.append(inline(" ".join(para)) + "\n")
            para.clear()

    def close_lists(level=0):
        while len(list_stack) > level:
            out.append(r"\end{" + list_stack.pop() + "}")

    while i < len(lines):
        ln = lines[i]
        s = ln.rstrip()
        if s.startswith("```"):
            flush_para()
            close_lists()
            j = i + 1
            block = []
            while j < len(lines) and not lines[j].startswith("```"):
                block.append(lines[j])
                j += 1
            out.append(r"\begin{Verbatim}[fontsize=\footnotesize]" + "\n" + "\n".join(block) + "\n" + r"\end{Verbatim}")
            i = j + 1
            continue
        if s.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[\s:|-]+\|?$", lines[i + 1].strip()):
            flush_para()
            close_lists()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(lines[i])
                i += 1
            out.append(table(rows))
            continue
        m = re.match(r"^!\[(.*)\]\((.*)\)$", s.strip())
        if m:
            flush_para()
            close_lists()
            out.append(r"\begin{figure}[H]\centering\includegraphics[width=0.95\linewidth]{" + m.group(2) + r"}\caption{" + inline(m.group(1)) + r"}\end{figure}")
            i += 1
            continue
        if s.startswith("#### "):
            flush_para(); close_lists()
            out.append(r"\paragraph{" + inline(strip_num(s[5:])) + "}")
        elif s.startswith("### "):
            flush_para(); close_lists()
            out.append(r"\subsection{" + inline(strip_num(s[4:])) + "}")
        elif s.startswith("## "):
            flush_para(); close_lists()
            out.append(r"\section{" + inline(strip_num(s[3:])) + "}")
        elif s.strip() == "---":
            flush_para(); close_lists()
            out.append(r"\clearpage")
        elif s.lstrip().startswith(("- ", "* ")) or re.match(r"^\s*\d+\. ", s):
            flush_para()
            indent = len(s) - len(s.lstrip())
            level = indent // 2 + 1
            kind = "enumerate" if re.match(r"^\s*\d+\. ", s) else "itemize"
            while len(list_stack) > level:
                out.append(r"\end{" + list_stack.pop() + "}")
            while len(list_stack) < level:
                list_stack.append(kind)
                out.append(r"\begin{" + kind + "}")
            text = re.sub(r"^\s*(\d+\.|-|\*) ", "", s)
            # continuation lines
            j = i + 1
            while j < len(lines) and lines[j].strip() and not re.match(r"^\s*(\d+\.|-|\*) ", lines[j]) and not lines[j].startswith(("#", "|", "```", "!")) and (len(lines[j]) - len(lines[j].lstrip())) > indent:
                text += " " + lines[j].strip()
                j += 1
            out.append(r"\item " + inline(text))
            i = j
            continue
        elif s.startswith("> "):
            flush_para(); close_lists()
            out.append(r"\begin{quote}\itshape " + inline(s[2:]) + r"\end{quote}")
        elif not s.strip():
            flush_para()
            close_lists()
        else:
            para.append(s.strip())
        i += 1
    flush_para()
    close_lists()
    out.append(r"\end{document}")
    return "\n".join(out)


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    with open(src, encoding="utf-8") as f:
        tex = convert(f.read())
    with open(dst, "w", encoding="utf-8") as f:
        f.write(tex)
    print("wrote", dst)
