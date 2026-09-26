#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v pandoc >/dev/null || { echo 'pandoc is required' >&2; exit 1; }
command -v xelatex >/dev/null || { echo 'xelatex is required' >&2; exit 1; }
pandoc TEORIA_Y_DEMOSTRACIONES.md --standalone \
 --from=markdown+tex_math_dollars+pipe_tables --lua-filter=pdf_filter.lua \
 --pdf-engine=xelatex --toc --toc-depth=2 \
 -V documentclass=article -V fontsize=11pt -V geometry:margin=0.82in \
 -V mainfont='DejaVu Serif' -V monofont='DejaVu Sans Mono' \
 -V colorlinks=true -V linkcolor=accent -V urlcolor=accent \
 -H pdf_header.tex -o TEORIA_Y_DEMOSTRACIONES.pdf
