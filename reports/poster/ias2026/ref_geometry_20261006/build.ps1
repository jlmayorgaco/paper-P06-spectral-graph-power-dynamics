# Geometry-only build: compile twice (remember-picture overlay), render 75 dpi (2700 x 3150 = 36 x 42 in) and compare with the reference.
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
foreach ($i in 1..2) { lualatex -interaction=nonstopmode -halt-on-error main.tex | Out-Null }
pdftoppm -r 75 -png -singlefile main.pdf render_75dpi
pdftoppm -scale-to-x 1161 -scale-to-y 1355 -png -singlefile main.pdf render_ref_scale
python compare.py
