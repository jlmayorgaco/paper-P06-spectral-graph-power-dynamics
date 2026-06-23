$ErrorActionPreference = "Stop"

New-Item -ItemType Directory -Force -Path .\build | Out-Null
Remove-Item -LiteralPath .\build\main.aux,.\build\main.bbl,.\build\main.blg,.\build\main.log,.\build\main.out,.\build\main.pdf -Force -ErrorAction SilentlyContinue
python .\figures\make_figures.py
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
if ($LASTEXITCODE -ne 0) { throw "pdflatex pass 1 failed" }
Copy-Item -LiteralPath .\refs.bib -Destination .\build\refs.bib -Force
Push-Location build
bibtex main
if ($LASTEXITCODE -ne 0) { Pop-Location; throw "bibtex failed" }
Pop-Location
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
if ($LASTEXITCODE -ne 0) { throw "pdflatex pass 2 failed" }
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
if ($LASTEXITCODE -ne 0) { throw "pdflatex pass 3 failed" }

Write-Host "Built build\main.pdf"
