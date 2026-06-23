$ErrorActionPreference = "Stop"

python .\scripts\generate_demo_and_figures.py
if ($LASTEXITCODE -ne 0) { throw "generate_demo_and_figures.py failed with exit code $LASTEXITCODE" }

python .\scripts\generate_enhanced_poster_assets.py
if ($LASTEXITCODE -ne 0) { throw "generate_enhanced_poster_assets.py failed with exit code $LASTEXITCODE" }

python .\scripts\generate_modal_substitutability_assets.py
if ($LASTEXITCODE -ne 0) { throw "generate_modal_substitutability_assets.py failed with exit code $LASTEXITCODE" }

python .\scripts\verify_inertia_lever_space.py
if ($LASTEXITCODE -ne 0) { throw "verify_inertia_lever_space.py failed with exit code $LASTEXITCODE" }

New-Item -ItemType Directory -Force -Path .\build | Out-Null
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build poster.tex
if ($LASTEXITCODE -ne 0) { throw "pdflatex failed with exit code $LASTEXITCODE" }

Copy-Item -LiteralPath .\build\poster.pdf -Destination .\poster_A0_ieee_session.pdf -Force

Write-Host "Built poster_A0_ieee_session.pdf"
