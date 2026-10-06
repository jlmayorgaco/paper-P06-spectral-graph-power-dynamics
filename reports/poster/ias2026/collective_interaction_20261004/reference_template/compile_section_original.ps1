param(
    [Parameter(Mandatory = $true)]
    [string]$Section,
    [switch]$StrictLayout
)

$ErrorActionPreference = 'Stop'
$posterDir = Split-Path -Parent $PSCommandPath
$repoRoot = (Resolve-Path (Join-Path $posterDir '../../..')).Path
$sectionDir = Join-Path $posterDir 'sections'
$buildDir = Join-Path $posterDir 'build'
$previewDir = Join-Path $buildDir 'sections'
$outputDir = Join-Path $repoRoot 'output/pdf'
$available = @(Get-ChildItem -LiteralPath $sectionDir -Filter '*.tex' -File | Sort-Object Name | ForEach-Object BaseName)

if ($Section -notin @('all', 'full') -and $Section -notin $available) {
    throw "Sección desconocida: $Section. Opciones: $($available -join ', '), full, all."
}

New-Item -ItemType Directory -Force -Path $previewDir | Out-Null
$env:TEXMFCACHE = Join-Path $buildDir 'texmf-cache'
$env:TEXMFVAR = $env:TEXMFCACHE
New-Item -ItemType Directory -Force -Path $env:TEXMFCACHE | Out-Null

Push-Location $posterDir
try {
    if ($Section -eq 'all') { $targets = $available }
    elseif ($Section -eq 'full') { $targets = @() }
    else { $targets = @($Section) }

    foreach ($name in $targets) {
        $source = Join-Path $sectionDir "$name.tex"
        $pdf = Join-Path $previewDir "$name.pdf"
        $prefix = Join-Path $previewDir $name
        $consoleLog = Join-Path $previewDir "$name.console.log"
        & lualatex -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$previewDir" "-jobname=$name" $source *> $consoleLog
        if ($LASTEXITCODE -ne 0) {
            Get-Content -LiteralPath $consoleLog -Tail 25
            throw "No se pudo compilar $name."
        }
        & pdftoppm -f 1 -l 1 -r 100 -png -singlefile $pdf $prefix
        if ($LASTEXITCODE -ne 0) { throw "No se pudo renderizar $name en PNG." }
        $layoutWarnings = @(Select-String -LiteralPath (Join-Path $previewDir "$name.log") -Pattern 'Overfull|Underfull|LaTeX Error|Missing character')
        if ($layoutWarnings.Count -gt 0) {
            if ($StrictLayout) { throw "$name tiene $($layoutWarnings.Count) avisos de composición; revisa $name.log." }
            Write-Warning "$name tiene $($layoutWarnings.Count) avisos de composición."
        }
        Write-Host "Vista individual: $pdf"
    }

    if ($Section -in @('all', 'full')) {
        New-Item -ItemType Directory -Force -Path $outputDir | Out-Null
        $consoleLog = Join-Path $buildDir 'main.console.log'
        & lualatex -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$buildDir" -jobname=main 'main.tex' *> $consoleLog
        if ($LASTEXITCODE -ne 0) {
            Get-Content -LiteralPath $consoleLog -Tail 25
            throw 'No se pudo compilar el póster completo.'
        }
        $finalPdf = Join-Path $outputDir 'IAS2026_Poster_Pulido.pdf'
        $finalPrefix = Join-Path $outputDir 'IAS2026_Poster_Pulido'
        Copy-Item -LiteralPath (Join-Path $buildDir 'main.pdf') -Destination $finalPdf -Force
        & pdftoppm -f 1 -l 1 -r 75 -png -singlefile $finalPdf $finalPrefix
        if ($LASTEXITCODE -ne 0) { throw 'No se pudo renderizar el póster completo en PNG.' }
        $layoutWarnings = @(Select-String -LiteralPath (Join-Path $buildDir 'main.log') -Pattern 'Overfull|Underfull|LaTeX Error|Missing character')
        if ($layoutWarnings.Count -gt 0) {
            if ($StrictLayout) { throw "El póster completo tiene $($layoutWarnings.Count) avisos de composición; revisa main.log." }
            Write-Warning "El póster completo tiene $($layoutWarnings.Count) avisos de composición."
        }
        Write-Host "Póster completo: $finalPdf"
    }
}
finally {
    Pop-Location
}
