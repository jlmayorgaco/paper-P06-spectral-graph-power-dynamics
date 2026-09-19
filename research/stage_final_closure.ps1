$ErrorActionPreference = 'Stop'

$repo = (Get-Location).Path
$sourceRoot = Join-Path $repo 'research\ias2026_bulletproof'
$targetRoot = Join-Path $repo 'research\ias2026_final_closure'

function Copy-FilteredTree {
    param(
        [string]$Source,
        [string]$Target,
        [switch]$SourceIsCode
    )
    $sourcePath = (Resolve-Path $Source).Path
    Get-ChildItem -LiteralPath $sourcePath -File -Recurse | ForEach-Object {
        $relative = $_.FullName.Substring($sourcePath.Length + 1)
        if ($relative -match '^bundle\\') { return }
        if ($relative -match '^paper\\') { return }
        if ($relative -match '__pycache__') { return }
        if ($relative -match '(FINAL_PAPER|FINAL_POSTER_DRAFT)\.pdf$') { return }
        if ($SourceIsCode) { $relative = Join-Path 'code' $relative }
        $destination = Join-Path $targetRoot $relative
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destination) | Out-Null
        Copy-Item -LiteralPath $_.FullName -Destination $destination -Force
    }
}

New-Item -ItemType Directory -Force -Path $targetRoot | Out-Null
Copy-FilteredTree -Source $sourceRoot -Target $targetRoot
Copy-FilteredTree -Source (Join-Path $sourceRoot 'src') -Target $targetRoot -SourceIsCode

function Copy-TreeAlias {
    param([string]$RelativeSource, [string]$RelativeTarget)
    $source = Join-Path $targetRoot $RelativeSource
    $target = Join-Path $targetRoot $RelativeTarget
    if (-not (Test-Path -LiteralPath $source -PathType Container)) { return }
    Get-ChildItem -LiteralPath $source -File -Recurse | ForEach-Object {
        $relative = $_.FullName.Substring($source.Length + 1)
        $destination = Join-Path $target $relative
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destination) | Out-Null
        Copy-Item -LiteralPath $_.FullName -Destination $destination -Force
    }
}

Copy-TreeAlias 'raw\gfl11' 'raw\canonical'
Copy-TreeAlias 'raw\true_same_model' 'raw\same_model'
Copy-TreeAlias 'raw\p2' 'raw\second_model\powerdynamics'
Copy-TreeAlias 'raw\p5' 'raw\second_model\simplegfldc'
Copy-TreeAlias 'raw\p4' 'raw\tds'
Copy-TreeAlias 'raw\p6' 'raw\holdout'

$inputMap = @(
    @('prereg\p6_genuine_holdout_inputs.csv', 'raw\holdout\p6_genuine_holdout_inputs.csv'),
    @('prereg\P6_GENUINE_HOLDOUT_PROTOCOL.md', 'raw\holdout\P6_GENUINE_HOLDOUT_PROTOCOL.md')
)
foreach ($item in $inputMap) {
    $source = Join-Path $targetRoot $item[0]
    $destination = Join-Path $targetRoot $item[1]
    if (Test-Path -LiteralPath $source -PathType Leaf) {
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destination) | Out-Null
        Copy-Item -LiteralPath $source -Destination $destination -Force
    }
}

New-Item -ItemType Directory -Force -Path (Join-Path $targetRoot 'raw\robust'), (Join-Path $targetRoot 'raw\scaling'), (Join-Path $targetRoot 'code\shared'), (Join-Path $targetRoot 'bundle') | Out-Null
$count = (Get-ChildItem -LiteralPath $targetRoot -File -Recurse).Count
Write-Output "staged_files=$count"
