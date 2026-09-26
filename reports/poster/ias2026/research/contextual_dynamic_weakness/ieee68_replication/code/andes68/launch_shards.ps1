param([string]$Phase, [int]$N, [string]$Extra = "")
$V = 'C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\.venv\xtool-andes-gfl'
$env:USERPROFILE = "$V\home"; $env:HOME = "$V\home"
$env:OPENBLAS_NUM_THREADS = '1'; $env:OMP_NUM_THREADS = '1'; $env:MKL_NUM_THREADS = '1'
for ($k = 0; $k -lt $N; $k++) {
  $args = @('run68.py', 'shard', $Phase, "$N", "$k")
  if ($Extra) { $args += $Extra }
  Start-Process -FilePath "$V\Scripts\python.exe" -ArgumentList $args -WorkingDirectory $PSScriptRoot -RedirectStandardOutput "$PSScriptRoot\..\..\logs\${Phase}_shard$k.out" -RedirectStandardError "$PSScriptRoot\..\..\logs\${Phase}_shard$k.err" -WindowStyle Hidden
}
