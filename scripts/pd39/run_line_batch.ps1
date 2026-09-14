param(
    [Parameter(Mandatory=$true)][int]$Start,
    [Parameter(Mandatory=$true)][int]$End,
    [Parameter(Mandatory=$true)][string]$Out,
    [Parameter(Mandatory=$true)][string]$SummaryOut
)
$env:PD39_LINE_START = [string]$Start
$env:PD39_LINE_END = [string]$End
$env:PD39_LINE_OUT = $Out
$env:PD39_LINE_SUMMARY_OUT = $SummaryOut
& 'C:\Users\walla\AppData\Local\Programs\Julia\julia.exe' --project=. scripts\pd39\17_line_repair.jl
