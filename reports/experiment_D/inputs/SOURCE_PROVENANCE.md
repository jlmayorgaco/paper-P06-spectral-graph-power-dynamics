# Experiment D input provenance

The copied CSV files in this directory come from the installed PowerDynamics.jl 5.0.0 IEEE-39 example data directory:

`C:\Users\walla\.julia\packages\PowerDynamics\VzOiZ\docs\examples\ieee39data\`

The package version and git tree recorded by Experiment A are `5.0.0` and `1a32897016af608de2d2134922017614051f87b6`. The SHA-256 values below identify the exact source files used when the copies were made. They are recorded for provenance only; the analytic pipeline reads the frozen copies under this directory.

| Source file | SHA-256 |
|---|---|
| `bus.csv` | `d732a896e95c013fa3906941a540f2d9161ca99a2c60e571a890174680cb6943` |
| `branch.csv` | `b1858adb3e409750b5db24bc5073dda1e0683a59ef6e5966840cd73a18293972` |
| `load.csv` | `635e38cb30ed983be884d2e1e88cbe4069d69097a98969648390da1464df0d30` |
| `machine.csv` | `8cf9a02574f650d2c459f536a02ae6bbfe16122e927dbbd54ca54ec2556616a1` |
| `avr.csv` | `46d9062bc0db558f863faf5706829d109b3ca817297edc43d492ae7af39457bf` |
| `gov.csv` | `2cd24f6c43c0fb4175b0b6e5f579239e23f6b7a01805e0c71eb88b1fb892d7ac` |

The independent component transcriptions were compared against these PowerDynamics source files:

| Source file | SHA-256 |
|---|---|
| `src/Library/Renewables/ComposableInverter.jl` | `37d2d22aaadd14d9f089912257a58494dd17651b0fa54dd121a319ff7a5577ae` |
| `src/Library/Machines/SauerPaiMachine.jl` | `f9a62d020eafe800d82b0bc7ad38d3f34f5d0edb01666e3ce0adfaa4a7fe4570` |
| `src/Library/Controls/AVRs.jl` | `44deec55e02ef90bfbc2cc998a39bb7d8b6b0eb82dec965f97b6094070d835ad` |
| `src/Library/Controls/Govs.jl` | `9db33df19a280ad8674660fa42adbc76ba7b3f0c35ba41593212c0fd167f5561` |
| `src/Library/Loads/StaticLoads.jl` | `5d894f5f1f334e3401a7daf296705685619a6afdcb57b996f0eabe9face1a703` |
| `src/Library/Branches/PiLine_fault.jl` | `c07ec8a9d2d35a937895a0b836a712a8e82a1554994f3c5c7fd72c3b7c3de103` |
| `docs/examples/ieee39_part1.jl` | `bd9d4c02ad10fcf0d1041c0f7baf1e0a38aed084c1eef39516a4071dab089d9d` |

These source files were read for transcription and provenance. Their package modules are not imported or called by `src/bnd_design/`.
