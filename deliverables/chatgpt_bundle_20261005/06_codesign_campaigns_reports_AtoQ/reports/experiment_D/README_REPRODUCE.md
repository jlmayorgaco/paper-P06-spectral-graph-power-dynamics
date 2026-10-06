# Experiment D reproduction

## Frozen outcome

ExpD certifies the saturated replacement optimum on the existing preregistered PLL domain: `rho in [0,1]`, `beta in [0.9,1.1]`, `Kp=beta*Kp0`, and `Ki=beta^2*Ki0`. Bus 33 reaches 632 MW at `rho*=1`, `Kp*=31.4159265359 rad/s`, and `Ki*=246.7401100272 rad/s^2`. The same endpoint is feasible at buses 30, 35, and 37. The exact objective certificate is the physical upper bound `P0*rho <= P0` plus a full-spectrum-feasible witness at `rho=1`; resultants and nonbinding pole-boundary branches are not claimed.

## Reproduction commands

Run the structural local-model audit:

```powershell
julia --project=. experiments/bnd_expD/run_expD_structure.jl
```

Freeze the analytical endpoint candidate using only frozen ExpC full-state matrices:

```powershell
julia --project=. experiments/bnd_expD/freeze_endpoint_candidates.jl
```

Only after the candidate and its SHA-256 sidecar are frozen, validate the four endpoints and evaluate the 132-cell PowerDynamics grid:

```powershell
julia --project=. experiments/bnd_expD/run_expD_validation_only.jl
```

Render the grid figure and final report artifacts:

```powershell
python experiments/bnd_expD/plot_validation_grid.py
python experiments/bnd_expD/finalize_expD_results.py
```

`run_expD_validation_only.jl` checks the candidate's frozen flag, status, and SHA-256 before importing PowerDynamics. The validation surface uses 11 values of rho and three preregistered PLL bandwidth scales at each of buses 30, 33, 35, and 37. It records unstable cells rather than removing them. Bus 37 has nine intermediate-share spectral violations at rho 0.7–0.9, despite a feasible rho=1 endpoint.

## Main artifacts

- `ANALYTIC_CANDIDATE.json` and `ANALYTIC_CANDIDATE.sha256`: frozen analytical witness and hash.
- `tables/TABLE_D09_analytic_optimum.csv`: selected optima.
- `tables/TABLE_D10_multimode_check.csv`: complete analytic endpoint spectra.
- `tables/TABLE_D11_PD_validation.csv`: fresh endpoint-to-PowerDynamics match.
- `tables/TABLE_D17_PD_rho_PLL_grid.csv`: all post-freeze grid cells.
- `tables/TABLE_D19_PD_grid_full_spectra.csv`: every PowerDynamics pole at all grid cells.
- `tables/TABLE_D18_cross_bus_optima.csv`: final four-bus summary.
- `figures/FIG_D01_PD_validation_grid.png`: spectral-margin surface.
- `REPORT_EXP_D.md`, `RESULTS_EXP_D.json`, and `CLAIM_LEDGER_EXP_D.md`: final result and scope.

The local PLL and port structural-audit run writes its preliminary artifacts with `_STRUCTURE` suffixes so that it cannot overwrite the frozen candidate, full-spectrum tables, validation results, or final report.

## Model limits

The current-share model scales terminal current and retains each device's internal equations for interior rho; zero-share devices are removed at the endpoints. The model has no current limiter. The report therefore makes a frozen-operating-point small-signal claim and documents the bus-37 interior stability hole; it does not claim nonlinear or current-limited deployment performance.
