# E6 runner: distributed-gain nonlinear delayed events

Runner: `src/nhop_events.jl` (usage: `julia --project=. src/nhop_events.jl CASE.toml LABEL [event indices]`). Read-only dependencies (`DelayedEvents.jl`, `ReducedDAE.jl`) are included, not modified.

## Design
- CASE.toml: `rho` (10), `tau` (scalar or 10), `Kp`, `KI` (10x10, row i = receiving PLL, column j = measured detector).
- Model built with the diagonal gains `kp = diag(Kp)`, `ki = diag(KI)`. The delayed injection is `dx += Bdiag*(ed - e(x)) + Boff*ed`, with `Bdiag = DE.injection(m)` and `Boff[:,j] = sum_{i!=j} e_omega_i Kp[i,j]/pll_tau_i + e_xi_i KI[i,j]`.
- State indices confirmed against `DE.injection`: omega = `m.gfidx[i][4]` (rhs `(xi+Kp e-omega)/pll_tau`), xi = `m.gfidx[i][5]` (rhs `Ki e`).
- Jacobian `J = Fx - Bdiag*C`; delayed terms are history (as in the existing runner).
- Solver: method of steps, Rodas5P, reltol = abstol = 1e-9, dtmax 0.01, 60 s, metrics via `R.metrics` (window 0.5 s), pass guards identical to `events.jl`. Steps have length min(tau_j); per-detector delays use one detector evaluation per distinct delay. Event 1..5 = `DE.CASES`.
- Outputs: `raw/events/events_LABEL.csv` (label, tau_ms, bus, delta, F, R, Vmin, Vmax, slack, runtime, pass, error; tau_ms = 1000*max tau) and `raw/events/traj/LABEL_bus<b>_<delta>_trajectory.csv`.
- Fix during development: the first draft assigned dual numbers into a Float64 vector (ForwardDiff time-gradient failure in Rodas5P); replaced by a type-promoting vector.

## Test 1: parity (Kp = 28.2743 I, KI = 246.7401 I, rho 0.875, tau 0.040 s)
Reference: `research_gold/checks/T9_FBK_base40.csv` (relative difference shown; run `raw/events/events_parity.csv`).

| event | quantity | gold | runner | rel. diff |
|---|---|---|---|---|
| 1 (bus 8, -100) | F | 0.4615490788454758 | 0.4615490788454758 | 0 |
| | R | 0.1382758495928848 | 0.1382758495928848 | 0 |
| | Vmin | 0.9874643533937851 | 0.9874643533937851 | 0 |
| | Vmax | 1.0777110517139958 | 1.0777110517139958 | 0 |
| | slack | 0.09790327434111055 | 0.09790327434111055 | 0 |
| 2 (bus 16, +100) | F | 0.4438596322629803 | 0.4438596322629803 | 0 |
| | R | 0.126035744299283 | 0.126035744299283 | 0 |
| | Vmin | 0.9721083659192574 | 0.9721083659192574 | 0 |
| | Vmax | 1.0612480273995026 | 1.0612480273995026 | 0 |
| | slack | 0.007426701452854956 | 0.007426701452854956 | 0 |

PASS (bitwise identical printed values, well within 1e-6).

## Test 2: sanity (Kp[1,2] = Kp[2,1] = 2.0, event 1)
F = 0.461498, R = 0.138290, Vmin = 0.987445, Vmax = 1.077664, slack = 0.097919, pass = true, runtime 223 s. Finite; differs from the parity run at the 1e-4 relative level (F: +2.5e-5 rel.), consistent with a small off-diagonal coupling. PASS.

Runtime: about 75-235 s per event with two runs in parallel.
