# Executed status and gates

| Stage | Status | Evidence | Boundary |
|---|---|---|---|
| M0 baseline | PASS, numerically reproduced | `M0_BASELINE_REPRODUCTION.csv`, `M0_EVENT_REPRODUCTION.csv` | 87.5% seed, all five frozen zero-delay events, 203 physical poles |
| M1 zero-delay search | BEST FOUND, not optimum | `M1_CANDIDATE_AUDIT.csv`, `M1_ZERO_DELAY_DESIGN.toml` | 88.4551408% and 623.7411616 MW retained SG; frozen interpolation plus one midpoint only |
| M2 exact action space | EXACT_IDENTITY algebraically; NUMERICALLY_VALIDATED in accepted samples | `M2_ACTION_SPACE_RANK.csv`, `M2_RECONSTRUCTION_VALIDATION.csv` | Rank at most 30 on a 282-dimensional descriptor; 20 random in-box samples are **not** event-validated. 20 out-of-gain-box samples were rejected by the model API, so the requested 20 feasible/20 infeasible physical validation set was not achieved. |
| M3 DDE oracle | NUMERICALLY_VALIDATED for two fixed designs at uniform 20/30/40 ms; not certified | `M3_DDE_ORACLE.csv`, `M3_ALL_REFINED_ROOTS.csv`, `M3_ORACLE_METHOD.md` | No arbitrary-precision/interval count and no heterogeneous-delay candidate oracle tested here. At 40 ms, 6 vs 12 roots right of the margin. |
| M4 finite-frequency graph damping | BLOCKED | `M0_ZERO_FREQUENCY_SCHUR_REPRODUCTION.csv`, theory §11 | The hidden Schur block at zero frequency has condition `3.095842e18` after gauge deflation. A physical finite-frequency angular-port reduction was not validated. |
| M5 conditional PI law | EXACT conditional identity; NUMERICALLY_VALIDATED target-root residual | `M5_CLOSED_FORM_PI_VALIDATION.csv` | Only 5/45 prescribed target roots have gains within frozen bounds; no nearest/rightmost-pole or reduced direct-lift prediction demonstrated. |
| M6–M8 analytical authority / predictor | BLOCKED | status-bearing `M6`–`M8` tables and theory | LP two-anchor theorem is reduced-model mathematics; its physical coefficients and an accurate full-model predictor are unavailable. M1 is not an analytical predictor. |
| M9–M14 delay frontier / graph placement | BLOCKED | status-bearing `M9`–`M13` tables | No delay-specific co-design and no full nonlinear DDE event solver. Fixed-design pole counts are not replacement capacities. |
| M15 comparison | PARTIAL | `M15_BASELINE_COMPARISON.csv` | Reproduced seed and best M1 point vs failing historical joint candidate; not a common optimum benchmark. |
| M16 nonlinear validation | PASS at zero delay only | `M16_FULL_VALIDATION.csv` | Five complete zero-delay events for the best found design; no method-of-steps delayed events. |
| M17 bounds | BLOCKED | `M17_OPTIMALITY_BOUNDS.csv` | Best zero-delay witness gives a feasible lower bound on maximum replacement only. No nontrivial certified upper bound or gap. |

The requested positive-delay maximum `P_GFL^max(tau)`, optimal `rho/Kp/Ki` and spatial-latency frontier are **not determined**. The zero-delay best found vector is not denoted an optimum. No previous frozen experiment was edited. No commit or push was made.
