# Experiment status and gates

| Gate | Result | Evidence |
|---|---|---|
| L0 baseline parity | PASS | `L0_REPRODUCTION.csv`; `L0_REPLACEMENT_ONLY_REPRODUCTION.csv` |
| L1 fixed-ρ gain action | PASS (exact identity, Float64 check) | `L1_ACTION_SPACE_VALIDATION.csv`: rank 10 in 48 evaluations; max relative error 4.72e-15 |
| L2 simple-root gradient | PASS for tested smooth cases | 9/9 signs, median relative error 1.75e-05; near-switch envelope nonsmooth in `L2_INTERMEDIATE_GRADIENT_VALIDATION.csv` |
| L3 two-start co-design | EXECUTED, best found; not converged | `L3_OPTIMIZATION_TRACE.csv`; `N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp` passes five events; prescribed <0.01 ms/KKT stopping gate not reached |
| L4 family switch | SUPPORTED NUMERICALLY | `L4_MODE_SWITCH_POINT.csv`; η≈0.880008, MAC≈0.0041 |
| Complete numerical contour around best threshold | PASS | `evaluations/N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp/ROOT_COUNTS.csv`; coverage: COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR |
| Positive-delay nonlinear DDE events | NOT RUN | No delayed nonlinear simulator in this experiment |
| Global optimum/interval root certificate | NOT AVAILABLE | No validated global bound or interval arithmetic |

The frequency and actuator event limits were frozen in `EXPERIMENT_MANIFEST.json`. The rejected full second step raised the spectral threshold but failed the bus16 +100 MW actuator guard (`event_screens/N_step1_multimode_lp_next_full/`). The winning design has only 2.88e-06 actuator-fraction slack above that guard, so it is sensitive to numerical/model uncertainty.
The next frozen-trust LP predictor, if present, is **unvalidated** and is excluded from the best-found result. Its predicted threshold is 41.05730843137167 ms; see `PENDING_CANDIDATES.md`.
