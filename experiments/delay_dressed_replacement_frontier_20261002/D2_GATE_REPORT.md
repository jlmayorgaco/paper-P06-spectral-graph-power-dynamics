# D2 eigenvalue-sensitivity validation

**Overall status: BLOCKED_PARAMETER_SET_INCOMPLETE.**

For the dominant locally tracked conjugate pair, the exact simple-root sensitivity formula was checked for τ, Kp and Ki at three frozen delay patterns and buses 30, 32 and 33. All 27 centered finite-difference comparisons passed the 1% relative-error target. Maximum relative errors were 6.74e-7 for τ, 2.27e-5 for Kp and 1.03e-4 for Ki. The corresponding maximum absolute errors are in `baseline_reproduction/TABLE_D02_DERIVATIVE_VALIDATION.csv`.

This is **SUPPORTED_LOCAL / NUMERICALLY_VALIDATED** only for the two tracked roots and stated parameter points. The ρ derivative was not derived or validated; the full sensitivity gate therefore remains blocked. DDE root coverage is also unresolved, so these sensitivities do not define a complete stability frontier.
