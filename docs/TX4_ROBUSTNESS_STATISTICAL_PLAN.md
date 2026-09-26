# TX4 Robustness Statistical Plan

## Estimands

The primary estimands are the probability of `H4_PRESENT`, the probability of
`EXACT_H4`, the probability of any `NONCOMPOSABLE` outcome, and the continuous
distribution of `g_star`, defined as the first sign-change boundary in the
H4 full-spectrum alpha along the g coordinate with all other coordinates fixed
at their sampled condition. Secondary estimands are `delta_H4`, `eta_H4`,
`ell_H4`, `q_H4`, return distance, and mode-scoped alpha.

## Sampling

The deterministic campaign maps structure and phase boundaries. Scrambled
Sobol samples quantify low-discrepancy coverage of the declared box. The Monte
Carlo sample uses independent uniforms because no calibrated physical weights
were found in the prior audit. Results are therefore bounded-box engineering
probabilities, not population probabilities.

## Sensitivity

Morris elementary effects rank all seven coordinates. A Saltelli design then
estimates first-order and total-order Sobol indices for `H4_PRESENT` and
`EXACT_H4`. Logistic models with preregistered main effects and pairwise
interactions estimate odds of those same binary endpoints; separation and
nonconvergence are reported rather than repaired by changing the model.

## Return and family robustness

Return robustness compares full-spectrum and mode-scoped decisions on paired
conditions. Mode-family robustness repeats the exact-H4 decision in the
0.3–1.5 Hz band and in the full spectrum. TDS robustness reports agreement,
trace failure, and sign consistency by condition stratum. All summaries retain
the portfolio identity and condition identifier so no portfolio-specific
refitting can be hidden by aggregation.
