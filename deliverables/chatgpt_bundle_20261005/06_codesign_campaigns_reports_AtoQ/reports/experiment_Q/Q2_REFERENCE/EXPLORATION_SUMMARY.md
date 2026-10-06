# Q2 continuation diagnostics — exploratory evidence only

This supplemental stage checks two useful points after the all-SG feasibility reference was frozen. It does not optimize the design and does not certify architecture infeasibility.

Wall time was not instrumented for the 80 boundary evaluations or the separate ExpG seed re-evaluation.

## One-sided local PLL RoCoF

Equation checked: immediately after the algebraic voltage response, the local PLL equation gives `RoCoF_PLL(0+) = Kp e_theta/(2*pi*tau_PLL)`. The gain-corner table has 80 evaluations: one tested generator at a time with all other generators all-SG, both `rho_GFL=1e-6` (GFL insertion) and `rho_GFL=0.999999`, `epsilon_SG=1e-6` (SG removal), and four Kp/Ki box corners. Every evaluated case exceeds the 0.5 Hz/s limit for the 100 MW step. Across the small-GFL-insertion boundary, values range from `1.0133` to `187.3196 Hz/s` (minimum at bus 39 and Kp's lower bound); across the SG-removal boundary they range from `9.4513` to `350.7038 Hz/s` (minimum at bus 37 and Kp's lower bound). The values are exactly independent of Ki across paired corners and scale by 16 between Kp bounds.

These are single-bus continuation boundary checks, not a support-wide result. The primary metric counts every present PLL state, so a tiny GFL share still creates a measured frequency output. The induced algebraic voltage change can differ when several generator buses are mixed; the result cannot prune all mixed supports.

## ExpG seed re-evaluation

The frozen ExpG candidate was re-evaluated with the ExpN exact analytical model, without modifying its frozen file and without a PowerDynamics call. At the ExpG coordinate point, retained SG is `3821.7712 MW`, alpha is `-0.13793 s^-1`, and the frequency peak is `0.09355 Hz`; however, the primary local PLL RoCoF is `49.1212 Hz/s` at bus 36 and the Float64 CARE check fails. It is therefore not a secure Q incumbent. This was one exact-model point evaluation; elapsed wall time was not instrumented.

## Continuation decision

These diagnostics permit prioritizing both GFL-insertion and SG-removal boundaries and reject the tested ExpG point under the Q frequency metric. They do not establish a local optimum, any support-wide infeasibility result, or a global lower bound. The support regions capable of beating the all-SG feasible reference remain **OPEN**.

Reproduce with `experiments/bnd_expQ/audit_one_sided_pll_rocof.jl` and `experiments/bnd_expQ/evaluate_expG_seed_q.jl`.
