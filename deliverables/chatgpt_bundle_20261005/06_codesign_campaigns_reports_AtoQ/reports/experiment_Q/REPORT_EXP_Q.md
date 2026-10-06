# Experiment Q — secure SG-to-GFL co-design

## Status

**PASS_SECURE_GAP_OPEN.** A frozen all-SG reference passes the declared 100 MW bus-16 frequency/RoCoF limits, the complete-spectrum margin, the P/Q equilibrium checks, and the normalized robustness requirement under a Float64 bounded-real CARE/LMI certificate. It retains all `5402.761089978847 MW` of original generator dispatch. This is a feasible reference point, not an optimized or locally KKT-certified design. The minimum secure retained-SG dispatch is not determined.

The rigorous objective bounds currently available are `0 MW ≤ J* ≤ 5402.761089978847 MW`; the `5402.761 MW` gap is open. The lower bound is only the nonnegativity of `Σ P_i ε_i`. The frequency-only knapsack value is not a valid bound because the exact steady-frequency response is nonaffine.

## Frozen baseline and frequency outputs

Q0 reproduced the ExpN/ExpP baseline and verified both frozen candidate hashes. The declared P/Q split and separate buses 31/39 loads passed to machine precision. The ExpN nominal point retains `1.124344477653 MW`, with analytical/PD alpha `-0.050000001105 / -0.050000000139 s⁻¹`. Its `ω=0` robustness witness is `3.74649e-11`, below `beta_req=1.69912e-6`. Its independent 100 MW bus-16 PD TDS peaks at `39.3334 Hz` and `1.07221 Hz/s` and does not settle in 60 s.

Q1A compared SG rotor frequency, GFL PLL frequency, and filtered bus-voltage frequency. The metric is `METRIC_DEPENDENT`: on the 100 MW 0.1 s pulse, the SG-only RoCoF peak passes 0.5 Hz/s while the local PLL and bus-frequency diagnostics exceed it. Q therefore freezes the prescribed conservative primary output as the maximum absolute local SG rotor or GFL PLL frequency state over buses 30–39; the bus-frequency estimator remains a diagnostic. It is an 11-sample cubic Savitzky–Golay differentiator at 100 Hz, with 0.1 s support and measured derivative bandwidth 11.7069 Hz.

## Linear model and steady response

The reduced load-to-frequency model reuses 69 immutable independent ExpN PD matrix/pole identity cases under the exact same model SHA. Their maximum reduced-A relative error is `3.20e-15`; maximum matched-pole error is `5.56e-8 s⁻¹`. Modal residues reconstruct direct matrix-exponential step and RoCoF outputs with maximum absolute errors `2.34e-9` and `1.63e-10` in the unit response.

Against independent ExpP PD traces for the same 15 s windows, the 1 MW and 10 MW primary frequency trajectories have relative RMSE `0.070%` and `0.632%`; peak errors are `0.118%` and `1.037%`. The 100 MW nonlinear trajectory departs substantially: its 60 s maximum is `39.3335 Hz`, while the linear trajectory over that window peaks at `21.1647 Hz`. Large-signal behavior is not inferred from the linear response.

Twenty centered gain perturbations give `max |dF∞/dK|=1.49e-12 Hz/(MW·gain-unit)`. Three controlled comparisons show Kp/Ki affect poles and transient quantities; the recorded beta values are sampled upper witnesses only. A 100-point exact-model test rejects `Kload+Σc_i ε_i` as an affine representation: inferred `Kload` ranges from `-3995.97` to `1438.82 MW/Hz`, relative range `43.93%`. Hence no frequency-only analytical global lower bound is claimed. The bus-30 `57.21 MW` knapsack result is marked `SCREENING_ONLY`.

## Secure feasibility reference

The pre-validation candidate is [Z_Q_ALL_SG_SECURE_REFERENCE.toml](Q2_REFERENCE/Z_Q_ALL_SG_SECURE_REFERENCE.toml), SHA-256 `f3286c45e32f0a69b8d0af62461e07b116121365c007161199741419fa80cb74`. It uses `rho=0` at every generator; Kp/Ki are inactive. It was frozen before the independent Q PowerDynamics/TDS validation.

- Complete physical spectrum: 113 poles, alpha `-0.098065409336 s⁻¹`; immutable ExpN all-SG PD identity reports max pole difference `2.31e-12 s⁻¹`.
- Robustness: Float64 bounded-real certificate at `beta_req`, beta lower bound `1.699120699918e-6`; CARE relative residual `1.86e-6`, `P_min=0.0024073`, CARE/LMI max eigenvalue `-0.08404`, closed-loop alpha `-0.04807 s⁻¹`. No directed rounding was used.
- Linear +100 MW bus-16 metrics: `F∞=0.03756 Hz`, `Fpeak=0.07581 Hz`, `Rpeak=0.12532 Hz/s`.
- Independent PD TDS: local SG frequency peak `0.075365 Hz` (bus 39), local SG RoCoF `0.123835 Hz/s` (bus 33), filtered bus-frequency peak `0.074642 Hz`; the event passes 0.5 Hz and 0.5 Hz/s. Trim residual `2.66e-13`; maximum P/Q error `0 pu`.

This reference establishes feasibility, not minimum dispatch. The ExpG coordinate point, re-evaluated using the ExpN PD-exact model, was rejected: its alpha `-0.13793 s⁻¹` and `Fpeak=0.09355 Hz` pass, but local PLL RoCoF is `49.12 Hz/s` at bus 36 and its CARE check fails. It was not frozen as a Q candidate and PD was not called on it.

The one-sided local PLL test covers both ends of each single-bus mixed branch. With `rho=1e-6` (GFL insertion) the minimum by bus is `1.013 Hz/s` at bus 39; with `rho=0.999999` (`epsilon=1e-6`, SG removal) it is `9.451 Hz/s` at bus 37. All 80 tested bus/boundary/gain-corner combinations exceed `0.5 Hz/s` for the 100 MW initial PLL RoCoF. The instantaneous value is proportional to Kp and independent of Ki. This is evidence at local branch boundaries, not a proof that all mixed interiors fail.

## Optimization and unfinished work

No Q active-set SQP/KKT continuation, robust disturbance continuation, support enumeration, branch certification, optimality multipliers, BND self-energy decomposition at a secure optimum, or GSP compression solve was run. The ExpG seed fails Q constraints and is not used as an incumbent. The only Q result is the frozen all-SG feasibility reference and the explicitly non-global gap above.

No `J*(ΔP)` curve or optimal `rho*/Kp*/Ki*` exists. Q shadow prices, a global architecture certificate, and poster figures that depend on a secure optimum are not claimed. All PLLs remain local; no GFM, droop, virtual inertia, controller architecture, event, gains, or dependencies were changed.

Two evidence-limited figures are included: [FIG_Q1](figures/FIG_Q1_stable_not_secure.png) compares ExpN nominal stability and 100 MW security metrics against the all-SG feasible reference; [FIG_Q4](figures/FIG_Q4_all_sg_linear_vs_PD.png) compares the linear all-SG event trajectory against independent PD TDS. No replacement frontier, authority map at an optimum, or BND/GSP optimum figure is presented because the corresponding design was not obtained.

## Reproducibility and runtime

Environment: Julia 1.11.9, PowerDynamics 5.0.0, NetworkDynamics 1.3.0; Project and Manifest hashes match the frozen ExpN environment. The Q all-SG reference run records total wall time `171.968 s`, including the numerical CARE test; its PD TDS runtime is `45.664 s`. Q0 used one independent PD rebuild and one 60 s event simulation; Q1A used four PD TDS runs; Q1B–Q1G used 69 inherited PD identity cases, 20 centered gain derivatives, 100 affine-structure points, and three dynamic gain comparisons. Q1 stage runner wall time was not instrumented. The separate ExpG seed test used one exact-model evaluation and one failed CARE calculation, with no PD call.

This branch is `research/expQ-secure-gfl-codesign`, based on parent `d0fecb3264aeb855ab6700fc6ef66120d4b6c33c`. No commit or push was made. ExpN/ExpP frozen candidates and dependencies remain unchanged.
