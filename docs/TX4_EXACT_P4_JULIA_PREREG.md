# TX4 Exact P4 / GFL11 Julia Reproduction - Preregistration

Preregistered before the exact P4/GFL11 run. No parameter fitting, retuning, or post hoc model edits are permitted after the first exact-case execution.

## Question

Can the exact frozen 11-state TX4 GFL and exact P4 surviving-SG policy be independently reimplemented in Julia and reproduce the frozen Python IEEE-39 V4 results?

## Frozen inputs

- Network: frozen IEEE-39 data and admittance snapshot already used by the independent validation harness.
- Candidate replacement buses: `30, 33, 35, 37`.
- Census: all 16 subsets of those four buses, including base and H4.
- GFL: exact 11-state controller, including the Q/V integrator state `xv`.
- GFL policy: `g=0.03625`, `kp_v=2`, `ki_v=20`, `leak=0.05`.
- SG policy: surviving AVR gain multiplied by `k=1.425`; surviving AVR time constant multiplied by `t=1.5`; heterogeneity exponent `h=1.0`.
- Reactive policy: matched P/Q replacement, so the AC equilibrium is shared across the census.
- Frequency: 60 Hz base; reported critical frequency is `abs(imag(lambda))/(2*pi)`.

## Acceptance gates

| Gate | Criterion |
|---|---|
| G1 | Every replaced GFL contributes exactly 11 dynamic states. |
| G2 | H4 has `nx=86`. |
| G3 | Python and Julia stable/unstable verdicts match for all 16 cases. |
| G4 | `abs(alpha_py-alpha_jl) <= 1e-4 s^-1` for H4. |
| G5 | `abs(f_py-f_jl) <= 1e-4 Hz` for H4. |
| G6 | Voltage-mode MAC after declared matching is at least 0.95 for H4. |
| G7 | Equilibrium voltage magnitude and angle differences are each at most `1e-6`. |
| G8 | Equation audit, state ordering, policy audit, and source-manifest checks report no hidden mismatch. |

## Primary case and proper-subset rule

The primary case is H4=`30+33+35+37`. All 15 proper subsets must be stable in both implementations. A case is stable when the maximum transverse real part is strictly negative after the common gauge tolerance `1e-4`.

## Failure handling

If any gate fails, the result is reported as Case B (cross-code implementation mismatch) unless the exact run is blocked before a result exists, in which case it is Case C. No tuning or second campaign is authorized by this preregistration.
