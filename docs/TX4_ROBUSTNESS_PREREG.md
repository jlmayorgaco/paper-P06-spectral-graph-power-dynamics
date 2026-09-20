# TX4 Robustness Preregistration

Date: 2026-09-20  
Branch: `research/tx4-final-robustness-statistics`  
Parent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`

## Frozen model and policy

The campaign evaluates the frozen matched-Q reduced semi-explicit IEEE-39 DAE
from the exact P4/GFL11 closure. The target family is buses 30, 33, 35, and 37;
all 16 subsets, including the empty set and H4, are evaluated for every
condition. Equilibrium is solved before central-difference DAE linearization;
the reduced state matrix is formed by the same index-one Schur reduction used
by the exact closure.

No result-dependent parameter refitting, portfolio-specific uncertainty set,
or closure-space radius is allowed.

## Frozen parameter domain

The campaign parameter vector is
`(g,k,t,h,epsilon,damping,inertia)`:

| parameter | meaning | closed interval | implementation |
|---|---|---:|---|
| `g` | global converter voltage-control gain | [0.020, 0.250] | converter `voltage_gain`; leak fixed at 0.05 |
| `k` | synchronous-machine `ka` multiplier | [0.75, 1.75] | fleet `machine_scaling["ka"]` |
| `t` | synchronous-machine `ta` multiplier | [0.75, 1.75] | fleet `machine_scaling["ta"]` |
| `h` | bus-37 converter gain heterogeneity multiplier | [0.50, 1.50] | bus-37 converter gain is `g*h` when replaced |
| `epsilon` | common active/reactive load multiplier | [0.90, 1.10] | frozen network loads scaled together |
| `damping` | additive synchronous damping | [0.00, 0.20] | `machine_damping` |
| `inertia` | surviving-machine inertia multiplier | [0.80, 1.20] | `machine_services["inertia"]` |

The intervals are bounded symmetric engineering assumptions because the prior
manifest declared no calibrated physical weights. They are not claimed to be a
certified physical uncertainty envelope.

## Seeds and sample sizes

- Deterministic one-dimensional sweeps: 201 equally spaced points per primary
  coordinate.
- Deterministic two-dimensional maps: 41 x 41 points for `g x k`, `g x t`,
  `g x h`, and the preregistered load-stress pairs `g x epsilon`, `k x epsilon`,
  `t x epsilon`.
- Scrambled Sobol QMC: 4096 points, seed 20260920.
- Monte Carlo: 5000 points, NumPy PCG64 seed 20260921, independent uniforms on
  each bounded interval.
- Morris screening: 40 trajectories, 8-level grid, seed 20260922; the six
  largest absolute elementary effects are retained for the Saltelli follow-up.
- Saltelli/Sobol follow-up: base size 1024, seed 20260923.
- Julia cross-code random spot-check: 32 conditions, seed 20260924, four
  preregistered portfolios per condition (empty, each single bus, H4, and one
  proper triple, selected deterministically from the target family).

## Modal definitions

- `alpha_all`: maximum real part of the complete reduced spectrum.
- `alpha_EM`: maximum real part among non-gauge modes with frequency in
  0.3–1.5 Hz and magnitude at least `1e-4`.
- `H4_PRESENT`: H4 has `alpha_all > 1e-8`.
- `EXACT_H4`: H4 is present and all 15 proper subsets have `alpha_all <= 1e-8`.
- `NONCOMPOSABLE`: any equilibrium/linearization failure, non-finite metric,
  or rating/power-flow infeasibility.
- `delta_H4`: H4 `alpha_all` minus the maximum proper-subset `alpha_all`.
- `eta_H4`: negative of the maximum proper-subset `alpha_all` (positive means
  the proper-subset family is stable under the threshold).
- `ell_H4`: minimum absolute real part among H4 non-gauge modes (a local modal
  margin, reported only as a numerical diagnostic).
- `q_H4`: minimum distance of a H4 non-gauge eigenvalue to the imaginary axis;
  it is a closure-distance diagnostic, not a physical robust radius.

The primary blocker threshold is `1e-8 s^-1`; a sensitivity band of
`1e-6 s^-1` is reported separately and is never silently substituted for the
primary verdict.

## Multiple testing and stopping

All declared endpoints are computed for every completed condition. Confidence
intervals use Wilson intervals for proportions and percentile bootstrap with
10,000 resamples for continuous summaries. Pairwise method comparisons use
paired McNemar tests with Holm correction across the two primary method labels
(full-spectrum and 0.3–1.5 Hz mode-scoped). No post-hoc stopping, boundary
selection, or interval retuning is permitted. Failures are retained as
`NONCOMPOSABLE`; they are not discarded.

## TDS and topology scope

The TDS section is a trace-level robustness audit over approximately 24
condition strata, up to 72 traces, using the existing reduced-D​​AE TDS
integrator and fixed kick/pulse definitions. It is not a new EMT or switching
simulation. No line/topology perturbation is introduced beyond the prior
audited continuous network contract.
