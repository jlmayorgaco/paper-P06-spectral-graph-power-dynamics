# TX4 Robustness Statistics Audit - Corrected Final Report

Date: 2026-09-20  
Branch: `research/tx4-robustness-statistics-audit`  
Frozen exact parent: `f64db0004026ceafdb08dd13b5e2ff59d6060742`

## Executive verdict

The retained ExtraTrees tier failed the prespecified strict surrogate gate
(false-minimal rows=18, missed-minimal rows=14, exact-H4 agreement=0.924883, H0 antichain agreement=0.701878). Primary QMC and MC blocker statistics therefore use the complete exact fallback: 4,096 QMC conditions and 5,000 MC conditions, all 16 portfolios, 145,536 exact portfolio evaluations.

The legacy H4_PRESENT flag had a complement/minimality bug: in the exact calibration it marked H4 unstable rather than requiring H4 to be inclusion-minimal. 156 of 426 exact calibration conditions were legacy-true but corrected-H0-false.

## Corrected definitions

For each condition, U0 contains portfolios with alpha_EM >= 0 in the frozen
0.3-1.5 Hz EM band. H0 is the inclusion-minimal antichain of U0. H4_UNSTABLE,
H4_PRESENT, EXACT_H4, NONCOMPOSABLE, and KAPPA are computed from U0/H0. The
full-set relation makes H4_PRESENT and EXACT_H4 numerically identical here:
any proper unstable blocker prevents H4 from being minimal.

## Primary exact results

| Endpoint | QMC coverage | MC assumed-box probability |
|---|---:|---:|
| H4_UNSTABLE | 0.577881 | 0.581400 |
| H4_PRESENT | 0.113770 | 0.115400 |
| EXACT_H4 | 0.113770 | 0.115400 |
| NONCOMPOSABLE | 0.427246 | 0.436600 |

QMC values are fixed-design coverage fractions, not probability claims. MC
values use independent uniform draws over the declared bounded engineering box;
Wilson intervals are in `TX4_MC_EXACT_BLOCKER_SUMMARY.csv`.

## Continuous diagnostics

The true-minimality margin is `min(alpha_EM(H4), -max proper alpha_EM)`. The
QMC median is -0.085269 s^-1 and the MC
median is -0.084810 s^-1, with the MC BCa
95% interval recorded in `TX4_MC_DELTA_H4_TRUE_MINIMALITY.csv`. The exact
128-condition eta stratum completed with zero failures and nonnegative
`min_i sigma_min(I+Q_(H4\i)(s_H))` values.

## Provenance and sensitivity

All 329,440 legacy master rows now carry `evaluation_source`: EXACT_DAE,
SURROGATE_EXTRATREES, DERIVED_FROM_EXACT, or UNKNOWN. The exact QMC/MC
fallback is labeled EXACT_DAE and is the only source used for primary blocker
statistics. Legacy Morris remains labeled approximate screening and the
existing Sobol table remains explicitly surrogate-based; neither is promoted
to exact causal sensitivity evidence. G-star values remain local exact-H4
phase-boundary diagnostics, not physical robust radii.

No U005/H005 threshold was found in the retained preregistration or TX4 code;
those endpoints are not analyzed and no threshold was invented.

## Verification outputs

- `results/TX4_QMC_ALL16_EXACT.parquet` and `results/TX4_MC_ALL16_EXACT.parquet`
- `results/TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv`
- `results/TX4_SURROGATE_MINIMALITY_VALIDATION.csv`
- `results/TX4_AUDIT_H4_FLAG_COMPARISON.csv`
- `results/TX4_AUDIT_ETA_H4.csv` and `results/TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv`
- `results/TX4_ROBUSTNESS_INVARIANT_CHECKS.csv`
- `results/TX4_AUDIT_FIGURES/F1_...png` through `F10_...png`

All final invariant checks must be PASS before the corrected bundle is called
complete.
