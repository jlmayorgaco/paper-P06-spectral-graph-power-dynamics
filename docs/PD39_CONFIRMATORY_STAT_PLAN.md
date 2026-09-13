# PD39 confirmatory statistical and numerical analysis plan

## Estimands

The 256 portfolios and 46 branches are finite exhaustive benchmark objects,
not random samples. Primary quantities are exact finite-benchmark values:
nominal alpha, `m_9`, blocker sets `H_0` and `H_0.05`, structured boundary
estimates under the frozen search protocol, repair effort, and deterministic
holdout coverage.

## Primary reporting

For every design and condition report equilibrium status, residual, voltage
feasibility, alpha, margin, spectrum status, and error provenance. Report
coverage, median, IQR, extrema, paired per-condition differences, and effect
sizes for the 24 designed holdout conditions. Do not turn coverage into a
population probability.

For node ranks (`n=8`), report exact score/rank tables, Spearman, Kendall,
top-3 overlap, and rank changes between portfolios. For links, report
Spearman, Kendall, top-5 and top-10 overlap, distinguishing all 46 branches
from the 11 outage cases that were not estimable in discovery. Use no p-value
as headline evidence for the eight-node ranking.

## Secondary inference

When useful, use paired permutation tests, sign tests, or a cluster bootstrap
over holdout conditions only as secondary descriptive context. Apply Holm
correction within any declared family of secondary tests. Do not treat the
exhaustive 256 portfolios as independent hypothesis tests.

## Radius inference

Each boundary must include a feasible/stable side and a violating side, the
final normalized bracket, exact perturbation, target eigenvalue, and
equilibrium residual. If no bracket is found in the frozen box, report
`not_found_in_box`; do not call it infinite or zero without evidence.

## Failure accounting

Separate graph disconnection, PF divergence, PF nonfinite state, dynamic
initialization failure, fixed-point failure, spectrum failure, true
instability, and low-margin robust failure. Never convert a failed solve into
an unstable label or assign a static penalty score.

