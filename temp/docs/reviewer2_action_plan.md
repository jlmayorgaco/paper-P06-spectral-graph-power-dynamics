# Reviewer-2 Action Plan: Framing, Novelty, and Critical Path

Date: 2026-06-08

This note records the action taken after the hard Reviewer-2 audit. It is meant
to be the internal source of truth for the next paper iteration.

## Decision

Reviewer 2 is right on the main point: the current paper should not be framed as
new mathematics or as an estimator paper before Experiment A is run. The
validated full-ANDES result is the controller-aware rational bridge and the
mechanism audit it enables. The second-order correction remains important, but
it is not yet the empirical headline.

The manuscript framing has therefore been shifted from:

> second-order rational graph-filter perturbation as the central result

to:

> static scalar damping bridges fail in IBR cases; the controller-aware rational
> self-energy bridge preserves controller-mediated damping and exposes
> network-family versus control-family mechanisms.

## What Is Not Original

Do not defend these as new:

- proportional modal formula `zeta_k = delta_k/(2 sqrt(nu_k))`,
- QEP formulation of second-order dynamics,
- second-order eigenvalue perturbation,
- Schur-complement condensation,
- rational self-energy as a general mathematical object,
- Beyn contour integration for nonlinear eigenvalue problems,
- generic damping-ratio sensitivity or inverter-placement by damping metrics.

These are established methods. They can be used, but they cannot be sold as the
paper's mathematical novelty.

## Literature Boundary Already Confirmed

Targeted web/literature lookup confirms the audit's caution:

- Contour integral and Beyn methods for NEPs are established, e.g.
  `https://arxiv.org/abs/1003.1580`,
  `https://epubs.siam.org/doi/10.1137/20M1389303`, and
  `https://www.sciencedirect.com/science/article/pii/S037704271500374X`.
- Power-system dynamic/frequency-dependent equivalents are established, e.g.
  `https://www.mdpi.com/1996-1073/15/4/1396` and related FDNE literature.

Implication: the novelty is not "we invented NEPs/self-energy/contour
integration." The defensible claim is a power-system formulation and audit:
controller-aware self-energy as the damping-margin bridge for IBR cases where
static nodal damping fails.

## What Remains Defensible

### Strongest verified contribution

The scalar bridge failure and rational bridge contrast:

- In the ANDES cases, physical `GENROU.D` is zero or insufficient to represent
  observed damping.
- A static `(L,M,D)` bridge therefore fails as a damping-margin object.
- The Schur NEP with controller self-energy retains the condensed control
  dynamics and recovers critical poles, including a Mix60 control-family mode
  missed by scalar and fixed-point reductions.

This is not an approximate validation of a new model against ANDES. It is an
exact Schur representation used to show completeness and to contrast the failure
of the static scalar bridge.

### Most interesting mechanism

The control-resonant reinforcement reversal:

- A line reinforcement can move a network-family mode toward a PLL-dominated
  condensed-control pole.
- The self-energy term can dominate the damping decrease.
- The mechanism disappears under the soft-GFL negative control.

Current limitation: it is modal-local and appears in a minority of audit records.
It does not yet produce a global-margin-collapse headline.

### Supporting theory

The second-order correction and rational graph-filter interpretation are
supporting theory:

- useful as an interpretable damping-margin decomposition,
- verified in synthetic and reduced PLL tests,
- not yet validated as a full-ANDES estimator,
- must be demoted if Experiment A shows diagonal all-mode screening already
  performs as well as the adaptive method.

## Manuscript Changes Made

The paper package has been reframed conservatively:

- title changed to `Controller-Aware Rational Self-Energy for Damping-Margin
  Analysis of IBR Grids`;
- abstract now opens with static bridge failure and rational self-energy;
- contribution order now puts the controller-aware bridge before the
  second-order estimator;
- "resonant Braess" has been renamed in the main claim path to
  "control-resonant reinforcement reversal";
- `paper_ieee_transactions/README.md` now states that the estimator is not the
  headline until Experiment A is run.

## Critical Path: Experiment A

Experiment A decides what paper this becomes.

Question:

> On the validated controller-aware bridge/full ANDES poles, does the adaptive
> estimator actually improve over the diagonal all-mode screen and registered
> baselines?

Required comparisons:

- full ANDES or exact Schur-NEP pole margin as reference;
- homogeneous/proportional formula baseline;
- Fiedler-only baseline;
- diagonal all-mode baseline;
- second-order correction;
- reduced-QEP fallback where applicable;
- adaptive workflow.

Required metrics:

- median error,
- p95 error,
- failure cases,
- diagnostic trigger precision/recall,
- whether the estimator changes the practical decision relative to the diagonal
  baseline.

Decision rule:

- If adaptive improves materially where the diagnostic triggers, an
  estimator-centered Transactions path remains possible.
- If diagonal all-mode is already enough, the estimator becomes an
  interpretation/fallback and the paper should be framed as a bridge/mechanism
  paper.
- If neither bridge-enabled estimator nor mechanism produces a strong empirical
  advantage, target a conference or thesis-method paper, not IEEE Transactions.

## Next Implementation Step

Do not polish more theory before Experiment A. Implement or complete:

`validation/ieee39_rational_filter/phaseA_estimator_over_nep.py`

The script should consume the Phase-0D bridge artifacts and produce:

- estimator comparison table,
- diagnostic-trigger table,
- pole/mode identity audit,
- `phaseA_status.json`,
- one figure showing error by method and case.

Until that exists, the manuscript should not claim final estimator or planning
validation.
