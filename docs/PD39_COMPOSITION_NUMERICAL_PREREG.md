# PD39 follow-up preregistration: penetration versus composition and numerical audit

Status: frozen before the follow-up numerical execution.

This is a short confirmatory follow-up to the frozen discovery campaign at
commit `61554336` and the confirmatory campaign at commit `902403cf`.  It does
not alter, overwrite, or reinterpret either campaign.

## Scope and frozen question

The primary question is whether the observed 7/8 to 8/8 cliff can be
explained by penetration and remaining synchronous resources alone, or whether
portfolio composition still changes the dynamic result after those quantities
are closely matched.

The same stock IEEE-39 PowerDynamics model, GFL model, dispatch convention,
load, network, nine discovery controller scenarios, equilibrium gate, and
non-gauge eigenvalue definition are used.  No new controller policy, load
shedding, or topology action is introduced.

## Frozen portfolio selection

Selection uses only the frozen static metadata in
`results/PD39_CLASSICAL_BASELINES.csv`; alpha, damping, eigenvectors, and any
other dynamic outcome are not used to select pairs.

The reported candidate scales are:

* candidate converted dispatch: 4620 MW;
* candidate machine rating: 6400 MVA;
* 10% matching limits: 462 MW and 640 MVA;
* remaining-inertia matching scale: 5% of the all-SG total-inertia scale
  reported by the frozen baseline table.

The selected set is:

* V8;
* all eight 7/8 predecessors;
* the following six 6/8 pairs, each selected because the two portfolios have
  equal converted MW and equal IBR MVA, with only the small dispatch-model
  inertia difference shown in the output table:

  1. `30;32;33;34;36;37` versus `30;33;34;35;36;37`;
  2. `30;32;33;34;36;38` versus `30;33;34;35;36;38`;
  3. `30;32;33;34;37;38` versus `30;33;34;35;37;38`;
  4. `30;32;33;36;37;38` versus `30;33;35;36;37;38`;
  5. `30;32;34;36;37;38` versus `30;34;35;36;37;38`;
  6. `32;33;34;36;37;38` versus `33;34;35;36;37;38`.

The four 7/8 matched pairs are disjoint and cover all eight predecessors:

1. `32;33;34;35;36;37;38` versus `30;32;34;35;36;37;38`;
2. `30;32;33;34;36;37;38` versus `30;33;34;35;36;37;38`;
3. `30;32;33;35;36;37;38` versus `30;32;33;34;35;36;37`;
4. `30;32;33;34;35;37;38` versus `30;32;33;34;35;36;38`.

The pair list is a prespecified descriptive comparison set.  The first and
fourth 7/8 pairs have the largest static mismatch within the disjoint cover;
their exact differences are retained rather than silently discarded.

All selected portfolios are evaluated under all nine frozen discovery
scenarios by joining to the existing discovery rows; those rows are retained
as the frozen numerical record and are not replaced.  A blanket rerun is not
part of this short campaign because it would only rebuild the same compiled
model and add no independent scenario.  The genuinely new numerical checks
below are rerun from fresh network constructions.

## Primary composition analysis

For every selected pair and scenario, report the signed and absolute alpha
difference, margin difference, static differences, and whether the two cases
have the same stability/robustness classification.  Summaries use exact
finite-set effect sizes: maxima, medians, and per-scenario paired differences.
No population inference is claimed from the 6/8 or 7/8 finite sets.

The V8 interpretation is explicitly split into:

* penetration/resource explanation: V8 versus the distribution of 6/8 and
  7/8 cases, together with converted MW, IBR MVA, remaining SG MW/MVA, and
  remaining inertia;
* composition explanation: within-cardinality matched pairs above.

## Frozen numerical audit

The numerical audit covers V8 and the two closest 7/8 pairs in the frozen
list, under nominal and the worst high-PLL discovery corner.  It contains:

1. solver tolerance sweep with `tol = abstol = reltol` in
   `{1e-8, 1e-10, 1e-12}` wherever the installed PowerDynamics/SciML API
   accepts the keyword;
2. reduced-matrix eigenspectrum recomputation through both `eigen` and
   `eigvals`;
3. generalized descriptor eigenvalue cross-check on `(A, M)`, retaining the
   finite non-gauge values and recording API failures rather than substituting
   a different model;
4. Float32, Float64, and BigFloat arithmetic recomputation of the already
   formed reduced Jacobian, explicitly labelled as an arithmetic audit rather
   than a full BigFloat equilibrium solve;
5. finite, full-equilibrium perturbations around V8 in the available control,
   operating-point, and line coordinates.  The frozen signed perturbations
   are +/-0.5% for PLL/filter/current-control/load/IBR coordinates and +/-1%
   for two representative line-impedance coordinates (the first and last
   branch in the frozen branch ordering).

Every row records equilibrium residual, fixed-point status, spectrum status,
condition estimates, solver settings, and any exception text.  A disagreement
is investigated and reported; no tolerance or solver result is silently
selected because it gives a preferred conclusion.

## Decision language

The follow-up may support, weaken, or refute a composition-dependent
interpretation.  It must not call V8 a numerical artifact unless the audit
identifies a reproducible pathology.  If the reduced Jacobian is ill
conditioned, that fact is reported separately from the physical interpretation.

The follow-up does not execute weak-node/link radii, co-design optimization,
or the 128-direction radius campaign.  Those remain explicitly out of scope
for this short audit.
