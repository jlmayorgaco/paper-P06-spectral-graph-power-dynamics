# PD39 255+1 — final mechanism and generalization validation preregistration

Status: frozen before new numerics.  
Base commit: `902403cfb0e8dc38323a5623caa51dd663957764`.  
Branch: `research/pd39-255plus1-mechanism-validation`.

This campaign validates the frozen 255+1 discovery mechanism.  It does not
rerun or overwrite the earlier confirmatory campaign, and it excludes weak
nodes, weak links, structured-radius searches, repair optimization, planners,
co-design, IEEE-68, and EMT.

## Frozen definitions

For portfolio `S` and condition `c`, `alpha(S,c)` is the real part of the
rightmost non-gauge eigenvalue after complete equilibrium initialization.  The
true-stability boundary is `alpha >= 0`; the engineering robustness target is
`alpha <= -0.05 s^-1`.  Thus `H0(c)` contains minimal true-stability blockers
and `H0.05(c)` contains minimal blockers of the 0.05 requirement.

The nine frozen discovery scenarios are read from the existing discovery CSV.
The 24 holdout conditions are read verbatim from
`results/PD39_HOLDOUT_CONDITIONS.csv`; they are not regenerated or screened.

## Gate A — numerical truth audit

The audit contains exactly the nine portfolios `V8` plus all eight 7/8
predecessors, crossed with all nine frozen discovery scenarios.  Each case is
evaluated at `tol=abstol=reltol` in `{1e-8, 1e-10, 1e-12}` when accepted by the
installed API.  The converged tight result is the `1e-12` row when available;
the `1e-10` row is the reproducibility reference if the strict initializer
rejects a residual at machine precision.

For every case and tolerance record equilibrium residual, fixed-point status,
rightmost alpha, critical eigenvalue, frequency, damping ratio, reduced
Jacobian condition, and smallest singular value.  Recompute the Jacobian
independently by central finite differences of the full network RHS and apply
the same documented DAE Schur reduction.  Compare `eigen`, `eigvals`, and the
descriptor generalized spectrum `(A,M)` when available.  API limitations are
recorded as unavailable, never replaced by a proxy.

Primary numerical truth pass: tight tolerance results agree within
`1e-6 s^-1`, the independent Jacobian agrees in alpha within `1e-4 s^-1`, and
the independent eigensolver/descriptor checks do not change the sign or the
0.05 classification of any case.  Any failure is retained case-by-case.

## Gate B — high-PLL nonlinear TDS

The four high-PLL corners are exactly the four frozen scenarios with
`pll_scale=1.2`.  The three portfolios are:

* V8;
* strongest 7/8: the 7/8 predecessor with largest frozen `m_9`;
* weakest 7/8: the 7/8 predecessor with smallest frozen `m_9`.

This produces exactly 12 simulations.  Each uses the preregistered common
disturbance: the largest positive active-power load bus in the stock base
case, +1% P and Q at constant power factor from 1.0 to 1.1 s, restored and
simulated to 20 s.  Integration settings and axes are identical.  Record
frequency spread, load-bus voltage, settling, estimated decay/growth rate,
and a comparison with linear alpha.

Pass means the nonlinear estimated sign and fragile/robust ordering are
consistent with the linear prediction for the three portfolios at the four
corners.  Exact equality of rates is not required.  A failed integration is a
reported failure, not a substituted disturbance.

## Gate C — complete holdout census

Evaluate every one of the 256 portfolios at every one of the 24 frozen
conditions: exactly `256 * 24 = 6144` portfolio-condition evaluations.
Every row records equilibrium status, voltage feasibility, alpha, margin,
true stability, and 0.05 robustness.  No case is screened out before the
equilibrium/spectrum attempt.

For each condition `c`, enumerate exact finite-set `H0(c)` and `H0.05(c)` by
testing every portfolio against its strict proper-subset relation.  Report
counts across the 24 conditions, the condition-wise identities, failures, and
coverage.  These are designed-condition census results, not population
probabilities.

## Gate D — penetration, composition, and mechanism

Use all eight 7/8 predecessors and V8, plus the aggregate-matched portfolio
comparisons defined from the frozen static table before outcomes are read.
Match cardinality, converted MW, IBR MVA, remaining SG MW/MVA, and remaining
inertia; report exact differences and paired alpha effects.  Do not attribute
the V8 cliff to one aggregate variable without the matched comparison.

Track the critical mode across the 7/8-to-V8 comparison using common
bus-voltage/network observables and modal assurance.  Classify the mechanism
only as physical electromechanical/control, converter-control oscillatory,
slow control, real aperiodic, DAE/index-related, or numerical/unresolved after
Gate A diagnostics.

A physical homotopy is attempted only after an API/model-semantics audit finds
a documented continuous interpolation between the SG and GFL component
models that preserves the intended load and dispatch.  If the model exposes
only discrete `replace_buses` compilation with no such interpolation, the
homotopy is `BLOCKED`; no artificial blend is constructed.

## Gate E — exact network-closure feasibility

Audit whether the repository and installed model define constructible exact
objects `K`, `D`, and `Q` for the requested closure identity.  The audit must
identify their dimensions, entries, units, and construction path from model
objects.  If all are constructible, test the determinant identity and the
`Q -> -1` boundary with symbolic or high-precision arithmetic.  If any object
or identity is undefined in the installed model, report `BLOCKED` with the
precise missing semantic/API element.  A graph Laplacian or arbitrary proxy is
not allowed to stand in for `K`, `D`, or `Q`.

## Frozen decision outputs

The final report must state pass/fail separately for Gates A–E, retain all
negative outcomes, and end with exactly one of `FINAL CASE: A`,
`FINAL CASE: B`, or `FINAL CASE: C`:

* A: mechanism and generalization validated;
* B: partially validated; claims must be narrowed;
* C: mechanism not validated or the required audit is blocked.

