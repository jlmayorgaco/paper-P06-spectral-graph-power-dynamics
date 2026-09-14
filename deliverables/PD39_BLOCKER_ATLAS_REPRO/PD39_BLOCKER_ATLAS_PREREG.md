# PD39 Physical Blocker Validation + First Compatibility Atlas

Preregistration frozen before new numerical results on branch
`research/pd39-physical-blockers-atlas-v1`.

## Decision and scope

The primary decision is whether low-order minimal H0 blockers in the C12
fresh holdout are physically meaningful dynamic instabilities. Only Gates A--C
and, conditional on a GO, the narrow atlas/continuation/theory/remediation
extensions are in scope. Weak-node/link ranking, structured-radius search,
planner comparison, broad optimization, IEEE-68, EMT, and new controller
searches are out of scope.

## Frozen blocker selection

From the complete committed holdout census, select the first two cardinality-2
minimal H0 blockers in `(cardinality, condition, portfolio)` order and the
first cardinality-3 blocker in the same order: `35;36` at C12, `37;38` at C12,
and `30;32;33` at C12. The fixed stable controls are `32;36`, `36;38`, and
`30;33;35`, respectively, selected by same cardinality, then minimum MW
difference, MVA difference, inertia difference, and lexicographic portfolio.
No selection may be changed after results are observed.

## Gate A: numerical/DAE physicality

For each selected blocker, every immediate predecessor, every proper subset if
cheap, and its matched control, solve equilibrium at the frozen condition C12
and at nominal condition. Use the native PowerDynamics reduction as E0 and a
separate dense/generalized eigensolver as E1. Use the descriptor generalized
problem whenever the mass matrix is accessible; otherwise write
`DESCRIPTOR=BLOCKED` with the exact API reason.

Use equilibrium tolerances `1e-8`, `1e-10`, and `1e-12` where supported, with no
silent relaxation. For the independent derivative use central finite
difference steps `1e-5`, `1e-6`, and `1e-7` relative to `1+abs(u[j])`, fixed
before inspection. Report convergence, residual, voltage extrema, P/Q balance,
load service, native/independent alpha, eigenvalue and frequency, Jacobian and
algebraic conditioning, singular values, eigenvector conditioning, and any
zero/gauge modes.

Mechanism rules are frozen: `EM_OSCILLATORY` if frequency >= 0.1 Hz and the
dominant participation is synchronous-machine angle/speed; 
`CONVERTER_OSCILLATORY` if frequency >= 0.1 Hz and converter/PLL/current
states dominate; `MIXED_OSCILLATORY` if both groups exceed 0.2 normalized
participation; `APERIODIC_DYNAMIC` if frequency < 0.1 Hz without an algebraic
singularity; `DAE_SINGULARITY` if algebraic conditioning or a structural/gauge
mode explains the boundary; `NUMERICALLY_UNRESOLVED` if independent methods
disagree in sign or the equilibrium is not reproducible. A large positive
alpha alone is never physicality proof.

Minimality separation is

`delta_H = min(alpha(H), -max(alpha(G) for proper G subset H))`.

The primary blocker-strength threshold is `delta_H >= 0.02 s^-1`; report
sensitivity at 0.01 and 0.05 without changing the threshold.

## Gate B: nonlinear TDS

Use the previous validated common disturbance: deterministic largest active-P
load bus 39, +1% P/Q constant-power-factor pulse from 1.0 to 1.1 s, restore,
integrate to 20 s, identical settings and axes for every case. For each
selected blocker, every immediate predecessor, and its control, run C12 and
nominal when distinct. Estimate growth/decay and dominant frequency using the
same fixed fitting window and record fit quality. Physical confirmation
requires feasible equilibrium, independent unstable sign agreement, no DAE or
unresolved classification, nonlinear predicted growth, decaying predecessors,
and oscillatory frequency error <= 0.10 Hz.

## Gate C: composition

Aggregate distance is normalized using the pre-existing candidate ranges and
the three coordinates converted MW, converted MVA, and remaining inertia. Use
the frozen gate: distance <= 0.05, absolute alpha difference >= 0.05 s^-1,
stable control alpha <= -0.02 s^-1, blocker alpha >= +0.02 s^-1. Exact MW/MVA
match with opposite sign is very strong evidence. If the gate is not met,
`COMPOSITION_STRONG=FALSE`; thresholds may not be loosened.

## Conditional atlas and boundaries

Proceed only on GO: at least two physically confirmed blockers and at least one
strong composition result. Use exactly two preselected coordinates: common PLL
gain scale x and filter scale y, both normalized in `[0.8,1.2]`. Evaluate an
11x11 grid for the selected blockers, their subsets, controls, V8, and 7/8
predecessors; full 256 portfolios are preferred only if the exact model solve
is feasible without changing scope. Refine only cells with a stability-sign
change until boundary uncertainty <=0.01 normalized units or solver limit.
At each cell compute H0, H0.05, kappa0 and kappa0.05, where kappa is the
smallest cardinality among blockers under the corresponding rule. Ground truth
is direct equilibrium plus spectrum, not a surrogate. Trace selected blocker
boundaries with predictor-corrector if practical and record mode MAC and
minimality along the boundary.

## Contextual return and theorem

Audit exact PowerDynamics extraction of device ports, local admittance changes,
network transfer, and K/D/Q. If unavailable, mark PD contextual return
`BLOCKED` and do not construct a proxy. Then use the existing frozen TX4
implementation, without modifying its data, for H4 and one additional TX4
blocker: test the exact Schur-return identity at the existing boundary and
record the nearest return eigenvalue to 1 and `det(I+Q_R)` status. The theorem
writeup must state only the block-partition/Schur-complement result under its
explicit regularity assumptions; no walk, cycle, or holonomy score is primary.

## Conditional remediation and diagnostics

If at least two physical blockers pass, evaluate only local first-order
boundary derivatives for the frozen controller/line/topology actions and solve
one first-order L1 LP. Perform one exact nonlinear validation. Do not iterate
constraint generation. Compare blocker membership with existing participation,
eigenvalue sensitivity, and valid static metrics; acknowledge exact recovery if
it occurs.

## Final decision

- CASE A STRONG: at least two physical blockers, one strong composition pair,
  coherent atlas, contextual-return verification in PD or TX4, and information
  beyond trivial aggregate/classical diagnostics.
- CASE B NARROW: physical evidence exists but atlas/composition/remediation has
  limited value.
- CASE C STOP: artifacts, contradictory TDS, no meaningful composition, or no
  practical information beyond conventional diagnostics.

All primary outputs, deviations, failures, and negative results are retained.
