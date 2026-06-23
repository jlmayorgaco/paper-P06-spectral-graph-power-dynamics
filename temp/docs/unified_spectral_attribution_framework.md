# Unified Spectral-Attribution Framework

Date: 2026-06-08

This note integrates the current theory into one calibrated framework.  It is
intended to prevent the paper from drifting between several attractive but
unequally validated claims.

## One Object

The common object is the retained-network nonlinear eigenvalue problem

```text
T(s) = s^2 M + s Sigma(s) + L_tilde,
Sigma(s) = D0 + C (sI - Acc)^(-1) B.
```

It is best interpreted as a spectral-attribution ledger:

- `L_tilde` is the network-geometry/topology lever.
- `M` is the inertia or virtual-mass lever.
- `Sigma(s)` is the controller self-energy lever.
- The damping margin is read from the oscillatory zeros of `T(s)`.

The framework asks: which lever moved the critical pole, what control comparison
rules out the simpler explanation, and does the reduced prediction agree with
the full NEP or full ANDES eigensolve?

This is not a claim that Schur complements, rational NEPs, contour integration,
Rayleigh sensitivities, or perturbation theory are new.  The possible novelty is
the power-system formulation: using a controller-aware self-energy as the
damping-margin bridge for IBR grids when static nodal damping is not where the
modal damping actually lives.

## Core Thesis

IBR damping-margin analysis should not be organized around a static scalar
`(L,M,D)` model unless the controller condensation is approximately constant at
the critical frequency.  In the IEEE39/ANDES cases tested here, that condition
fails: `GENROU.D = 0` does not contain the measured damping, while controller,
exciter, stabilizer, PLL, inverter, and auxiliary states do.

Therefore the physically correct reduced object is the rational self-energy
bridge above.  The diagonal damping region, second-order correction, inertia
sensitivity, reinforcement sensitivity, and cross-term expansion are lower-order
views or perturbations of this same object.

## Contribution Ledger

| Claim | What is concrete | Status | Originality boundary |
| --- | --- | --- | --- |
| Controller-aware rational bridge | The Schur NEP with `Sigma(s)` reproduces the full ANDES critical poles in no-PSS, base, and Mix60 cases, including the Mix60 PLL-family mode missed by scalar and fixed-point bridges. | Verified full-ANDES bridge. | The math tools are known; the contribution is the IBR damping-margin bridge, the negative scalar-bridge audit, and the family attribution. |
| Network/control family taxonomy | Each NEP zero is classified by condensed-state participation `pi_c`; Mix60 criticality is control-family (`pi_c=0.853`). | Verified in diagnostic ANDES cases. | Likely original as a practical damping-margin taxonomy for controller-condensed IBR grids, but needs final citation audit. |
| Damping-region zero-order screen | Modes are placed in `(nu,delta)` and compared with a damping-margin region. | Theoretical/diagnostic. | Not standalone novelty; proportional modal damping and consensus-region ideas have prior art. |
| Second-order margin correction | Closed-form local correction explains off-diagonal modal damping and improves reduced-model median errors; fallback is needed in tails. | Verified in surrogate QEP and reduced PLL model; full-ANDES estimator validation pending. | Perturbation theory is known; the contribution is its use as an interpretable IBR damping-margin correction and rational graph-filter view. |
| Control-resonant reinforcement reversal | In Mix60/no-PSS, 40 strict records satisfy the resonant criteria: a network-family mode moves toward a PLL-dominated condensed-control pole and the self-energy term dominates the damping decrease. Soft-GFL removes the mechanism. | Confirmed modal-local mechanism over validated bridge. | Not a universal Braess paradox and not dominant global-margin collapse. Novelty is the bridge-enabled network-control mechanism audit. |
| Hidden control-resonance margin | A tracked pole can be stable today but close to a condensed-control pole and highly sensitive through `s Sigma'(s)`. Define distance-to-control-pole, self-energy fraction, and action margin to a damping threshold. | Definition and diagnostic route; supported by the resonant audit logic but not yet validated as a predictor. | Potentially new if it forecasts fragile stable cases that scalar eigenvalue or static `(L,M,D)` screens miss. Needs a registered prediction test. |
| Modal substitutability of virtual inertia | The action space decomposes into existing lines, new lines, and shunt/virtual-inertia components. The shunt-defect projection on the critical mode decides whether topology can substitute for local inertia. | Verified in the six-bus poster gate; first-order formulas match finite differences at about `1e-4`. Multi-graph/full-ANDES validation pending. | The edge/diagonal decomposition is classical; the candidate novelty is the robustness-lever taxonomy and the modal classifier for inertia-vs-topology choice. |
| Weak links / damping-aware reinforcement | Frequency-only line strengthening can disagree with damping-margin effects in the surrogate. | Surrogate planning diagnostic; full-ANDES topology planning pending. | QEP sensitivities exist; this is an application unless full-ANDES regret beats registered baselines. |
| Inertia sign reversal | With fixed `D_i/M_i`, 5447/16229 reduced-model perturbations reduce both the diagonal screen and full-QEP damping margin. | Reduced-model sign diagnostic. | Potentially useful, but not yet a full-ANDES IBR claim. Must be distinguished from existing "more inertia can hurt" results in other metrics. |
| Joint inertia-control cross term | The local second-order expansion contains an irreducible `E_D x E_L` term. | Mathematically local; Phase-E0 ANDES found 0/13 measurable singleton GENROU-PLL cross cases. | Do not headline. It becomes empirical only if a physically justified virtual-inertia/control case makes the cross term measurable and improves over additive prediction. |
| Conversion and planning rankings | The framework can define post-action damping-margin ranking and regret. | Pending. | Similar placement and damping-index work exists; this becomes original only if it beats registered controls on full-ANDES/NEP validation. |

## What Is Actually New Enough To Build The Paper Around

The strongest currently supported paper is:

```text
Controller-aware rational self-energy as the damping-margin bridge for IBR grids,
with network/control family attribution and a bridge-enabled control-resonant
reinforcement mechanism.
```

That claim has full-ANDES support at the bridge level and a mechanism audit that
uses the bridge in a way a scalar model cannot.

The second-order estimator, inertia reversal, and planning maps are important
supporting theory, but they are not yet the empirical headline because their
final comparisons over the validated bridge are still pending.

The joint inertia-control cross term is not the current headline.  The local
formula is real, but the registered ANDES E0 gate did not find a measurable
singleton GENROU-PLL effect.  This negative result should stay in the manuscript
because it prevents overclaiming the no-separability narrative.

## How The Pieces Fit

1. Start with full ANDES eigensolve as ground truth.
2. Condense non-retained controller states by Schur complement to obtain
   `Sigma(s)`.
3. Solve `T(s)=0` by contour integration.  This validates the measurement
   instrument and separates network-family from control-family modes.
4. For network-family modes, use the diagonal damping-region screen as the
   zero-order graph picture.
5. Add the second-order correction to explain how non-proportional damping moves
   the critical pole away from the diagonal prediction.
6. Use exact NEP sensitivity to decompose line actions into retained-network
   and self-energy effects.
7. Treat inertia sensitivity and cross terms as local perturbation diagnostics,
   not finite-conversion estimators, unless full-ANDES gates pass.
8. Use planning rankings only after comparing regret against registered
   baselines.

## Matrix, Geometric, And Topological Push

The next layer of originality should not be "another centrality index."  It
should use the NEP matrix geometry to define vulnerability and graph
reconfiguration in pole space.

### 1. Hidden Control-Resonance Margin

For a NEP zero `s_k` and a condensed-control pole `mu_j in spec(Acc)`, define

```text
d_ctrl(k,j) = |s_k - mu_j| / (|s_k| + |mu_j|).
```

Define the self-energy fraction

```text
phi_Sigma(k) =
  | y_k^* s_k Sigma'(s_k) x_k |
  / ( | y_k^* T_s(s_k) x_k | + eps ).
```

A pole has hidden control-resonance exposure when `d_ctrl` is small and
`phi_Sigma` is large.  This is not the same as being unstable.  It means the
system is stable at the present operating point, but the damping-margin
sensitivity is already controlled by nearby controller dynamics.

For an admissible action `p` and a target damping ratio `zeta_star`, define the
local action margin

```text
m_hidden(k,p) =
  (zeta_k - zeta_star) / [ - d zeta_k / dp ]_+.
```

The system hidden margin is the minimum over tracked modes and admissible
actions.  This gives a concrete way to say "the system looks stable, but it is
only a small line/control/inertia perturbation away from crossing the damping
threshold."

Claim status: candidate diagnostic.  It becomes a real contribution only if a
registered test shows that small `m_hidden` predicts future damping-margin loss
better than scalar eigenvalue distance, damping ratio alone, or `(L,M,D)` scores.

### 2. Resolvent-Weighted Weak Nodes And Weak Links

Weak elements should be defined by damping-margin sensitivity, not by graph
centrality:

```text
weak node for damping:   W_i^D     = - d zeta_c / d D_i
weak node for inertia:   W_i^M     = - d zeta_c / d M_i
weak node for control:   W_i^kappa = - d zeta_c / d kappa_i
weak link:               W_e       = - d zeta_c / d w_e.
```

A resonant weak link is stronger:

```text
W_e > 0,
d d_ctrl / d w_e < 0,
phi_Sigma is large.
```

That means strengthening the line reduces damping, pulls the mode toward a
control pole, and the self-energy derivative dominates the pole motion.  This is
topological, but not in the classical centrality sense: the topology is weak
because its perturbation vector enters a dangerous direction in the rational
NEP.

Candidate novelty: "resolvent-weighted topological vulnerability."  The exact
phrase should wait for literature review, but the idea is distinct from Fiedler
participation, impedance distance, and static line centrality because it depends
on `spec(Acc)` and the frequency-dependent self-energy.

### 2b. Line-Shunt Lever Taxonomy And Modal Substitutability

A separate structural result lives in the conservative network-inertia layer.
For the complete graph edge basis `b_ij=e_i-e_j`,

```text
L_all = span{ b_ij b_ij^T : i < j },
S_shunt = { diag(g) : g in R^n }.
```

As vector spaces,

```text
Sym(n) = L_all direct-sum S_shunt.
```

If existing and candidate-new lines are separated,

```text
L_all = L_existing direct-sum L_new.
```

An infinitesimal local inertia action changes the inertial Laplacian as

```text
d L_tilde / d M_i = -(E_i L_tilde + L_tilde E_i) / (2 M_i).
```

Its off-diagonal support is inherited from existing lines, while its row-sum
defect is diagonal/shunt-like.  Structurally,

```text
d L_tilde / d M_i in L_existing direct-sum S_shunt,
```

not in `L_new`.  This does not make inertia equivalent to line reinforcement:
inertia and line reinforcement can move the critical modal stiffness in opposite
directions.  It says where the action lives in the lever space.

For a critical mode `q_c`, define the shunt-defect vector

```text
r_i = 0.5 L e_i,
```

and the modal substitutability score

```text
s_i = 1 - | r_i^T q_c | / || r_i ||.
```

High `s_i` means the inertia action's irreplaceable shunt component is nearly
orthogonal to the critical mode, so topology can be tried first.  Low `s_i`
means the shunt component aligns with the critical mode, so virtual inertia or a
synchronous condenser is the honest lever.

Current evidence: the poster package verifies the sign-correct first-order
formulas on a six-bus system and obtains `s=[0.54,0.14,0.58,0.36,0.90,0.50]`
with an `8.9x` discrimination ratio.  This is promising but not yet a
full-paper claim.  It needs robustness over random graphs and blind validation
against full ANDES damping outcomes.

### 3. Pole-Motion Cone

For action variables

```text
u = [ Delta w, Delta M, Delta D, Delta kappa ],
```

NEP sensitivity gives

```text
Delta s_k ~= J_k u,
J_k = [ ds_k/dp_1, ..., ds_k/dp_r ].
```

In the complex plane, feasible actions form a cone of possible pole motions. The
damping-ratio gradient defines a half-space:

```text
Delta zeta_k ~= g_k^T u.
```

This gives a geometric interpretation:

- safe actions move the pole into the damping-improving half-space;
- dangerous actions move it into the damping-decreasing half-space;
- hidden-resonant actions additionally move it toward a control pole.

This is probably the cleanest "new geometry" in the project.  It turns graph
configuration into a pole-space cone problem instead of a scalar ranking problem.

### 4. Minimum-Cost Graph Reconfiguration

The local planning problem becomes

```text
maximize     t
subject to   zeta_k + g_k^T u >= t,        for protected modes k
             d_ctrl(k,j) + h_kj^T u >= d_min,
             c^T |u| <= B,
             u in engineering limits.
```

The first constraint protects damping.  The second protects hidden
control-resonance distance.  The third is the budget.

This is the route from diagnosis to optimization:

1. compute NEP poles and families;
2. compute weak node/link sensitivities;
3. build the pole-motion cone;
4. solve the local reconfiguration problem;
5. recompute the full NEP or ANDES eigensolve after each selected action.

Claim status: not validated yet.  This is the natural Transactions-level
extension if Experiment A shows that the NEP-based sensitivities predict full
ANDES pole motion reliably.

## Best Candidate For "New New"

The strongest not-yet-fully-exploited idea is:

```text
Hidden control-resonance margin as a topology/control vulnerability certificate:
a stable grid can be close to a damping-margin loss not because the present
critical pole is poorly damped, but because a network-family pole lies in the
resolvent shadow of a controller pole and small admissible graph actions move it
through the self-energy geometry.
```

This is more original than another weak-node score because it is a different
question:

- classical stability asks whether the present poles are stable;
- damping-margin analysis asks whether the present damping ratio is acceptable;
- hidden control-resonance margin asks whether admissible topology/control
  changes can quickly turn a stable network-family pole into a controller-driven
  damping problem.

Required gate:

1. Select cases with a range of hidden margins.
2. Apply registered small line/control/inertia perturbations.
3. Check whether low hidden margin predicts observed damping loss.
4. Compare against baselines: present `zeta`, distance to imaginary axis,
   scalar `(L,M,D)` damping screen, Fiedler participation, and frequency-only
   line sensitivity.
5. Report negatives.  If hidden margin does not predict better than these
   baselines, keep it as interpretation rather than contribution.

## Venue-Dependent Claim Strength

For a master thesis chapter or methodological paper, the current package is
coherent: it contains a validated bridge, derivations, reduced-model estimator
evidence, and one bridge-enabled mechanism audit.

For a conference paper, the cleanest version should focus on the bridge and the
control-resonant reinforcement mechanism, while reporting estimator/planning
gates as future work.

For an IEEE Transactions submission, the next required gate is Experiment A:
diagonal, second-order, reduced-QEP, and adaptive estimator errors must be
compared against the validated NEP/full-ANDES poles over calibrated IBR
operating ranges.  Planning rankings need regret against frequency-only,
Fiedler, damping-index, and full-eigensolve baselines.

## One-Sentence Contribution Statement

This paper introduces a controller-aware spectral-attribution framework for IBR
damping margins, in which the rational self-energy `Sigma(s)` replaces static
nodal damping as the bridge from full controller dynamics to graph-modal
analysis; the framework is validated against full ANDES poles and used to
separate network-family and control-family mechanisms, including a
control-resonant reinforcement reversal that a scalar `(L,M,D)` bridge cannot
see.
