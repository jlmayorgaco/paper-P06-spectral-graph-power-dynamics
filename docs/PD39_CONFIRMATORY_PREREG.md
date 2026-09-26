# PD39 confirmatory campaign preregistration

**Frozen before confirmatory numerical outcomes**  
**Branch:** `research/pd39-robust-transition-confirmatory`  
**Discovery parent:** `61554336`  
**Discovery audit:** `docs/PD39_CONFIRMATORY_DISCOVERY_AUDIT.md`

This document freezes the confirmatory extension. Discovery data are retained
as discovery data and are not overwritten, relabelled, or used as an
additional holdout. Any departure from this protocol is recorded in
`PD39_CONFIRMATORY_DEVIATIONS.md` before the affected result is interpreted.

## 1. Frozen model and intervention universe

- Network: the maintained PowerDynamics IEEE-39 tutorial network, 39 buses
  and 46 branches.
- Dynamic SG model: the stock Sauer-Pai / AVR Type I / TGOV1 composition.
- GFL model: the documented `PowerDynamics.Library.ComposableInverter.SimpleGFLDC`.
- Candidate SG buses: `{30, 32, 33, 34, 35, 36, 37, 38}`.
- Slack bus: 31; it is not replaceable.
- A replacement preserves the original bus power-flow model and dispatch.
- No artificial load shedding is allowed. Every accepted design must solve
  the requested load, preserve active/reactive balance, remain connected, and
  satisfy the frozen voltage bounds of the stock model audit.
- The eight candidate active dispatches and their sum are read from the
  official machine/bus tables at runtime and copied into every design table.

## 2. Three margins

For a qualified equilibrium, let `alpha(S,c)` be the real part of the
rightmost non-gauge eigenvalue. Eigenvalues with `abs(lambda) <= 1e-8` are
gauge/numerical modes and are excluded.

```text
M1: m_alpha(S,c) = -alpha(S,c)

M2: m_9(S) = min over the nine frozen discovery scenarios c of m_alpha(S,c)

M3: rho_0(S)   = inf ||delta||_infinity such that alpha(S,delta) = 0
    rho_0.05(S)= inf ||delta||_infinity such that alpha(S,delta) = -0.05 s^-1
```

`m_9` is a finite nine-scenario margin, never a structured radius. `rho_0`
is distance to true instability and `rho_0.05` is distance to violation of
the engineering target. The reported radii are the smallest boundaries found
under this frozen search protocol, not globally certified optima.

For nominally unstable designs, a radius search may start on the violating
side; a reported boundary is accepted only if the search also finds a
qualified point on the safe side and a qualified point on the violating side.

## 3. Confirmatory stress box

The primary normalized coordinates are physically grouped as follows.

### 3.1 Control coordinates that actually exist

For each active GFL template, the common-policy coordinates are:

1. PLL bandwidth, changing the implemented `PLL_Kp` and `PLL_Ki` consistently;
2. filter reactance `Xf` (with `Rf` fixed as in the stock model);
3. current-controller bandwidth, changing the implemented current-loop gains
   consistently.

The primary range for each is `[-0.10, +0.10]` relative to its nominal
constructed value. The secondary stress range is `[-0.25, +0.25]`. No
unimplemented outer-loop, limiter, DC-link, AVR, or governor parameter is
invented or tuned.

### 3.2 Operating-point coordinates

- total load P multiplier: `[-0.05,+0.05]`;
- total load Q multiplier: `[-0.05,+0.05]`;
- aggregate replaced-GFL active injection multiplier:
  `[-0.05,+0.05]`.

Every point is re-solved to equilibrium. Active balance is preserved by the
deterministic proportional redispatch rule over surviving synchronous
machines with available positive dispatch; the rule and resulting dispatch
are stored. No load is shed. Reactive balance is solved by the network power
flow and is rejected if the resulting equilibrium is infeasible.

### 3.3 Network coordinates

There is one coordinate per official branch. For branch `e`, the primary
coordinate changes its series impedance as

```text
R_e(delta) = R_e0 * (1 + delta_e)
X_e(delta) = X_e0 * (1 + delta_e),  delta_e in [-0.10,+0.10].
```

Branch shunts and topology are unchanged for the continuous radius box. The
secondary range is `[-0.25,+0.25]`. The reinforcement convention used later
is `gamma_e >= 1`, with `R_e=R_e0/gamma_e` and `X_e=X_e0/gamma_e` and shunts
unchanged.

The complete primary box therefore contains the three control coordinates,
three operating-point coordinates, and 46 branch coordinates. The coordinate
normalization is always actual change divided by the stated allowed scale.

## 4. Fresh holdout conditions

Generate exactly 24 maximin Latin-hypercube points over the complete primary
box, without screening by stability or margin. Use deterministic seed
`39024`, a fixed implementation, and publish the generated coordinates as
`C01` through `C24` in row order. These are designed deterministic test
conditions, not a probability sample and not probabilistic coverage.

The holdout reports, for every proposed design:

- equilibrium success and voltage feasibility;
- 24/24 coverage;
- worst and median `alpha`;
- IQR and extrema of `m_alpha` where qualified;
- worst `rho` estimate where the radius protocol is run;
- paired per-condition differences between planners.

## 5. Modal analysis and numerical audit

For V8 and all eight 7-of-8 predecessors, archive the complete spectrum,
rightmost complex eigenvalue, frequency, damping ratio when oscillatory,
common-coordinate mode shape, participation, equilibrium residual, reduced
Jacobian condition estimates, smallest singular value, critical-eigenvalue
conditioning, duplicate/pinned poles, and active/disabled/integrator state
metadata when exposed by the API.

Mode comparison uses a common bus-voltage/network-observable representation,
not device-state indices whose dimensions change with replacements. A mode is
classified as electromechanical, converter-control oscillatory, slow control,
real aperiodic, DAE/index/singularity-related, or numerical/unresolved only
after these diagnostics pass. An unresolved mode cannot support a headline
mechanistic claim.

The 7/8-to-8/8 table reports, for each missing bus,
`Delta alpha = alpha(V8)-alpha(S_-i)`, damping and frequency changes, and a
common-observable modal assurance coefficient. The interpretation is
explicitly one of: same mode moves right; mode switch; real control pole;
algebraic/index pathology; or initialization/model singularity.

## 6. Structured-radius search protocol

The primary search is deterministic:

1. calculate local eigenvalue sensitivities/automatic derivatives where the
   installed API supports them;
2. construct 128 maximin-LHS directions with seed `39025`, including the
   positive and negative coordinate axes;
3. for each direction, perform nonlinear equilibrium re-solution and exact
   spectrum evaluation while expanding a bracket within the frozen box;
4. retain only brackets with a qualified safe side and a qualified violating
   side;
5. refine the boundary by bisection to normalized infinity-norm tolerance
   `1e-3` and record the final bracket and equilibrium residual;
6. audit selected radii with 2048 deterministic random/maximin directions
   using seed `39026`.

The reported record contains `rho`, exact `delta*`, active coordinates,
boundary eigenvalue, frequency, target, side labels, bracket tolerance, and
equilibrium residual. The L2 norm of the same `delta*` is a secondary metric.
No claim of global optimality is made.

## 7. Frozen selection for hidden-fragility analysis

Selection uses discovery data only and is executed before confirmatory radius
outcomes are read:

- all eight 7-of-8 portfolios;
- V8;
- all-SG reference;
- bus-37 singleton;
- one deterministic representative nearest each 25%, 50%, and 75% target,
  choosing highest discovery `m_9` among ties and then bit-mask order;
- all matched pairs satisfying: both nominally stable, equal cardinality,
  converted-MW difference at most 10% of candidate MW, and
  `abs(alpha_A-alpha_B) <= 0.02 s^-1`.

For every eligible pair, compare `rho_0` and `rho_0.05`. Primary evidence for
hidden fragility is a ratio of at least 3 in either direction; strong evidence
is at least 5. If no pair meets the frozen pairing rule, H3 is not supported.

## 8. Static baselines

For every mathematically meaningful portfolio/site, calculate SCR, generalized
SCR only if the model assumptions and data are actually available, Thevenin
strength, electrical distance, effective resistance, minimum voltage,
voltage sensitivity, remaining SG inertia/MVA/MW, and IBR penetration by MW
and MVA. An unavailable metric is recorded as unavailable; a proxy is never
labelled gSCR.

Mode-based baselines include participation and conventional eigenvalue
sensitivity. Report exact ranks, Spearman, Kendall, top-k overlap, and rank
changes; for the eight-node universe effect sizes and rank tables are primary,
not p-values.

## 9. Corrected transition design problem

Let `P_cand` be the sum of active dispatch of the eight candidate SGs. Targets
are exactly 25%, 50%, 75%, and 100% of `P_cand`. Each target is solved under
the following intervention families:

1. placement only;
2. placement plus controller tuning;
3. placement plus line reinforcement;
4. placement plus controller and line reinforcement;
5. placement plus controller, line, and at most one admissible topology action.

At 100%, placement is fixed to V8. The primary robust design must satisfy,
in every discovery scenario: replaced MW at least target, full load served,
balance, connected network, voltage bounds, and `alpha <= -0.05 s^-1`.
It is then evaluated for `rho_0`, `rho_0.05`, and all 24 holdout conditions.

Controller-only V8 repair uses only the three common GFL coordinates above and
the primary ±10% bounds. Line-only repair tests every branch with
`gamma_e in [1,1.25]`; the secondary test permits `[1,1.50]`. A single line
is solved first; two-line combinations are tested only if no single-line
solution works, using deterministic increasing pair order.

Topology-only repair tests zero-switch and one-switch cases only. A switch is
admissible only when the graph stays connected, PF and dynamic initialization
solve, load remains served, voltages remain inside the frozen bounds, and no
known thermal rating is violated. The discovery outage list has no persisted
admissibility labels; therefore the confirmatory diagnostic creates those
labels before topology optimization.

Joint designs report a Pareto set with controller effort, reinforcement effort,
number of topology actions, worst-case alpha, `rho_0`, and `rho_0.05`. These
efforts are not converted into fake common dollars.

## 10. Planner comparison

For each target, compare the following deterministic planners on the same
24 holdout conditions:

- **P0 penetration only:** smallest bit-mask portfolio meeting the target and
  equilibrium/voltage/connectedness checks.
- **P1 nominal small-signal:** smallest bit-mask portfolio meeting the target
  and nominal `alpha <= -0.05 s^-1`.
- **P2 static strength:** among target-feasible portfolios, choose the one
  minimizing the sum of the discovery-frozen static node weakness proxy;
  ties use bit-mask order. If a stronger defensible static metric is
  implementable without violating the no-fake-gSCR rule, it is documented and
  selected before outcomes are evaluated.
- **P3 robust dynamic:** require all nine discovery scenarios to meet the
  0.05 margin and maximize the frozen structured-radius objective; report
  intervention effort and all failures if no feasible design exists.

Compare feasible coverage, worst alpha, `rho_0`, `rho_0.05`, voltage range,
and remedial effort. If P1 or P2 performs as well as P3, practical advantage
is reported as weak.

## 11. TDS preregistration

The common primary disturbance is determined mechanically from the stock base
case: select the bus with the largest positive active-power load, breaking ties
by the smallest bus number. At `t=1.0 s`, apply a +1% P and Q load pulse at
constant power factor for 100 ms, restore the original load at `t=1.1 s`, and
simulate to `t=20 s`. The secondary disturbance is the same event at +5%.

Use identical solver settings and plot axes for all cases. The primary cases
are:

1. best 7-of-8 by discovery robust margin;
2. worst 7-of-8;
3. V8;
4. V8 plus best controller-only repair, if one exists;
5. V8 plus best line-only repair, if one exists;
6. selected Pareto joint repair;
7. optional all-SG reference.

Observables are bus-frequency spread, the relevant voltage observable, modal
envelope, settling time, and estimated decay/growth rate. The primary TDS
claim is sign and relative-damping consistency with linear alpha, not exact
rate equality.

If the installed PowerDynamics API cannot faithfully implement the specified
temporary load event and restoration, TDS stops. The API limitation and a
replacement physical disturbance must be documented and preregistered before
any TDS run; no silent substitution is allowed.

## 12. Figures, tables, and kill criteria

The campaign produces the requested F1–F13 figures and T1–T7 tables as
vector PDF/SVG plus PNG where the underlying quantity is estimable. Missing or
non-estimable quantities are shown as such, never replaced by arbitrary bad
scores. The final review explicitly audits novelty, correctness, model
adequacy, evidence, usefulness, IAS impact, and TPWRS readiness on a 0–10
scale.

The story is downgraded if V8 is pathology, trivially predicted by ordinary
metrics, radius is nearly deterministic from nominal alpha, classical ranks
match closely, robust co-design has no holdout benefit, repairs fail holdout,
or TDS consistently contradicts the linear ordering. Negative outcomes are
preserved.

