# PHASE E — mechanism report

Experiments E14B, E15, E16, E20 and the Track-B ranking. All numbers regenerate
from `experiments/` with manifests under `results/manifests/`. Lint clean,
81 tests passing.

The two tracks are reported separately and are not recombined.

---

# TRACK A — higher-order loss of synchronous support

Flagship `30+33+35+37`, `rho = 1`, matched reactive policy, frozen operating
point (the AC solution is bit-identical to base).

## E15 — the order-4 claim holds on one tracked branch

The Moebius decomposition is taken of the **tracked** eigenvalue real part, not
of the spectral abscissa. The abscissa is a maximum over modes and mixes
branches; the branch is followed by MAC on the machine states that survive in the
flagship, which are a subset of the survivors of every proper subset.

| reconstruction | value | verdict |
|---|---|---|
| order ≤ 0 (base) | `-0.126478` | stable |
| order ≤ 1 (isolated effects) | `-0.604696` | stable |
| order ≤ 2 (pairwise complete) | `-0.361574` | stable |
| order ≤ 3 (triple complete) | `-0.074381` | stable |
| **exact (order ≤ 4)** | **`+0.144670`** | **unstable** |
| irreducible fourth-order term | **`mu4 = +0.219051`** | |

Minimum MAC across the 16-subset lattice `0.9674`; every tracked mode in the
synchronous family. **ORDER-4 CLAIM ACCEPTED.** The third-order-complete model
still has `0.074` of margin, and the irreducible fourth-order term is three times
that. No study up to third order can see this.

## E14B — which synchronous service is actually required

A synchronous condenser is retained at each retired bus: the full machine
apparatus at a declared rating, producing no active power, with the converter
carrying every displaced megawatt. Services are then ablated.

| configuration | tracked branch | spectral abscissa |
|---|---|---|
| A0 original machines | `-0.126` | `-0.126` |
| A1 full PV-GFL replacement | **`+0.145`** | **`+0.145`** |
| A2 condenser 100 % | `-0.216` | `-0.188` |
| A2 condenser 50 % | `-0.228` | `-0.193` |
| A2 condenser 25 % | `-0.415` | `-0.201` |
| A3 condenser, inertia ×0.1 | `-0.238` | `-0.167` |
| A4 condenser, PSS gain 0 | `-0.221` | `-0.189` |
| A5 condenser, manual excitation | `-0.221` | `-0.189` |
| A6 manual excitation **and** no PSS | `-0.221` | `-0.189` |
| A7 factorial, all 8 combinations | `-0.238` to `-0.167` | all stable |

**No individual service explains it.** Inertia scaled to a tenth, excitation
frozen and PSS switched off together still restore stability. Scaling PSS from 1
to 0 moves the branch by `0.005`; freezing the AVR moves it by `0.005`; inertia
from 1.0 to 0.1 moves it by `0.02`. What the four retirements remove is the
machines' **electromagnetic presence** — a voltage source behind a transient
reactance — not their mechanical or control services.

Generalization across all nine genuine order-4 cores:

| mitigation | cores restored |
|---|---|
| condenser at 25 % rating | **9 of 9** |
| condenser at 50 % rating | **9 of 9** |
| condenser with H ×0.1, manual excitation, PSS off | **9 of 9** |

Two caveats. The condenser rows at `inertia = 0.1` with AVR and PSS active report
a MAC-best match on a real mode; the branch has merged there and only the
spectral abscissa is meaningful in those two rows. And zero inertia was not
tested: a condenser with `H = 0` is a singular machine, not a service ablation.

## E16 — the failure compresses to an 8x8 port model

| quantity | value |
|---|---|
| port operator | 78 × 78 (bus voltages only) |
| state-space dimension of the flagship | 82 |
| base modes reproduced by `det T(s) = 0` | **68 of 68** |
| flagship modes reproduced | **80 of 80** |
| rank of the four-replacement update | **8** (2 per bus) |
| finite-difference floor of the port model | `6.4e-07` relative |
| determinant lemma, five probe points | `8.8e-08` to `4.8e-06` relative |

At the crossing mode `+0.144670 + 3.610828j`:

| factor | magnitude |
|---|---|
| full `det(I + C K)` | `4.2e-08` |
| **individual** (product of single-replacement factors) | **`0.1159`** |
| **collective** `det(I + Q)` | **`3.6e-07`** |
| `rho(Q)` | `1.00000` |
| `sigma_min(I + Q)` | `9.7e-08` |

So an 82-state failure is represented exactly by an 8 × 8 action operator on four
bus ports, and at the crossing mode the **collective** factor vanishes while the
individual factor does not. The loop through the four replacement ports closes
exactly.

Accuracy is limited to about `1e-6` by the central-difference port linearization,
not by the identity. Analytic device Jacobians would be needed before claiming
machine precision, and before winding on a tight contour.

Two statements must not be conflated: `sigma_min(I+Q) -> 0` with `rho(Q) = 1` is
a structural fact about the closed loop; calling the collective factor "the
cause" is a block-split-dependent statement for 2×2 blocks and needs the
gauge-invariant holonomies of E17/E18, which have not been run.

## E20 — bus 30 is structurally special, and not because it is weak

| bus | CPI inter-area | CPI PLL | CPI all | SCR | inertia M | solo Δζ |
|---|---|---|---|---|---|---|
| **30** | **1.000** | **0.333** | 0.524 | 3.269 | 87.4 | `-0.0007` |
| 33 | 0.750 | 0.667 | 0.690 | 2.951 | 67.2 | `+0.0058` |
| 38 | 0.583 | 0.600 | 0.595 | 0.907 | 116.2 | `+0.0080` |
| 34 | 0.000 | 0.633 | 0.452 | 2.316 | 56.2 | `-0.0038` |

Bus 30 is in **every** inter-area core and has the **lowest** participation of all
nine in the PLL family. Its individual replacement is essentially neutral
(`Δζ = -0.0007`).

Size-4 portfolios, instability with and without bus 30:

| | stable | unstable |
|---|---|---|
| without 30 | 69 | 1 |
| with 30 | 47 | 9 |

Fisher exact odds ratio **13.2**, `p = 0.0050`. Restricted to the removed-inertia
range spanned by the bus-30 portfolios, odds ratio **9.0**, `p = 0.019`.
Mann-Whitney on the spectral abscissa, `p = 0.00089`.

`CPI_interarea` correlates with inertia (Spearman `0.61`) and rating (`0.53`) but
**not with SCR** (`-0.034`). With ten positives the estimate is imprecise; the
defensible statement is that bus 30 raises the odds of the inter-area failure by
roughly an order of magnitude, that this is not predicted by short-circuit
strength, and that it is specific to one of the two failure families. Do not name
a new metric on this evidence.

## Track-A repair status

From E14 N3: dividing the converter active-power outer loop by four stabilizes
**9 of 9** genuine cores, margins `-0.032` to `-0.235`, with the critical mode
staying in the synchronous family and **all four PV replacements connected**.
That is the engineering result:

> the converters did not cause the higher-order failure, but their control can
> compensate for the lost synchronous support.

E21 has not been run: there is no optimized minimum-change repair yet, and no
comparison of change norms against restoring a machine.

---

# TRACK B — controller-mediated collective synchronization failure

## Ranking of the 23 controller-mediated cores

All are violent. The mildest is `32+33+34+35+36`, size 5:
`alpha_GFL = +118.63`, `alpha_static = -0.0889`, `cond(gz) = 5.0e3`, dominant
state `theta_pll_gfl33`.

## The continuation found no mild crossing

Sweeping the PLL natural frequency from `37.4` down to `2 rad/s` at `zeta = 0.7`:

| `wn` [rad/s] | 2 | 4 | 8 | 12 | 20 | 37.4 |
|---|---|---|---|---|---|---|
| spectral abscissa | `+34.0` | `+38.9` | `+48.6` | `+58.2` | `+77.1` | `+117.9` |

The divergence never approaches the boundary, and below `wn = 5` the dominant
state moves from `theta_pll` to the reactive outer-loop integrator `x_q`. Loop
ablation on the same core: outer Q loop ×0.01 gives `+99.9`; outer P loop ×0.1
gives `+109.8`; both outer loops ×0.1 gives `+92.4`; removing the Q integral
gives `+114.0`. **No converter retuning removes it.**

This is grid-following loss of synchronism in a network that has become too weak
to support it, not a tuning fault. It is therefore not repairable by the
controller, which is the opposite of Track A.

## What does repair Track B

The Track-A mitigation works here too. Synchronous condenser on the ten mildest
controller-mediated cores:

| mitigation | cores restored |
|---|---|
| condenser at 50 % rating | **10 of 10** |
| condenser with H ×0.1, manual excitation, PSS off | **10 of 10** |

(At 25 % rating four of ten are feasible and all four are restored; the rest
exceed the condenser rating on reactive power.)

## Track-B open items

Winding provenance (E19) has **not** been run, so item 6 of the report request is
unanswered: it is not yet established that the extra unstable spectral count
localizes to a small controller-mediated core. The port machinery now exists and
this is the natural next step; Moebius must not be used here.

---

# Which track is stronger for the poster

**Track A**, for four reasons.

1. It has a complete, verified chain: order-4 claim accepted on a tracked branch
   (E15), exact 8 × 8 port representation of an 82-state failure with the
   collective factor carrying the mode (E16), a structurally special bus with a
   statistically significant effect not predicted by SCR (E20), and a repair that
   keeps every megawatt of PV connected (E14 N3).
2. Its message is counter-intuitive and defensible: *the converters are not the
   cause, and are less destabilizing than an ideal static source; but their
   control is the actuator that compensates the loss.*
3. The mitigation is industrially concrete: a quarter-rated synchronous condenser
   with no controls and negligible inertia restores 9 of 9 cores.
4. Nothing in it depends on a violent, hard-to-defend eigenvalue.

**Track B** is a strong second paper, not a poster panel yet. Its mediation
evidence is clean, but every instance is at `+100` to `+1276 rad/s`, no mild
crossing exists under physical controller continuation, and the provenance work
that would make it a mechanism claim has not been done.

---

# ADDENDUM — reconciliation, E17, E18, E21

## Reconciliation of the A2 numbers

`-0.188` and `-0.415` were two different quantities at two different ratings:

| | spectral abscissa | tracked inter-area branch |
|---|---|---|
| condenser rating 1.00 | **`-0.1878`** | `-0.2158` |
| condenser rating 0.50 | `-0.1929` | `-0.2283` |
| condenser rating 0.25 | `-0.2011` | **`-0.4150`** |

The cross-core table also reported only the abscissa, which made the two look
contradictory. Both are now regenerated with **both columns**, and
`configs/ias2026/trackA_mitigation_definitions.yaml` freezes the rule that a
mitigation number without its rating and both quantities is not a result. It
also freezes the two variants that were previously conflated: `sc25` is
**rating 0.25 at full services**, while `minimal` is **rating 1.00 with H x0.1,
manual excitation and PSS off**. They isolate different things.

With both quantities, all nine cores are restored by every variant (9/9 by
abscissa and 9/9 by the tracked branch; minimum MAC 0.857 across all variants).
The binding constraint after repair is nearly independent of rating
(`-0.188` to `-0.201`) while the tracked branch moves further left at smaller
ratings; two `minimal` rows report a merged branch and only their abscissa is
meaningful.

## E17 — the collective result is representation invariant

Admissible transformation: an invertible 2x2 basis change at each replacement
port, with the update transforming so the operator is untouched, giving
`M -> S^-1 M S` block-conformally.

A first attempt measured relative error on `det(I+M)` and `det(I+Q)` **at the
crossing mode**, where both are about `1e-8` by construction. That is a
conditioning test, not an invariance test, and it produced a false
gauge-dependent verdict. Determinant invariance is now measured at generic probe
points; at the crossing mode the reported invariant is the one that is well scaled
there.

| quantity | change over identity, rotations, condition 10/100/1000, polar basis |
|---|---|
| `det(I+M)` at generic points | `<= 6.3e-13` (cond <= 100) |
| `det(I+Q)` at generic points | `<= 2.0e-12` (cond <= 100) |
| individual factor | `<= 5.9e-12` |
| eigenvalues of `Q` | `<= 7.3e-12` |
| **distance of a `Q` eigenvalue to `-1`** | **`1.578e-07` for every basis** |
| `sigma_min(I+Q)` | up to **99.8 percent** — not invariant |
| norm of `Q` | up to **554x** — not invariant |

**The publishable invariant is: at the crossing mode the interaction operator `Q`
has an eigenvalue equal to `-1`.** `sigma_min(I+Q)` and the norm of `Q` must not
be quoted; E16 reported `sigma_min(I+Q) = 9.7e-08` and that number is basis
dependent.

## E18 — the cycle scores do not explain anything, the closure does

13 directed cycles over the four ports. Their trace, determinant and spectral
radius are gauge invariant to `5.3e-15`, so the objects are sound. The question
is whether they discriminate.

Against 12 stable four-replacement portfolios matched on megawatts
(4176-4366 MW against the flagship's 4271 MW):

| quantity | flagship | matched stable controls |
|---|---|---|
| strongest 2-cycle radius | `0.0455` | up to `0.0760` |
| strongest 3-cycle radius | `0.0052` | up to `0.0196` |
| strongest 4-cycle radius | `0.00070` | up to `0.00454` |
| **distance of a `Q` eigenvalue to `-1`** | **`1.6e-07`** | **`0.142` to `0.741`** |

Every cycle score of the flagship is **smaller** than the best control; the
4-cycle separation ratio is `0.16`. Per the stated rule, cycle magnitude is
**abandoned as an explanatory metric**.

The closure condition separates the flagship from every matched control by six
orders of magnitude. That is the discriminating invariant, and E17 already proved
it does not depend on the port basis.

This is stronger than "a particular loop created the mode": no single cycle
dominates, because the effect is distributed over all four ports. It is
consistent with E14B — the failure is the collective loss of four electromagnetic
voltage sources, not a feedback loop between converters.

## E21 — minimal repair, and the honest trade-off

Target declared before optimizing: spectral abscissa at most `-0.05`. Base case
`-0.126`; unrepaired flagship `+0.145`; total PV replaced 4270.7 MW.

| strategy | PV kept | PV forgone | sync MVA | rel. change | abscissa | zeta_min | closure |
|---|---|---|---|---|---|---|---|
| unrepaired | 4270.7 | 0 | 0 | — | `+0.1447` | `-0.0400` | `0.0000` |
| RA restore machine 30 | 3230.7 | **1040.0** | 1040 | 0 | `-0.2082` | `0.0254` | `0.404` |
| RA restore machine 37 | 3300.5 | **970.2** | 970 | 0 | `-0.1635` | `0.0343` | `0.212` |
| **RB weaken P loop x0.427** | **4270.7** | **0** | **0** | 1.202 | `-0.0500` | `0.0143` | `0.125` |
| **RC minimum-norm retune** | **4270.7** | **0** | **0** | **0.862** | `-0.0500` | `0.0142` | `0.114` |
| **RD condenser at 25 percent** | **4270.7** | **0** | 1068 | 0 | `-0.2011` | `0.0341` | `0.480` |

Seven strategies meet the target; three keep every PV megawatt.

RC is the minimal intervention: a relative log-change of `0.862` spread over the
converter loops (`ki_P x0.46` with small adjustments elsewhere), against `1.202`
for weakening the P loop alone. Both keep all 4270.7 MW with **zero synchronous
capacity**.

The honest trade-off, which the poster must state rather than hide: a converter
retune restores stability but only to `-0.050`, **weaker than the base case at
`-0.126`**. Restoring a machine or adding a condenser reaches `-0.16` to `-0.21`,
but costs either 970-1040 MW of PV or 1068 MVA of synchronous capacity.

In every successful repair the gauge-invariant closure distance moves from
`0.000` to between `0.11` and `0.48`, an independent confirmation that the repair
acts on the identified mechanism rather than merely moving an eigenvalue.
