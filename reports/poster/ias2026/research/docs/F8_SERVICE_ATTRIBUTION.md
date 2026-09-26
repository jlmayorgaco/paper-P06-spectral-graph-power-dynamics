# F8 — Which synchronous services reshape the incompatibility hypergraph

Post-freeze work (after tag `IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE`); F7
results were read, never modified. Code: `experiments/F8_service_attribution.py`,
`F8_report.py`, `F8_tables.py`, `F8B_mean_vs_heterogeneity.py`,
`F8C_service_closure.py`, `F8D_avr_blend.py`; tests
`tests/test_service_interventions.py`. Data: `results/F8/`, `results/F8B/`,
`results/F8C/`. Same model scope as F7 (IEEE-39, harmonized first-order AVR,
matched dispatch, core 30/33/35/37, `Gamma` = RHP x 0.3–1.5 Hz).

**Verdict.**

1. The service whose restoration returns the fleet to full composability is
   **electromagnetic presence at the retired buses** — a voltage behind
   transient reactance. A condenser with 1 % of the retired machine's kinetic
   energy, frozen EMFs, field held, no PSS, no damping and no reactive output
   empties `H` at four of five unsafe points at 25 % rating and at all five at
   full rating. The effect is graded: growing the condenser rating removes the
   hyperedges one at a time, reversing the F7 contraction.
2. Inertia is **not** the restoring service. Synthetic inertia without
   electromagnetic presence leaves `H` non-empty at 4 of 5 unsafe points, and
   **surviving-fleet inertia is destabilizing** in both directions (x0.5 empties
   `H` at 4 of 5 points; x2 makes the base unstable or `kappa` smaller).
3. The surviving fleet's **excitation dynamics are necessary for
   non-composability but cannot be switched off**: manual excitation makes the
   base grid aperiodically unstable. Slowed continuously, the AVR moves the same
   fleet through `kappa = 1 -> 2 -> 3 -> 4 -> inf` and reaches composability
   before the base fails.
4. PSS, damping and reactive share mostly move margins; they change `H` only at
   specific points, and several change the witness coalition itself.
5. Port closure marks every one of 1 511 service-induced imaginary-axis
   boundaries (closure distance at most `9.7e-6`); the underlying physical cause
   does not matter to the reduced descriptor.
6. Excitation **heterogeneity is causal** at matched fleet mean under every
   definition of the mean, but the **mean dominates** the margin, and the sign
   of the heterogeneity effect depends on which parameter is spread and on the
   policy.

## 1. Design (declared before running)

Points, from the frozen F7A map (interior of the largest `H` region of each
`kappa`, distance-transform maximum), plus a point inside the instability tongue
0.01 above its tip:

| point | `g` | `k` | frozen `H` | `kappa` |
|---|---|---|---|---|
| `P_inf` | 1.000 | 0.500 | `EMPTY` | inf |
| `P4` | 0.036 | 1.425 | `30+33+35+37` | 4 |
| `P3` | 0.000 | 1.631 | `30+33+35 | 30+33+37` | 3 |
| `P2` | 0.110 | 1.850 | `30+33 | 33+35 | 33+37 | 30+35+37` | 2 |
| `P1` | 0.133 | 2.178 | `30 | 33 | 35 | 37` | 1 |
| `P_fold` | 0.103 | 1.260 | `30+33+35+37` | 4 |

Interventions are physical model configurations, applied consistently to all 16
subsets, never Jacobian scalings:

- **Class R — restore services of the retired machines.** A synchronous condenser
  (no active power) at each retired bus. Full `2^7` factorial: electromagnetic
  presence (rating 0.25 / 1.0 of the retired machine), S1 inertia (1 % / 100 %
  of its kinetic energy), S2 damping (0 / 2 pu), S3 electromagnetic transients
  (EMFs frozen = classical / dynamic), S4 AVR (field held / native), S5 PSS
  (off / native), S6 reactive support (condenser carries 0 / all bus Q). Plus the
  plain replacement and converters with synthetic inertia equal to the retired
  machine's `H` (inertia **without** electromagnetic presence).
- **Class A — the surviving fleet.** Inertia x0.5/1/2, D 0/2, EMF transients,
  AVR, PSS; `3 x 2^4` factorial.
- **S7 governor**: not represented in the harmonized model (TGOV1N dropped,
  reheat constant outside the band, `Dt = 0`); excluded, not ablated.

Every switch keeps the operating point (tests), and reduces to a documented
configuration at its ends (AVR blend 0 = manual excitation; flux blend 0 =
classical machine; synthetic inertia 0 = original spectrum plus one decoupled
pole).

1 068 (point, intervention) lattices were solved on the direct path. The full
table (`intervention, old H, new H, old kappa, new kappa, d alpha_IA, d m_cl`) is
`results/F8/F8_intervention_table.csv`; the key rows are
`results/F8/F8_key_table.md`.

## 2. The eight questions

**1. Which individual service can change kappa?** Counted over all factorial
pairs that differ in one service only (`results/F8/F8_flip_summary_*.csv`):

| service | where it changes `kappa` |
|---|---|
| electromagnetic presence (condenser vs none) | every unsafe point (`kappa` -> inf at P2–P4, P_fold; at P1 only at full rating) |
| condenser rating 0.25 -> 1.0 | P1 (64 of 64 pairs), P2 (20 of 64) |
| condenser inertia | P2 (20 of 64); never at P3, P4, P_fold |
| condenser EMF transients | P2 (20 of 64) |
| condenser damping, reactive share | P2 (4 and 8 of 64) |
| condenser AVR | P2 (4 of 64) |
| condenser PSS | never |
| surviving-fleet inertia | P1–P4, P_fold |
| surviving-fleet EMF transients, damping | P1–P4, P_fold (damping: P2, P3, P_fold) |
| surviving-fleet PSS | P3, P_fold |

**2. Which services only move margins?** Condenser PSS (0 `H` changes in 384
pairs), and at P3, P4, P_fold and P_inf **every** condenser service other than
presence itself: all 128 condenser configurations give `H = EMPTY` there, and the
services only move `alpha_IA` and `m_cl`.

**3. Genuine synergy.** Pairs whose joint change moves `H` while neither single
change does (`results/F8/F8_synergy_*.csv`): condenser EMF transients x AVR,
AVR x reactive share, EMF transients x reactive share, damping x EMF transients
(P1, P2 only; 22 instances), and surviving-fleet inertia x EMF transients
(9 instances). The first is structural — an AVR acts only through the EMF
dynamics — and is reported as such, not as a discovery.

**4. Is AVR necessary, or one coordinate?** Both, in a precise sense.

- Manual excitation of the surviving fleet (field held) makes the **base**
  aperiodically unstable (a real mode at `+0.23 s^-1`, loss of synchronizing
  torque) with or without damping. So "remove the AVR" is not a viable
  configuration of this network.
- A surviving fleet without EMF transients (classical machines) with damping
  `D = 2` is composable (`H = EMPTY`) at all six points; with `D = 0` its base is
  unstable.
- Slowing the surviving AVR continuously (`efd' -> beta efd'`, `F8D_avr_blend.py`)
  moves the same fleet through the regions before the base fails:

  | point | sequence as `beta` falls from 1 | `H` empty at | base fails at |
  |---|---|---|---|
  | P1 | `kappa 1 -> 2 -> 3 -> 4 -> inf` | 0.53 | 0.20 |
  | P2 | `2 -> 3 -> 4 -> inf` | 0.62 | 0.20 |
  | P3 | six (D = 0) to seven (D = 2) reorganizations, non-monotone | 0.24 | 0.17 |
  | P4, P_fold | `4 -> inf` | 0.85 | 0.17 |

  AVR speed is therefore a coordinate that orders the regions, with a window
  (`beta` roughly 0.2–0.5) in which every point is viable **and** composable. It
  is not a switch that creates the phenomenon out of nothing.

**5. Does electromagnetic presence without inertia recover composability?** Yes.
The minimal condenser (1 % inertia, EMFs frozen, field held, no PSS/D/Q)
empties `H` at P2, P3, P4 and P_fold at 25 % rating and at P1 at full rating.
Along its rating (F8C), the last hyperedge disappears at

| point | rating at which `H` becomes empty |
|---|---|
| P_fold | 0.25 % of the retired machine |
| P4 | 2.5 % |
| P3 | 9.1 % |
| P2 | 19.5 % |
| P1 | between 25 % and 100 % |

and on the way the hypergraph **expands back**: at P2 the hyperedges vanish in the
order `33+37, 33+35, 30+35+37, 30+33, 33+35+37, 30+33+35, 30+33+37, 30+33+35+37`,
each at an imaginary-axis crossing with closure distance about `1e-8`. It is
*presence* — a stiff voltage source behind transient reactance — not dynamics:
freezing the EMFs changes nothing.

**6. Does inertia alone recover composability?** No. Synthetic inertia on the
converters (the retired machine's `H`, no electromagnetic presence) empties `H`
only at the marginal P_fold; P4 and P1 keep their hypergraphs, P2 keeps
`kappa = 2` with a smaller `H`, and P3 **gains** a hyperedge. Added to a minimal
condenser, full kinetic energy is neutral at P3/P4/P_fold and **harmful** at P2
(`EMPTY -> 30+33+35`). Surviving-fleet inertia is destabilizing in both
directions (x0.5 empties `H` at P1, P2, P4, P_fold; x2 makes the base unstable at
P1–P4 and P_fold worse, `kappa 4 -> 1`).

**7. Does PSS alter the hypergraph?** Rarely. Condenser PSS: never.
Surviving-fleet PSS off: margin-only at P1, P2, P4; adds two triples at P3
(`kappa` stays 3); empties `H` at P_fold (the power-input PSS is locally
destabilizing there).

**8. Which service changes the minimal witness coalition itself?** Changes of `H`
between two non-empty hypergraphs with different `kappa`-witnesses: condenser
inertia (32 pairs), EMF transients (22), damping (14), AVR (10) and reactive
share (10) at P1; damping and reactive share at P2 (4 each); surviving-fleet
inertia at P3/P4/P_fold, damping at P1–P3 and PSS at P1, P3 and P_fold.

## 3. Both directions

| claim | restore at the unsafe side | remove at the safe side | verdict |
|---|---|---|---|
| electromagnetic presence at retired buses is protective | minimal condenser: P2, P3, P4, P_fold -> `EMPTY`; P1 at full rating | removing the condenser from any safe condenser configuration returns the frozen unsafe `H` (same pairs); shrinking the rating reverses the hyperedge sequence at the same thresholds | **both directions** |
| converter voltage support is a substitute for it | F7: raising `g` at P4 -> `EMPTY` | F7: lowering `g` at P_inf (`g = 1 -> 0`, `k = 0.5`) -> `kappa = 4` | **both directions**, non-monotone in between (F7 tongue) |
| surviving-fleet inertia is destabilizing | x0.5 at P1, P2, P4, P_fold -> `EMPTY` | x2 at P_fold: `kappa 4 -> 1`; at P_inf margin worsens (`+0.020`) but stays safe | consistent in both directions |
| synthetic inertia is protective | only P_fold | — | **not supported** |
| surviving-fleet AVR speed drives non-composability | slowing it at P1–P4, P_fold -> `EMPTY` before the base fails | restoring native speed returns each frozen `H` | both directions along the viable window |

## 4. F8B — excitation mean versus heterogeneity

Matched cases `TE_i(m, h) = c(m, h) (TE_i/G)^h`, with `c` chosen so that the fleet
mean equals `m` under three definitions (geometric, arithmetic, harmonic = mean
exciter bandwidth), 33 means x 21 heterogeneity amplitudes, at `g = 0, 0.05, 0.2`;
likewise for the gain `KA`. 8 316 cases; 120 of 120 random cases agree with the
direct path.

| question | answer |
|---|---|
| does heterogeneity change `H` at fixed mean? | **yes**, on 17–33 of 33 mean levels, under every mean definition and every `g` |
| does the mean change `H` at fixed heterogeneity? | yes, on 19–21 of 21 heterogeneity levels |
| share of flagship-margin variance (TE): mean / heterogeneity / interaction | arithmetic 78–91 % / 3–11 % / 2–19 %; geometric 47–75 % / 2–21 % / 7–51 %; harmonic 17–55 % / 4–35 % / 12–79 % |
| sign of the TE-heterogeneity effect | with voltage support (`g = 0.05, 0.2`) `kappa` never decreases as heterogeneity grows (protective); at fixed Q it goes both ways depending on the mean level |
| sign of the KA-heterogeneity effect | `kappa` only decreases as the gain spread grows (destabilizing) |

**Heterogeneity is causal** in the interventional sense — changing it at a matched
mean changes `H` — but it is **not dominant**: the fleet mean carries most of the
margin variance, and "mean versus heterogeneity" is itself definition-dependent
(the interaction share is large when the geometric or harmonic mean is held).
N21 stands and is sharpened: neither "heterogeneity is required" nor
"heterogeneity is irrelevant" is supportable.

## 5. F8C — port closure at service-induced boundaries

For every factorial pair with different `H` (both ends feasible), the service was
moved continuously and every flipping subset bisected on the direct path
(296 paths, `results/F8C/F8C_events.csv`):

| | |
|---|---|
| located subset boundaries | 1 531 |
| imaginary-axis crossings | 1 511 |
| upper-band-edge re-classifications | 20 (condenser rating at P1/P2) |
| port-visible imaginary-axis crossings | **1 511 of 1 511** |
| multi-port witnesses with closure distance below `1e-4` | **1 193 of 1 193** (max `9.7e-6`) |
| single-port witnesses, individual factor | 318, at most `2.7e-6` |
| port-invisible transitions | **0** |

Port closure remains the right reduced descriptor whatever physical service moves
the boundary: presence, inertia, transients, AVR, damping, reactive share,
synthetic inertia.

## 6. What this changes

- O20 is sharpened: the service removed by retirement is electromagnetic
  presence, and it is dose-dependent.
- A new, IEEE-39-specific observation: surviving-fleet inertia is destabilizing
  for this mechanism. It is **not** a general statement about inertia.
- The phenomenon is an interaction between the loss of local voltage stiffness
  at the retired buses and the excitation dynamics of the remaining fleet; the
  converter's voltage regulation partially substitutes for the lost stiffness,
  non-monotonically.
