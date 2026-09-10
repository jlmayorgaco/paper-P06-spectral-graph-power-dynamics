# Gate 1 — Whole-RHP composability as the primary object

Journal gate after the post-F7 freeze (`IAS2026_TRACKA_F8_F12_POST_F7_FREEZE`).
Frozen results are read, never modified. Code: `experiments/G1_whole_rhp.py`,
`experiments/G1_f8_rhp_followup.py`, `experiments/G1c_f8d_rhp.py`. Data:
`results/G1/`. Manifests:
`outputs/ias2026_journal_gates/`.

    Gamma_RHP = { Re s > 0 }                        primary: planning safety
    Gamma_IA  = { Re s > 0, 0.3 <= f <= 1.5 Hz }    secondary: one declared mechanism

Both exclude the angle reference's double zero (`|s| > 1e-3`). The
discrepancy classes and their thresholds were fixed in the script header before
the run:

| class | definition |
|---|---|
| APERIODIC | an RHP eigenvalue with `f < 1e-3 Hz` (real) |
| BAND_EDGE_ARTIFACT | an oscillatory RHP eigenvalue within 0.05 Hz outside a band edge |
| OSCILLATORY | any other oscillatory RHP eigenvalue outside the band |
| NO_DISCREPANCY | `H_RHP = H_IA` |

A discrepancy is classified from the out-of-band RHP modes of the hyperedges in
`H_RHP \ H_IA`.

**Planning rule adopted from here on.** Every planning claim, meaning "this
portfolio is safe" or "this service restores composability", uses `H_RHP`. The
band object `H_IA` may be used only to study the declared inter-area mechanism.

## 1. IEEE-39: the frozen F7 maps

| map | points | `H_RHP = H_IA` | `kappa` equal | discrepancies | distinct `H_RHP` / `H_IA` |
|---|---|---|---|---|---|
| F7A (`g`, `k`) | 120 390 | **100 %** | 100 % | none | 30 / 30 |
| F7B (`g`, `t`) | 79 028 | **96.46 %** | 99.93 % | 2 800 BAND_EDGE_ARTIFACT | 36 / 33 |
| F7C (`g`, `h`) | 145 811 | **100 %** | 100 % | none | 16 / 16 |

- There are no aperiodic discrepancies and no out-of-band oscillatory
  discrepancies.
- The 2 800 F7B points are one mechanism: the same inter-area branch crosses
  the band edge while still in the RHP. The band misses those points, so the
  planning object there is `H_RHP`. It adds three hypergraphs not seen in the
  band (36 against 33).
- All 42 436 points with `H_RHP = EMPTY` are fully small-signal stable for every
  portfolio.
- On IEEE-39 the declared mechanism is therefore the only small-signal
  instability of the 16 portfolios in these slices, except near a band edge.

**Boundaries.** 128 RHP boundaries on pure-policy edges of F7B were located.
All 128 are oscillatory imaginary-axis crossings, and all are port-visible.

## 2. Kundur: the frozen F12 grids

| slice | nodes | `H_RHP = H_IA` | `kappa` equal | discrepancies | `H_RHP` values | `H_RHP = EMPTY` |
|---|---|---|---|---|---|---|
| K12A (`g`, `k`) | 3 721 | **8.2 %** | 84.4 % | APERIODIC 3 116; BAND_EDGE 129; OSCILLATORY 101; mixtures 69 | `2\|3+4`, `2+3\|2+4\|3+4`, `2\|3\|4`, `2\|4` | **0** |
| K12B (`g`, `t`) | 2 196 | **0.3 %** | 81.0 % | APERIODIC 1 930; BAND_EDGE 142; OSCILLATORY 16; mixtures 102 | the same four, plus `2` | **0** |

- Under both objects `H_RHP` changes with the reactive policy on every
  pure-policy line (61/61 and 36/36).
- `kappa_RHP` takes only the values 1 and 2.
- No policy in either slice makes the Kundur fleet composable.
- The band hypergraph reproduces the frozen F12 labels at 100 %, but for
  planning it is the wrong object on this grid.

**Boundaries.**

- 667 RHP boundaries were located: 473 oscillatory imaginary-axis crossings
  (single witnesses 2, 3, 4; 0.18–0.73 Hz) and 194 aperiodic crossings through
  the origin (witnesses `3+4` and `2+3+4`).
- All 473 oscillatory crossings are port-visible. The largest single-port
  individual factor is `7.8e-7`.
- The aperiodic crossings cannot be localized by the port closure. At `s = 0`
  the reduced port operator contains the reference zero, so no frequency-domain
  closure test is evaluable there. These crossings are reported and classified,
  not port-certified.

**Port closure over both benchmarks.** 601 of 601 oscillatory RHP boundaries
are port-visible, with a largest multi-port closure distance of `1.0e-6`.

## 3. F8 planning claims under Gamma_RHP

Recomputed over all 1 068 frozen (point, intervention) lattices:
`H_RHP = H_IA(frozen)` in **85.4 %** of cases. All 156 differences are
condenser configurations (class R) whose condenser has **zero damping**:

| condenser | fraction of configurations where `H_RHP != H_IA` |
|---|---|
| `D = 0`, frozen EMFs (classical) | 68.8 % |
| `D = 0`, dynamic EMFs | 12.5 % (only 1 % inertia at full rating) |
| `D = 2` (either EMF model) | **0 %** |

The plain replacement, converter synthetic inertia and every surviving-fleet
(class A) configuration are unchanged.

**Mode anatomy** (`results/G1/G1_f8_rhp_modes.csv`):

- All 228 out-of-band modes behind the differences are **the condenser's own
  electromechanical (swing) mode**: dominant states `delta` and `omega` of the
  condenser, condenser participation at least 98 %.
- They lie at 7.7–14.3 Hz, a high frequency because the condenser carries 1 %
  of the retired kinetic energy.
- Growth rates reach `+0.50 s^-1`.
- The other 40 modes are the same condenser mode just above 1.5 Hz (BAND_EDGE,
  `Re <= 0.0014`).
- An undamped low-inertia classical condenser is thus unstable on its own
  swing mode. The band study could not see this.

**Consequence for O72.** The "minimal condenser" of F8 — 1 % inertia, frozen
EMFs, field held, and no PSS, damping or reactive output — is **not** a safe
planning configuration. Under `Gamma_RHP`:

- It makes `H_RHP` non-empty at every point, including the safe `P_inf`.
- This holds at every rating from 0.1 % to 100 % (`results/G1/G1_em_presence_path_rhp.csv`).

**Restated planning claim.** Electromagnetic presence at the retired buses
restores whole-RHP composability once the condenser's own swing mode is damped
by either of two single physical services: mechanical/damper damping
(`D = 2 pu`) or rotor-flux dynamics.

- 80 condenser configurations give `H_RHP = EMPTY` at every unsafe point
  P2–P_fold, and every one of them has `D = 2` or dynamic EMFs. Having one of
  those two services is **necessary but not sufficient**. Of the 96
  configurations that have one, 16 stay non-empty somewhere:
  - 12 are band-level service effects already in F8, at P2 only: full condenser
    inertia with damping, and AVR plus reactive share with dynamic EMFs;
    `H_RHP = H_IA` there.
  - 4 are the undamped condenser swing mode again: full rating, 1 % inertia,
    dynamic EMFs, `D = 0`, reactive share on, giving `30+33` at P2–P_fold.
- The electromagnetic-presence thresholds, bisected on the rating with the
  per-unit services held fixed as in F8C, are:

| point | band threshold (F8C, classical `D = 0`) | RHP threshold, `+ D = 2` | RHP threshold, `+ dynamic EMFs` | hypergraph just below |
|---|---|---|---|---|
| P_fold | 0.25 % | 0.242 % | 0.251 % | `30+33+35+37` |
| P4 | 2.5 % | 2.48 % | 2.56 % | `30+33+35+37` |
| P3 | 9.1 % | 9.13 % | 10.5 % | `30+33+35+37` |
| P2 | 19.5 % | 18.8 % | 19.1 % | `30+33+35+37` |
| P1 | between 25 and 100 % | 44.8 % | 44.8 % | `30+33+37` |

- Time-domain check (G2, nonlinear phasor model): at P4 the `D = 2`
  condenser's inter-area `alpha` is +0.0264 at 2 % and −0.0294 at 3 %, which
  places the threshold at 2.47 % (bisected: 2.48 %). The undamped condenser's
  8.65 Hz swing mode grows at +0.034 s^-1, as predicted.
- The inter-area dose-response survives unchanged.
- What changes is the definition of a minimal viable condenser: the condenser
  must damp its own swing mode.
- Per factor level, the fraction of RHP-safe configurations is:

| factor | low level | high level |
|---|---|---|
| damping | 49 % (0) | 89 % (2) |
| EMF dynamics | 54 % | 84 % |
| inertia | 61 % (1 %) | 77 % (100 %) |
| PSS | no effect | no effect |
| reactive share | 73 % (0) | 66 % (1) |

The following F8 verdicts are unchanged under `Gamma_RHP`; their
configurations have identical `H_RHP` and `H_IA`:

- inertia alone does not restore composability;
- surviving-fleet inertia is destabilizing;
- manual excitation makes the base unstable.

**The F8D AVR-speed path** (`experiments/G1c_f8d_rhp.py`,
`results/G1/G1_f8d_rhp.csv`): 480 cases, with `H_RHP = H_IA` in 97.3 % of
them. Slowing the surviving AVRs still orders the regions. Near the slow end,
though, the inter-area family drops **below** the band, and single
replacements become RHP-unstable at 0.21–0.27 Hz, which the band cannot see.

| point | band: `H` empty from `beta` = | RHP: `H` empty from `beta` = | RHP-composable window before the base fails |
|---|---|---|---|
| P1 (D 0 / 2) | 0.53 / 0.73 | 0.53 / 0.73 | unchanged |
| P2 | 0.62 / 0.73 | 0.62 / 0.73 | unchanged |
| P3 | 0.24 | **never** | **none**: `30\|33\|35` at `beta = 0.24–0.28`, `30\|33\|35\|37` at 0.20 |
| P4 | 0.85 | 0.85 | narrowed: `H_RHP = 30` for `beta <= 0.24` (D 0) and `<= 0.28` (D 2) |
| P_fold, P_inf | as band | as band | unchanged |

F8's statement that there is "a window (`beta` about 0.2–0.5) in which every
point is viable and composable" therefore does **not** hold under whole-RHP
safety. P3 has no composable `beta`. The restated form is: AVR speed orders
the regions, and it reaches RHP-composability before the base fails at P1, P2,
P4 and P_fold, but not at P3.

## 4. Answer to the gate

| benchmark | `H_RHP` vs `H_IA` | discrepancy type | planning consequence |
|---|---|---|---|
| IEEE-39, F7A/F7C | identical | none | band results are planning results |
| IEEE-39, F7B | 96.5 % identical | band-edge artifact only (same inter-area branch) | use `H_RHP`; 3 extra hypergraphs |
| IEEE-39, F8 services | 85.4 % identical | out-of-band oscillatory (undamped condenser swing mode, 7.7–14.3 Hz) | O72 restated: condenser needs damping or flux dynamics; thresholds unchanged |
| IEEE-39, F8D AVR-speed path | 97.3 % identical | inter-area family below the band (0.21–0.27 Hz) at slow AVR | P3 never RHP-composable; the composable window at P4 narrows |
| Kundur | 0.3–8 % identical | mostly aperiodic (pairs of grid-following replacements lose synchronism), also out-of-band oscillatory and band-edge | `H_RHP` never empty; no policy makes the fleet composable |

`H_Gamma = EMPTY` does **not** imply small-signal stability unless `Gamma`
covers every instability mechanism of interest. Kundur is the counterexample
(`H_IA = EMPTY` at 556 K12A and 374 K12B nodes; `H_RHP` is never empty). So is
the undamped condenser on IEEE-39.
