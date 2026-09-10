# v2B — boundary campaign: two passes, one non-executable test, and a defect in the observable

Protocol `configs/ias2026/trackA_v2B.yaml`, seed 20260911, experiment `E28`.
400 draws, **293 accepted**, 107 rejected, **every rejection for
non-dispatchability** — a resource boundary, never counted as an instability. No
v1 or v2A sample is reused. v1 stays frozen and refuted; the v2A global criterion
stays unsatisfied.

Audit: 0 of 293 accepted samples has a generator on a reactive limit; minimum MAC
is 0.966; worst `cond(T_0(j omega))` over the band is `8.91e+03`, so the port
operator is well conditioned everywhere the margin was measured.

The three verdicts are **independent**. There is no combined criterion, and H5B
is not required to rescue anything.

| | verdict | one line |
|---|---|---|
| **H4B** response surface | **PASS** | the load-by-availability interaction is real and both zero-level curves lie inside the dispatchable region |
| **H2B** conditional repair | **PASS** | 96 of 96 unstable samples are stabilized by voltage regulation |
| **H5B** port margin | **FAIL** — non-executable | the preregistered near-boundary stratum is empty, 0 of 293 |

And one finding that was not a hypothesis at all, reported in full below: **the
preregistered observable `alpha_Q` is defective**, and that defect is what
emptied the H5B stratum.

## H4B — response surface: PASS

Ordinary least squares for `delta_alpha_Q` over the dispatchable support, with
the analogous logistic model for the indicator of `alpha_Q > 0`. Estimators are
implemented in `ibr_cycles.uncertainty.regression` and checked in
`tests/test_regression.py` against the closed-form simple regression, the exact
noise-free fit, the logistic score equation and the observed information.

| term | coefficient | 95 % interval | p |
|---|---|---|---|
| const | `+116.585` | `[+91.056, +142.114]` | `4.5e-17` |
| `L` | `-285.349` | `[-335.951, -234.747]` | `6.8e-24` |
| `A` | `+41.086` | `[+23.452, +58.720]` | `7.4e-06` |
| **`L*A`** | **`-49.129`** | **`[-65.172, -33.087]`** | **`5.9e-09`** |
| `L^2` | `+174.113` | `[+146.543, +201.684]` | `1.7e-28` |
| `A^2` | `+2.506` | `[-4.701, +9.713]` | `0.496` |

`R^2 = 0.490`; logistic pseudo-`R^2 = 0.686`. The preregistered PASS condition is
the interaction coefficient `b3` differing from zero at the 5 % level **and**
both zero-level curves existing inside the dispatchable region. Both hold.

The logistic interaction is **not** significant, `p = 0.682`. So the interaction
is established for the size of the damping shift, not for the binary outcome.
That is reported as it stands and not smoothed over.

### The three boundaries, kept apart

| availability | dispatchable up to | `delta_alpha_Q = 0` at load | 95 % band | `alpha_Q = 0` at load | 95 % band |
|---|---|---|---|---|---|
| 0.751 | 1.007 | 0.9629 | [0.9543, 0.9701] | 0.9570 | [0.9505, 0.9643] |
| 0.813 | 1.008 | 0.9746 | [0.9708, 0.9784] | 0.9701 | [0.9670, 0.9739] |
| 0.875 | 1.013 | 0.9867 | [0.9832, 0.9904] | 0.9829 | [0.9798, 0.9867] |
| 0.938 | 1.038 | 0.9991 | [0.9963, 1.0025] | 0.9953 | [0.9918, 0.9987] |
| 1.000 | 1.029 | 1.0118 | [1.0053, 1.0177] | 1.0059 | [0.9973, 1.0146] |

Both curves lie strictly below the dispatchable edge at every availability, so
they are located inside physically feasible support and nothing is drawn where
the fleet cannot serve the load.

**Higher availability is not worse.** Both boundaries move to *higher* load as
availability rises: at 0.751 availability the effect changes sign at load 0.963,
at full availability at load 1.012. The negative `L*A` coefficient sits on top of
a positive `A` coefficient, and their combination over the sampled range means
more photovoltaic resource buys headroom in load, with diminishing return. Any
wording that reads "more inverter-based resource is worse" is contradicted by
this table.

Degradation fraction, tracked observable:

| load \ availability | 0.75–0.85 | 0.85–0.92 | 0.92–0.96 | 0.96–1.00 |
|---|---|---|---|---|
| 0.90–0.95 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.95–0.98 | 0.400 | 0.050 | 0.000 | 0.000 |
| 0.98–1.00 | 1.000 | 0.700 | 0.300 | 0.200 |
| 1.00–1.02 | 1.000 | 1.000 | 0.833 | 0.650 |
| 1.02–1.05 | — | — | 1.000 | 1.000 |

Cell counts thin out toward the top right because that is where dispatchability
runs out: the non-dispatchable share is 0.000, 0.000, 0.075, 0.388 and 0.875 by
load stratum. Empty cells are absent support, not stable operating points.

### A correction inside the reporting code

The fitted `delta_alpha_Q` surface is convex in load, so it has a second root on
the far side of its vertex. At high availability that root lands on the first
grid point, which is the fit reaching the edge of its own data rather than a
located boundary. `zero_curve` now returns the **upcrossing** — the crossing from
improvement to degradation as load rises, which is what the protocol asks for.
The verdict is unchanged either way, since it rests on `b3` and on the curves
lying inside the dispatchable region, and both held before and after. No
threshold was touched.

## H2B — conditional repair: PASS

Among the **96** samples where the reactive-setpoint replacement leaves the
tracked branch unstable, voltage regulation gives `alpha_V < 0` in **96 of 96**,
fraction `1.000`, bootstrap interval `[1.000, 1.000]`. The preregistered PASS
condition is a lower bound above 0.90.

This is the third independent confirmation of the same conditional statement:
v1 post-hoc observation P2, v2A preregistered H2 at 63 of 63, and now v2B at
96 of 96 on fresh seeds.

The diagnostic below strengthens it further. Voltage regulation does not merely
stabilize the tracked branch: over all 293 accepted samples the **worst** mode in
the whole 0.3–1.5 Hz band under voltage control has real part at most `-0.0682`.
The band is stable everywhere, not just the branch that was being followed.

The v2A wording stands unchanged and remains the required wording. Voltage
regulation restores the **damping**, at a **higher frequency**: median band
frequency 0.642 Hz in the reference, 0.576 Hz under setpoint replacement,
0.727 Hz under voltage control. It does not restore the original mode.

## H5B — port margin: FAIL, and the reason is not the port margin

The preregistered primary test compares `m_4` between samples with
`abs(alpha_Q) <= 0.05` and samples with `abs(alpha_Q) >= 0.20`.

**The near-boundary stratum contains 0 of 293 samples.** The smallest
`abs(alpha_Q)` observed is `0.0980`. The Mann-Whitney statistic is undefined, so
the primary condition is not met and the frozen protocol returns **FAIL**.

This is a **non-executable test, not an empirical refutation of the port
mechanism**, and it must be quoted that way. Unlike v2A's H5, the observable
`m_4` is not degenerate here: it ranges over `0.0067` to `0.3570` across the
campaign, at condition number below `9e+03`. The port margin measured something.
The stratum that was supposed to test it was empty.

## The defect: the replacement splits the branch, and the tracker chatters

`alpha_Q` never comes within 0.05 of zero because the two sides of the transition
sit in **disjoint frequency clusters**: 197 samples at 0.651–0.708 Hz with
`alpha_Q` in `[-0.790, -0.259]`, and 96 samples at 0.570–0.585 Hz with `alpha_Q`
in `[+0.098, +0.435]`. No overlap, and a gap in damping of width 0.36.

A deterministic load sweep at fixed availability, run outside the campaign,
explains it. In the **base** case the band holds one inter-area branch at
0.644 Hz. After replacement the band holds **two** branches of
indistinguishable shape — 0.566 Hz and 0.666 Hz, MAC 0.972 and 0.982 to the same
base mode. As load rises the 0.666 Hz branch becomes steadily better damped,
`-0.475` to `-0.890`, while the 0.575 Hz branch crosses zero smoothly at load
`0.9765`. Their MAC values cross around load 0.988, within about 0.001 of each
other, so the argmax-MAC selector flips back and forth:

| load | worst mode in the band | branch the tracker reports |
|---|---|---|
| 0.9775 | `+0.0120` | `-0.6158` at 0.662 Hz |
| 0.9800 | `+0.0348` | `+0.0348` at 0.575 Hz |
| 0.9825 | `+0.0571` | `-0.6544` at 0.662 Hz |
| 0.9850 | `+0.0819` | `+0.0819` at 0.576 Hz |
| 0.9875 | `+0.1075` | `-0.6998` at 0.662 Hz |

`alpha_Q` is therefore a **chattering selection between two branches**, not a
continuous curve. The physical stability boundary is smooth and well defined; the
gap around zero was manufactured by the selector.

Diagnostics recorded for every accepted sample, feeding no verdict, using the
worst mode in the band among those above the same MAC floor of 0.80:

| diagnostic | value |
|---|---|
| band modes of comparable shape, base case | 1.41 mean |
| band modes of comparable shape, after replacement | **2.00 mean** |
| samples whose band carries two such branches | **293 of 293** |
| tracked branch and band envelope disagree in **sign** | **43 of 293** |
| smallest `abs(alpha_Q)`, tracked | `0.0980` |
| smallest `abs(alpha_Q)`, band envelope | **`0.0013`** |
| unstable by band envelope | 139, against 96 tracked |
| of those, whole band stable under voltage control | **139 of 139** |

So the tracked observable **under-reports instability in 43 samples, 14.7 % of
the campaign**, and places the stability boundary optimistically high by 0.002 to
0.022 in load, the error growing with availability:

| availability | `alpha_Q = 0` at load, tracked | same, band envelope |
|---|---|---|
| 0.751 | 0.9570 | 0.9553 |
| 0.875 | 0.9829 | 0.9753 |
| 1.000 | 1.0059 | 0.9842 |

### What this does and does not touch

- **H2B is unaffected and strengthened.** Under either observable every unstable
  sample is repaired, 96 of 96 tracked and 139 of 139 by band envelope.
- **H4B's interaction coefficient survives in sign and in significance**, but the
  `alpha_Q = 0` curve must **not** be quoted as the physical stability boundary.
  It is the boundary of a chattering indicator and it is biased optimistic. The
  `delta_alpha_Q = 0` curve carries the same bias, since it is built on the same
  `alpha_Q`.
- **H5B's FAIL is fully explained** by the defect and carries no information
  about the port mechanism.

### Descriptive, explicitly not a test, explicitly not a retroactive pass

Repeating the H5B comparison on the band envelope instead of the tracked branch —
a **non-preregistered observable**, so this changes no verdict and passes
nothing:

| | n | median `m_4` |
|---|---|---|
| `abs(alpha)` ≤ 0.05 | 31 | `0.0221` |
| `abs(alpha)` ≥ 0.20 | 164 | `0.1710` |

one-sided Mann-Whitney `p = 5.7e-19`; median `abs(f_port - f_band)` near the
boundary `0.0047` Hz, against the 0.05 Hz secondary threshold; Spearman between
`m_4` and `abs(alpha)` is `+0.938`.

These numbers are recorded **only** to establish that a corrected campaign is
worth running. **H5 is not passed, not partially passed, and not rescued.** A
hypothesis tested on an observable chosen after seeing the data is not tested at
all. The correct disposition is a new preregistration on fresh seeds.

## Proposed v2C, awaiting authorization

Not authorized and not run. One change only: replace the argmax-MAC single-branch
observable with the **band envelope**, the worst real part among band modes above
the MAC floor, which is well defined whether the branch splits or not and which
reduces to the current observable when it does not. Everything else — strata,
nuisance variables, feasibility gates, `m_4` and its band — stays as frozen here,
with a new seed and no reuse of any v2B sample.

The v2B verdicts above stand as recorded regardless of whether v2C runs.

## Status

- v1: frozen, refuted, unchanged.
- v2A: frozen, global criterion not satisfied, nothing passed retroactively.
- **v2B: frozen.** H4B PASS, H2B PASS, H5B FAIL as non-executable. The
  observable defect is recorded as a limitation of v2B itself, not as a reason to
  revise any of the three verdicts.
- Poster: **untouched**, as required until v2B is frozen.
