# v2C — the corrected observable: three passes, and the order-4 claim survives

Protocol `configs/ias2026/v2c_modal_family_protocol.yaml`, frozen before the seed
was used. Seed 20260912. **400 draws, 295 accepted, 105 rejected, every rejection
for non-dispatchability, 0 tracking failures.**

One methodological change from v2B and nothing else: the stability endpoint is
the **envelope of the inter-area modal family**, not one selected descendant.

    alpha_IA(theta) = max over lambda in C_IA(theta) of Re(lambda)

v1, v2A and v2B stay frozen and untouched. H5B stays **non-executable**, not
falsified. The v2B `alpha_Q = 0` and `delta_alpha_Q = 0` curves stay **invalid**
as physical boundaries; nothing below rehabilitates them. H4B and H2B keep their
v2B wording and status and are not restated here.

| file | contents |
|---|---|
| `configs/ias2026/v2c_modal_family_protocol.yaml` | the frozen protocol |
| `src/ibr_cycles/dynamics/modal_family.py` | family definition and subspace metrics |
| `experiments/_v2c_common.py` | the frozen definitions, shared by all v2C scripts |
| `experiments/E29_family_continuation.py` | deterministic pre-registration evidence |
| `experiments/E28_v2c_modal_family.py` | the 400-draw campaign |
| `experiments/E30_v2c_numerics.py` | numerical validation of the port margin |
| `experiments/E31_v2c_order4_audit.py` | re-audit of the order-4 claim |

Tables follow the repository convention of carrying their own experiment prefix,
so the campaign writes `E28_v2c_modal_family_*.csv` and the two audits write
`E30_v2c_numerics_audit.csv` and `E31_v2c_order4_audit_*.csv`.

## 1–3. Sampling, tracking, and the shape of the family

| | |
|---|---|
| drawn / accepted | 400 / **295** |
| rejected | 105, **all** for non-dispatchability |
| tracking failures | **0** |
| family tracked in every accepted sample | **yes, 295 of 295** |
| number of descendants | **2 in 295 of 295** |
| per-member overlap to the reference | 0.9530 to 0.9895 |
| member independence, smallest singular value | min 0.0591, median 0.0958 |

The two descendants are **nearly collinear in shape** — an independence of 0.059
is an angle of about 4.6 degrees — while being **far apart in damping**, median
separation 0.548 in real part. That combination is precisely why a single-vector
assurance criterion cannot choose between them and why the choice matters so
much: the selector picks almost at random between two branches whose stability
verdicts differ by half a unit of damping.

Pre-registration evidence, `E29`, deterministic and drawing no random sample: at
all 175 dispatchable continuation points the family has exactly 2 members, and
the **subspace** it spans is continuous between consecutive steps at a worst
principal-angle cosine of **0.999952**, median 0.999989, while the argmax
selector changes sign 19 times over the same continuation.

Following the protocol, this is **not** called a mode splitting. The supported
wording is: *the base inter-area modal family develops two closely related
descendants.*

## 4. H4C — response surface on `alpha_IA`: **PASS**

Independent replication. H4B was not assumed to transfer.

| term | coefficient | 95 % interval | p |
|---|---|---|---|
| const | `+27.321` | `[+20.417, +34.226]` | `1.5e-13` |
| `L` | `-73.502` | `[-87.496, -59.509]` | `2.2e-21` |
| `A` | `+11.559` | `[+6.806, +16.312]` | `3.0e-06` |
| **`L*A`** | **`-13.710`** | **`[-18.124, -9.295]`** | **`3.7e-09`** |
| `L^2` | `+48.007` | `[+40.304, +55.711]` | `6.0e-28` |
| `A^2` | `+0.547` | `[-1.447, +2.541]` | `0.591` |

`R^2 = 0.902`, against 0.490 for the same design on the defective observable.
Removing the chattering roughly halves the unexplained variance, which is what a
tracking artefact looks like when it is removed.

Residual diagnostics: mean `-2.3e-11`, standard deviation `0.0719`, skew `+0.08`,
kurtosis `2.70`, Spearman of `|residual|` against fitted `+0.245`. Symmetric and
near-Gaussian, with mild heteroscedasticity that the confidence intervals do not
account for; the interaction p-value has nine orders of magnitude of margin, so
this does not touch the verdict.

**The load-by-availability interaction survives the correction.** `b3` keeps its
sign and its significance, `3.7e-09` against `5.9e-09` in v2B.

The logistic interaction is **not** significant, `p = 0.838`, with the same sign
as the continuous one. The protocol declares this an acceptable outcome and it is
reported as it stands. Binary significance was not forced.

Degradation fraction, `delta_alpha_IA > 0`:

| load \ availability | 0.75–0.85 | 0.85–0.92 | 0.92–0.96 | 0.96–1.00 |
|---|---|---|---|---|
| 0.90–0.95 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.95–0.98 | 0.900 | 0.650 | 0.500 | 0.500 |
| 0.98–1.00 | 1.000 | 1.000 | 1.000 | 1.000 |
| 1.00–1.02 | 1.000 | 1.000 | 1.000 | 1.000 |
| 1.02–1.05 | — | — | 1.000 | 1.000 |

Instability fraction, `alpha_IA > 0`:

| load \ availability | 0.75–0.85 | 0.85–0.92 | 0.92–0.96 | 0.96–1.00 |
|---|---|---|---|---|
| 0.90–0.95 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.95–0.98 | 0.450 | 0.250 | 0.050 | 0.000 |
| 0.98–1.00 | 1.000 | 0.950 | 0.850 | 0.700 |
| 1.00–1.02 | 1.000 | 1.000 | 1.000 | 1.000 |
| 1.02–1.05 | — | — | 1.000 | 1.000 |

Empty cells are absent dispatchable support, never stable operating points.

## 5. The physical boundary

This is the first campaign entitled to state one, because it is the first with a
well-defined endpoint. Upcrossing convention, drawn only over dispatchable
support, with percentile bootstrap bands at full support, 400 of 400 resamples.

| availability | dispatchable to | `delta_alpha_IA = 0` | 95 % band | **`alpha_IA = 0`** | 95 % band | worst descendant frequency at the boundary |
|---|---|---|---|---|---|---|
| 0.753 | 1.011 | 0.9459 | [0.9411, 0.9507] | **0.9602** | [0.9517, 0.9695] | 0.580 Hz |
| 0.814 | 1.011 | 0.9513 | [0.9489, 0.9534] | **0.9664** | [0.9626, 0.9709] | 0.575 Hz |
| 0.876 | 1.015 | 0.9568 | [0.9544, 0.9592] | **0.9732** | [0.9691, 0.9770] | 0.575 Hz |
| 0.938 | 1.028 | 0.9630 | [0.9606, 0.9650] | **0.9804** | [0.9773, 0.9832] | 0.572 Hz |
| 0.999 | 1.037 | 0.9695 | [0.9657, 0.9726] | **0.9880** | [0.9832, 0.9941] | 0.570 Hz |

Both boundaries lie strictly inside the dispatchable region at every
availability, and both move to **higher** load as availability rises. Higher
photovoltaic availability is not worse on this benchmark; it buys load headroom.

The deterministic continuation, estimated independently and never used as a
target, put the crossing near load 0.9765 at availability 0.92. The campaign puts
it at 0.9732 at availability 0.876 and 0.9804 at 0.938, which brackets it. The
two estimates agree without either being fitted to the other.

## 6. H2C — repair: **PASS**

The frozen voltage-control policy, unchanged.

| criterion | unstable | repaired | fraction | exact 95 % interval |
|---|---|---|---|---|
| `alpha_IA` | 140 | **140** | 1.0000 | Clopper-Pearson [0.9740, 1.0000] |
| worst mode over the whole band, no tracking at all | 140 | **140** | 1.0000 | [0.9740, 1.0000] |

Least negative repaired `alpha_IA`: `-0.0683`. Worst full-band real part under
voltage control **over all 295 accepted samples**: `-0.0683`. The two criteria
select the same 140 samples and agree completely, which means the repair result
does not depend on any tracking decision.

Fourth independent confirmation of the same conditional statement, after v1's
post-hoc P2, v2A's 63 of 63 and v2B's 96 of 96.

The wording is unchanged and remains required: voltage regulation restores the
**damping** at a **higher frequency**. Median worst-descendant frequency 0.642 Hz
in the reference, 0.575 Hz after replacement, 0.727 Hz under voltage control. It
does not restore the original mode.

## 7. H5C — port margin: **PASS**

Near-boundary stratum by fixed count, `k = 50`, frozen before execution, so it
cannot be empty. The two groups were checked to be disjoint before the test was
allowed to report.

| | n | `abs(alpha_IA)` | median `m_4` |
|---|---|---|---|
| near | 50 | ≤ 0.0809 | **0.0292** |
| far | 50 | ≥ 0.3118 | **0.1851** |

- primary 1, one-sided Mann-Whitney: `p = 3.5e-18` — **met**
- primary 2, Spearman `m_4` against `abs(alpha_IA)`: `+0.9183`, `p = 6.6e-120` — **met**
- preregistered convenience form, `log(m_4)` on `abs(alpha_IA)`: slope `+6.91`,
  `p = 3.0e-142`. This is **not** a claimed power law and no functional form was
  introduced after seeing the data
- secondary, `Delta_f = f_port - f_IA,worst` in the near group: median
  `abs(Delta_f) = 0.0050 Hz`, 95th percentile `0.0098 Hz`, against a 0.05 Hz
  threshold — **met**. Over all accepted samples the 95th percentile is 0.1622 Hz,
  as expected far from the boundary where the band minimum need not sit on the mode

Robustness check, option A, `abs(alpha_IA) <= 0.05`: 29 samples, `p = 8.6e-14`.
Declared in the protocol as a robustness check and not part of the verdict.

So the port closure tightens as the inter-area family approaches the axis, and
the frequency at which it tightens is the frequency of the descendant that is
crossing. That is the statement H5 was written to test in v2A and could not test
in either v2A or v2B.

## 8. Numerical validation — the margin is not a solve artifact

21 cases: 10 generic, 10 nearest the boundary, and the deterministic nominal
boundary. Extended precision was declared unavailable in advance — `mpmath` is
not installed and the pinned interpreter has no `pip` — so the perturbation
ensemble and the multi-algorithm comparison were preregistered as the substitute.

| check | worst over 21 cases |
|---|---|
| backward residual of the `T_0` solve | `3.0e-17` |
| LU vs QR vs SVD spread in `m_4` | `8.3e-14` |
| independent reimplementation of the closure algebra vs production | `0.0` |
| one step of iterative refinement | `8.1e-15` |
| perturbation ensemble at relative `1e-13`, worst induced shift | `3.1e-11` |
| `cond(T_0(j omega))` | `3.2e+03` |

The smallest margin audited is `1.115e-02`, which exceeds the worst numerical
shift by a factor of **3.6e+08**. `m_4` becoming small near closure is a property
of the operator, not of the solve.

**One real limitation, and it runs against the finding rather than for it.** The
frozen 61-point grid is adequate far from the boundary but loose near it. Golden
section refinement around the grid minimum gives:

| `alpha_IA` | `m_4`, 61-point grid | 121-point grid | locally refined |
|---|---|---|---|
| `+0.00032` | 0.03146 | 0.00194 | **0.00017** |
| `-0.00271` | 0.01115 | 0.01115 | **0.00131** |
| `-0.00439` | 0.01456 | 0.01456 | **0.00219** |
| `+0.01465` | 0.01472 | 0.01472 | **0.00724** |
| `-0.07359` | 0.04067 | 0.03232 | **0.03226** |
| `+0.35134` | 0.30127 | 0.30127 | **0.30088** |

Ratio of refined to grid margin: median 0.993 in the generic cohort, median 0.180
and minimum 0.005 in the near-boundary cohort. Refined, the Spearman against
`abs(alpha_IA)` rises from `+0.890` to `+0.988`, and at `alpha_IA = +0.00032` the
margin is `1.7e-04`, that is, the closure reaches `-1`.

The preregistered grid metric is therefore a **conservative upper bound whose
looseness grows exactly where the effect is**, so its bias works against H5C.
The verdict stands on the frozen metric and the metric was **not** swapped. The
refinement is reported as a numerical property, and it carries the caveat that
golden section assumes unimodality on the bracket; the monotone agreement between
the 61-point, 121-point and refined values is the evidence that it holds here.

## 9. Order-4 audit — the Track-A claim **survives**

All 16 subsets of the frozen flagship 30+33+35+37, re-run with `alpha_IA` and the
comparison basis pinned to the 12 flagship survivors so the threshold means the
same thing at every node. Family size 2 at all 16 subsets, minimum member overlap
0.9405, zero tracking failures.

| truncation | family envelope | | forward argmax selector | |
|---|---|---|---|---|
| order ≤ 0 | `-0.126478` | stable | `-0.126478` | stable |
| order ≤ 1 | `-0.207821` | stable | `-0.207821` | stable |
| order ≤ 2 | `-0.335225` | stable | `-0.335225` | stable |
| order ≤ 3 | `-0.288384` | stable | `-0.555177` | stable |
| **order ≤ 4, exact** | **`+0.144670`** | **UNSTABLE** | `-0.734356` | stable |
| irreducible `mu_4` | **`+0.433054`** | | `-0.179179` | |

Reconstruction residual `2.8e-17`. **The order-4 claim survives a strictly
stronger modal definition**, which is a stronger result than E15 had.

Branch-independent corroboration, depending on no tracking decision: the count of
right-half-plane modes inside the band is **0 for all fifteen proper subsets and
1 for the full portfolio**. Every single safe, every pair safe, every triple safe,
the portfolio unstable — as an integer count of unstable modes, not as a tracked
scalar.

The audit also shows how fragile the old observable was. A forward argmax
selector anchored on the base inter-area mode returns `-0.734` at the flagship, a
**false negative**: it follows the well-damped 0.66 Hz descendant while the
0.575 Hz one crosses. E15 escaped this only by anchoring backwards on the
flagship's own critical mode, which cannot fail at the flagship by construction
and therefore could not have revealed the problem. The family envelope gives
`+0.14467` at the flagship, agreeing with E15, and is well defined at every
subset in either direction.

## 10. Claim status

**VALIDATED — supported by v2C**

- `V1` the load-by-availability interaction on `delta_alpha_IA`, `b3 = -13.71`
  `[-18.12, -9.30]`, `p = 3.7e-09`, `R^2 = 0.902`
- `V2` the physical inter-area stability boundary `alpha_IA = 0`, from load 0.960
  at availability 0.75 to load 0.988 at full availability, with bootstrap bands,
  inside dispatchable support, at a worst-descendant frequency of 0.570–0.580 Hz
- `V3` higher photovoltaic availability moves both boundaries to higher load
- `V4` the frozen repair stabilizes 140 of 140 unstable samples on `alpha_IA`
  and, independently, 140 of 140 on the whole band, exact interval [0.974, 1.000]
- `V5` the port closure margin tightens as the family approaches the axis,
  Spearman `+0.918`, and `omega_port` agrees with the crossing descendant to a
  median of 0.005 Hz
- `V6` the base inter-area family develops two closely related descendants,
  2 of 2 at all 295 samples and all 16 lattice subsets, nearly collinear in shape
  and separated by a median of 0.548 in damping

**SURVIVES CORRECTION — old claims that remain true under `alpha_IA`**

- `S1` the order-4 claim: orders 1, 2 and 3 stable, exact order 4 unstable,
  `mu_4 = +0.433`, plus band right-half-plane counts 0 for all proper subsets
  and 1 for the portfolio
- `S2` conditional repair, now confirmed a fourth time on fresh seeds
- `S3` the existence of a load-by-availability interaction and of a frontier
- `S4` voltage regulation restores damping at a higher frequency and does not
  restore the original mode

**INVALIDATED OBSERVABLE**

- `I1` every `alpha_Q` physical-boundary claim from v2B: the `alpha_Q = 0` and
  `delta_alpha_Q = 0` curves and their bands. Not rehabilitated by v2C
- `I2` any use of a forward argmax single-branch selector as a stability
  endpoint on this system

**NON-EXECUTABLE**

- `N1` H5B, and H5B only. The stratum was empty; the hypothesis was never tested.
  H5C is a new test on a new observable and fresh seeds, not a re-scoring of H5B

**FAILED — genuine falsifications, unchanged**

- H2b, modal restoration on the complex branch distance, refuted in v2A
- the v1 unconditional claim, refuted at the held-out gate, still withdrawn

## 11. Poster safety

Safe to state, with the operating condition attached: the interaction, the
`alpha_IA = 0` boundary with its band, the repair result with its exact interval,
the order-4 result with the branch-independent mode count, and the port-closure
correlation.

Not safe, and not to be written: anything from the v2B `alpha_Q` boundaries;
"higher photovoltaic availability is worse"; "voltage control restores the
original mode"; "mode splitting"; a power law for `m_4`; and any statement that
H5B was refuted.

**No poster figure has been touched.**

## 12. Status and stop

- v1, v2A, v2B: frozen, unchanged.
- **v2C: complete.** H4C PASS, H2C PASS, H5C PASS, order-4 survives, numerics
  clean.
- Per the protocol stop rule, **no further Monte Carlo has been launched.**
