# v2A — stratified campaign: three of four hypotheses pass, and two real negatives

Protocol `configs/ias2026/trackA_v2.yaml`, seed 20260910. 360 draws, **264
accepted**, 96 rejected, every rejection for "not dispatchable". No v1 sample is
reused. v1 remains frozen and refuted.

Dispatcher rebuilt before sampling: resource-capping availability, iterative
saturation with redistribution among units that still have headroom, the slack
taking only the final residual, and a reactive-limit audit. Eight tests cover
power conservation, active limits, availability semantics, the slack not acting
as a default balancer, and reproducibility.

**Limit audit: 0 of 264 accepted samples has a generator on a reactive limit and
the slack is within limits in all 264.** The high-load points are genuinely
active-limited, not hiding a reactive infeasibility.

## The map

Fraction of samples where the replacement **degrades** the tracked branch:

| load \ availability | 0.70–0.80 | 0.80–0.90 | 0.90–0.95 | 0.95–1.00 |
|---|---|---|---|---|
| 0.85–0.95 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.95–1.00 | 0.750 | 0.345 | 0.200 | 0.067 |
| **1.00–1.05** | — | **1.000** | **0.846** | **0.833** |

Median `delta_alpha_Q`:

| load \ availability | 0.70–0.80 | 0.80–0.90 | 0.90–0.95 | 0.95–1.00 |
|---|---|---|---|---|
| 0.85–0.95 | `-0.136` | `-0.131` | `-0.137` | `-0.148` |
| 0.95–1.00 | `+0.294` | `-0.367` | `-0.342` | `-0.405` |
| 1.00–1.05 | — | `+0.351` | `+0.392` | `+0.399` |

Samples per cell: 30 / 30 / 30 / 30, 20 / 29 / 30 / 30, 0 / 4 / 13 / 18.

## Preregistered hypotheses

| | result | evidence |
|---|---|---|
| **H1** operating dependence | **PASS** | top cell 0.833 [0.667, 1.000] against bottom cell 0.000 [0.000, 0.000], non-overlapping |
| **H2** conditional mitigation | **PASS** | among the 63 samples with `alpha_Q > 0`, voltage control gives `alpha_V < 0` in **63 of 63**, 1.000 [1.000, 1.000] |
| **H2b** modal restoration | **FAIL** | `D_V < D_Q` in 0.144 [0.102, 0.186] |
| **H4** frontier | **PASS** | 3 cells above one half, 6 below |
| **H5** closure tracks instability | **INVALID, non-informative by construction** | median `d_4` is `0.0000` in every band; the correlation was never evaluable |

Declared success was H1 and H2 and H4 and H5. **The preregistered global PASS
criterion is NOT SATISFIED and nothing is passed retroactively.** H1, H2 and H4
pass on fresh samples under a corrected dispatcher; H2b is refuted; H5 was not
evaluable with the frozen definition.

H5 was **not falsified empirically**. The observable was mis-specified: it is an
identity, not a measurement. That invalidates the chosen observable, it does not
bear on the port mechanism either way.

H2 is the important one: it is the v1 post-hoc observation P2, restated as a
preregistered conditional hypothesis and confirmed on 63 new samples, 63 of 63,
with the unstable samples sitting at median load 0.998 and median availability
0.904.

## Negative 1 — voltage control restores the damping but moves the mode

| quantity | median |
|---|---|
| `D_Q = abs(lambda_Q - lambda_base)` | `0.401` |
| `D_V = abs(lambda_V - lambda_base)` | **`0.522`** |
| damping distance, setpoint | `0.298` |
| damping distance, voltage control | **`0.021`** |
| branch frequency, reference | 0.6449 Hz |
| branch frequency, setpoint | 0.6706 Hz |
| branch frequency, voltage control | **0.7281 Hz** |

Voltage control brings the **damping** back to essentially the synchronous
reference — the damping distance falls from `0.298` to `0.021`, and the auxiliary
damping-only measure holds in 264 of 264 — while pushing the **frequency**
further away, from `0.047` Hz off the reference to `0.083` Hz.

So the controller does not put the branch back where it was. It produces a
differently placed, well damped mode: it stiffens the mode as well as damping it.
The real-part-only measure would have reported perfect restoration and hidden
this entirely. Requiring the complex branch position was the right call.

The honest wording is therefore **not** "voltage control restores the
synchronous behaviour of the branch". It is "voltage control restores the damping
of the branch to the synchronous reference while raising its frequency".

## H5 — invalid by construction, and that is my error

Median `d_4` by band of `alpha_Q`:

| band | n | median `d_4` |
|---|---|---|
| `alpha_Q < -0.2` | 201 | `0.0000` |
| `0 < alpha_Q < 0.2` | 19 | `0.0000` |
| `alpha_Q > 0.2` | 44 | `0.0000` |

`d_4` is zero everywhere, so the correlation test was never going to mean
anything. The reason is algebraic, not numerical: `d_4` is evaluated **at**
`lambda_Q`, and `lambda_Q` is by construction an eigenvalue of the replaced
system, so `det(I + M(lambda_Q)) = 0` identically and some `mu` equals `-1` at
every operating point. Evaluating the closure at the system's own eigenvalue can
only ever return zero.

E23 worked because it evaluated the closure at a **fixed reference point** — the
flagship crossing eigenvalue — while varying the portfolio, so each subset's
operator was probed away from its own spectrum. The descriptor is a
portfolio-order coordinate measured at a common probe, and it was specified here
as if it were an operating-condition coordinate measured at a moving probe.

My first proposed correction — evaluate at each sample's own `lambda_base` — was
**rejected on review and correctly so**. `Q_4` is built on the resolvent of the
base operator, and `det T_0(lambda_base) = 0` by definition of an eigenvalue, so
that probe is singular or badly conditioned. It would have replaced a degenerate
definition with a different degenerate definition.

The correction adopted for v2B is an **imaginary-axis port closure margin**,
minimised over the frozen inter-area band, which sits away from the spectrum of
both systems and carries a direct physical meaning at an oscillatory boundary.
See `configs/ias2026/trackA_v2B.yaml`.

## Status

- v1: frozen, refuted, unchanged.
- v2A: frozen. H1, H2 and H4 pass on fresh samples; H2b is refuted; H5 is invalid
  by construction. The preregistered global PASS criterion is NOT satisfied and
  nothing is passed retroactively.
- v2B: authorized as a NEW preregistered follow-up, justified by H1, H2 and H4
  reproducing on fresh data. It is not an attempt to make the v2A criterion pass.
- Poster: untouched.
