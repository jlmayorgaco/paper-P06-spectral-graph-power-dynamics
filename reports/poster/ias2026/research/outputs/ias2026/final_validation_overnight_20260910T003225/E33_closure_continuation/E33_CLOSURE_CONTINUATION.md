# E33 — closure continuation against the full-system boundary: **PASS**

This is the candidate main mathematical figure for a Transactions submission.

Deterministic continuation, 427 points over load 0.900–1.050 in steps of 0.0025
and photovoltaic availability 0.70/0.75/0.80/0.85/0.90/0.95/1.00. **290 points
dispatchable and evaluated**, the rest outside the fleet's envelope. Frozen v2C
modal-family observable, frozen 61-point closure metric.

## The closure sees the boundary the eigenproblem sees

| | value |
|---|---|
| Spearman, frozen `m₄` against `abs(α_IA)` | **+0.861**, p = 1.5e-86 |
| Spearman, locally refined `m₄` | +0.861 |
| `log m₄` on `abs(α_IA)`, preregistered convenience form | slope +7.20 |
| median `abs(Δf)` = `f_port − f_IA`, all points | **0.0075 Hz** |
| median `abs(Δf)`, within 0.05 of the boundary | **0.0049 Hz** |
| 95th percentile `abs(Δf)` near the boundary | **0.0093 Hz** |
| smallest `m₄` near the boundary, frozen metric | 0.00488 |
| smallest `m₄` near the boundary, refined | **0.00116** |
| closest approach of the interaction eigenvalue `μ` to −1 | **0.00488** |
| worst `cond(T₀(jω))` | 8.8e+03 |
| worst backward residual of the `T₀` solve | **3.6e-17** |

## The boundary located two independent ways

For each availability, the load at which `α_IA` crosses zero, against the load at
which `m₄` reaches its minimum:

| availability | `α_IA = 0` at load | `m₄` minimum, frozen | error | `m₄` minimum, refined | error |
|---|---|---|---|---|---|
| 0.70 | 0.9475 | 0.9475 | **0.0000** | 0.9500 | 0.0025 |
| 0.75 | 0.9550 | 0.9575 | 0.0025 | 0.9575 | 0.0025 |
| 0.80 | 0.9625 | 0.9650 | 0.0025 | 0.9650 | 0.0025 |
| 0.85 | 0.9700 | 0.9700 | **0.0000** | 0.9700 | **0.0000** |
| 0.90 | 0.9750 | 0.9775 | 0.0025 | 0.9750 | **0.0000** |
| 0.95 | 0.9800 | 0.9825 | 0.0025 | 0.9800 | **0.0000** |
| 1.00 | 0.9825 | 0.9825 | **0.0000** | 0.9825 | **0.0000** |

**Median boundary error: one grid step, 0.0025 in load, on the frozen metric; and
exactly zero on the refined metric.** The port-space closure minimum locates the
full-system stability boundary to a quarter of one percent of load.

The refinement is reported as a **secondary numerical validation only**. The
frozen 61-point metric was not replaced, and it is the metric the verdict rests
on. E40 confirms that the frequency grid, not floating point, is the dominant
error in the frozen metric: grid variation 2.3e-2 against a worst numerical shift
of 2.8e-10.

## Figures

`E33_closure_tracks_boundary.png` — three aligned panels along the availability
0.90 slice: `α_IA` crossing zero, `m₄` collapsing on a log axis, and `f_IA`
against `f_port`, with the crossing marked on all three.

`E33_mu_approaches_minus_one.png` — the interaction eigenvalue nearest `−1` in
the complex plane over all 290 points, coloured by `α_IA`, with the continuation
traces drawn; and the `m₄` against `abs(α_IA)` scatter on log axes.

## What may be claimed

That the gauge-invariant port closure `d_cl = min_μ |μ + 1|`, evaluated on the
imaginary axis over the frozen inter-area band, reaches its minimum at the
operating point where the inter-area modal family crosses the imaginary axis, to
within one continuation step; and that the frequency at which it does so is the
frequency of the crossing mode, to a median of 0.005 Hz.

That is a statement about an operator on eight action coordinates predicting the
stability boundary of a 110-state system.

## Files

`E33_closure_continuation.csv` (+ `.parquet`), `E33_boundaries.csv`,
`E33_closure_tracks_boundary.png`, `E33_mu_approaches_minus_one.png`,
`E33_figure_source.csv`, `manifest.json`.
