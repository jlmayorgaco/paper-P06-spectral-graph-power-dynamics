# E36 — held-out controller Monte Carlo: **PASS**, and it is decisive

Seed **20260914**, never used before. 1500 draws, **1500 accepted, 0 rejected**.
The operating point is the frozen nominal one; only the converter is uncertain,
so controller uncertainty is separated from operating-point uncertainty rather
than mixed with it.

Coordinates are physical, not raw gains — sampling `kp` and `ki` independently
would wander outside the region where the loops mean anything:

| coordinate | range |
|---|---|
| PLL natural frequency | nominal ± 20 % |
| PLL damping ratio | 0.50 to 1.00 |
| outer active bandwidth | ± 25 % |
| outer reactive bandwidth | ± 25 % |
| measurement bandwidth | ± 25 % |
| current-loop bandwidth | ± 15 % |

Every draw carries the whole 16-subset lattice — 24 000 case solutions in total.

## Result

| | successes | estimate | 95 % exact interval |
|---|---|---|---|
| all proper subsets stable | **1500 / 1500** | 1.000 | [0.9975, 1.000] |
| full portfolio unstable | **1500 / 1500** | 1.000 | [0.9975, 1.000] |
| **genuine order-4 crossing** | **1500 / 1500** | **1.000** | **[0.9975, 1.000]** |
| inter-area family retained | **1500 / 1500** | 1.000 | [0.9975, 1.000] |
| condenser repair succeeds | **1500 / 1500** | 1.000 | [0.9975, 1.000] |

`α_IA` of the flagship ranges from **+0.0799 to +0.2090** across the whole
controller box — unstable at every single draw, never marginal. Median `m₄` is
0.0919. Minimum member overlap 0.800; family size 2 or 3.

## What this settles

**The Track-A effect is not a property of one controller tuning.** Over a
six-dimensional physical controller box, every draw reproduces the full pattern:
each proper subset stable, the portfolio unstable, and the reconstruction through
order 3 stable while the exact portfolio is not. Not 95 % of draws — all of them.

This is the strongest available evidence for the claim that Track A is **not
controller-caused**. If the failure were a controller interaction, a ±20 % change
in PLL bandwidth and ±25 % changes in three outer loops would move it; instead
the damping wanders between +0.08 and +0.21 and the structure never changes.

It also means the **converter is not where the fix has to come from**: the
condenser repair, which does not touch the converter at all, succeeds in every
one of the 1500 draws.

Note the contrast with E35, where the *operating point* was uncertain: there the
portfolio was unstable at only 24 % of dispatchable points. Controller
uncertainty does not move this effect; operating condition does. That is the
central engineering message and the two Monte Carlos separate it cleanly.

## Files

`E36_MC_controller.csv` (+ `.parquet`), `E36_MC_controller_summary.csv`,
`manifest.json`.
