# E37 — machine and synchronous-support Monte Carlo: **MIXED**, and it is the sensitive axis

Seed **20260915**, never used before. 1000 draws, **all accepted**, family
retained **1000 of 1000**, family size 2 in every draw, minimum overlap 0.886.

The operating point and the converter are the frozen nominal ones; the
**synchronous machines** are uncertain — inertia ±20 %, transient reactances
±10 %, AVR gain and time constant ±15 %, PSS gain and time constants ±20 %.
Machine damping `D` is **zero throughout the source case**, so a ±25 %
perturbation of it is identically no perturbation; it is reported as not
applicable rather than sampled and presented as though it had been varied.

## Result

| | estimate | 95 % exact interval |
|---|---|---|
| all proper subsets stable | **0.721** | [0.692, 0.749] |
| full portfolio unstable | **0.710** | [0.681, 0.738] |
| **genuine order-4 crossing** | **0.418** | [0.387, 0.449] |
| order-4 **given** the portfolio is unstable | **0.589** | [0.552, 0.625] |
| inter-area family retained | **1.000** | [0.996, 1.000] |
| condenser repair succeeds given unstable | **0.966** | [0.950, 0.978] |

`alpha_IA` of the flagship ranges from **−0.266 to +0.575** — stable in 29 % of
draws, unstable in 71 %.

## This is the axis that moves the result

Set against E36, the contrast is the whole point:

| uncertainty varied | genuine order-4 crossing |
|---|---|
| **converter tuning** (E36, 1500 draws) | **1500 of 1500, 100 %** |
| **machine parameters** (E37, 1000 draws) | **418 of 1000, 42 %** |
| **operating point** (E35, 1004 accepted) | 209 of 1004, 21 % |

The effect is completely insensitive to how the converters are tuned and strongly
sensitive to the synchronous machines' own parameters. That is **consistent with
the physical claim** — this is a loss-of-synchronous-support phenomenon, so the
machines' dynamics are exactly what should matter — and it simultaneously
**limits the generality** of the order-4 structure: with realistic uncertainty in
inertia, reactance and excitation data, the irreducible fourth-order pattern
holds in fewer than half the draws.

**Not explained by a degenerate base case.** The base system itself is unstable in
only **2 of 1000** draws (`alpha_base` up to `+0.024`), and conditioning on a
stable base changes nothing: 0.7224 against 0.7210 for proper-subset stability,
0.4188 against 0.4180 for the order-4 crossing.

Base-case damping is most sensitive to the AVR time constant (Spearman −0.51),
inertia (+0.49), transient reactance (+0.48) and AVR gain (+0.45).

**The condenser repair is nearly but not entirely robust here**: it fails in 24 of
710 unstable draws, 3.4 %.

## What may be claimed

That the order-4 structure is a property of *this machine data*. A Transactions
paper must either state the machine parameters as given and the result as
conditional on them, or characterise the region of machine-parameter space in
which the structure holds. It may not be presented as robust to machine
uncertainty.

## Files

`E37_MC_machine.csv` (+ `.parquet`), `E37_MC_machine_summary.csv`,
`E37_MC_machine_conditional.csv`, `manifest.json`.
