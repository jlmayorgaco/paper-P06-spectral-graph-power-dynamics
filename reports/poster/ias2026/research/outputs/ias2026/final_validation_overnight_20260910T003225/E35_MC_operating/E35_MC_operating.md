# E35 — fresh held-out operating-point Monte Carlo: **G3 MIXED**, with one real negative

Seed **20260913**, never used in v1, v2A, v2B or v2C. 2000 draws, **1004
accepted**, 996 rejected — **every rejection for non-dispatchability**, counted as
neither stable nor unstable. **Zero tracking failures.** Every accepted sample
carries the whole 16-subset lattice, so the order-4 statement is tested per
operating point rather than once at the nominal point.

Sampling: global P load 0.85–1.15, global Q load 0.90–1.10, per-bus load scatter
±5 %, dispatch redistribution ±10 %, PV availability 0.70–1.00. The high
rejection rate is the dispatchable envelope doing its job: the fleet cannot serve
load much above 1.05 at any availability, which v2B and v2C already established.

## Preregistered probabilities, exact Clopper-Pearson intervals

| | estimate | 95 % interval |
|---|---|---|
| **P1** every proper subset stable | **0.947** | [0.932, 0.960] |
| **P2** full portfolio unstable | **0.239** | [0.213, 0.267] |
| **P3** order ≤ 3 stable while exact unstable | **0.208** | [0.183, 0.235] |
| **P4** inter-area family retained | **1.000** | [0.996, 1.000] |
| **P5** closure smaller near the boundary | Mann-Whitney **p = 7.3e-38** | median `m₄` 0.0214 near vs 0.1873 far |
| **P6** RC converter repair succeeds | **0.488** | [0.423, 0.553] |
| **P7** synchronous condenser succeeds | **1.000** | [0.985, 1.000] |

Conditional on the portfolio actually being unstable:

| | estimate | 95 % interval |
|---|---|---|
| **genuine order-4 crossing given unstable** | **0.871** | [0.822, 0.911] |
| every proper subset stable given unstable | 0.779 | [0.721, 0.830] |

## What this says

**The mechanism is real but conditional.** A genuine irreducible fourth-order
crossing occurs at about **21 %** of dispatchable operating points. That is not a
majority, and the unconditional P3 criterion is not met, which is why the gate
reads MIXED rather than PASS. But **when this portfolio fails, it fails as a
fourth-order effect 87 % of the time** — the lower-order reconstruction stays
stable while the exact portfolio does not. Both numbers must be quoted together;
either alone misleads.

**The family observable is completely stable.** 1004 of 1004 samples retained the
inter-area family, size 2 in every one, minimum member overlap 0.841. No sample
needed a tracking-failure rejection. That is the v2C correction working across a
much wider envelope than v2C itself sampled.

**The closure margin separates the boundary decisively.** Median `m₄` is 0.0214
within 0.05 of the boundary against 0.1873 away from it, `p = 7.3e-38`.

## The negative result: the frozen converter retune is not robust

RC was optimized once, at the nominal operating point, in E21. Across the
held-out envelope it restores stability in **117 of 240** unstable samples. Its
authority is bounded, and bounded in an interpretable way:

| severity of the instability | RC succeeds |
|---|---|
| `α_IA` ≤ 0.122 (mild) | **80 of 80, 100 %** |
| 0.122 < `α_IA` ≤ 0.242 | 37 of 80, 46 % |
| `α_IA` > 0.242 (severe) | **0 of 80, 0 %** |

Failed repairs leave `α` between `+0.0004` and `+0.2558`. Success also falls with
load: 80 % in the lowest load tercile, 14 % in the highest.

The **25 % synchronous condenser succeeds in all 240**, leaving `α_IA` between
`−0.340` and `−0.221` — a large and uniform margin.

**What may and may not be claimed.** "A converter retune repairs stability while
retaining all PV" is true at the discovery point and for mild instabilities, and
is **not** a robust claim across the operating envelope. What is robust is that
*some* repair always worked: either repair succeeded in 240 of 240, because the
condenser always did.

**What was not tested:** an *adaptive* retune, re-optimized per operating point.
E34 shows a retune exists that meets every declared margin at the nominal point;
whether one exists at every operating point is a different question and is
recorded in `TPWRS_REMAINING_GAPS.md`.

## Gate status

- **G3 MIXED.** The corrected family order-4 effect does not fail in held-out
  Monte Carlo; it occurs at 21 % of operating points unconditionally and 87 %
  conditional on failure.
- **G4 PASS.** The closure separates boundary from interior on fresh data,
  `p = 7.3e-38`.
- **G5 MIXED.** The *fixed* repair does not work only on discovery operating
  points — it works on 49 % of held-out unstable samples and on all mild ones —
  but it is not robust. The condenser is.

## Files

`E35_MC_operating_2000.csv` (+ `.parquet`), `E35_MC_operating_summary.csv`,
`manifest.json`.
