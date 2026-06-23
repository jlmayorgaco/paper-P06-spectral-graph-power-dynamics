# Gate 1: Coordination Value of Planning Levers

This gate tests whether coordinated packages of topology, damping, and inertia
outperform the best single lever under equal normalized cost and DC flow
feasibility.  The ground truth is the complete eigensolve of the tested
second-order dynamic model.

- Commit: `46f48b23d2ca168678e843f32d844eae5f4d93e5`
- Random seed: `20260609`
- Costs: `c_w=1.0`, `c_D=2.0`, `c_M=8.0`
- Equal-cost budget: `1.0`
- Combined-package search: top-6 screened candidates per lever
- Cases passing base filters: `62` (60 synthetic, 2 IEEE-topology)
- Median topology-action feasibility ratio: `1.000`

## Aggregate Table

| metric | pair | eps | n | median SynNorm | % SynNorm < -0.10 | median CV | % CV < 0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| S_zeta | w+D | 0.05 | 62 | 0.007668 | 0.0% | 0.0173 | 31.1% |
| S_zeta | w+D | 0.10 | 62 | 0.02207 | 0.0% | 0.0173 | 31.1% |
| S_zeta | w+D | 0.20 | 62 | 0.04991 | 4.8% | 0.0173 | 31.1% |
| S_zeta | w+M | 0.05 | 62 | -0.002966 | 11.3% | -0.001815 | 50.8% |
| S_zeta | w+M | 0.10 | 62 | -0.00512 | 22.6% | -0.001815 | 50.8% |
| S_zeta | w+M | 0.20 | 62 | 0.02065 | 17.7% | -0.001815 | 50.8% |
| S_zeta | D+M | 0.05 | 62 | 0.004641 | 0.0% | 0.03641 | 21.0% |
| S_zeta | D+M | 0.10 | 62 | 0.01641 | 6.5% | 0.03641 | 21.0% |
| S_zeta | D+M | 0.20 | 62 | 0.03434 | 4.8% | 0.03641 | 21.0% |
| S_g | w+D | 0.05 | 62 | 0.0007638 | 0.0% | -0.05294 | 72.1% |
| S_g | w+D | 0.10 | 62 | 0.04141 | 0.0% | -0.05294 | 72.1% |
| S_g | w+D | 0.20 | 62 | 0.1427 | 0.0% | -0.05294 | 72.1% |
| S_g | w+M | 0.05 | 62 | 0.0003818 | 12.9% | -0.07363 | 67.2% |
| S_g | w+M | 0.10 | 62 | 0.2555 | 11.3% | -0.07363 | 67.2% |
| S_g | w+M | 0.20 | 62 | 0.5423 | 11.3% | -0.07363 | 67.2% |
| S_g | D+M | 0.05 | 62 | 0.02298 | 1.6% | -0.006759 | 51.6% |
| S_g | D+M | 0.10 | 62 | 0.08476 | 0.0% | -0.006759 | 51.6% |
| S_g | D+M | 0.20 | 62 | 0.1519 | 1.6% | -0.006759 | 51.6% |

## Metric Agreement

- Pair-case comparisons: `184`
- Same CV sign for modal and resolvent margins: `50.0%`
- Same single-lever ranking for modal and resolvent margins: `28.3%`

## Verdict

- Status: **COORDINATION_VALUE_OBSERVED**
- Implication: Equal-cost mixed packages beat the best single lever in a majority of cases for at least one metric/pair.
- Rows with CV majority support: `4`
- Rows with nonlinear-synergy majority support: `0`

## Scope

- This is not a full ANDES IBR validation. IEEE cases are used as topology/flow templates.
- The resolvent margin uses identity protection scalings because no protection matrices were provided.
- The equal-cost package search is favorable but screened: it evaluates combinations among the top single-lever candidates.