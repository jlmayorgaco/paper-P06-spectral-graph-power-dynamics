# Phase Braess-NEP Report

Gate/verdict: **PASS_RESONANT_FOUND**
Random seed: `20260608`
Git commit: `46f48b23d2ca168678e843f32d844eae5f4d93e5`

## Method
Each line reinforcement reduces the physical ANDES `r` and `x` fields by `1/(1+eps)` and reruns power flow, full eigensolve, and the Phase-0D Beyn Schur-NEP contour solver.  The finite-difference reference is the NEP/ANDES pole movement, not a scalar fixed-point model.

The first-order Schur-NEP sensitivity is decomposed as:
- A: retained omega-row / delta-column stiffness proxy;
- B: residual retained-network/modal-rotation proxy;
- C: condensed-control self-energy contribution from the Schur term;
with the control-frequency signature reported as the fraction of the denominator coming from `A_ec (sI-A_cc)^(-2) A_ce`.

A resonant Braess classification requires a network-family damping decrease, decreasing distance to a matched `A_cc` pole, dominant negative C contribution, relative pole distance below 10%, and a nearest `A_cc` pole whose eigenvector is PLL-dominated.  Otherwise the reversal is not claimed resonant.

## Summary
- Cases tested: mix60, mix60_soft_gfl, mix60_no_pss
- Line-mode perturbations completed: 460
- Margin reversals: 348
- Network-family mode reversals: 232
- Resonant reversals: 40
- Smooth analytic-vs-NEP median relative error: 0.000807729516004942

## Clearest Case
- Case/line/mode-rank: `mix60_no_pss` `Line_42` (20-34), mode 2
- Classification: **resonant**
- Margin zeta: 0.123566 -> 0.123558 (Delta=-7.985e-06)
- Network zeta: 0.133884 -> 0.133816 (FD derivative=-6.860e-03)
- Distance to matched Acc pole: 3.33227 -> 3.3279 (Delta=-4.370e-03)
- Relative Acc distance: 3.683e-02 -> 3.679e-02; nearest Acc PLL-dominated=True
- Base nearest Acc top states: ae PLL1 G9:0.5613; ae PLL1 G5:0.5097; ae PLL1 G7:0.4108; Qsen_y REGF1 G4:0.2273; Qsen_y REGF1 G6:0.1969; S1_y REECA1 1:0.1706
- A dzeta=-1.091e-09, B dzeta=7.974e-09, C dzeta=-6.874e-03
- Analytic total dzeta=-6.874e-03, FD network dzeta=-6.860e-03, relative error=2.015e-03
- Control denominator fraction=9.594e-01

## Verdict
At least one line reinforcement produced a network-family damping decrease, moved the tracked mode closer to a matched condensed-control pole, and was dominated by the self-energy term.  The resonant Braess mechanism is observed in this benchmark.

Artifacts:
- CSV: `outputs\phase_braess_nep\phase_braess_nep_line_decomposition.csv`
- JSON: `outputs\phase_braess_nep\phase_braess_nep_status.json`
- Figure: `outputs\phase_braess_nep\fig_braess_nep_s_plane.png`