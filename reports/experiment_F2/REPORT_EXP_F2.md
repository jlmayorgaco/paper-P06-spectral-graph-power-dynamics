# Experiment F2 — Actual graph-modal GFL characteristic

**F2_MODAL_STRUCTURE: `NON_MODAL`**

## Exact characteristic structure

The analytic SimpleGFLDC model contributes nine differential states per port. At common parameters its network interconnection is the exact descriptor determinant `det(P(s))=0`, with `P(s) = [I⊗(sI₉−A(Kp,Ki))  −I⊗B; I⊗C  Yport,rect + I⊗D]` under the frozen component-current sign convention. It has 90 differential and 20 algebraic variables for ten generator ports; its characteristic polynomial degree is at most 90 in `s`. Kp and Ki enter the local A matrix affinely through one rank-one update each; the exact gain dependencies are reported in TABLE_F01. The controlled homogeneous limit repeats the bus-33 analytic GFL state and parameters at all ten ports. Fixed injections balance this common state against the frozen passive network and have zero derivative, so they do not alter the characteristic pencil.

A scalar law `F(s;ν,Kp,Ki)=0` requires a common node basis for both the conductance and susceptance parts of the port admittance. In the F0 conductance eigenbasis, the conductance off-diagonal residual is `5.511393448606591e-17` and the reactive residual is `0.2913839107215721`; the normalized commutator is `0.08255205726800578`. The homogeneous descriptor residual and total port-transfer residual are `2.491091916863042e-6` and `0.2807917283445948`. The small raw descriptor ratio is dominated by its diagonal `sI` state terms; the gate is applied to the condensed network/port coupling so that those diagonal terms do not dilute the coupling measure.

## Heterogeneous IEEE-39 analytic model

The detailed analytic endpoint uses every generator port, the frozen network plus archived impedance-load linearizations, independent local nine-state SimpleGFLDC Jacobians, and the same F0 basis. At the frozen ExpC critical pole `eig94` (`s=-0.35790147141068607 + 3.9211749225740378im`), its port-transfer and descriptor off-diagonal residuals are `0.2836827293809271` and `0.061415305507763626`. Exact 2-channel Schur self-energy Gamma by graph mode is stored in TABLE_F03; this diagnoses coupling, it does not create an independent modal characteristic.

This assembly imports no PowerDynamics module. Its saved matrices and all mode blocks are reproducible from the archived analytical inputs.

## PLL-only limit

With `e≈−Vθ`, the implemented PLL equations retain the first-order frequency state and give `τs³+s²+VKp s+VKi=0`. Only in the ideal frequency-loop limit `τ→0` does this reduce to `s²+VKp s+VKi=0`. No graph eigenvalue appears in that isolated PLL limit; any graph dependence has to come from the coupled network equations.

## STOP gate

**NON_MODAL.** The homogeneous detailed model exceeds the preregistered 5% modal-coupling tolerance, so no exact/approximately separated `F(s;ν,Kp,Ki)` is available. F3/F4 spectral gain-law work is stopped. Continue with the exact heterogeneous self-energy/pathway formulation; do not claim `K=h(Lc)`.
