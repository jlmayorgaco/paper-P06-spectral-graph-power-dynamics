# Q1 — PLL delay contribution has rank at most one

## Source equation and delayed channel

The frozen repository PLL is the `ComposableInverter.PLL_LPF` used by `PhysicalGFLDC` in `experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl`. In the repository's state notation, its two PI-channel states obey

\[
\dot\theta_i=\Delta\omega_i,\qquad
\dot{\Delta\omega}_i=\frac{\Delta\omega_{I,i}+K_{p,i}e_i-\Delta\omega_i}{\tau_{\rm LPF}},\qquad
\dot{\Delta\omega}_{I,i}=K_{i,i}e_i,
\]

where the actual scalar phase-detector error is

\[
e_i(x)=-\sin\theta_i\,u_{r,i}(x)+\cos\theta_i\,u_{i,i}(x).
\]

The experiment delays only this measured scalar `e_i`; it does not delay current control, voltage injection, or setpoints. Linearizing the delayed scalar about the frozen equilibrium gives one row functional \(c_i^T\delta x=e_{i,x}\delta x\). The same delayed scalar multiplies both PI branches. Its state-space injection is

\[
b_i=K_{p,i}b_{p,i}+K_{i,i}b_{I,i},\quad
b_{p,i}=\tau_{\rm LPF}^{-1}e_{\Delta\omega_i},\quad
b_{I,i}=e_{\Delta\omega_{I,i}}.
\]

Therefore the delayed Jacobian contribution is exactly

\[
A_{\tau_i}=b_i c_i^T,
\]

and \(\operatorname{rank} A_{\tau_i}\le1\). This conclusion depends on one shared scalar measurement path. It would not hold if the proportional and integral branches used distinct delayed measurements or if another delayed PLL measurement channel were present.

After the fixed algebraic network closure and removal of the rotational gauge coordinate, \(c_i\) is the derivative of the repository's actual PLL error through the voltage solution. The numeric audit constructs that derivative with automatic differentiation and checks the outer-product reconstruction and singular spectrum device by device.

**Claim class:** source-equation identity, conditional on the stated shared-error delay contract; numerical residual is reported separately in `TABLE_Q01_PLL_RANK_RESIDUAL.csv`.
