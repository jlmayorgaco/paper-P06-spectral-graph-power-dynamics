# Final summary — ExpH

**Status: PARTIALLY_SUPPORTED.** The local Jordan/Feshbach theory is supported by the analytical model, but the complete robust co-design is not validated independently.

- Gauge quotient: full chain length 2; one simple physical zero remains after quotient.
- Root law: \(\lambda_{phys}=-A_i(K)\epsilon_i+O(\epsilon_i^2)\); measured exponent 1.01088.
- Gain dependence: complex-step derivative median error (1.28\times10^{-6}), p95 (1.78\times10^{-6}); mixed-authority rank 1.
- Structure: exact one-scalar PLL Woodbury form; graph basis captures 99.24%; direct and self-energy terms exhibit strong cancellation.
- Frozen analytical candidate: bus 36, \(\rho=0.98\), 11.2 MW retained, 5391.561 MW converted; \(\alpha=-0.064796\,s^{-1}\); direct robust \(\beta=10^{-4}\), margin 0.22515.
- Search limits: fixed single-support epsilon grid only. No joint KKT, LICQ/SOSC, local, or global certificate. The 100 MW frequency excursion is active and not fed back into design.
- Independent PowerDynamics: equilibrium residual (1.40\times10^{-12}), 102 finite poles, 3 gauge poles, \(\alpha_{PD}=-0.001326\,s^{-1}\). This fails the required \(-0.05\,s^{-1}\) margin and differs from analysis by 0.06347 s⁻¹.
- Nonlinear TDS: two bus-16 pulse amplitudes pass response scaling; this does not resolve the spectral mismatch.
- ExpG analytical comparison: +3810.571 MW converted at stronger requested \(\beta\); the comparison does not transfer to the PD validation because its margin fails.
- ExpE: 1.189211882 MW retained is only the user's nominal reference; no final candidate artifact was available.

Frozen candidate SHA-256: `2967b772d0330e71e528298940e02812b0642c0a71610f2a72bd1c16eedaeac3`. No push.
