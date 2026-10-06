# ExpQ2 baseline regression gate

**Status: PASS_BASELINE.** This gate reuses the just-reproduced ExpQ Q0 PowerDynamics build/TDS results and the frozen all-SG PowerDynamics TDS reference, then independently checks their expected values and the immutable ExpN/ExpP candidate hashes. It makes no design-model or candidate changes.

- Passed checks: 14/14.
- ExpN: support 38, retained 1.124344477653 MW, rho38 0.998645368099213, alpha analytic/PD -0.050000001105/-0.050000000139 s⁻¹.
- P/Q max error 1.481e-15 pu; ExpN 100 MW bus-16 peak/RoCoF 39.333401 Hz / 1.072208 Hz/s; pointwise beta witness 3.7465e-11 < 1.6991e-06.
- All-SG PD TDS peak/RoCoF 0.075365 Hz / 0.123835 Hz/s, alpha -0.098065409 s⁻¹.
- Evidence: `TABLE_Q2_BASELINE_REGRESSION.csv`; no optimizer has run.
