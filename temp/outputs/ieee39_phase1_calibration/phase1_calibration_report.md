# Phase-1 Controller Calibration Audit

This audit sweeps documented controller profiles for the derived IEEE39 IBR cases. The reference is the full ANDES eigensolve for each modified case.

- Profiles tested: 10
- Case/profile runs: 50
- Stable case/profile runs: 9

## Best profile by case
- ieee39_ibr_gfl20: controller-only very_soft_all (pos=1, max_real=1.011); diagnostic no_pss_soft_all (pos=0, max_real=4e-10); critical top=e1d GENROU 10
- ieee39_ibr_gfm20: controller-only baseline (pos=6, max_real=4.291); diagnostic no_exciter_pss (pos=2, max_real=4.286); critical top=omega GENROU 2
- ieee39_ibr_mix20: controller-only slow_pll (pos=3, max_real=1.066); diagnostic no_exciter_pss (pos=0, max_real=8.298e-08); critical top=WO_x BusFreq 4
- ieee39_ibr_mix40: controller-only slow_pll (pos=2, max_real=1.065); diagnostic no_exciter_pss (pos=0, max_real=8.867e-08); critical top=PI_xi PLL1 G7
- ieee39_ibr_mix60: controller-only soft_gfl (pos=0, max_real=4e-10); diagnostic no_pss (pos=0, max_real=2e-09); critical top=omega GENROU 2

## Interpretation
- A profile that makes a case stable is a calibration candidate, not final tuning.
- High-frequency critical modes are classified by participation, not by frequency alone.
- Estimator validation remains blocked until a calibrated full-ANDES IBR case is matched by a reduced `(L,M,D)`/control-coupling surrogate.
