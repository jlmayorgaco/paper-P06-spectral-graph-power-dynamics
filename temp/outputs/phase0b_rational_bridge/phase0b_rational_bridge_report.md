# Phase-0B Controller-Aware Rational Bridge

Gate status: **BLOCKED**

This audit uses the Schur complement `A_eff(s)=A_ee + A_ec (sI-A_cc)^(-1) A_ce`.
Success is non-circular: poles are obtained by fixed-point/root refinement from independent initial poles, not by evaluating the bridge at the known ANDES pole.

## Cases
### no_pss
- Accepted: True
- States: x_e=20, x_c=110
- Converged bridge poles: 9 from 54 independent starts
- Critical ANDES pole: -1.34635 + j8.61104, f=1.37 Hz, zeta=0.1545
- Critical bridge pole: -1.34635 + j8.61104, f=1.37 Hz, zeta=0.1545
- Errors: frequency=8.1e-14%, damping=6.83e-13%, sign=True
- Scalar D norm: 0
- Controller-aware D_eff real Frobenius norm: 4.053
- Controller-aware D_eff imaginary Frobenius norm: 2.178
- Partition CSV: `outputs\phase0b_rational_bridge\no_pss_partition.csv`
- Validation CSV: `outputs\phase0b_rational_bridge\no_pss_mode_validation.csv`

### base
- Accepted: True
- States: x_e=20, x_c=140
- Converged bridge poles: 9 from 54 independent starts
- Critical ANDES pole: -1.34601 + j8.61148, f=1.371 Hz, zeta=0.1544
- Critical bridge pole: -1.34601 + j8.61148, f=1.371 Hz, zeta=0.1544
- Errors: frequency=1.3e-13%, damping=1.26e-12%, sign=True
- Scalar D norm: 0
- Controller-aware D_eff real Frobenius norm: 4.053
- Controller-aware D_eff imaginary Frobenius norm: 2.178
- Partition CSV: `outputs\phase0b_rational_bridge\base_partition.csv`
- Validation CSV: `outputs\phase0b_rational_bridge\base_mode_validation.csv`

### mix60
- Accepted: False
- States: x_e=11, x_c=167
- Converged bridge poles: 5 from 40 independent starts
- Critical ANDES pole: -10.7514 + j86.3419, f=13.74 Hz, zeta=0.1236
- Critical bridge pole: -1.52556 + j9.35844, f=1.489 Hz, zeta=0.1609
- Errors: frequency=89.2%, damping=30.2%, sign=True
- Scalar D norm: 0
- Controller-aware D_eff real Frobenius norm: 2.216
- Controller-aware D_eff imaginary Frobenius norm: 0.82
- Partition CSV: `outputs\phase0b_rational_bridge\mix60_partition.csv`
- Validation CSV: `outputs\phase0b_rational_bridge\mix60_mode_validation.csv`

## Interpretation
The controller-aware Schur bridge did not pass the non-circular gate for all required cases. Do not advance to Phase 2-4 as full-ANDES evidence; inspect the validation tables and participation factors to identify which modes are not captured by the selected electromechanical partition.
