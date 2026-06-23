# Phase-0D Rational NEP Contour Solver

Gate status: **PASS**

This audit solves the exact first-order Schur nonlinear eigenproblem
`S_e(s)=sI-A_ee-A_ec(sI-A_cc)^(-1)A_ce` with Beyn contour integration.
No fixed-point iteration is used for acceptance.  The full ANDES eigensolve
is the reference.

## Non-Circularity Audit
- Beyn zeros are first extracted from the contour moments `A0` and `A1` and the reduced eigenproblem.
- The local polish starts only from each Beyn Ritz value and solves for the smallest eigenvalue of `S_e(s)` to vanish.
- The polish function receives no ANDES pole, mode index, frequency, or damping target. ANDES eigenvalues are used only after extraction for validation and count diagnostics.
- Therefore a zero can be produced from the Jacobian/Schur NEP even if the corresponding ANDES pole is not supplied to the solver.

## Contours
- `network_margin`: Re=[-4.0, 0.8], Im=[3.5, 12.0] rad/s (0.557-1.91 Hz), 512 quadrature nodes.
- `pll_control`: Re=[-25.0, 0.8], Im=[60.0, 95.0] rad/s (9.55-15.1 Hz), 512 quadrature nodes.

## Required Gate Checks
### no_pss
- Accepted: **True**
- Argument-count checks pass: **True**
- States: retained=20, condensed=110
- Beyn zeros accepted: 9
- Sigma/D_eff at s=0: scalar D norm=0, real Fro=9.968, imag Fro=0, solve=solve, cond=8.94e+04.
- Sigma/D_eff at matched critical pole: real Fro=4.053, imag Fro=2.178. Phase-0B critical-pole reference: real Fro=4.053, imag Fro=2.178.
- ANDES critical pole: -1.34635 + j8.61104, f=1.3705 Hz, zeta=0.15447.
- NEP-Beyn matched pole: -1.34635 + j8.61104, f=1.3705 Hz, zeta=0.15447, family=I_network, pi_c=0.680.
- Critical errors: frequency=4.861e-14%, damping=6.289e-13%, relative residual=5.467e-17.
- Nearest Acc pole to matched NEP zero: -0.44778 + j0.584957, distance=8.08.
- Top condensed-state reconstruction: LAG_y TGOV1N TGOV1 1:3.536; LAG_y TGOV1N TGOV1 8:2.216; LAG_y TGOV1N TGOV1 9:1.144; LAG_y TGOV1N TGOV1 3:1.122; LA_y IEEEX1 9:1.025; LA_y IEEEX1 7:0.7354
- Contour counts:
  - network_margin: extracted=9, rank=14, argument-count=9.000+j3.51e-06, Acc poles inside=0, argument+Acc=9.000, ANDES poles inside=9, count error=6.42e-05, pass=True.
  - pll_control: extracted=0, rank=20, argument-count=-0.000+j1.56e-06, Acc poles inside=0, argument+Acc=-0.000, ANDES poles inside=0, count error=2.19e-06, pass=True.
- Zeros table: `outputs\phase0d_nep_contour_solver\no_pss_zeros.csv`
- Contour diagnostics: `outputs\phase0d_nep_contour_solver\no_pss_contour_counts.csv`

### base
- Accepted: **True**
- Argument-count checks pass: **True**
- States: retained=20, condensed=140
- Beyn zeros accepted: 9
- Sigma/D_eff at s=0: scalar D norm=0, real Fro=9.968, imag Fro=0, solve=pinv, cond=7.71e+18.
- Sigma/D_eff at matched critical pole: real Fro=4.053, imag Fro=2.178. Phase-0B critical-pole reference: real Fro=4.053, imag Fro=2.178.
- ANDES critical pole: -1.34601 + j8.61148, f=1.3706 Hz, zeta=0.15443.
- NEP-Beyn matched pole: -1.34601 + j8.61148, f=1.3706 Hz, zeta=0.15443, family=I_network, pi_c=0.653.
- Critical errors: frequency=1.62e-13%, damping=1.06e-12%, relative residual=2.981e-17.
- Nearest Acc pole to matched NEP zero: -0.495739 + j0.71038, distance=7.95.
- Top condensed-state reconstruction: LAG_y TGOV1N TGOV1 1:3.536; LAG_y TGOV1N TGOV1 8:2.215; LAG_y TGOV1N TGOV1 9:1.145; LAG_y TGOV1N TGOV1 3:1.123; LA_y IEEEX1 9:0.9703; LA_y IEEEX1 7:0.7134
- Contour counts:
  - network_margin: extracted=9, rank=14, argument-count=9.000+j4.14e-06, Acc poles inside=0, argument+Acc=9.000, ANDES poles inside=9, count error=6.54e-05, pass=True.
  - pll_control: extracted=0, rank=20, argument-count=-0.000+j1.56e-06, Acc poles inside=0, argument+Acc=-0.000, ANDES poles inside=0, count error=2.19e-06, pass=True.
- Zeros table: `outputs\phase0d_nep_contour_solver\base_zeros.csv`
- Contour diagnostics: `outputs\phase0d_nep_contour_solver\base_contour_counts.csv`

### mix60
- Accepted: **True**
- Argument-count checks pass: **True**
- States: retained=11, condensed=167
- Beyn zeros accepted: 7
- Sigma/D_eff at s=0: scalar D norm=0, real Fro=4.94, imag Fro=0, solve=pinv, cond=3.85e+20.
- Sigma/D_eff at matched critical pole: real Fro=0.05444, imag Fro=0.5387. Phase-0B critical-pole reference: real Fro=2.216, imag Fro=0.82.
- ANDES critical pole: -10.7514 + j86.3419, f=13.742 Hz, zeta=0.12357.
- NEP-Beyn matched pole: -10.7514 + j86.3419, f=13.742 Hz, zeta=0.12357, family=II_control, pi_c=0.853.
- Critical errors: frequency=8.475e-11%, damping=2.479e-09%, relative residual=8.606e-11.
- Nearest Acc pole to matched NEP zero: -10.5593 + j86.246, distance=0.215.
- Top condensed-state reconstruction: ae PLL1 G7:14.18; ae PLL1 G5:13.83; Qsen_y REGF1 G6:7.926; Qsig_y REGF1 G6:3.478; am PLL1 G7:3.266; am PLL1 G5:3.185
- Contour counts:
  - network_margin: extracted=4, rank=9, argument-count=4.000+j-6.17e-06, Acc poles inside=0, argument+Acc=4.000, ANDES poles inside=4, count error=2.87e-05, pass=True.
  - pll_control: extracted=3, rank=11, argument-count=-0.000+j1.14e-06, Acc poles inside=3, argument+Acc=3.000, ANDES poles inside=3, count error=1.89e-06, pass=True.
- Zeros table: `outputs\phase0d_nep_contour_solver\mix60_zeros.csv`
- Contour diagnostics: `outputs\phase0d_nep_contour_solver\mix60_contour_counts.csv`

## Interpretation
Phase-0D passed for all three cases. The rational Schur NEP contour solver can be used as the bridge for later validation phases.
