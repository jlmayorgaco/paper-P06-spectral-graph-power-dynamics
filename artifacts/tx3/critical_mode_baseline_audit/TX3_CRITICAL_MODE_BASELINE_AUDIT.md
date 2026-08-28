# TX3 Critical-Mode Baseline Audit

## Scope

This is a descriptive post-mortem at `kappa_grid=1`. It uses the eight frozen E05D seeds, all 15 frozen
E05B dynamic coalitions, and the original A1--A8 endpoints. It creates no claim gate and cannot change
`C3a=REJECTED`, `C3b-A=REJECTED`, `C3b-B=REJECTED`, or `C3c=UNRESOLVED`.

## Actual margin-setting mode

All eight seeds select the same physical family under the frozen rule: {'synchronous_electromechanical': 8}. Frequencies
span 1.21973--1.22906 Hz and damping spans
0.041654--0.045275. Synchronous-machine plus governor/
excitation participation is 0.9862--0.9898;
PLL participation is 0.0025--0.0031; converter-control
participation is 0.0076--0.0107.
The descriptor-weighted endogeneity audit gives the same family-level interpretation. Classification uses a
fixed 0.60 dominance threshold on normalized absolute left/right participation; all state-level values and the
left/right vectors are released.

## Connected effects on that mode

Every one of the 200 unique seed/subset tracks and 120 candidate/seed connected cases is
retained. Tracking success is 100.000%; the minimum successful BMAC is
0.999969. The global rho_comp range is
-9.63916e-06--0.000139286; median absolute rho_comp is
6.38782e-06. The largest absolute damping externality is
5.14674e-06. Hard/screen composition failures are
0/0. These values are descriptive and do
not retroactively activate any historical gate.

| Coalition | Order | Valid | Median Re Delta lambda | Median Delta zeta | Median rho_comp | Min rho | Max rho |
|---|---:|---:|---:|---:|---:|---:|---:|
| A7-A8 | 2 | 8/8 | 3.9846e-05 | -4.2179e-06 | 1.1699e-04 | 6.6925e-05 | 1.3929e-04 |
| A3-A6 | 2 | 8/8 | 1.3796e-05 | -1.7363e-06 | 4.0921e-05 | 3.1238e-05 | 4.9539e-05 |
| A2-A8 | 2 | 8/8 | 3.9989e-06 | -5.0509e-07 | 1.1685e-05 | 5.9138e-06 | 1.4665e-05 |
| A2-A7 | 2 | 8/8 | 9.6226e-06 | -1.2218e-06 | 2.8241e-05 | 1.8286e-05 | 4.1976e-05 |
| A3-A8 | 2 | 8/8 | 4.9970e-06 | -6.1182e-07 | 1.4952e-05 | 1.2006e-05 | 1.7577e-05 |
| A1-A4 | 2 | 8/8 | 1.6738e-05 | -2.3508e-06 | 5.0240e-05 | 4.3315e-05 | 6.5633e-05 |
| A2-A5 | 2 | 8/8 | 1.4122e-05 | -1.8667e-06 | 4.2932e-05 | 3.9374e-05 | 4.9120e-05 |
| A5-A8 | 2 | 8/8 | -2.2747e-06 | 3.0841e-07 | -6.7495e-06 | -9.6392e-06 | -5.2587e-06 |
| A6-A8 | 2 | 8/8 | -9.5244e-07 | 1.2872e-07 | -2.8777e-06 | -4.8147e-06 | -2.3745e-06 |
| A5-A6 | 2 | 8/8 | -1.4437e-07 | 1.8971e-08 | -4.2617e-07 | -5.4835e-07 | -2.7234e-07 |
| A2-A7-A8 | 3 | 8/8 | 3.1655e-07 | -4.0219e-08 | 9.2716e-07 | 2.3847e-07 | 1.8434e-06 |
| A2-A5-A8 | 3 | 8/8 | 1.0647e-06 | -1.4085e-07 | 3.2219e-06 | 2.9972e-06 | 3.8138e-06 |
| A5-A7-A8 | 3 | 8/8 | -2.1555e-07 | 3.0790e-08 | -6.5123e-07 | -8.7577e-07 | -5.0173e-07 |
| A2-A5-A7 | 3 | 8/8 | 1.3593e-06 | -1.7964e-07 | 4.1078e-06 | 3.6080e-06 | 4.7476e-06 |
| A3-A6-A8 | 3 | 8/8 | 3.3338e-07 | -4.3716e-08 | 9.9750e-07 | 7.2349e-07 | 1.2007e-06 |

## Why the externality and margin stories differ

The E05C development locks, which selected the largest real pole externality from the E05B population, have
damping 0.4939--0.9917
and are dominated by PLL/current-command/electrical-control states. The true margin-setting family has roughly
0.0436 damping and is synchronous-electromechanical. Thus the poles carrying the
demonstrated controller externalities are not the poles setting system damping margin in this benchmark.

The audit supports an explanation, not a new hypothesis test: finite controller interactions can be real and
reproducible without materially consuming the decay margin of a different, synchronous-dominated critical mode.
No cycle/SCC causality, mechanism surgery, E06, or new ParaEMT claim is authorized.
