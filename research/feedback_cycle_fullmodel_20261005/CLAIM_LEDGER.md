# CLAIM LEDGER

(label | claim | evidence file) — filled as results arrive.

| label | claim | evidence |
|---|---|---|
| POWERDYNAMICS_VALIDATED | Julia Jacobians equal the exact-delay Python reader (rel. 1e-16, rho 0.125-0.875); baseline roots, counts (0 beyond margin at 40 ms; 8 unstable/10 beyond at 44 ms) reproduced | derived/TABLE02_baseline_parity.csv |
| NONLINEAR_TDS_VALIDATED | Baseline 40 ms passes the five frozen events (max dF 0.471 Hz) | derived/TABLE12_nonlinear_events.csv |
| NUMERICALLY_VERIFIED | All-SG has no mode above 1.5 Hz; with GFL a 4.5-5.0 Hz PLL family (participation 0.90-0.99) exists for every rho > 0 | TABLE03, TABLE04a, TABLE_C01 |
| NUMERICALLY_VERIFIED | Margin crossings of the PLL family at 41.57, 41.68, 42.27, 43.26 ms; counts 0/4/6/8 unstable roots at 41/42/43/44 ms (floating-point, not certified) | TABLE_D01b, TABLE_D01c |
| EXACT_IDENTITY | det Delta = det(sI-A0) prod L_i det(I+Q), 24/24 random points, max error 1.1e-12 | TABLE_E00 |
| NUMERICALLY_VERIFIED | At critical roots eig(Q) = -1 with |L_i| in [0.011, 0.28]; closest pair cycles 35-36, 33-34, 36-37 (|1-p| 0.005-0.018) | TABLE_F02 |
| EXPLORATORY | SCCs at 0.1 max|Q| are 7-10 buses (not localised); frozen pairs = top-3 physical-coupling pairs | TABLE_F03, TABLE_G02 |
| REFUTED | P4: Kp authority per log-gain is 1.7-2.0x Ki at critical roots | TABLE08 |
| REFUTED | S0 = {30,33,36,37} (Ki/Kp/joint) and all four-site supports restore the margin at 44 ms | TABLE09 |
| NEGATIVE_RESULT | Feedback-core order needs 10 sites; sensitivity order 9; physical order 6; random 8-10 | TABLE09c |
| NUMERICALLY_VERIFIED | Uniform Kp -14.4% (or all-ten nodal) restores the margin; 5/5 events | TABLE09, TABLE12 |
| EXPLORATORY | Forced response: resonant ratio 2.2 (40 ms), 12 (sparse repaired); 44 ms case unstable | TABLE11 |
| NEGATIVE_RESULT | No RHP zeros in the PLL return channels; load-to-angle RHP zeros already in all-SG; not independently confirmed | TABLE_L01 |
| UNRESOLVED | Whether a successful support smaller than 6 sites exists (search was rule-ordered, not exhaustive) | TABLE09c |
