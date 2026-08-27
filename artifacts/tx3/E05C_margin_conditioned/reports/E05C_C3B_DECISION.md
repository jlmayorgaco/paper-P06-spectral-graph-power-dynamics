# E05C C3b decision

Frozen population: 15 coalitions (10 pairs, 5 triples), each with one development-locked family. New validation:
8 seeds, gamma_end 1.05225--1.08312, 100% valid tracking.

| Coalition | Family | Valid | Median rho tau=0 | Median rho tau=1 | Maximum rho | C3b-A | C3b-B |
|---|---|---:|---:|---:|---:|---|---|
| A7-A8 | R012 | 8/8 | 0.000564407 | 0.000601267 | 0.000621027 | REJECTED | REJECTED |
| A3-A6 | R012 | 8/8 | 0.000787077 | 0.000617077 | 0.000954389 | REJECTED | REJECTED |
| A2-A8 | R012 | 8/8 | 0.00186431 | 0.00197909 | 0.00211758 | REJECTED | REJECTED |
| A2-A7 | R012 | 8/8 | 0.00154141 | 0.00160632 | 0.0016372 | REJECTED | REJECTED |
| A3-A8 | R011 | 8/8 | 0.000791472 | 0.000832142 | 0.000931017 | REJECTED | REJECTED |
| A1-A4 | R017 | 8/8 | 0.000226414 | 0.000208533 | 0.000241365 | REJECTED | REJECTED |
| A2-A5 | R008 | 8/8 | -0.000286291 | -0.000294496 | -0.000254755 | REJECTED | REJECTED |
| A5-A8 | R008 | 8/8 | -0.000240235 | -0.000266062 | -0.000219433 | REJECTED | REJECTED |
| A6-A8 | R008 | 8/8 | -4.98529e-05 | -4.84682e-05 | -4.70431e-05 | REJECTED | REJECTED |
| A5-A6 | R017 | 8/8 | -0.00256774 | -0.00239706 | -0.00217231 | REJECTED | REJECTED |
| A2-A7-A8 | R012 | 8/8 | 0.000178194 | 0.000193083 | 0.000208951 | REJECTED | REJECTED |
| A2-A5-A8 | R017 | 8/8 | 2.82437e-05 | 2.91068e-05 | 3.17191e-05 | REJECTED | REJECTED |
| A5-A7-A8 | R008 | 8/8 | -2.54671e-05 | -2.86629e-05 | -2.28701e-05 | REJECTED | REJECTED |
| A2-A5-A7 | R008 | 8/8 | -1.63942e-05 | -1.60651e-05 | -1.39225e-05 | REJECTED | REJECTED |
| A3-A6-A8 | R012 | 8/8 | 7.09634e-05 | 5.07837e-05 | 9.3005e-05 | REJECTED | REJECTED |

1. Candidate count: 15. 2. Every frozen family is listed above. 3. Seed count: 8. 4. Every baseline endpoint is
in `e05c_seed_stress_limits.parquet`. 5. Valid tracks: 5760/5760. 6. Destabilizing real
shifts occur in ['A7-A8', 'A3-A6', 'A2-A8', 'A2-A7', 'A3-A8', 'A1-A4', 'A2-A7-A8', 'A2-A5-A8', 'A3-A6-A8'].
7--8. Nominal/final rho are tabulated. 9. No candidate reaches the frozen amplification gate. 10. C3b-A=REJECTED.
11. Hard failures=0. 12. Screen failures=0. 13. Points with rho_comp>=0.25=0.
14. Points with Delta zeta<=-0.0025=0. 15. C3b-B=REJECTED.
16. Neither pairs nor triples approach 0.10. 17. Conditioning is reported separately and does not explain a
material effect. 18. Strongest result: A2-A8 at C06, tau=1.0,
rho_comp=0.00211758. 19. Strongest negative evidence: zero material damping points and zero
composition failures despite complete valid tracking. 20. E06 consideration authorized: **false**.
