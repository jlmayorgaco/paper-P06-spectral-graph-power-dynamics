# Poster claim ledger — Experiment Q

Claims below are limited to the computations completed in this experiment.

- **Q0 — zero-delay feasibility witness (full nonlinear numerical result):** Uniform `rho=0.875` with nominal gains passed the full 203-pole ODE check and all five frozen PowerDynamics events, giving 4,727.416 MW GFL (87.5%) and 675.345 MW retained SG. This is a feasible point, not a maximum.
- **Q1 — PLL delayed-channel rank (source-level identity):** Under the repository PLL equations, proportional and integral directions share the same scalar delayed phase error, yielding a rank-at-most-one delayed contribution per GFL.
- **Q2 — low-rank DDE closure (exact algebraic identity plus numerical check):** The determinant-lemma closure reproduced tested operators with maximum residual `3.13e-18`; log-determinant magnitude and phase discrepancies were at most `7.18e-11` and `8.01e-11` rad on the evaluated points.
- **Q4 — conditional PI target equation (exact conditional identity):** The prescribed target-root equation was satisfied across 45 checked cases; only 5/45 gain pairs met both frozen gain bounds. A nearest/rightmost-pole boundary is not established.
- **Q5 — SG terminal action (reduced descriptor theorem and numerical check):** For the fixed interior descriptor, changing a retained share produces an affine rank-at-most-two KCL-row update. None of 120 tested conditional quadratic roots was real and within the physical retention interval.
- **Q6 — joint action space (reduced descriptor theorem and numerical check):** The fixed-coordinate update dimension is at most 30 (10 gain directions and 20 KCL-row directions); the tested factorization residual was at most `9.07e-15`.

**Do not claim:** a positive-delay SAFE/UNSAFE classification; a delay-aware feasible optimum; maximum replacement; unique or optimal gains; a validated delay frontier; current-limiter/DC-energy safety; or heterogeneous-delay spatial effects. The 20/40 ms winding values are unverified diagnostics, not root counts.
