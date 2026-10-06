# Claim ledger — Experiment D

| ID | Claim | Status | Evidence / scope |
|---|---|---|---|
| D-C01 | Local PLL gain dependence is affine and low rank. | EXACT | Kp and Ki each rank one; joint rank two in the independently coded nine-state model. |
| D-C02 | Local rated SG/GFL ports form an affine rho mixture at the frozen equilibrium. | LOCAL_EXACT | Device-port identity and cross-bus local-port audit; no current-limiter claim. |
| D-C03 | A full-system determinant-lemma closure was used for the optimizer. | NOT_USED | The attained physical upper bound proves the objective without pole-branch elimination; no closure is claimed. |
| D-C04 | Bus-33 replacement maximum is globally certified on the preregistered domain. | PASS | rho is bounded above by 1 and the all-mode-feasible rho=1, beta=1 witness attains that bound. |
| D-C05 | The bus-33 optimum matches detailed PowerDynamics. | PASS | Exact stored rightmost-pole match after the candidate SHA gate. |
| D-C06 | The optimum repeats at buses 30, 35, and 37. | PASS | All full-replacement endpoints pass the all-mode margin and fresh PowerDynamics validation. |
| D-C07 | The post-freeze surface confirms feasibility over the whole design box. | PARTIAL | 132/132 cells evaluated; 123 meet the margin, with nine bus-37 interior violations. The optimum endpoint is confirmed. |
| D-C08 | Every rho boundary branch or stationary point was eliminated algebraically. | NOT_CLAIMED | This is unnecessary for the saturated global maximum and was not computed. |
| D-C09 | The result supports nonlinear current-limited deployment. | NOT_ESTABLISHED | The current-sharing model has no current limiter; only equilibrium and small-signal spectra are covered. |
