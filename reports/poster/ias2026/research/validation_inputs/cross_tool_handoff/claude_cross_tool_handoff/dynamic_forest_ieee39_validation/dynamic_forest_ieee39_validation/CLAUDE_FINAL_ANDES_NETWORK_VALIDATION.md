# FINAL independent ANDES validation — network diagnostic layer only

Purpose: independently validate the physical-line diagnostic without pretending that ANDES implements the paper's GFL model.

Use the already reconciled equation-equivalent IEEE-39 model available in both implementations: same live network/Ybus/shunts/taps and dispatch, same electromechanical machine equations, and the equation-equivalent first-order SEXS/AVR configuration from F1. Keep PSS off unless the PSS implementations are equation-equivalent. Use a static-injection replacement/removal case if the exact GFL is unavailable. Do not tune parameters to preserve P4.

Preregister these 12 branches selected by RNG seed 20260911 before the holdout result was evaluated:

    [3, 9, 15, 17, 22, 26, 31, 32, 35, 39, 42, 44]

Define gamma_e as a multiplier on the complete branch two-port admittance and use gamma=1 +/- 0.002.

For each branch in the internal implementation:
1. solve PF/equilibrium at gamma+/-;
2. track the same inter-area family;
3. compute full-DAE finite-difference d lambda / d gamma_e;
4. rebuild the reduced port operator at each perturbed equilibrium;
5. compute ds*/d gamma_e = -(p^H Q_gamma q)/(p^H Q_s q).

For exactly the same cases in ANDES:
1. solve PF;
2. run eigenanalysis with the equation-equivalent dynamics;
3. track the same inter-area family;
4. compute d lambda_ANDES / d gamma_e.

Frozen report metrics:
- PF voltage mismatch internal vs ANDES;
- base mode frequency/damping mismatch;
- Spearman correlation of Re(dlambda/dgamma) internal vs ANDES;
- sign agreement;
- median/max relative complex derivative error;
- top-5 stabilizing branch overlap;
- internal port derivative vs internal DAE;
- internal port derivative vs ANDES.

Freeze strong-pass criteria before execution:
- Spearman >= 0.90;
- sign agreement >= 11/12;
- top-5 overlap >= 4/5;
- median relative derivative error <= 10%;
- no strong branch (|d alpha/d gamma| above median) has opposite sign.

Allowed claim if it passes:
'On an equation-equivalent IEEE-39 electromechanical configuration, the reduced network-port derivative independently reproduces the physical branch-reinforcement sensitivity obtained in ANDES.'

Prohibited claim:
'ANDES independently validates the GFL portfolio hypergraph.'

Outputs: FINAL_ANDES_NETWORK_SENSITIVITY.csv, FINAL_ANDES_NETWORK_VALIDATION.md, one scatter plot, source-data CSV, manifest. Then STOP.
