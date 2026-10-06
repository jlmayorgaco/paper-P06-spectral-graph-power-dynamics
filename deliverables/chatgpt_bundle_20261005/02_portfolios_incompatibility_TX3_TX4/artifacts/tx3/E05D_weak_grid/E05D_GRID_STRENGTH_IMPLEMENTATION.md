# E05D conventional grid-strength implementation

The diagnostic follows the conventional NERC definition `SCR = S_sc / S_rated`. On system base,
`S_sc,pu = V_POI^2 / |Z_th|`, hence `SCR_i = V_i^2/(|Z_th,i| S_rated,i,pu)`. The passive positive-sequence
Ybus retains lines, transformers, line charging, and fixed shunts; constant-power loads and the three GFL
current sources are opened, while synchronous source buses 30--35 and slack 39 are ideal voltage sources
(incrementally grounded). `Z_th` is the corresponding diagonal of the inverse reduced Ybus. Each GFL rating
is 1000 MVA on the 100-MVA system base. This is a diagnostic, not a claim gate, and conventional SCR can be
optimistic when nearby IBRs interact. Reference: NERC, *Short-Circuit Modeling and System Strength* (2018),
https://www.nerc.com/globalassets/programs/rapa/ra/short_circuit_whitepaper_final_1_26_18.pdf.

Tests in `tests/unit/test_e05d_weak_grid.py` verify the one-source analytical Thevenin formula, inverse
impedance scaling, and the exact 34/12 topology classification.
