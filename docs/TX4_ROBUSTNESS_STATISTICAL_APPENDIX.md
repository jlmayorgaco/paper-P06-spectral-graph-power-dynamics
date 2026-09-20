# TX4 Robustness Statistical Appendix

## A. Master-table dimensions

The master table contains 329,440 portfolio rows over 20,590 conditions. Exact H4 rows: 20,590; exact all-portfolio calibration rows: 6,816; surrogate proper-subset rows: 302,460.

## B. QMC summary

campaign      endpoint  n_conditions  success  proportion  ci95_low  ci95_high                                                                             interpretation
     QMC    H4_PRESENT          4096   2367.0    0.577881  0.562689   0.592926 bounded-box engineering probability; proper subsets outside calibration are surrogate-tier
     QMC      EXACT_H4          4096    301.0    0.073486  0.065889   0.081883 bounded-box engineering probability; proper subsets outside calibration are surrogate-tier
     QMC NONCOMPOSABLE          4096      0.0    0.000000  0.000000   0.000937 bounded-box engineering probability; proper subsets outside calibration are surrogate-tier
     QMC delta_H4_mean          4096      NaN    0.028325 -0.095992   0.225307                                                             empirical 2.5/97.5 percentiles

## C. MC summary

campaign      endpoint  n_conditions  success  proportion  ci95_low  ci95_high                                                                             interpretation
      MC    H4_PRESENT          5000   2907.0    0.581400  0.567668   0.595007 bounded-box engineering probability; proper subsets outside calibration are surrogate-tier
      MC      EXACT_H4          5000    356.0    0.071200  0.064396   0.078662 bounded-box engineering probability; proper subsets outside calibration are surrogate-tier
      MC NONCOMPOSABLE          5000      0.0    0.000000  0.000000   0.000768 bounded-box engineering probability; proper subsets outside calibration are surrogate-tier
      MC delta_H4_mean          5000      NaN    0.029846 -0.092035   0.218205                                                             empirical 2.5/97.5 percentiles

## D. Method comparison

campaign    n  both_stable  full_only  mode_only  neither     mcnemar_p                 status
      MC 5000          823          0       1270     2907  0.000000e+00 PAIRED_MODE_COMPARISON
     QMC 4096          743          0        986     2367 3.058118e-297 PAIRED_MODE_COMPARISON

## E. Sobol summary

  endpoint parameter        S1       ST  N_base             status                                                                  deviation
H4_PRESENT         g  0.004659 0.039040    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
H4_PRESENT         k  0.510143 0.699922    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
H4_PRESENT         t  0.024651 0.054318    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
H4_PRESENT         h  0.003227 0.022243    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
H4_PRESENT   epsilon  0.492282 0.634389    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
H4_PRESENT   damping  0.001464 0.001925    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
H4_PRESENT   inertia  0.010765 0.025557    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4         g  0.083624 0.263662    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4         k -0.040083 0.869212    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4         t -0.002415 0.183221    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4         h  0.028033 0.040984    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4   epsilon -0.221062 0.747133    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4   damping  0.000464 0.021851    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations
  EXACT_H4   inertia  0.004287 0.071698    1024 SURROGATE_SALTELLI indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations

## F. Logistic summary

  endpoint status parameter  coefficient  odds_ratio    n                                         deviation
H4_PRESENT    FIT         g    -1.270247    0.280762 5000 predictive logistic fit; no causal interpretation
H4_PRESENT    FIT         k     5.262715  193.004856 5000 predictive logistic fit; no causal interpretation
H4_PRESENT    FIT         t    -4.610627    0.009946 5000 predictive logistic fit; no causal interpretation
H4_PRESENT    FIT         h    -0.334342    0.715809 5000 predictive logistic fit; no causal interpretation
H4_PRESENT    FIT   epsilon    -5.730459    0.003246 5000 predictive logistic fit; no causal interpretation
H4_PRESENT    FIT   damping    -0.235769    0.789964 5000 predictive logistic fit; no causal interpretation
H4_PRESENT    FIT   inertia     1.901780    6.697806 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT         g    -0.376294    0.686401 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT         k    -0.176499    0.838199 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT         t     0.780460    2.182476 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT         h    -0.197558    0.820733 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT   epsilon    -0.455609    0.634062 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT   damping    -0.019448    0.980740 5000 predictive logistic fit; no causal interpretation
  EXACT_H4    FIT   inertia     0.617434    1.854165 5000 predictive logistic fit; no causal interpretation

## G. File inventory

Master CSV, compressed CSV, Parquet, 1D sweeps, 2D phase Parquet, QMC/MC summaries, Morris, Sobol, logistic, g-star, method, return, mode, TDS, Julia, claim matrix, metadata, and headline JSON are retained under results/.

## H. Audit limitations

No physical uncertainty weights, physical robust radius, fresh random Julia parity, fresh same-model PowerDynamics network parity, nonlinear TDS, EMT, current-limit, DC-link, protection, or hardware certification is claimed.
