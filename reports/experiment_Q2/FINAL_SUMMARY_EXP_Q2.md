# Final summary — Experiment Q2

```text
EXP_Q2_STATUS: FAIL_FREQUENCY_METRIC

BASELINE_PASS: YES (14/14)
MODEL_SHA: e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a
PQ_CONTRACT_PASS: YES (max trim error 5.73e-13 pu in F1)

GRID_FREQUENCY_METRIC_STATUS: PASS_CONTINUOUS; RoCoF application METRIC_DEPENDENT
GRID_FREQUENCY_DEFINITION: 1/(2*pi) d(unwrapped arg(Vbus))/dt at buses 30:39; passive, same channels in every architecture
GRID_FREQUENCY_FILTER_OR_ESTIMATOR: no frozen project bandwidth; 100 Hz sampling, 0/0.04/0.10/0.20/0.40 s sensitivity; 0.10 s is diagnostic only
RHO_ZERO_CONTINUITY_MAX_ERROR: 5.618558e-11 at rho=1e-8 (condition-aware tolerance 2.997007e-5)

EXPN_GRID_FREQ_100MW_PEAK: 21.254654 Hz (0.10 s diagnostic, analytic)
EXPN_GRID_ROCOF_100MW_PEAK: unbounded continuous unfiltered (phase jump 0.001266211 rad); smooth finite-state term 12334.913 Hz/s; sampled/filtered estimates 2.018155 / 1.399921 / 0.803533 / 0.783376 / 0.742607 Hz/s
EXPN_GRID_FREQ_100MW_STEADY: 23.104449 Hz (analytic)

PLL_DC_SUPPORT_RESULT: PLL_DC_SUPPORT_INDEPENDENT
MAX_ABS_DFINF_DKP: 7.656891e-14 Hz/MW per gain unit
MAX_ABS_DFINF_DKI: 1.245934e-14 Hz/MW per gain unit

STEADY_FREQ_STRUCTURE: NONAFFINE
AFFINE_MODEL_MAX_ERROR: 1.780236 relative; 3298.2703 MW/Hz absolute on 86 finite holdouts

ROBUST_D0_STATUS: NOT_RUN (secure problem stopped at unresolved Rmax measurement)
ROBUST_D0_RETAINED_SG_MW: NOT_RUN
ROBUST_D0_SUPPORT: NOT_RUN
ROBUST_D0_RHO: NOT_RUN
ROBUST_D0_KP: NOT_RUN
ROBUST_D0_KI: NOT_RUN
ROBUST_D0_ALPHA: NOT_RUN
ROBUST_D0_BETA_LOWER: NOT_RUN

CONTINUATION_TO_100MW_PASS: NO
CONTINUATION_POINTS: 0
ACTIVE_SET_SWITCHES: 0
SUPPORT_SWITCHES: 0

SECURE_100MW_STATUS: NOT_SOLVED — R_peak constraint has no frozen finite measurement definition
SECURE_100MW_RETAINED_SG_MW: NOT_RUN
SECURE_100MW_GFL_MW: NOT_RUN
SECURE_100MW_GFL_FRACTION: NOT_RUN
SECURE_100MW_SUPPORT: NOT_RUN
SECURE_100MW_RHO: NOT_RUN
SECURE_100MW_KP: NOT_RUN
SECURE_100MW_KI: NOT_RUN

SECURE_100MW_ALPHA_ANALYTIC: no secure candidate; ExpN reference -0.050000001105 s^-1
SECURE_100MW_BETA_LOWER: no secure candidate; ExpN pointwise witness 3.746493e-11 < beta_req (failure)
SECURE_100MW_GRID_FREQ_PEAK: NOT_RUN for candidate; ExpN reference 21.254654 Hz analytic
SECURE_100MW_GRID_FREQ_STEADY: NOT_RUN for candidate; ExpN reference 23.104449 Hz analytic
SECURE_100MW_GRID_ROCOF: METRIC_DEPENDENT / unbounded ideal unfiltered step

KKT_PRIMAL: N/A
KKT_STATIONARITY: N/A
KKT_COMPLEMENTARITY: N/A
KKT_DUAL_FEASIBILITY: N/A
LICQ: N/A
SOSC: N/A

PD_ALPHA: -0.050000000139 s^-1 (existing ExpN independent validation)
PD_ALPHA_ERROR: 9.666e-10 s^-1
PD_GRID_FREQ_PEAK: NOT_RUN at 100 MW under the new estimator
PD_GRID_FREQ_STEADY: NOT_RUN at 100 MW under the new estimator
PD_GRID_ROCOF: NOT_RUN at 100 MW under the new estimator
PD_100MW_BUS16_PASS: NO under previous ExpN COI/PLL metric; new bus-frequency RoCoF remains undefined without bandwidth
PD_OUT_OF_SAMPLE_STATUS: NOT_RUN (no frozen secure candidate)

MODAL_AUTHORITY_BEST_BUS: NOT_RUN
FREQUENCY_AUTHORITY_BEST_BUS: NOT_RUN
SAME_BUS: NOT_RUN

BND_ACTIVE_MODE: ExpN rightmost pole -0.050000001105 s^-1
BND_DIRECT_TERM: NOT_RUN
BND_SELF_ENERGY_TERM: NOT_RUN
BND_TOTAL_SENSITIVITY: NOT_RUN
BND_RECONSTRUCTION_ERROR: NOT_RUN

GLOBAL_LOWER_BOUND_MW: NOT_AVAILABLE
GLOBAL_UPPER_BOUND_MW: NOT_AVAILABLE
GLOBAL_GAP_MW: NOT_AVAILABLE
GLOBAL_CERTIFIED: NO

MAIN_THEORETICAL_RESULT: PLL Kp/Ki change transient dynamics but sampled dF_inf/dK is at numerical zero; ideal bus-phase RoCoF of an instantaneous step is distributional.
MAIN_POWER_SYSTEM_RESULT: ExpN remains spectrally stable but fails grid-frequency peak and steady-offset limits; its finite RoCoF estimate also exceeds 0.5 Hz/s at each tested window.
MAIN_BND_RESULT: rightmost mode dominates peak and steady offset; fast modes dominate smooth RoCoF.
MAIN_LIMITATION: no frozen measurement bandwidth or load-step rise-time semantics for the 0.5 Hz/s constraint.
POSTER_MAIN_NUMBER: ExpN filtered grid-frequency peak 21.254654 Hz at 100 MW; ideal unfiltered RoCoF unbounded.
POSTER_MAIN_CLAIM: stable does not imply frequency secure, and a finite RoCoF claim requires a frozen measurement/event bandwidth.
PUSH: NO
```
