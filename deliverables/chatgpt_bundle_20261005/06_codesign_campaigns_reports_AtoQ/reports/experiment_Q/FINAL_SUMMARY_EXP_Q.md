EXP_Q_STATUS: PASS_SECURE_GAP_OPEN

BASELINE_REPRODUCED: YES
MODEL_SHA: e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a
PQ_CONTRACT_PASS: YES

FREQUENCY_METRIC_STATUS: METRIC_DEPENDENT
PRIMARY_FREQUENCY_OUTPUT: maximum absolute local SG rotor / GFL PLL state frequency at generator buses 30:39
BUS_FREQUENCY_OUTPUT: diagnostic only; SG derivative with cubic 11-sample Savitzky–Golay filter at 100 Hz
PLL_FREQUENCY_OUTPUT: Δω_PLL/(2π) Hz
SG_FREQUENCY_OUTPUT: 60(ω_SG−1) Hz

PLL_STEADY_SUPPORT_RESULT: PLL shapes dynamics but does not create sustained primary power response
MAX_ABS_DFINF_DKP: 4.9489e-13 Hz/(MW·gain-unit)
MAX_ABS_DFINF_DKI: 1.4935e-12 Hz/(MW·gain-unit)

STEADY_FREQ_STRUCTURE: NONAFFINE
AFFINE_STIFFNESS_MAX_ERROR: relative inferred K_load range 43.9293%
K_LOAD: sample mean 123.7166 MW/Hz; inferred range [-3995.97, 1438.82] MW/Hz
C_I: TGOV1 coefficients in TABLE_Q1G_frequency_authority.csv; bus 39 is 0 (no TGOV1)

FREQ_ANALYTIC_LOWER_BOUND_100MW: NOT AVAILABLE; affine law failed
FREQ_ANALYTIC_SUPPORT: SCREENING_ONLY bus 30
FREQ_ANALYTIC_FRACTIONAL_BUS: 30, epsilon=0.2288503 in SCREENING_ONLY allocation

ROBUST_D0_STATUS: NOT OPTIMIZED; all-SG feasible reference numerically passes bounded-real CARE/LMI
ROBUST_D0_RETAINED_SG_MW: 5402.761089978847
ROBUST_D0_RHO: [0,0,0,0,0,0,0,0,0,0]
ROBUST_D0_KP: inactive; nominal values stored in frozen reference
ROBUST_D0_KI: inactive; nominal values stored in frozen reference
ROBUST_D0_ALPHA: -0.09806540933624713
ROBUST_D0_BETA_LOWER: 1.6991206999182038e-6 (Float64 bounded-real CARE/LMI; no directed rounding)

CONTINUATION_COMPLETED_TO_100MW: NO KKT continuation; isolated +100 MW feasibility test completed
CONTINUATION_POINTS: 1 all-SG frozen validation point
ACTIVE_SET_SWITCHES: none
SUPPORT_SWITCHES: none optimized; one-sided PLL boundary tests recorded for ten buses

SECURE_100MW_STATUS: FEASIBLE ALL-SG REFERENCE; NOT AN OPTIMUM
SECURE_100MW_RETAINED_SG_MW: 5402.761089978847
SECURE_100MW_GFL_MW: 0
SECURE_100MW_GFL_FRACTION: 0
SECURE_100MW_SUPPORT: {30,31,32,33,34,35,36,37,38,39}
SECURE_100MW_RHO: [0,0,0,0,0,0,0,0,0,0]
SECURE_100MW_KP: inactive
SECURE_100MW_KI: inactive

SECURE_100MW_ALPHA: -0.09806540933624713
SECURE_100MW_BETA_LOWER: 1.6991206999182038e-6
SECURE_100MW_ROCOF: 0.12383507394416611 Hz/s (PD TDS; bus 33)
SECURE_100MW_FREQ_PEAK: 0.07536489665023405 Hz (PD TDS; bus 39)
SECURE_100MW_FREQ_STEADY: 0.03755683117282334 Hz (linear model)
SECURE_100MW_SETTLING: not recorded

KKT_PRIMAL: NOT RUN
KKT_STATIONARITY: NOT RUN
KKT_COMPLEMENTARITY: NOT RUN
LICQ: NOT RUN
SOSC: NOT RUN

PD_ALPHA: -0.09806540933624897 (inherited all-SG identity regression)
PD_ANALYTIC_ALPHA_ERROR: 1.83e-15
PD_TDS_100MW: PASS at bus 16
PD_MAX_ROCOF: 0.12383507394416611 Hz/s
PD_MAX_FREQ_DEVIATION: 0.07536489665023405 Hz
PD_SETTLING: not computed

MODAL_DOMINANT_BUS: not determined for a secure optimum; ExpP rightmost branch is 99.62% GFL-state participation
FREQUENCY_DOMINANT_BUS: all-SG reference bus 39; one-sided GFL-insertion minimum is bus 39 (1.013 Hz/s), SG-removal minimum is bus 37 (9.451 Hz/s)
SAME_LOCATION: NOT APPLICABLE

BND_ACTIVE_MODE: no secure optimum; decomposition not run
BND_DIRECT_TERM: NOT RUN
BND_SELF_ENERGY_TERM: NOT RUN
BND_RECONSTRUCTION_ERROR: NOT RUN
BND_CANCELLATION_RATIO: ExpP modal-response ratios 1.974 (RoCoF), 2.526 (steady frequency)

GLOBAL_LOWER_BOUND_MW: 0
GLOBAL_UPPER_BOUND_MW: 5402.761089978847
GLOBAL_GAP_MW: 5402.761089978847
GLOBAL_CERTIFIED: NO

GSP_FREE_GAIN_PARAMETER_COUNT: 20
GSP_BEST_REDUCED_PARAMETER_COUNT: NOT RUN
GSP_REPLACEMENT_LOSS_MW: NOT RUN
GSP_SECURE_VALIDATION: NOT RUN

MAIN_THEORETICAL_RESULT: SimpleGFLDC PLL gains change dynamics but have no sustained active-power/frequency support; the exact full-model DC stiffness is nonaffine.
MAIN_POWER_SYSTEM_RESULT: A numerically robust all-SG model passes the 100 MW bus-16 event; a high-replacement ExpG seed fails local PLL RoCoF in the Q metric.
MAIN_BND_RESULT: the ExpP rightmost mode does not dominate RoCoF; modal contributions cancel materially.
MAIN_GSP_RESULT: NOT RUN because no secure optimum was obtained.
MAIN_LIMITATION: No Q active-set SQP/KKT/support search was completed. The all-SG candidate is a feasible upper bound only; the 5402.761 MW global gap remains open.
BEST_POSTER_CLAIM: Stable ExpN replacement can still fail sustained frequency/RoCoF security; local PLL RoCoF can reject a high-replacement candidate even when alpha and frequency peak pass.
PUSH: NO
