# Final summary — ExpP

```text
EXP_P_STATUS: FAIL_TDS
LAST_COMPLETED_STAGE: P6
EXPN_MODEL_REUSED: YES; frozen PD-exact model, model SHA e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a
EXPN_REGRESSION: PASS; analytic alpha error 0, PD alpha error 9.66594e-10 1/s
PQ_CONTRACT_PASSED: YES; P0 component P/Q maximum errors 1.13687e-15 / 1.48138e-15 pu
MODEL_SHA: e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a
DESIGN_PD_CALLS: 0 before freeze; P0 used 2 independent PD comparisons, P5 used 1 independent PD spectral build plus TDS
PLL_SINGLE_DEVICE_RANK: 1; maximum residual 1.00822e-17
PLL_PAIR_AFFINE_IDENTITY: PASS conditionally per device; same-device cross term residual 2.13177e-8, reference pencil condition up to 6.00787e11; cross-device determinant terms are nonzero
REAL_BOUNDARY_HANDLED_SEPARATELY: YES; omega=0 uses a real line
SINGLE_SG_QUADRATIC_IDENTITY: PASS; rank-2 operator residual 6.93642e-18; 2-port rank-4 residual 6.09522e-18
EXPN_RETENTION_ROOT_REPRODUCED_MW: 1.1243444848621267 MW versus frozen 1.124344477653483 MW; difference 7.20864e-9 MW
CONDITIONAL_ROOT_ALL_POLES_PASS: selected bus-38 root PASS; every root was checked and roots limited by other poles were rejected; 10 bus architectures, 12 real roots recorded
NOMINAL_OBJECTIVE_TYPE: minimize dispatch-weighted retained SG MW, weights from initialized generator P components; not installed capacity
NOMINAL_RHO: [1,1,1,1,1,1,1,1,0.9986453680992127,1] for buses 30:39
NOMINAL_KP: [7.853981633974483 x 10]
NOMINAL_KI: [986.9604401089358 x 10]
NOMINAL_RETAINED_SG_MW: 1.124344477653483
NOMINAL_ALPHA: -0.0500000011050941 1/s analytic; -0.050000000138500054 1/s PD
NOMINAL_KKT_SCOPE: LOCAL_KKT_CERTIFIED on fixed support {38}, one epsilon plus 20 gains, historical guard 1e-9 1/s; not global
NOMINAL_KKT_RESIDUALS: primal 0; stationarity 0; complementarity 8.19547e-10; LICQ true, active rank 1; SOSC strict gain bounds/vacuous critical cone
NOMINAL_NUMERICAL_GUARD: frozen candidate uses historical 1e-9 1/s; proposed new-candidate guard is 1e-6 1/s. Conditional root moves cost +7.79732e-6 MW but was not jointly corrected or PD validated.
ROBUST_UNCERTAINTY_DEFINITION: frozen ExpG normalized static full-block, scale 1; scenario SHA 5c047290a6fa119c7aefc19a894da433452d27027cbcfa3befe549cd08ec242e; uncertainty SHA 70492604e2b0cada4418591605f72a8dfbae2440f8afb499af9ff4b3e298ffc1
ROBUST_BETA_REQUIRED: 1.6991206999182038e-6 (normalized)
HINF_GAMMA_LOWER: nominal observed gamma lower bound 2.6691625907696247e10; reciprocal pointwise beta upper bound 3.746493388818478e-11
HINF_GAMMA_UPPER: NONE VALID; attempted 588539.7077 rejected (CARE relative residual 4.23175e-4; max residual eigenvalue +0.0661052)
ROBUST_BETA_CERTIFIED: ExpN incumbent fails by pointwise witness; one conditional fixed-gain point has a provisional Float64 Lipschitz interval lower bound above beta_req
ROBUST_CERTIFICATE_TYPE: nominal beta_star upper bound 3.746493388818478e-11; conditional point beta_star lower bound 1.6991528445276632e-6 and sampled upper bound 1.7329480509001635e-6; interval uses no outward rounding, so not a formal validated-arithmetic certificate
ROBUST_CONDITIONAL_POINT: epsilon_38=0.0013695147207744968; rho=[1,1,1,1,1,1,1,1,0.9986304852792255,1]; Kp=[7.853981633974483 x 10]; Ki=[986.9604401089358 x 10]; retained SG=1.136697218242845 MW; alpha=-0.05169088501717586 1/s; 101 physical poles pass
ROBUST_CONDITIONAL_POINT_STATUS: scalar fixed-gain continuation only; beta lower-bound margin 3.214460945939892e-11; not robust KKT, not frozen, not independently PowerDynamics validated
ROBUST_KP: no robust optimum
ROBUST_KI: no robust optimum
ROBUST_RETAINED_SG_MW: no robust optimum
TRANSIENT_PROFILE_AND_LOCATIONS: six 0.1-s small pulses at buses 8/16/29; declared event is sustained +100 MW load step at bus 16 from t=1 s, Q fixed, output retained-SG COI frequency at bus 38
TRANSIENT_CONSTRAINTS_ENFORCED: NO in P3/P4 optimizer; independent post-freeze event fails
MAX_ROCOF: 1.072208443371192 Hz/s for the declared step; inherited limit 0.5 Hz/s
MAX_FREQUENCY_EXCURSION: 39.33340125578227 Hz absolute deviation; inherited limit 0.5 Hz
TRANSIENT_TAIL_STATUS: NOT_SETTLED_WITHIN_60S
GSP_PARAMETER_COUNT: 2 common gains; 4/6/8/10 coefficients for feature prefixes/full rank-5 basis; 20 free-gain reference
GSP_REPLACEMENT_LOSS_MW: 0 at the existing uniform incumbent only
GSP_CERTIFIED_ERROR_OR_EMPIRICAL: exact interpolation of incumbent gains (fit residual 0); no GSP optimization or self-energy truncation/remainder certificate
GLOBAL_LOWER_BOUND: 0 MW, TRIVIAL_PHYSICAL_LOWER_BOUND
GLOBAL_UPPER_BOUND: 1.124344477653483 MW for nominal complete-spectrum-only problem; none for the combined robust/transient problem
GLOBAL_GAP_MW: 1.124344477653483 MW nominal; combined gap unavailable
GLOBAL_BOUND_TYPE: valid but trivial nominal interval, not within 0.01 MW tolerance
OPEN_ARCHITECTURE_REGIONS: 1023/1024 support masks lack a joint KKT search; incumbent disconnected branches, robust/transient feasible regions, and full GSP co-design remain open
PD_FULL_SPECTRUM_VALIDATION: PASS; 101 poles at candidate, alpha error 9.66594e-10 1/s, max matched pole error 4.27729e-8 1/s, equilibrium residual 5.39530e-12
TDS_SMALL_SIGNAL_VALIDATION: PASS; 6 pulses; maximum direct linear/nonlinear frequency/RoCoF relative errors 0.00208787 / 0.00319901
TDS_DESIGN_EVENT_VALIDATION: FAIL; 100 MW sustained bus-16 step violates frequency and RoCoF limits
NONLINEAR_SECURITY_THEOREM: NOT_CLAIMED
PRIMARY_DESIGN_ITERATIONS: one fixed-support SQP iteration, zero accepted descent steps
PRIMARY_DESIGN_TIME: 34.598 s reported for P3; counters do not include separate sensitivity/QP operations
PRIMARY_DESIGN_EIGENSOLVES: 2 spectrum evaluations recorded by P3 wrapper; supporting function-call counts not fully instrumented
MAIN_VERIFIED_RESULT: ExpN P/Q and complete-spectrum identity reproduced; rank-one/rank-two/rank-four identities and conditional retention roots verified; frozen nominal point independently matches PD
MAIN_UNRESOLVED_RESULT: declared step fails; robust/transient joint optimizer and independent validation absent; conditional robust point has only provisional Float64 interval evidence without outward rounding; new 1e-6 guard not KKT/PD validated; global gap open; ExpP withheld gradient and fault-injection mutation suites not run
PUSH: NO
```
