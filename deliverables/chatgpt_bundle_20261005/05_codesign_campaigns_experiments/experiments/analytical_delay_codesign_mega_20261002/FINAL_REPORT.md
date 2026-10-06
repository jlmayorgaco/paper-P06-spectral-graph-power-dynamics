# Final gated report

The intended analytical delay-aware maximum replacement method has **not** passed its final gates. The experiment recovered a defensible zero-delay witness and an exact low-rank action identity, and it corrected a serious DDE root-count diagnostic from earlier work. It did not establish a delayed maximum frontier or an optimal tuning vector.

The best fully validated zero-delay point **found on the frozen search path** is 4779.019928 MW GFL (88.4551408%) with 623.741162 MW retained SG. The equilibrium residual is `1.52e-11`; 203 physical poles have rightmost real part `−0.0801498 s⁻¹`. All five frozen 61-second nonlinear events pass. The worst frequency excursion is `0.485513 Hz` in `bus16_minus100`; the smallest SG actuator slack is `0.00330036` in `bus16_plus100`, above the frozen `0.002` minimum. The next higher tested point (88.7735210%) fails actuator slack at bus 16 +100 MW. The historic 90.0470421% joint point fails both actuator and frequency holdouts. This brackets one searched path, not the global problem.

For fixed designs, the exact linear DDE characteristic has an unexpected high-frequency failure. The rank-10 return-ratio trace integral and independently accumulated determinant phase both count zero roots to the right of `−0.05 s⁻¹` at 20 and 30 ms uniform PLL delay. At 40 ms they count 6 for the 87.5% seed and 12 for the 88.455% design. Three and six distinct positive-imaginary roots, respectively, were refined on the exact nonlinear characteristic; their conjugates exhaust those counts numerically. The unstable roots lie around 5.18–5.40 Hz, while the zero-delay rightmost low-frequency mode is about 0.015 Hz. Thus tracking only that low-frequency pair misses this PLL cluster. These calculations are Float64 and do not include delayed nonlinear event trajectories or an interval certificate.

The exact descriptor update for fixed-model SG→GFL/PLL design factors through at most 30 action columns. Across 21 accepted designs and three frozen heterogeneous delay fields, maximum reconstruction error is `3.23e-16`; determinant-lemma log-magnitude error is `1.20e-12`. The 20 random in-bounds designs were not screened by events and therefore cannot be called feasible. Twenty deliberately out-of-box gain draws were rejected by the source model API before factorization. This is a protocol shortfall, openly recorded. The conditional PI target-root law also reproduces its prescribed roots to `1.02e-16`, but only 5/45 gain pairs lie within the frozen bounds; it is not a nearest-pole or global design law.

The previously proposed zero-frequency graph damping Schur expansion remains invalid for the current partition: the hidden block has condition `3.095842e18` after rotational-mode deflation. A physically defensible finite-frequency angular-port damping model, controller-authority calibration, and reduced-to-full predictor have not been validated. The reduced LP two-anchor theorem is mathematically correct under its stated two-constraint structure, but no full nonlinear partial-retention claim follows from it.

The strongest next experiment is narrow: implement an actual nonlinear PLL measurement-delay method of steps, validate the five frozen events at 20–40 ms, then optimize only in a delay region where the exact characteristic count is zero. A second path is to turn the 30-dimensional action factor into a measured computational advantage and a useful tuning predictor. Without one of those bridges, this is research evidence rather than a competitive maximum-replacement poster.

```text
MEGA_EXPERIMENT_STATUS: PARTIAL_GATED; no delay-specific maximum or global optimum

BASELINE_REPRODUCED: YES; 87.5% GFL, five events, 203 physical poles
BEST_ZERO_DELAY_VALIDATED_GFL_PERCENT: 88.45514078456176% (best found on frozen path)
BEST_ZERO_DELAY_RETAINED_SG_MW: 623.7411615845335

EXACT_ACTION_SPACE: rank <= 30 descriptor factorization, exact algebraic identity for declared model
FULL_MODEL_DIMENSION: 282 descriptor; 203 physical differential modes
ACTION_SPACE_DIMENSION: <= 30 (numerical difference rank 30 in tested samples)
MAX_RECONSTRUCTION_ERROR: 3.222843298282425e-16 relative, 21 accepted designs x 3 delay fields

DDE_ORACLE: NUMERICALLY_VALIDATED at tested fixed uniform-delay designs; NOT interval certified
ROOT_COUNT_METHOD: balanced 10x10 return-ratio phase + trace integral + exact-characteristic root refinement at 40 ms
POSITIVE_DELAY_CASES_CERTIFIED: 0 directed-rounding certificates; 6 numerical cases evaluated
INDETERMINATE_CASES: arbitrary heterogeneous patterns and untested designs/delays

SINGLE_MODE_CLOSED_FORM: exact conditional PI target-root equation only; reduced design law unvalidated
CLOSED_FORM_DIRECT_LIFT_ERROR: NOT MEASURED
TWO_MODE_VARIANT_REQUIRED: UNDETERMINED
MINIMAL_MODAL_CORE_SIZE: UNKNOWN
ROUCHE_CERTIFICATE: NONE

CONTROLLER_AUTHORITY: symbolically derived; not calibrated against full-model optimum
REDUCED_RETENTION_LP: derived, physical coefficient validation blocked
PARTIAL_RETENTION_THEOREM: PROVED for reduced LP with two aggregate active constraints only
REDUCED_OPTIMAL_PARTIAL_SGS: NOT COMPUTED

ANALYTICAL_PREDICTOR: NOT VALIDATED
FULL_MODEL_CORRECTOR: NOT IMPLEMENTED; M1 is a frozen interpolation search
FULL_MODEL_ITERATIONS_REQUIRED: UNKNOWN
FINAL_PREDICTION_ERROR_PP: NOT MEASURED

MAX_VALIDATED_GFL_0MS: 4779.019928394313 MW / 88.45514078456176% (best found, not maximum)
MAX_VALIDATED_GFL_20MS: NOT ESTABLISHED; no nonlinear DDE events
MAX_VALIDATED_GFL_40MS: NOT ESTABLISHED; tested fixed designs are spectrally unsafe
MAX_VALIDATED_GFL_50MS: NOT ESTABLISHED; not evaluated

UNIFORM_DELAY_EFFECT_MW: NOT DETERMINED as a replacement frontier
UNIFORM_DELAY_EFFECT_PP: NOT DETERMINED

HETEROGENEOUS_DELAY_EFFECT: NOT DETERMINED
SAME_MULTISET_RANGE_MW: NOT DETERMINED
SAME_MULTISET_RANGE_PP: NOT DETERMINED

GRAPH_DAMPING: zero-frequency Schur invalid for tested partition; finite-frequency model not validated
COMMUTATOR_PREDICTOR: NOT VALIDATED
DELAY_SENSITIVITY_PREDICTOR: NOT VALIDATED
BEST_SPATIAL_PREDICTOR: UNKNOWN
BEST_SPATIAL_SPEARMAN: NOT MEASURED

DELAY_SHADOW_PRICE_RANGE_MW_PER_MS: NOT ESTABLISHED

BEST_FULLY_VALIDATED_DESIGN: eta_0375_bisect at tau=0, best found only
rho*: NOT CERTIFIED; best-found vector in M1_ZERO_DELAY_DESIGN.toml
Kp*: NOT CERTIFIED; best-found vector in M1_ZERO_DELAY_DESIGN.toml
Ki*: NOT CERTIFIED; best-found vector in M1_ZERO_DELAY_DESIGN.toml
retained_SG_MW: 623.7411615845335
GFL_MW: 4779.019928394313
GFL_percent: 88.45514078456176
critical_alpha: -0.08014983362991729 s^-1 at tau=0
critical_frequency: 0.014964234743335707 Hz at tau=0
worst_event: bus16_minus100 by frequency; bus16_plus100 by actuator margin
max_frequency_deviation: 0.48551314387927513 Hz
max_RoCoF: 0.21672650008249233 Hz/s

OPTIMALITY_BOUND: feasible lower witness only at tau=0; no nontrivial upper bound
R_F: 4779.019928394313 MW at tau=0 only
R_U: NOT ESTABLISHED
GAP: NOT ESTABLISHED

STRONGEST_THEOREM: exact rank-<=30 descriptor action factorization; reduced-LP <=2-partial-SG theorem under its assumptions
STRONGEST_NUMERICAL_RESULT: 40 ms exact-characteristic count 6 vs 12 roots, matched by 3 vs 6 refined conjugate pairs
STRONGEST_NEGATIVE_RESULT: zero-delay retuning increases replacement by 51.604 MW yet tested 40 ms PLL delay creates twice as many margin-violating roots; zero-frequency Schur damping invalid

RECOMMENDED_POSTER_STORY: No award-ready maximum-replacement claim. Provisional action-space/PLL-cluster research lead only.
POSTER_READY: NO

IF NO:
EXACT_BLOCKER: no nonlinear delayed five-event validation, no delay-specific co-design frontier or analytical predictor, no certified upper replacement bound; spatial placement untested

IF YES:
HEADLINE_CLAIM: NOT APPLICABLE
HEADLINE_EQUATION: NOT APPLICABLE
HEADLINE_NUMBER: NOT APPLICABLE
HEADLINE_FIGURE: NOT APPLICABLE
```
