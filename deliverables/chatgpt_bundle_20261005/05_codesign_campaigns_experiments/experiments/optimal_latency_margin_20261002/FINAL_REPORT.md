# Final targeted-experiment report

The preregistered experiment ends with a **decisive negative result for the proposed replacement–latency poster claim**, not a fabricated optimized frontier. The fast PLL failure is real in the tested linear DDE: full numerical contour counts change from zero violating roots to two or more near the limiting delay. However, the replacement-only, fixed-gain path has no material latency-margin loss. The gap between the two previously stored designs is almost entirely associated with their different PLL gain vectors. Under the experiment's <2 ms stop rule, further gain optimization and frontier inversion were not run.

The key comparisons (all delays are first locally tracked margin crossings, with full-contour safe/unsafe brackets) are:

| Nodal design | Gains | GFL replacement | Local crossing | Full-contour bracket |
|---|---|---:|---:|---:|
| Stored seed | Seed nominal | 87.500000% | 39.380219 ms | [39.375, 39.4375] ms |
| Same fixed-gain ρ path endpoint | Seed nominal | 88.455141% | 39.383659 ms | [39.375, 39.453125] ms |
| Stored zero-delay-tuned design | Zero-delay-tuned | 88.455141% | 37.387054 ms | [37.375, 37.4375] ms |

All fixed-gain intermediate levels are in `T04_PRECISE_CROSSINGS.csv`. The 1.993164 ms stored-design difference is below the 2 ms materiality threshold, and it confounds ρ with gains. The replacement-only local change is **+0.003433 ms** of latency margin, 580 times smaller than that threshold; the full-contour brackets alone establish that all five crossing regions lie within 0.078125 ms. The same-high-ρ in-bounds gain swap raises the local margin by **1.996605 ms**. That is an observed comparison, not an optimized delay margin, and delayed nonlinear-event feasibility was not tested.

At 40 ms, numerical exact-characteristic root counts find 6 margin-violating roots for the seed and 12 for the zero-delay-tuned high-ρ design; at 20–36 ms both have zero. The high design already has 4 at 38 ms. Active frequencies at first crossing are 5.26249 and 5.64012 Hz. Ten local delay sensitivities agree with centered finite differences to at most 7.79×10⁻⁶ relative error. The two active root right-vector MAC is 0.00182, evidencing a gain-dependent family switch; fixed-gain endpoints have MAC 0.99997. Fast branches were followed 40→20 ms; their origin at zero delay is unresolved.

Contour integration and root refinement are Float64 numerical evidence. The phrase *exact characteristic* means the physical linear DDE's exponential delays were retained; it does not mean interval-certified roots, certified global first delay, or nonlinear safety. Near ±0.001 ms of fixed-gain endpoint crossings, quadrature was unresolved (`T04_ENDPOINT_BOUNDARY_VERIFICATION.csv`); ±0.01 ms count checks passed (`T04_ENDPOINT_BOUNDARY_VERIFICATION_0P01MS.csv`). The failed finer calls remain explicitly INDETERMINATE. The previous campaign's rank≤30 action identity remains the fallback analytical result.

EXPERIMENT_STATUS: COMPLETED_WITH_WEAK_EFFECT_STOP_RULE

DELAY_TRANSITION_REPRODUCED: YES; numerical exact-characteristic contour, 20–40 ms
ROOT_COUNT_RELIABLE: NUMERICALLY_VALIDATED ON REPORTED CONTOURS; trace/phase agree; Float64, no interval certificate; ±0.001 ms spot checks INDETERMINATE

LOWER_REPLACEMENT_PERCENT: 87.500000
LOWER_REPLACEMENT_TAU_CRIT_MS: 39.3802185, first locally tracked crossing, not certified global infimum

HIGHER_REPLACEMENT_PERCENT: 88.4551408
HIGHER_REPLACEMENT_TAU_CRIT_MS: 37.3870544 with stored zero-delay-tuned gains; 39.3836594 with nominal seed gains

FIXED_GAIN_FRONTIER_MONOTONIC: LOCAL ROOT CROSSINGS INCREASE SLIGHTLY AT ALL FIVE TESTED LEVELS; tiny sign not separately certified by contour brackets
FIXED_GAIN_DELAY_LOSS_MS: −0.0034332 (local crossing; negative denotes margin gain); common full-contour bracket width 0.078125 ms

GAIN_OPTIMIZATION_SUCCESS: NOT RUN — preregistered weak-effect stop rule
OPTIMAL_GAIN_DELAY_MARGIN_RECOVERY_MS: NOT ESTABLISHED; observed in-bounds gain-vector swap at high ρ gives +1.9966049 ms, not optimum

TAU_CRIT_STAR_RHO:

| ρ (%) | Optimized τcrit*(ρ) |
|---:|---|
| 87.500000, 87.750000, 88.000000, 88.250000, 88.455141 | NOT ESTABLISHED — weak-effect stop rule |

RHO_MAX_SPECTRAL_TAU:

| Uniform delay (ms) | Optimized ρmax^spectral(τ) |
|---:|---|
| 0, 10, 20, 25, 30, 35, 40 | NOT ESTABLISHED — optimized frontier absent |

CRITICAL_MODE_FAMILY: Fast PLL-dominated cluster; seed family 3, tuned high-ρ family 5; near-zero right-vector MAC across those stored designs
CRITICAL_MODE_FREQUENCY_RANGE_HZ: 5.25994–5.64012 at reported local crossings
FAST_PLL_CLUSTER_CONFIRMED: YES, for the modeled linear DDE and tested delays
MODE_SWITCHES: YES between the two stored gain vectors; NO along the fixed-nominal-gain ρ path in the tested range

HETEROGENEOUS_DELAY_TESTED: NO — weak-effect stop rule
SAME_MULTISET_SPREAD_MS: NOT ESTABLISHED
SAME_MULTISET_REPLACEMENT_SPREAD: NOT ESTABLISHED
BEST_SPATIAL_PREDICTOR: NOT ESTABLISHED
SPATIAL_SPEARMAN: NOT ESTABLISHED

EXACT_ACTION_SPACE_DIMENSION: rank ≤30 descriptor action columns; delayed PLL return-ratio count uses 10 channels; identity established in preceding campaign
FULL_SYSTEM_DIMENSION: 282 descriptor variables, 203 physical modes
MAX_ACTION_SPACE_RECONSTRUCTION_ERROR: 3.222843298282425×10⁻¹⁶ relative, from preceding campaign; not rerun here

STRONGEST_ANALYTICAL_RESULT: The prior exact rank≤30 descriptor action identity, plus the exact simple-root local DDE sensitivity formula used here
STRONGEST_NUMERICAL_RESULT: Full-contour 0→2 transition around 39.3802/37.3871 ms and a gain-dependent fast critical-family switch
STRONGEST_NEGATIVE_RESULT: At fixed gains, 0.95514 percentage points more GFL replacement produce only +0.00343 ms local delay-margin change, not the hypothesized material loss

BEST_POSTER_HEADLINE: No validated replacement–latency headline; fallback is exact action space plus gain-dependent fast PLL mode
BEST_POSTER_EQUATION: Δ(s)=sE−A₀−ΣᵢAτᵢe^(−sτᵢ), with rank≤30 exact design-action factorization for the declared model
BEST_POSTER_NUMBER: 1.99316 ms combined stored-design gap versus −0.00343 ms replacement-only fixed-gain delay loss
BEST_POSTER_FIGURE: F02_FIXED_GAIN_RHO_TAU_FRONTIER.png, supported by F01 and F05

POSTER_READY: NO

IF NO:
EXACT_MISSING_RESULT: A material, fully evaluated optimized τcrit*(ρ) and its inverse ρmax^spectral(τ); positive-delay nonlinear events would additionally be needed for a dynamic-security/replacement claim. The preregistered materiality gate stopped the optimization rather than silently supplying these values.
