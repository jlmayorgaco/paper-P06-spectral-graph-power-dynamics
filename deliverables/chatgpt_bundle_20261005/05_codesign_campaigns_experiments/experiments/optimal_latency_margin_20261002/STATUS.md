# Gate and claim status

| Stage | Status | Evidence / reason |
|---|---|---|
| Frozen inputs and definitions | PASS | `FROZEN_PROTOCOL.json`, source SHA-256 hashes; five fixed-gain vectors frozen before evaluating them. |
| 30–40 ms transition | NUMERICALLY_VALIDATED | `T01`: seed 0 violating roots at 36/38 ms, 6 at 40 ms; high design 0 at 36 ms, 4 at 38 ms, 12 at 40 ms. |
| Rightmost-root coverage | NUMERICALLY_VALIDATED ON TESTED CONTOURS | At all 13 safe transition-grid points, no roots were counted to the right of candidate +0.001 s⁻¹ (`T01_RIGHTMOST_EXCLUSION.csv`). Unsafe counts match refined conjugate fast roots. Float64, no interval certificate. |
| Fixed-design crossing | SUPPORTED_LOCAL | `T02`: seed 39.3802185 ms, high zero-delay-tuned design 37.3870544 ms. The full contour switches 0→2 at ±0.01 ms around each (`T02_BOUNDARY_VERIFICATION.csv`). |
| Critical fast-family continuation | SUPPORTED_LOCAL | Nine positive-imaginary roots tracked from 40 down to 20 ms with neighboring right-vector MAC ≥0.999993. Origin as τ→0 remains unresolved. |
| Exact local delay sensitivity | NUMERICALLY_VALIDATED | Ten selected bus/design derivatives match centered FD; maximum relative error 7.79×10⁻⁶ (`T03_TAU_DERIVATIVE_VALIDATION.csv`). Kp/Ki/ρ columns are centered relinearized finite differences, not analytical derivatives. |
| Fixed-gain replacement curve | SUPPORTED_LOCAL | Five local roots cross 39.38023–39.38366 ms; each has full-contour [39.375, 39.453125] ms bracket. |
| Materiality | WEAK_EFFECT / STOP | Combined stored-design difference 1.99316 ms <2 ms. Replacement-only fixed-gain local difference is −0.003433 ms (slight margin gain); even common contour bracket width is 0.078125 ms. |
| Gain optimization and τcrit*(ρ) | BLOCKED_WEAK_EFFECT_STOP_RULE | Protocol forbids spending large compute after weak fixed-gain effect. No optimized gain vectors or global gain-optimality claim. |
| Inverse ρmax^spectral(τ) | BLOCKED_WEAK_EFFECT_STOP_RULE | Requires validated optimized frontier; no extrapolation. |
| Heterogeneous delay / spatial predictor | BLOCKED_WEAK_EFFECT_STOP_RULE | Triggered only by meaningful uniform-delay result. |
| Delayed nonlinear event safety | NOT TESTED | Optional; no reliable method-of-steps event campaign executed here. |
| Poster-ready replacement–latency frontier | NO | Missing optimized τcrit*(ρ), inverse frontier, and positive-delay event evidence. |

The failed ±0.001 ms full-contour spot checks are preserved in `T04_ENDPOINT_BOUNDARY_VERIFICATION.csv` as **INDETERMINATE**, never reclassified SAFE. The contour quadrature became unresolved near the boundary. Local exact-characteristic root refinement resolves the crossing to about 0.000061 ms, but does not turn that tiny slope into a global/root-count certificate.
