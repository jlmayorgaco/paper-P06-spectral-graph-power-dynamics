# TX4 Robustness Statistics Audit - Final Report

## 1. Executive verdict

Case B. The exact reduced-DAE instability benchmark remains strong, while the corrected minimal-blocker story is materially smaller and composite blockers are common. The surrogate tier failed its strict gate, so exact QMC/MC all-16 recomputation was required and completed.

## 2. What was wrong in the previous campaign

The previous H4_PRESENT flag measured H4 alpha_EM instability. It did not enumerate the H0 minimal-blocker antichain. NONCOMPOSABLE represented failure states rather than minimal blockers of cardinality at least two. eta_H4 used a proper-subset alpha margin rather than the registered local/collective sigma-minimum definition. Proper-subset rows outside calibration were surrogate-derived.

## 3. Flag-definition audit

U0={S: alpha_EM(S,c) >= 0} and H0 is its inclusion-minimal antichain. EXACT_H4 implies H4_PRESENT implies H4_UNSTABLE and implies NONCOMPOSABLE. Six synthetic tests pass.

## 4. H4 exact instability

Exact QMC coverage is 0.577881 (2,367/4,096). Exact MC probability is 0.581400 (2,907/5,000), Wilson 95% CI [0.567668, 0.595007]. These are instability endpoints, not minimality endpoints.

## 5. Exact-versus-surrogate provenance

The master has 329,440 rows: 26,980 exact legacy rows and 302,460 surrogate rows. Every row carries EXACT_DAE, SURROGATE_EXTRATREES, DERIVED_FROM_EXACT, or UNKNOWN provenance. Primary QMC/MC blocker results use only the new exact all-16 files.

## 6. Surrogate minimality validation

The reconstructed five-fold ExtraTrees audit found 18 false-minimal and 14 missed-minimal conditions; exact-H4 agreement was 0.924883 and H0-antichain exact match was 0.701878. The strict gate failed.

## 7. Exact QMC recomputation

All 4,096 retained QMC coordinates were evaluated for all 16 portfolios: 65,536 exact rows.

## 8. Exact MC recomputation

All 5,000 retained MC coordinates were evaluated for all 16 portfolios: 80,000 exact rows. No samples or seeds were regenerated.

## 9. Correct H4_PRESENT

Exact QMC coverage is 0.113770; exact MC probability is 0.115400 with Wilson CI [0.106838, 0.124553].

## 10. Correct EXACT_H4

Exact QMC coverage is 0.113770; exact MC probability is 0.115400 with Wilson CI [0.106838, 0.124553]. For full-set H4, H4_PRESENT and EXACT_H4 coincide because any proper blocker prevents H4 from being minimal.

## 11. Correct NONCOMPOSABLE

Exact QMC coverage is 0.427246; exact MC probability is 0.436600 with Wilson CI [0.422907, 0.450391].

## 12. KAPPA

The KAPPA PMF is reported for 0, 1, 2, 3, 4, and NULL in the exact QMC/MC summary tables.

## 13. delta_H4

The true-minimality margin is min(alpha_EM(H4), -max proper alpha_EM). QMC median is -0.085269 s^-1; MC median is -0.084810 s^-1 with BCa 95% CI [-0.088820, -0.081605].

## 14. Corrected eta_H4

The exact 128-row stratified eta audit gives median eta 0.111537 and all nonnegative values. It uses min_i sigma_min(I+Q_(H4\i)(s_H)); no physical-radius interpretation is permitted.

## 15. Collective/local validation

The same exact stratum has minimum local sigma 0.294663 and median collective sigma 0.111537, with zero failures.

## 16. QMC statistical interpretation

QMC values are fixed-design coverage fractions over the retained scrambled Sobol coordinates, not probability claims.

## 17. MC statistical interpretation

MC values are assumed independent-uniform bounded-box probabilities with Wilson intervals; no calibrated physical weights exist.

## 18. Sensitivity-provenance audit

Morris remains approximate screening and Sobol remains surrogate-based. Neither is promoted to exact sensitivity evidence.

## 19. g_star provenance

g_star is an exact-H4 local phase-boundary diagnostic, not a physical robust radius.

## 20. Negative results

No U005/H005 threshold was present in the retained preregistration/code, so no such endpoint was invented. The audit does not establish EMT, switching, current-limit, DC-link, protection, hardware, universal IEEE-39, or physical uncertainty claims.

## 21. Corrected poster implications

Use the exact wording in the handoff: distinguish instability from minimality and report composite blockers. The ten corrected F1-F10 figures are the poster-ready audit evidence.

## 22. Corrected paper implications

Replace legacy H4 presence language with H0 definitions, exact provenance, KAPPA, NONCOMPOSABLE, and exact QMC/MC results. Preserve the bounded-box and model-scope limitations.

## 23. FINAL CASE A/B/C

CASE B: exact H4 instability remains strong, but corrected minimality is smaller and composite blockers are common.

All 32 invariant checks PASS. No push was performed.
