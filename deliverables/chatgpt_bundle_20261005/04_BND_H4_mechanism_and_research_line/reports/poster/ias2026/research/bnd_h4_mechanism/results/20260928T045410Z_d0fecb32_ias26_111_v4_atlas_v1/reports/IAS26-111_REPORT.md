# IAS26-111 - Exact V4 blocker atlas

- Run: `20260928T045410Z_d0fecb32_ias26_111_v4_atlas_v1`
- Source: `20260927T144629Z_d0fecb32_ias26_060_operating_v1`
- Status: **PASS_WITH_LIMITATIONS**
- Type: post-hoc descriptive reanalysis; no new power-flow, DAE, or eigenvalue solves.
- Frozen threshold: `tau_dec = 1e-08 s^-1`.
- Valid scenarios: 967; physically invalid IDs retained without replacement: 33.
- Outcome-dependent filtering: none.

## Reconciliation gate

The exact required partition is 412 no-blocker + 190 H4-only + 365 alternative-witness = 967 valid scenarios. The recomputation produced 412 + 190 + 365 = 967. The scenario-level blocker-mask mismatch count against frozen `SCENARIO_METRICS.csv` is 0; the `(x_s,y_s)` mismatch count is 0.

## Results

| Exact minimal blocker set | Scenarios | Share of valid cohort |
|---|---:|---:|
| `NONE` | 412 | 42.61% |
| `H4_ONLY` | 190 | 19.65% |
| `{30,33,35} | {30,33,37}` | 151 | 15.62% |
| `{30,33,35} | {30,33,37} | {33,35,37}` | 137 | 14.17% |
| `{30,33,35}` | 77 | 7.96% |

The minimum blocker order counts are reported in `../derived/BLOCKER_ORDER_COUNTS.csv`; the complete row-level result is in `../derived/SCENARIO_ATLAS.csv`. For each valid scenario, `x_s=max_{R proper subset H4} alpha_perp(R)`, `y_s=alpha_perp(H4)`, and `m_comp=min(y_s,-x_s)`. The latter is a signed coordinate margin in the `(x,y)` plane, not a probability.

H4 was unstable in 555 valid scenarios. It was minimal in 190 and nonminimal in 365 alternative-witness scenarios. Every alternative witness is represented in `COUNTEREXAMPLES.csv`.

## Scope boundary

This atlas establishes descriptive witness switching only within the frozen synthetic ensemble and the 16 V4 portfolios. It does not establish mechanism invariance, bus-33 causality, V9/IEEE-39-wide generality, or naturalistic failure probabilities. No poster source, poster figure, or historical run was modified. See `../claims/SAFE_CLAIMS.md` and `../claims/FAILED_OR_UNSUPPORTED_CLAIMS.md`.
