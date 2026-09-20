# Second-model policy discovery

status: `FRESH_SECOND_MODEL_POLICY_DISCOVERY`
model: `PowerDynamics.ComposableInverter.SimpleGFLDC`
policy coordinates: PLL bandwidth multiplier and current-loop bandwidth multiplier
grid: 0.125, 0.5, 1.0, 2.0, 8.0
representative portfolios: 30+33+35 (proper stable reference), 30+33+35+37 (flagship blocker candidate)

| coordinate | scale | initialized | stable | unstable | mixed |
|---|---:|---:|---:|---:|---|
| pll | 0.125 | 2 | 2 | 0 | false |
| pll | 0.5 | 2 | 2 | 0 | false |
| pll | 1.0 | 2 | 2 | 0 | false |
| pll | 2.0 | 2 | 2 | 0 | false |
| pll | 8.0 | 2 | 2 | 0 | false |
| current | 0.125 | 2 | 2 | 0 | false |
| current | 0.5 | 2 | 2 | 0 | false |
| current | 1.0 | 2 | 2 | 0 | false |
| current | 2.0 | 2 | 2 | 0 | false |
| current | 8.0 | 2 | 2 | 0 | false |

No mixed stable/unstable policy region was found on this frozen second-model grid. This negative search result is recorded rather than converted into a mixed holdout.

All classifications use gauge-aware transverse spectra, a common matched-dispatch system-base operating point, and the corrected rating scope: matched_rating=false and rating_equivalence_established=false.
