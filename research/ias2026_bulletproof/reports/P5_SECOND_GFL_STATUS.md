# P5 — official SimpleGFLDC IEEE-39 census

component_status: PASS
component_result: raw/p5/p5_simplegfldc_result.csv
status: NEGATIVE_HOLDOUT
result_label: FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT
evidence_class: FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT
model: PowerDynamics.ComposableInverter.SimpleGFLDC
portfolios: 16 of 16 attempted
initialized: 16
stable_transverse: 16
matched_dispatch_and_rating: true
target_buses: {30,33,35,37}
gauge_filter: remove eigenvalues with |lambda| < 1e-8 before transverse stability label
no_retune: true
census_csv: raw/p5/p5_simplegfldc_ieee39_portfolios.csv

All 16 portfolios initialized and were stable after the preregistered gauge-aware transverse filter. This is a valid negative result for the official SimpleGFLDC second dynamic model. It is not a same-model cross-code parity gate and does not reproduce the frozen no-governor H4 blocker.
