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
matched_dispatch_PQ_system_base: true
matched_rating: false
rating_equivalence_established: false
dynamic_rating_basis: global_100_MVA_system_base
candidate_MVA: reported_from_retired_machine_ratings_only
no_retune: true
census_csv: raw/p5/p5_simplegfldc_ieee39_portfolios.csv

This is a second dynamic model negative result with matched scheduled P/Q on the common 100 MVA system base. Candidate MVA is reported from retired machine ratings, but dynamic rating equivalence is not established. It is not a same-model cross-code parity gate and does not reproduce the frozen no-governor H4 blocker.
