# P2 — PowerDynamics V4 portfolio census

status: NEGATIVE_HOLDOUT
result_label: FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT
evidence_class: FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT
portfolios: 16 of 16 attempted
initialized: 16
target_buses: {30,33,35,37}
matched_dispatch_and_rating: true
dynamic_model: official PowerDynamics IEEE-39 machines with official AVR/governor devices; GFL11 replacement harness
same_model_cross_code_gate: NOT_TESTED
network_tolerance: 1.0e-8
census_csv: raw/p2/p2_powerdynamics_portfolios.csv
hasse_csv: raw/p2/p2_hasse_edges.csv
mode_shapes_csv: raw/p2/p2_mode_shapes.csv

All requested portfolios initialized and had stable gauge-aware transverse spectra. This is a fresh alternative dynamic-model negative holdout. It does not reproduce or refute the frozen custom no-governor H4 model and is not TRUE_SAME_MODEL_CROSS_CODE_PASS.
