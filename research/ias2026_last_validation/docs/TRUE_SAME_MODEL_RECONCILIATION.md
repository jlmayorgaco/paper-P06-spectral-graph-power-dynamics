# True same-model Julia/Python reconciliation

status: SAME_MODEL_CROSS_CODE_VALIDATED
result_label: SAME_MODEL_CROSS_CODE_VALIDATED
evidence_class: TRUE_FROZEN_CUSTOM_MODEL_CROSS_CODE

The exact frozen custom SynchronousMachine, first-order AVR, two-state IEEEST PSS, D=0, constant Pm, no-governor, constant-power-load, and frozen IEEE-39 network conventions were implemented in Julia.

Base and full V4 equilibria were reconciled with the frozen Python oracle. State and algebraic errors are below 1e-6; critical real-part and frequency errors are below 1e-4; and critical eigenvector MAC is above 0.95. Gauge-like modes with absolute eigenvalue below 1e-4 are excluded from A_perp stability labels.

Evidence: `raw/reconciliation/same_model_reconciliation.json`, `raw/reconciliation/same_model_mode_family.csv`, and the Julia/Python state, algebraic, spectrum, Jacobian, and eigenvector files.

The all-target V4 case is the frozen H4-style portfolio for this validation run; the complete 16-portfolio enumeration remains the separate alternative-model census.
