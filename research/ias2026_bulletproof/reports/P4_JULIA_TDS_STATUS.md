# P4 — Julia common small-disturbance TDS

status: PASS_WITHOUT_MODAL_MATCH
evidence_class: FRESH_OFFICIAL_POWERDYNAMICS_TDS
cases: base=none, one=30, full=30+33+35+37
disturbance: common 1e-6 perturbation of the first network state; the highest-variance state trace is selected for decay/frequency estimation; identical solver/tolerances
eigensolve_comparison: raw/p4/p4_pd_tds.csv
traces: raw/p4/tds_<portfolio>.csv when the solve completed
This is official-model TDS execution evidence. All three trajectories solved and selected-trace decay diagnostics were recorded. The selected trace did not provide a reliable oscillation crossing, so measured frequency is `NaN` and no frequency agreement is claimed; the decay diagnostic is not treated as the critical eigendecay. It is not a true same-model frozen-no-governor TDS gate and does not promote the H4 blocker.
