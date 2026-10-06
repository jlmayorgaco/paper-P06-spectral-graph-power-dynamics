# Experiment F7 — Frozen-design PowerDynamics validation

**F7_STATUS: BLOCKED_NO_FROZEN_CANDIDATE**

F7 requires a hashed F3–F6 analytical design. F2 stopped the graph spectral
program as `NON_MODAL`, so there is no candidate to validate and no
PowerDynamics result is reported. The validation hash guard is available in
`src/bnd_validation/GraphPLLPowerDynamicsValidation.jl`; it does not import
PowerDynamics. The SimpleGFLDC component has no hard current limit, so
`HARD_CURRENT_LIMIT = NOT_MODELED` and no saturation-safety claim is made.
