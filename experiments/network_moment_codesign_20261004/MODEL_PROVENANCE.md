# Model provenance and scientific boundary

This experiment reads the existing full physical-supply IEEE39 SG/GFL DAE,
then uses its exact index-one algebraic elimination. Source hashes, git HEAD,
dirty state and package versions are in BASELINE_MANIFEST.json and
JULIA_VERSIONS.csv. Neither historical source nor any frozen result is changed.

- Julia physical model: `experiments/nonlinear_codesign_20261001/ReducedDAE.jl`.
- Nonlinear pure detector-delay integrator:
  `experiments/graph_gsp_codesign_20261003/DelayedEvents.jl`.
- Parametric Python matrices and roots: read-only sources from
  `experiments/graph_gsp_codesign_20261003/model/` and
  `experiments/interaction_decision_20261004/model.py`.
- Fixed replacement/gain starting point:
  `experiments/interaction_decision_20261004/designs/corrected.toml`.
- Complete-region spectral oracle:
  `experiments/analytical_delay_codesign_mega_20261002/m3_a_trace_integral.jl`.
- New full-state disturbance/output export: `export_inputs.jl` and `model/`.

The fixed ten-bus replacement weights the SG Norton contribution by1-rho and
GFL current contribution by rho, with the original per-device equations and
initial dispatch unchanged. The objective denominator is initialized generator
dispatch5402.761089978847MW, not net bus injection. No rho is optimized here.
Ki is fixed for controller design; selected Ki perturbations are analytical
identity tests only. The exogenous delay changes are never decision variables.

Delay enters the actual PLL phase-detector error only. The filtered type-II
PLL is theta_dot=omega, tf*omega_dot=xi+Kp*e(t-tau)-omega,
xi_dot=Ki*e(t-tau), tf=1/(300*2*pi). Current loops, P/Q references and other
device signals are not delayed. No Pade approximation supplies ground truth.

The all-PLL hidden system has174 states; all ten PLL triples are retained in
the equivalent port equation. Bus-frequency output is the0.5s voltage-phase
window over all39 buses, with direct algebraic event feedthrough included.
Inputs are changes of constant-impedance load parameters;100MW denotes the
declared parameter conversion at initial voltage, not an exact time-varying
constant-power withdrawal. Five events and all thresholds are inherited.

For literal Taylor coefficients the rounded exported matrices have a small
rotational defect. `moments.py` explicitly applies G0 -> G0(I-11T/10) and the
corresponding phase-output correction T0*1=1. Raw matrices remain untouched.
Finite-frequency204D parity is checked on the raw export separately.
`verify_moment_balls.py` proves regularity/rank/weight signs for this explicit
symmetry-restored binary model. It does not turn independent rounded entries
into an exact uncertain physical model or certify robust operational safety.

The mathematical port reduction and coefficient laws are exact identities
of the full linearized model under their assumptions. The time-integral
interpretation additionally requires stable nongauge dynamics. The nonlinear
event tests are separate numerical evidence, not an exact nonlinear moment law.
