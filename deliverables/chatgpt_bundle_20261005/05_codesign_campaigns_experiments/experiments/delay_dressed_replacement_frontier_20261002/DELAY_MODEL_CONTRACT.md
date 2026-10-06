# Delay model contract (pre-result freeze)

## Frozen no-delay model

Use the checked-in ten-port PowerDynamics IEEE-39 design lineage, with candidate generator ports at buses 30:39, 100 MVA and 60 Hz. Keep the ExpN dispatch-weighted component semantics: each initialized device P/Q is split at the frozen operating voltage, SG rating is scaled by epsilon_i=1-rho_i while its H is held fixed, and the GFL terminal-current port is weighted by rho_i. The equilibrium uses initialized device dispatch, never net-bus injection.

The full physical PLL uses the detector signal
e_i(t) = -sin(theta_i(t))*u_r,i(t) + cos(theta_i(t))*u_i(t).
The stock PLL also has a frequency-output first-order low-pass with tau_lpf=1/(2*pi*300) s. This filter remains in the model and is not the exogenous measurement/computation delay studied here.

## Delayed channel

For each installed GFL, delay only the actual phase-detector error before both PI paths:
dot(DeltaOmega_i) = (xi_i + Kp_i*e_i(t-tau_i) - DeltaOmega_i)/tau_lpf
dot(xi_i) = Ki_i*e_i(t-tau_i)
dot(theta_i) = DeltaOmega_i

All CC1, filter, DC-link, network-voltage, active/reactive setpoint, SG, AVR, and governor equations remain current-time equations. No other signal is delayed. At tau_i=0, the vector field is algebraically identical to the existing ReducedDAE PLL equations. Initial history for t<=0 is the constant pre-event equilibrium state, with e_i=0 at lock.

tau_i is fixed and exogenous in every optimization. Variables are rho_i, Kp_i and Ki_i only. No delay is fitted, tuned, optimized, or changed in response to a result.

## Gain and operating limits

Use existing ExpN bounds unchanged: 0.25*Kp0 <= Kp_i <= 4*Kp0 and 0.25*Ki0 <= Ki_i <= 4*Ki0, where Kp0=2*pi*5 and Ki0=(2*pi*5)^2/4. Record every bound-active gain. The support endpoints rho_i=0 or 1 change state dimension and must be handled as separate supports; do not differentiate through structural state deletion.

The historical event and output contract is frozen from reports/analytic_iteration_20261001/refined/protocol.toml:
- design: +100 MW load step at bus 8;
- five external cases: -100 MW at bus 8; +/-100 MW at buses 16 and 29;
- frequency peak <=0.5 Hz and 0.5-s estimator window;
- filtered/windowed RoCoF <=0.5 Hz/s with the same 0.5-s metric;
- modal decay margin >=0.05 s^-1;
- voltage 0.9-1.1 pu;
- minimum normalized SG governor/AVR actuator slack 0.002.
These are the campaign's declared limits, not claims of grid-code compliance. Hard GFL current and DC-energy limits have not been certified and cannot be claimed.

## Delay families frozen before delayed results

Uniform delay coordinates are exactly [0, 2, 5, 10, 20, 30, 40, 50] ms. They are screening coordinates, not a claim about typical or measured converter latency.

For spatial placement, freeze the ten-value multiset tau_j = (j-1)*50/9 ms for j=1,...,10 (mean 25 ms; same values, variance and maximum for every permutation). The graph-smooth, graph-rough, low-GSP-frequency, high-GSP-frequency, and 100 random permutations will be generated deterministically and hashed before any delayed pole, optimization, or time-domain result is evaluated. The random seed is 20261002. Pattern-generation rules and all resulting vectors belong in FROZEN_DELAY_PATTERNS.json. If the exact physical-angle operator fails graph-Laplacian tests, graph placement labels are exploratory and no Laplacian/GSP claim may use it; the separately defined lossless approximation must be named as such.

## Baseline and event gates

D0 requires the no-delay ReducedDAE and independently compiled PowerDynamics trim, complete finite ODE spectrum, rightmost physical eigenvalue and nonlinear event to match declared tolerances. The no-delay delayed-channel implementation must reproduce the same vector field and poles at tau=0.

D1 requires nonlinear-eigenvalue residual checks, tracked root identity through delay continuation, and a numerical root-count/coverage check over a declared rightmost region. A set of Newton roots seeded only by the zero-delay poles is not complete DDE spectrum. No Pade model is primary truth. If root coverage is unresolved, record BLOCKED_EXACT_DDE_SPECTRUM and do not optimize or label a capacity maximum.

Every candidate called fully feasible must pass the small-signal/DDE constraints, the design event, and all five frozen external events. Incomplete integration is not proof of instability. A sampled delay grid is not a certified upper bound.
