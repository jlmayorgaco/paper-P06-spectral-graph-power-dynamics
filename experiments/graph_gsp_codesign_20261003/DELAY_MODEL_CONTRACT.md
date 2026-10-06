# Physical model contract

Use the existing ReducedDAE model with dc_convention=:physical_supply. Its
algebraic network is eliminated exactly for this model; SG, governor, AVR,
converter current, DC and PLL differential states are retained (204 states).
No original model/source/frozen result is modified.

For generator buses 30--39, initialized full-dispatch P0 is shared as
P_SG=(1-rho)P0 and P_GFL=rho P0, with the corresponding repository reactive and
rating sharing. Network current contributions scale by these fractions. Each
device keeps its initialized full-rating internal trim. Effective physical SG
inertia/rating contribution scales with (1-rho). The denominator is initialized
generator dispatch, including the original slack-generator dispatch; it is not
net bus injection. This experiment stays strictly inside 0<rho<1, so derivatives
do not cross a change of device-state support.

The actual PLL has theta_dot=omega, omega_dot=(xi+Kp*e-omega)/tau_lpf,
xi_dot=Ki*e, where e=-sin(theta)*u_real+cos(theta)*u_imag.
Only e is replaced by its delayed value e(t-tau), in both indicated equations.
The original tau_lpf=1/(300*2*pi) seconds remains a distinct physical filter.
All ten pure delays equal 0.04 seconds, fixed during design. Gains are positive;
coordinates are rho, log(Kp), log(Ki). No power or current-control signal is delayed.

Five frozen admittance-load events are translated from event t=1 / end t=61
to event t=0 / end t=60 with equilibrium history for t<0. The actual frozen
perturbations are the repository's admittance-load changes parameterized by
initialized MW, not ideal constant-power steps. Use both polarities at 16 and
29 and the negative case at 8. The error history is zero before the event;
the algebraic voltage jump changes the measured error at t=0, reaching the PLL
at t=tau. This discontinuity requires left-limit integration at the preceding
interval endpoint and right-limit history in the following interval.

Frequency is obtained from unwrapped bus-voltage phase over a backward 0.5-s
window; RoCoF is the corresponding second backward phase difference divided
by 2*pi*(0.5 s)^2. Pre-event phase is retained. They are windowed bus metrics,
not instantaneous inertial RoCoF. Monitor all 39 buses for the full 60 s.

Acceptance uses |Delta f|<=0.5 Hz, |RoCoF|<=0.5 Hz/s, V in [0.9,1.1] pu,
normalized SG governor/AVR actuator slack>=0.002, and numerical absence of
non-gauge poles with real part>-0.05/s. The search uses slightly tighter guards.
No converter current-limiter or DC-energy safety is asserted. A finite event set
and finite horizon do not establish safety against every disturbance.
