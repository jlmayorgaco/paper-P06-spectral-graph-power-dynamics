# Selected-design validation, fixed before Julia execution

The selected numerical pair is damping_candidate_6.json by the recorded rule.
Its two exact interval root boxes already prove violation of the inherited
-0.05/s requirement. Four continuous common-contour counts establish one root
in every baseline/singleton/joint model. The weak-coupling walk-series criterion
fails on that contour and is not used as a certificate.

Run independent Julia equilibrium/ForwardDiff parity and numerical all-root
contour checks on anchor, both singletons, joint, complex_pair and corrected.
Expected results are predictions only: joint should have a conjugate pair to
the right of -.05; the other designs may pass or fail. Save every result.
The full-spectrum contour oracle is floating point, not an interval certificate.
Delays remain exact exponentials, fixed at .04s; gauge is quotiented explicitly.

Only after the corrected design passes the all-root numerical spectral check,
run exactly the five historical load events, 60s, dtmax=.01s, tolerance1e-9,
Rodas5P method of steps, frequency and RoCoF with .5s windows at all39 buses.
Limits: F<=.5Hz, R<=.5Hz/s, V in[.9,1.1], retained SG actuator slack>=.002.
No current-limiter or robust uncertainty claim. The corrected design was not
selected using these event outcomes. Do not retune it if an event fails.

The correction distributes an added modal-decay budget equally across the
two sites while their rho and prescribed +2% Kp remain fixed. Its target is
-.06/s, providing .01/s guard beyond the inherited requirement. The solved
scalar fixed point is not a global/local maximum of replacement or a minimum
gain-effort problem. Sampled contraction slopes are not uniform certificates.

Retain all exploratory solver failures. New evidence does not revise the
historical negative or baseline results. Fix implementation bugs in new
versions/validation records; do not interpret unresolved cases as feasible.
