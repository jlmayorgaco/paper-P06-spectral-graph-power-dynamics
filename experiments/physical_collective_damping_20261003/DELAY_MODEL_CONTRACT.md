# Actual delay channel and retained physical dynamics

For candidate i, the source model states include theta_i, omega_PLL,i and xi_i.
Define the actual PLL detector, after exact algebraic voltage elimination,

    e_i(x) = -sin(theta_i) v_r,i(x) + cos(theta_i) v_i,i(x).

The implemented delayed equations are

    dot(theta_i) = omega_PLL,i,
    T_PLL dot(omega_PLL,i) = xi_i + Kp_i e_i(x(t-tau_i)) - omega_PLL,i,
    dot(xi_i) = Ki_i e_i(x(t-tau_i)).

T_PLL=1/(2 pi 300) s is an existing physical low-pass state time constant. It is
not the pure delay tau_i, and was not changed. The whole scalar error is delayed,
including its past angle reference. This is distinct from delaying terminal
voltage alone and applying a present-time angle transform. The latter is not
the declared model and was not simulated.

No delay is inserted in P/Q setpoints, voltage injection, current control or
DC control. The network algebraic voltage at each history time is solved from
the state at that time. This campaign has no load-switching events, so its
network parameters are constant throughout every history and simulation.

The nonlinear system is the physical_supply source convention; gains and rho
are the three previously frozen design vectors. Positive pure delays are fixed
uniform 20/30/40/45 ms for spectral diagnostics; nonlinear tests use 40/45 ms and
the zero-delay reference. Zero-delay simulations invoke the undelayed RHS directly.
Heterogeneous delay is not tested in this campaign and no claim is made for it.

The ten SGs retain positive ratings. Their internal per-unit equations stay
unchanged, while network injection and physical mechanical inertia use the
retained rating (1-rho_i) S_i. Initialized dispatch P_i^0, not net bus injection,
is used for reported SG and GFL MW.
