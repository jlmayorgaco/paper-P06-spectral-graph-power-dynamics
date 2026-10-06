# Causal algebraic audit: fixed negative-evidence gate

This audit preserves an already explored, exact algebraic calculation. It is
not a preregistered discovery campaign or a global optimization result.

The hypothesis tested is whether the instantaneous, gain-independent response
to the existing five load events excludes full SG-to-GFL replacement under the
existing voltage and 0.5-second phase-window frequency/RoCoF limits.

The frozen events are (bus, delta MW) = (8,-100), (16,+100), (16,-100),
(29,+100), (29,-100). Delta is the original constant-impedance load-event
parameter, not an imposed voltage-independent power injection. The frozen
common-rho grid is [0, 0.875, 0.95, 0.99, 0.999, 1]. Every combination is
evaluated. No event or grid point is selected retrospectively.

Use ReducedDAE.jl with physical_supply DC convention, fixed initialized
dispatch and trim, baseline Kp=0.9*(2*pi*5), Ki=(2*pi*5)^2/4, and all 39 buses.
At t=0+, dynamic states retain their pre-event trim; solve only the algebraic
voltage equation. The window is W=0.5 seconds, with zero pre-event phase
deviation. Limits are F<=0.5 Hz, R<=0.5 Hz/s, and voltage in [0.9,1.1] pu.

Check the full-network algebraic jump identity against the original Kron
implementation, pre-event trim invariance across rho, and direct invariance
under a second allowed PLL gain choice. Read the already saved baseline
trajectories only to summarize samples from t=0 through 0.04 seconds. Run no
new trajectories and do not infer continuous-time extrema from those samples.

Report every algebraic case, the existing-trajectory summaries, source hashes,
Julia version, runtime, and numerical parity checks. Files outside this audit
directory are read-only. The all-rho=1 result decides only whether this
instantaneous necessary test can exclude that particular endpoint. No
global-rho bound, full-spectrum result, nonlinear feasibility claim, or
inference about an optimal replacement limit is permitted.
