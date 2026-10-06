# M1 zero-delay search contract (frozen before evaluating candidates)

The only fully physical five-event feasibility witness available at the start is the uniform 87.5% design. The prior analytical-iteration joint design reports 90.047% but fails two of these holdout events. It is used only to define a search direction, not as feasible evidence.

For interpolation parameter `eta` in `[0,1]`, set `rho(eta)=rho_seed+eta(rho_joint-rho_seed)`. Interpolate `Kp` and `Ki` geometrically, equivalently linearly in log-gain coordinates. This keeps gains positive and within the unchanged frozen bounds. The predeclared coarse evaluation grid is `eta = 0.25, 0.50, 0.75, 1.00`, with the already reproduced seed as `eta=0`. These are search coordinates, not physical delay values. No other gain bounds, events, model states, or security thresholds change.

Each candidate must pass equilibrium, complete finite spectrum, and all five full PowerDynamics events under the M0 contract. A failed or incomplete event makes that candidate infeasible. The best fully validated point among executed candidates is retained. There is no monotonicity assumption, global-optimum claim, local-optimum claim, or certified upper bound from this one-dimensional search.

After the coarse grid, a single bisection refinement between the highest feasible evaluated `eta` and the nearest higher evaluated infeasible `eta` may be run to improve the feasible witness. This refinement is fixed now, before seeing candidate results. It does not establish a one-dimensional global maximum or full-space stationarity.
