# Closure status — executed evidence

The frozen parent Z, N, and previous-best points were reconciled from stored
roots and five-event records. The stored 41 ms candidate passed. At the same
88.4551407846% GFL replacement (4779.019928 MW GFL and 623.741162 MW SG),
the largest design passing numerical full-characteristic root coverage and
all five frozen **zero-delay** nonlinear events is `full20_step15_medium`:

| Quantity | Measured value |
|---|---:|
| Z security-delay margin | 37.387187500 ms |
| N security-delay margin | 39.383359375 ms |
| Parent previous best | 40.844903557 ms |
| Best validated here | **43.797191610 ms** |
| Recovery versus Z | **6.410004110 ms / 17.144922%** |
| Zero-delay rightmost physical pole | −0.080148351 s⁻¹ |
| Worst frozen-event frequency excursion | 0.490018849 Hz |
| Worst frozen-event RoCoF | 0.205474354 Hz/s |
| Minimum normalized SG actuator slack | 0.002020306 (frozen limit 0.002) |

Five delayed modal crossings occur between 43.797191610 and
43.798385186 ms. Their near balance is numerical evidence for a multimode
minimax design, not a proof of an optimum. The local LP makes all five modal
inequalities and its actuator surrogate active, with KKT residual
`4.44e−15`; those multipliers do not certify the complete nonlinear problem.
The final design still improves materially over the preceding accepted point,
so no declared convergence criterion has been met.

The 203-state linear DDE method of steps decays at 0.90 and 0.98 times the
security-delay threshold and grows at 1.02 times it. There is **no full
nonlinear positive-delay DDE event validation**. The threshold uses
`Re(s)=−0.05 s⁻¹`, not the absolute `Re(s)=0` stability boundary. Numerical
contour coverage is not an interval-certified global root proof.

The bus-16 event has almost no actuator robustness at the tested step-13
design: 100.01 MW passes, 100.02 MW fails, and 101 MW fails. These are
one-event tests on that specific design, not a universal limit on the
architecture. The `Kp33−1%` screen at step 12 passes. Both `±1%` all-ZIP-load
runs reequilibrated but the chosen nonlinear solver did not complete; their
physical stability is indeterminate.

Multistart, a controlled generic-optimizer comparison, an optimized
nominal-versus-latency Pareto frontier, full-model KKT convergence, and a
robust design surviving the adverse event were **not established**. The
requested table files mark unexecuted work explicitly. This package provides
a strong feasible lower bound and a falsifiable multimode explanation, but
does not justify a maximum, global-optimum, or robust-safety headline.
