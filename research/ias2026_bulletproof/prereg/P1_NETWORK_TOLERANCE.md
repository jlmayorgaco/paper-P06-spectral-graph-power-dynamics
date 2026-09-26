# P1 network gate — preregistered thresholds

The documented current-source topology is evaluated with the installed
PowerDynamics 5.0.0 API. The numerical initialization tolerances are
`tol=1e-6` and `nwtol=1e-6`; these are solver settings, not post-hoc acceptance
criteria. The acceptance criteria are fixed before the sensitivity run:

- at least 11 retained dynamic states in the initialized network state;
- full-network dynamic residual at or below `1e-6`;
- repeated initialization state delta at or below `1e-8`;
- repeated small-signal spectrum delta at or below `1e-6`;
- all three weak-shunt sensitivity magnitudes must satisfy the same gate.

No residual threshold was relaxed after the failed direct-loopback attempt.
The direct-loopback topology is retained only as a negative failure-mode record;
the executed gate uses GFL current source → LoopbackConnection → dynamic
shunt/network bus → PiLine → `Library.VδConstraint(V=1, δ=0)` slack.
