# Multimode restoration after the single-mode test failed

The frozen one-root corrector succeeded at uniform +0.1ms, missed another
unstable catalog root at +0.5ms, and stopped on an infeasible QP at +1ms.
`single_mode_attempt/` preserves its code, log, designs and tables.
That is a negative result for the one-critical-mode assumption.

The revised algorithm uses all eleven inherited catalog branches and checks
branch separation, with the same exact moment equality, design bounds,
delays and events. Its QP has a nonnegative common restoration slack and
penalty1000; each step is at most5% of the original Kp. Up to60 iterations.
Every catalog root must end below -0.0600001/s. A finite catalog is still not
a complete-spectrum certificate. The full contour is the independent gate.
The final comparison is exploratory because this algorithm was revised
after observing the single-mode failure; no superiority/generalization claim.
