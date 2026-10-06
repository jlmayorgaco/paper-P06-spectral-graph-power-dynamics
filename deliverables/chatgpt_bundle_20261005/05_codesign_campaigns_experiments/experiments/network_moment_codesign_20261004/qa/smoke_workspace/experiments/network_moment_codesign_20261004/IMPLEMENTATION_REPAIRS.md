# Repairs before scientific evaluation

The first Python run subtracted one from the already zero-based ports.csv
indices. Its hidden block was singular and the run stopped before producing
any moment or design result. `moments_v0_failed.py`, the original input lock,
and the v1 lock preserve this correction. No experiment definition changed.

The first multimode run stopped when coarse root continuation converged two
seeds to the same pole. The collision was not accepted as a valid catalog.
The preserved `codesign_multimode_v0.py/log` and initial multimode lock record
that failure. Version1 uses analytic eigenvalue prediction, adaptive path
steps, and a minimum root separation check. Design constraints are unchanged.

The next run reached the frozen60-iteration limit at uniform+1ms and stopped
without saving the trial. Version2 saves every iteration and returns that
last candidate for the same independent feasibility checks, explicitly
tagged MAX_ITERATIONS_NO_CONVERGENCE_CLAIM. It does not increase the budget,
relax any constraint, or declare the candidate feasible from solver status.
