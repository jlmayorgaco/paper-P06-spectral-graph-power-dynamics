# Deviations

No post-hoc retuning has been performed.

## D001 — runtime package-version mismatch

The frozen closure record reports Python 3.13.14 with NumPy 2.5.2, SciPy 1.18.1, and pandas 3.0.5. The current runtime reports NumPy 2.3.5, SciPy 1.16.3, and pandas 2.3.3. The audit therefore treats frozen artifacts as retrospective evidence and records runtime checks separately. No expected value is changed to fit the current runtime.

## D002 — isolated worktree

The active user worktree had 981 staged/modified entries. The campaign was executed in a linked worktree from the verified freeze so those changes were not touched. This is a provenance-preserving execution choice, not a scientific model change.

## D003 — uncompleted dependent gates

The official PowerDynamics package was available and its IEEE-39 tutorial
equilibrium gate passed. Its model class is not the frozen L0 GFL model, so the
same-model parity, second-converter, physical-robustness, new blind-holdout,
and new Julia-TDS gates were stopped or left `NOT_TESTED`. These are explicit
scope outcomes, not substitutions for the missing experiments.
