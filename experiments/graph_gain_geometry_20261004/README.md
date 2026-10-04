# Gain geometry and causal-floor audit — preserved negative results

This exploratory campaign asked whether the current network admits a useful
all-gain exclusion of a critical pole, and whether initial algebraic frequency
metrics impose a synchronous-support floor. It did not establish a maximum
SG-to-GFL replacement or an all-gain stability barrier.

- `geometry.py`: exact fixed-frequency gain-parallelogram conditions lifted
  to a rank-one Hermitian matrix, then relaxed to an SDP.
- `sectors.py`: stronger real-sector lift with reality and RLT constraints.
- `radius_screen.py`: boundary point screens around preceding finite decisions.
- `causal_audit/`: exact initial algebraic event map and small-time checks.

Both full historical-gain-box SDP screens exclude 30 of50 sampled
rho/frequency points numerically. The band3--6Hz remains inconclusive. These
are point exclusions; they are neither a continuous contour certificate nor
a sampled certificate of global feasibility/optimality. Six actual roots
pass the primal consistency checks and are not falsely excluded. Some solver
statuses are `optimal_inaccurate`; recomputed Hermitian slacks are preserved.

The causal initial-event metric passes even at100% GFL under the existing
0.5s window. It supplies no useful retained-SG floor here. That instantaneous
algebraic check does not establish the stability of a100% GFL trajectory.

The structured-uncertainty lift follows established real-mu/IQC ideas. No new
general robust-control theorem is claimed. The scientific outcome is that
these particular relaxations do not close the requested all-gain bound.
The following network-moment campaign is a separate question, not a renamed
success of this experiment.

Reproduce from the repository root with Python3.13, numpy/scipy, pandas and
cvxpy/CLARABEL: run `geometry.py`, `sectors.py`, then `radius_screen.py`.
The causal audit has its own Julia protocol and input hashes. Frozen input
locks bind the original definitions; rerunning is not required to read results.

Status: **NEGATIVE_RESULT** for a useful support floor or continuous all-gain
barrier; **NUMERICALLY_VALIDATED** for the finite point screens only.
