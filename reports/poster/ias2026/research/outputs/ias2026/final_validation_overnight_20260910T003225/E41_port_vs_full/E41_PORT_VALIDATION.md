# E41 — port model against the full-state model: **PASS**

If the port representation did not reproduce the full model's engineering
conclusion across the whole lattice, the theory would be an elegant description
of something other than this system. It does, to ten significant figures.

18 cases: base, all 15 proper subsets, the flagship, the RC converter repair and
the 25 % synchronous condenser.

| check | result |
|---|---|
| max eigenvalue error, port NEP against full DAE | **7.3e-10** |
| median eigenvalue error | 1.3e-10 |
| max damping error | 6.4e-10 |
| max frequency error | 5.6e-11 Hz |
| `σ_min(T(λ))/σ_max(T(λ))` at the full model's own mode, worst | **2.7e-13** |
| inter-area family modes visible at the ports | **2 of 2 in every case** |
| band modes invisible at the ports | **0** |
| band right-half-plane count agrees | **every case**, flagship 1 = 1 |
| dimension | 110 dynamic states → 78 port coordinates → **8×8 action operator** |

Two independent things are established here.

**Visibility.** Every mode of the tracked inter-area family is a zero of
`det T(s)`: the relative smallest singular value of `T` evaluated at the full
model's own eigenvalue never exceeds `2.7e-13`. No inter-area mode is confined to
a device and invisible at its terminals. The port operator is not missing part of
the physics that matters here.

**Accuracy.** Solved as a nonlinear eigenvalue problem in its own right — secant
iteration on the smallest eigenvalue of `T(s)`, no explicit inverse of `T` formed
anywhere — the port model's own answer differs from the full model's by at most
`7.3e-10`. That is roundoff, not approximation.

**The engineering conclusion transfers exactly.** The port model puts the
flagship's single unstable band mode in the right half plane and puts none there
for any of the fifteen proper subsets, for the converter repair or for the
condenser. The classification a reviewer cares about is identical in the reduced
representation.

This is what justifies the theory over brute-force simulation: the same verdict
is reached from an operator on 8 action coordinates that a 110-state eigenvalue
problem reaches, and the operator carries an invariant — the closure distance —
that the state-space spectrum does not expose.

## Files

`E41_port_vs_full.csv`, `E41_port_vs_full.png`, `E41_figure_source.csv`,
`manifest.json`.
