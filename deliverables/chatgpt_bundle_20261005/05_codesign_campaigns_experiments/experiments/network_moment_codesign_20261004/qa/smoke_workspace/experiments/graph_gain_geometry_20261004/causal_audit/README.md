# Instantaneous causal gate: negative evidence

Run from the repository root:

```powershell
julia --project=. --startup-file=no experiments/graph_gain_geometry_20261004/causal_audit/run_audit.jl
```

This script performs 30 exact-model floating-point algebraic evaluations and
reads five existing baseline trajectories. It generates two CSV tables and a
source/output hash manifest in this directory. It runs no time trajectories.
The protocol preserves the exploratory status and fixed event/rho grid.

## Result and boundary

All tested instantaneous cases pass. At common rho=1, all five frozen events
pass: worst frequency is approximately 0.012190 Hz and worst RoCoF is
0.024380 Hz/s, against limits 0.5 and 0.5. Voltages lie approximately in
[0.966167,1.084829] pu, within [0.9,1.1]. Exact generated values and numerical
parity diagnostics are in TABLE_01_INSTANTANEOUS_ALGEBRAIC.csv.

Consequently, this instantaneous necessary condition cannot exclude complete
replacement at this particular common-rho endpoint. It supplies no positive
minimum retained-SG requirement. This is not a global-rho feasibility result:
no exhaustive rho search, nonlinear trajectory acceptance, spectral stability,
monotonicity theorem, or optimal replacement claim follows. Indeed, the phase
jump is nonmonotone on the tested common-rho grid.

TABLE_02_EXISTING_BASELINE_FIRST_40MS.csv summarizes the five saved samples
at 0, 0.01, 0.02, 0.03, and 0.04 seconds. These are sampled extrema of existing
baseline data, not continuous-time extrema or new integration results.

## Exact algebraic identity

In full real/imaginary voltage coordinates let D_i be the frozen SG Norton
matrix embedded at bus i, v* the pre-event voltage, and

    G_rho,e = Y_0 + DeltaY_e + sum_i (1-rho_i) D_i.
    DeltaY_e = -delta/(100 Vset_e^2) E_e E_e^T.
    v_e(0+) = v* - inv(G_rho,e) DeltaY_e v*.

The script verifies this full-network identity against ReducedDAE's original
Kron voltage solve, trim invariance across rho, and unchanged voltage under a
second allowable PLL gain choice. The event delta is a constant-impedance
parameter, so actual disturbed power depends on voltage. This calculation
tests algebraic SG Norton strength, not a direct inertia inequality.

The exact rho derivative within a nonsingular region is

    d v_e(0+) / d rho_i = inv(G_rho,e) D_i [v_e(0+) - v*].

## Causality and metric interpretation

For positive minimum PLL detector delay and equilibrium prehistory, all delayed
phase errors remain zero before that delay. DelayedEvents replaces the current
detector term by the delayed one, so subtracting B(K)e(x) from the ordinary RHS
removes every Kp/Ki dependence on this interval. Initially omega=xi=0; PLL
angles stay fixed. Current, DC, and SG controls still respond. This independence
applies only to the frozen PLL gain variables and controller architecture.

The frozen frequency/RoCoF metrics are 0.5-second phase windows with zero
pre-event phase deviation. At zero time, with phase jump dtheta,

    F(0+) = dtheta/(2*pi*W),
    R(0+) = dtheta/(2*pi*W^2),   W=0.5 s.

The same relation R=F/W holds for times before both the window and minimum
delay. Raw differentiation of the algebraic phase jump would create an impulse;
that quantity is not the declared windowed metric. A constant common phase
rotation cancels in arg(v/v*). A time-varying reference must retain its
accumulated phase to preserve the declared metric.

A future 40-ms gain-independent trajectory check could strengthen this
necessary test, but has not been run here. Any rigorously bounded violation
would reject all allowable PLL gains at that specified rho; passing this test
would remain only a necessary condition.
