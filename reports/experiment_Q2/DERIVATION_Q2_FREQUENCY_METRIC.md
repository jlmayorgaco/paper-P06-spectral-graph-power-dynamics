# Q2 frequency measurement derivation and gating

## Passive output

For a complex nominal bus voltage `V0=ur+j ui`, the first variation of its phase is

```text
δθ = (-ui δur + ur δui) / |V0|².
```

The architecture-invariant electrical-frequency output is

```text
δf_grid = (1/(2π)) d(δθ)/dt
        = C_grid x + D_grid ΔP_load
```

after eliminating algebraic voltages with the regular nodal block `Gy` and applying the validated rotational gauge quotient. `C_grid` and `D_grid` are constructed for the same bus-voltage channels (30–39) in all-SG, mixed, and all-GFL models; no controller-internal PLL or surviving-machine weighting enters the output. The implementation is `src/bnd_expQ2/GridFrequency.jl` using the fixed-support reduction in `src/bnd_expQ/LinearSecurity.jl`.

For a step, the algebraic voltage change at the event is

```text
δv(0+) = -Gy⁻¹ y_load ΔP.
```

It produces

```text
δθ(0+) = (1/(2π))[-ui, ur] δv(0+) / |V0|².
```

At the ExpN point this phase jump is `0.0012662112 rad` at bus 38 for +100 MW. Therefore an ideal continuous `dθ/dt` includes a Dirac impulse, and its unfiltered RoCoF has no finite supremum. Finite differences at 100 Hz and Savitzky–Golay derivatives with different windows are distinct measurement operators; they yield distinct RoCoF peaks. No inherited file freezes a grid-frequency filter or a load-step rise time for the 0.5 Hz/s limit.

## Reduced step response

For a fixed architecture, with `Aq,Bq` the gauge-quotiented reduced matrices and the passive output `(Cq,Dq)`, the finite part of the step response is

```text
δf_grid(t) = Cq Aq⁻¹(exp(Aq t)-I) Bq ΔP + Dq ΔP,
δf_grid(∞) = (-Cq Aq⁻¹Bq + Dq) ΔP.
```

The derivative for `t>0` is `Cq exp(Aq t)Bq ΔP`; it excludes the event impulse and is labeled the smooth finite-state RoCoF. This distinction is why a numerical finite-difference sample must not be advertised as the unfiltered continuous peak.

## Gain/DC check

At synchronized steady state the local GFL PLL integrator enforces frequency synchronization while the active-power current/DC equilibrium sets dispatched P. The PLL gains occur in the dynamic equation and do not alter the steady active-power share under the frozen trim contract. The 30-point centered test found maximum absolute DC sensitivities `7.66e-14` for Kp and `1.25e-14` for Ki in the sampled units. They do shape finite-time poles and filtered RoCoF.

## Affine stiffness falsification

The tested scalar `Keff=-1/mean(H0)` is a diagnostic reduction of the exact multi-output `H0=-C A⁻¹B+D`. Across 200 deterministic mixed/support designs, a 100-fit/86-holdout affine model in retained-SG fractions reached 178% maximum relative holdout error. It is classified NONAFFINE and is not a necessary-condition certificate or a pruning bound.
