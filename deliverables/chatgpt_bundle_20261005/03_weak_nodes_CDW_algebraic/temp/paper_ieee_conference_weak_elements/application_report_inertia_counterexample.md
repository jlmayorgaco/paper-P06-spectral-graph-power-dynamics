# Application Witness: Inertia Can Hurt at the Most Participating Generator

This note records a compact IEEE 39-bus/ANDES witness for the claim that a
traditional inertia-placement intuition can fail.

## Case

- Source audit: `outputs/phase_inertia_bridge_gate`
- Case: `base` (`andes/cases/ieee39/ieee39_full.xlsx`)
- Perturbation: increase documented `GENROU.M` by 1%
- Full ANDES eigensolve is the reference
- Controller-aware bridge: Schur NEP, validated against full ANDES
- Base critical pole: `s = -1.346013 + j 8.611485`
- Base damping margin: `zeta_min = 0.154429`
- Base frequency: `1.37056 Hz`
- Critical family: `I_network`

## Why a traditional rule would pick bus 30

The critical mode has its largest state participation at `GENROU_1`:

- `omega GENROU 1: 0.3682`
- `delta GENROU 1: 0.3422`
- next largest generator state is much smaller:
  `omega GENROU 8: 0.09063`

A common engineering heuristic would therefore add inertia at the most
participating machine, i.e. `GENROU_1` at bus 30.

## What actually happens

| Action | Bus | M change | Post zeta | Delta zeta | NEP local derivative | Full ANDES derivative |
|---|---:|---:|---:|---:|---:|---:|
| Add inertia at most participating generator (`GENROU_1`) | 30 | 8.400 -> 8.484 | 0.153794 | -6.3557e-4 | -6.4704e-2 | -6.3557e-2 |
| Add inertia at best audited generator (`GENROU_8`) | 37 | 4.860 -> 4.9086 | 0.154848 | +4.1836e-4 | positive | +4.1836e-2 |

Thus the naive participating-generator action **reduces** damping, while a
different generator increases damping.  The controller-aware NEP predicts the
negative sign at bus 30 and matches the full ANDES finite difference within about
1.8% in magnitude.

## Interpretation

This is not a claim that inertia is bad.  It is a counterexample to the
monotone rule "add inertia at the critical participating generator."  The sign of
the damping-margin sensitivity depends on the controller-aware eigenvalue
sensitivity

```math
\frac{\partial s_k}{\partial M_i}
=-\frac{y_k^*\,S_{M_i}(s_k)\,x_k}{y_k^*\,S_s(s_k)\,x_k},
```

not only on modal participation.  The denominator contains the control
self-energy derivative, so a local mass change can rotate the critical network
mode into a less damped direction even when it is applied at the dominant
machine.

## Claim discipline

- Verified: documented `GENROU.M` perturbation in full ANDES IEEE 39-bus.
- Not claimed: adding a brand-new generator with its own AVR/governor and power
  dispatch.  That is a finite network/model change and must be audited
  separately.
- Not claimed: GFM virtual inertia, because the tested `REGF1` sheets do not
  expose an unambiguous virtual-inertia parameter.
