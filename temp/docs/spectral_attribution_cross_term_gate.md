# Spectral Attribution Framework and E0 Cross-Term Gate

Date: 2026-06-08

This note records the current state of the "spectral attribution" reframing.
It accepts the strategic point that the strongest conceptual object is not an
individual perturbation formula, but the common operator

```text
T(s) = s^2 M + s Sigma(s) + L_tilde,
Sigma(s) = D0 + C (sI - Acc)^(-1) B.
```

The contribution should be framed as an attribution framework only if it
produces at least one measurable effect that cannot be obtained by treating
network geometry, inertia, and control as separate levers.

## Candidate Irreducible Effect

The candidate is the inertia-control cross term.  For a local perturbation

```text
M(rho) = M0 + rho Delta M,
D(rho) = D0 + rho Delta D,
```

the off-diagonal modal channel is

```text
H1(s0) = s0 Gamma1 + Lambda1,
```

where `Gamma1` is the damping/control channel and `Lambda1` is the
stiffness/inertia channel in the unperturbed inertial-Laplacian basis.  The
modal-coupling part of the second-order pole shift is

```text
s2_coup =
  1/(2s0 + delta_c) *
  sum_{l != c}
    [(s0 E_D[c,l] + E_L[c,l]) (s0 E_D[l,c] + E_L[l,c])]
    / (s0^2 + delta_l s0 + nu_l).
```

The non-additive term is

```text
s0 ( E_D[c,l] E_L[l,c] + E_L[c,l] E_D[l,c] ).
```

This term is irreducible in the local expansion: it vanishes if either the
inertia channel or the control/damping channel is absent, and it is not obtained
by adding an inertia-only pole shift to a control-only pole shift.

## Registered Empirical Gate E0

The implementation is:

```text
validation/ieee39_rational_filter/phase_cross_inertia_control.py
```

The registered test uses full ANDES eigensolve poles as the reference.  The Schur
partition is used only to report control participation of the tracked mode.

For a half step `h`, estimate:

```text
dM      = [s(M+h, C) - s0] / h
dC      = [s(M, C+h) - s0] / h
dCross  = [s(M+h,C+h) - s(M+h,C) - s(M,C+h) + s0] / h^2
```

Then test the predictions at `eps = 2h`:

```text
s_add   = s0 + eps (dM + dC)
s_cross = s_add + eps^2 dCross.
```

The cross term is marked measurable only if:

- the cross contribution is at least 5% of the true joint shift;
- the cross prediction improves the complex-pole error or damping-ratio error;
- the cross prediction improves over the additive prediction.

## Current Result

Initial ANDES test:

```text
python validation/ieee39_rational_filter/phase_cross_inertia_control.py --mode network --eps 0.01
```

Result:

- status: `NOT_MEASURABLE`
- base mode: `-0.656641 + j4.88482`, `0.777443 Hz`, `zeta=0.133226`, `pi_c=0.5479`
- `|true joint shift| = 1.556843e-02`
- `|cross shift| = 3.485490e-06`
- cross fraction: `2.23882e-04`
- additive complex error: `2.804225e-05`
- cross complex error: `3.152763e-05`
- additive zeta error: `4.679823e-07`
- cross zeta error: `5.318994e-07`

Singleton sweep:

```text
python validation/ieee39_rational_filter/phase_cross_inertia_control.py \
  --mode network --eps 0.01 --sweep-singletons \
  --out outputs/phase_cross_inertia_control_sweep_network
```

Result:

- 13 registered tests: all-active plus all active `GENROU.M x PLL1.Kp/Ki`
  singleton pairs;
- measurable cross cases: `0/13`;
- strongest cross fraction was the all-active case: `2.23882e-04`.

Critical/control-family tests were also negative:

- cross fraction: `3.83427e-09`;
- complex-error improvement: approximately `1.3e-06`;
- zeta-error improvement: negative.

## Interpretation

The cross term is mathematically real in the local expansion, but it is not
empirically measurable in the current Mix60/no-PSS ANDES gate with registered
1% GENROU-inertia and PLL-gain perturbations.  Therefore it should not be used
as the paper's empirical headline at this point.

The spectral-attribution framework remains useful as an organizer of network,
inertia, and control levers, but the current full-ANDES evidence supports the
controller-aware self-energy bridge and control-resonant reinforcement audit
more strongly than the inertia-control cross term.

## Next Honest Options

1. Search for a physically justified case where virtual inertia and PLL/droop
   retuning are actually coupled in the same controller, rather than perturbing
   retained `GENROU.M` and `PLL1` gains separately.
2. If no measurable case appears, keep the cross term as a theoretical local
   interaction and do not title the paper around it.
3. Continue with Experiment A, because it still decides whether the estimator
   itself has practical advantage over the diagonal all-mode screen.
