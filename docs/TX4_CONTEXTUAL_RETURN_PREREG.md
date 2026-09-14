# TX4 contextual-return final-closure preregistration

Version: `TX4-CR-v1`  
Frozen before new campaign numerics: yes  
Parent: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`  
Branch: `research/tx4-contextual-return-final`

## Scope and stop rules

This is a narrow closure campaign on the frozen TX4 full-order IEEE-39
phasor-domain DAE. Existing TX4 outputs are read only. New outputs use the
`TX4_CONTEXTUAL_RETURN_` prefix. The campaign stops at a gate if an exact test
is impossible; `BLOCKED` is never converted into a positive claim.

Out of scope: PowerDynamics/PD39, weak nodes, weak links, structured radii,
repair optimization, planners, co-design, IEEE-68, EMT, Shapley/cumulant/
cycle/holonomy claims, and any new converter model.

## Frozen witness and coordinates

Use the exact TX4 full-order model and direct equilibrium re-solution path.

```
H4 = {30, 33, 35, 37}
P4 = (g=0.03625, k=1.425, t=1.5, h=1)
```

The g-only path fixes `(k,t,h)=(1.425,1.5,1)` and varies only `g`. The
archived boundary anchor is `g*=0.20768140450381395`; it is independently
recomputed. The in-band window is the frozen TX4 window `0.3--1.5 Hz`.
Positive-frequency selection is used. The transverse spectral abscissa is the
largest real part after removing the rotational/reference direction. A true
instability is `alpha_perp>0`; the root tolerance is `|alpha_perp|<=1e-8 s^-1`
or the solver-limited equivalent recorded in the results.

## Exact closure and contextual return

For the common fixed port basis, let `M(s)` be the exact network-mediated
action matrix, `L_i(s)=I+M_ii(s)`, and `C_S(s)=I+Q_SS(s)`. Test the identity

```
det(sI-A) = det(sI-A_hh) * product_i det(L_i(s)) * det(C_S(s)).
```

For `i in H4`, `R=H4\\{i}`, define

```
G_R = (I + Q_RR)^(-1)
R_i|R = Q_iR G_R Q_Ri.
```

The contextual-return identity is

```
det(I+Q_H4) = det(I+Q_RR) * det(I+R_i|R).
```

The theorem is a standard finite-dimensional Schur-complement consequence
under exact full-order realization, compatible block partition, finite local
factors, invertible `I+Q_RR`, and a simple tracked root. It is not presented
as a new theorem. Return eigenvalues are tracked in the fixed port basis;
singular values are reported in that basis and are not claimed basis-invariant.

## Numerical truth audit

For H4, its four one-device removals, and four fixed proper-subset controls
(BASE, singleton 30, pair 30+33, and triple 30+33+35), record equilibrium
residual, dimensions, reduced-Jacobian condition estimate, smallest tested
singular values, and the critical eigenvalue. Use tolerance factors
`{0.5,1,2}` for the frozen finite-difference step and solve tolerances
`{1e-8,1e-9,1e-10}` where supported. Construct an independent central-difference
Jacobian of the DAE residual blocks, use `scipy.linalg.eig` in addition to
`numpy.linalg.eig`, and perform a descriptor/generalized check if a valid
`(E,A)` pair is exposed. Missing descriptor data is `NOT AVAILABLE`.

Pass requires no verdict reversal, tracked-mode agreement within `1e-6 s^-1`,
and residuals within the declared tolerance. Otherwise report FAIL or BLOCKED.

## Boundary and g sweep

Recompute the H4 crossing with a fixed bracket `[0.03625,0.25]` using
safeguarded bisection/Newton. At the frozen grid

```
0.03625, 0.05000, 0.07500, 0.10000, 0.12500, 0.15000,
0.17500, 0.19000, 0.20000, 0.20500, 0.20750, 0.20768,
0.2076814045, 0.20770, 0.21000, 0.22500, 0.25000
```

solve equilibrium, compute `alpha_perp`, track the full critical mode, evaluate
`C_H`, all contextual returns, and local factors at the matched critical
frequency. Refine the return crossing with the same fixed bracket and no
post-outcome grid changes. Primary agreement target:
`|g_eig*-g_return*|<=1e-4`, or solver-limited equivalent.

The local-versus-collective gate passes only if local factors remain
nonsingular throughout and the collective factor reaches singularity at the
H4 crossing. If a local factor fails first, the collective interpretation fails.

## Minimality, remediation, and derivative

At the recomputed root evaluate all 15 proper subsets at the same complex root:
`alpha_perp`, `sigma_min(I+Q_S)`, nearest Q eigenvalue to -1, and label. H4
minimality passes only if all 15 are stable and H4 is unstable on the P4 side.

Use only the existing g coordinate. The stable after-point is frozen before
TDS as `g_after=0.25`, selected by `min(0.25,g*+0.05)` from the archived
anchor. Report the P4-to-boundary first-order prediction, Newton iterations,
final boundary, and `alpha_perp(g_after)`; do not claim economic optimality.

For `a=g`, use

```
dR/da = Q_iR,a G_R Q_Ri + Q_iR G_R Q_Ri,a
        - Q_iR G_R Q_R,a G_R Q_Ri
dmu/da = y^H (dR/da) x / (y^H x).
```

Use biorthogonal eigenvectors for the tracked simple return eigenvalue and
central finite differences with `dg=1e-5` after full equilibrium re-solution.
Agreement threshold is absolute error `<=1e-4` or the recorded numerical floor.

## Nonlinear phasor TDS

Reuse the frozen G2 IEEE-39 D2 disturbance exactly: bus-20 active-load pulse,
`+2%`, duration `0.2 s`, two solver segments, frozen G2 BDF settings,
documented guards, and existing observables. Run H4 at P4 and H4 at
`g_after=0.25`; if computationally cheap also run the four P4 triples. Use
the same disturbance, horizon, axes, and acceptance rules. This is nonlinear
phasor-domain TDS, not EMT. A contradiction weakens the claim and is retained.

## Independent phasor, conventional screen, and literature

Reuse the authoritative E31 ANDES validation without rebuilding a model.
Report only its declared phasor scope; do not upgrade it to custom-GFL EMT.

Compare singleton, pair, and triple screens; aggregate MW/MVA; g-invariant
static quantities; full H4 eigenanalysis; and exact closure. Full eigenanalysis
detects H4 when H4 is evaluated; closure supplies the collective explanation.

Perform a focused 2020--2026 primary-source literature audit covering gSCR/grid
strength, impedance/gain-phase stability, stability manifolds and operating-
point sets, root-cause methods, placement, and multi-converter coordination.
No “first” claim is permitted without direct support.

## Claim levels and final case

- Level A: exact determinant/contextual-return identities under stated assumptions.
- Level B: TX4 numerical verification of H4 minimality, g boundary, unity
  crossing, local regularity, and remediation.
- Level C: only the existing independent phasor evidence supported by ANDES.
- Level D: EMT portfolio behavior and universal cross-model/network
  generalization are not validated.

Case A requires all primary gates and TDS consistency. Case B means the narrow
mechanism is supported but an independence, generalization, or TDS gate is
mixed/limited. Case C means numerical truth, local/collective interpretation,
or the root mechanism fails/is blocked; stop the headline.
