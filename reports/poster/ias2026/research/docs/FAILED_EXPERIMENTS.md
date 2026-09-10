# FAILED AND NEGATIVE RESULTS

Kept because they constrain what may be claimed.

## F1 — The IEEE-9 triple is not an irreducible third-order effect

Expected, from the toy case, that the genuine third-order Moebius term would be
what crosses the axis. On IEEE-9 it is not: the pair term BC contributes
`+0.018913` while the irreducible triple term is `-0.000806`, i.e. an order of
magnitude smaller and of the opposite sign. A pairwise-complete model already
predicts instability (`+0.004023`).

The `kappa = 3` structure is real; the mechanism behind it is pairwise
accumulation, not an irreducible triple. Recorded in `RESULTS.md` R3 and emitted
by `E01` as `mechanism_verdict`.

**Consequence for the poster:** "every single safe, every pair safe, the
portfolio unstable" is supported on IEEE-9. "Because of a genuine three-way
interaction" is **not**, and must not be written.

## F2 — Uniform contour sampling silently returns winding zero

First IEEE-9 winding run reported `0` for every portfolio including the unstable
one, with a clean integer residual and no warning. Cause: the contour passed
`0.005` from the created mode while sampling at `2.0`. Phase-driven adaptive
refinement did not fire because the excursion cancels between samples.

Fixed by tying the step to the known clearance (`max_segment`). Any winding
reported before this fix is void.

## F3 — No pole/residue analysis of Sigma for the natural controller partition

The hidden block of the IEEE-9 GFL is defective (`cond(V) = 1.2e19`, a fourfold
eigenvalue at the origin from cascaded integrators). Specification section 1.1
asks for a per-controller-pole residue table; that object does not exist for this
partition. Not worked around, not hidden: `expand_self_energy` returns
`DEFECTIVE_HIDDEN_BLOCK`.

## F4 — The specification eigenvalue envelope could not be reproduced exactly

The four reference scripts it cites are absent from this repository and the
textual description fixes neither the dq sign convention nor the data to six
decimals. An independent implementation from the stated equations agrees on all
eight portfolios to `<= 5.3e-4`, including the sign change at the triple, which
is reported as corroboration and not as reproduction.

## F5 — First IEEE-9 implementation was unstable at base

A real pole at `+7.63 rad/s` appeared in the base case. Cause: with the d axis
aligned to the terminal voltage, injected reactive power is `Q = v_q i_d - v_d i_q`,
so `dQ/di_q = -v_d < 0` and the q-axis PI acts as positive feedback. The
equilibrium was consistent either way, so only the linearization exposed it. The
inversion is now explicit and commented in `ieee9_gfl.py`.

## F6 — The toy orientation coincidence does not generalize

`tr(H_ABC) = tr(H_ACB)` holds to `1e-17` on the toy while `K` is demonstrably
non-symmetric there. No theorem forces it. Recorded as a case property in
`test_toy_orientations_coincide` and explicitly not generalized.

## F7 — The argmax-MAC single-branch tracker fails once the branch splits (E28)

Expected `alpha_Q` to move continuously through zero as the operating point
crosses the stability boundary. It does not. Under replacement the single base
inter-area branch at 0.644 Hz appears as **two** branches of indistinguishable
shape, 0.566 Hz and 0.666 Hz, at MAC 0.972 and 0.982 to the same base mode, in
293 of 293 accepted samples. Their MAC values cross near load 0.988 within about
0.001 of each other, so the selector flips between them: at loads 0.9800, 0.9825,
0.9850 and 0.9875 it reports `+0.035`, `-0.654`, `+0.082`, `-0.700` while the
worst mode in the band rises smoothly through `+0.035`, `+0.057`, `+0.082`,
`+0.108`.

Consequences, measured: the tracked observable disagrees in **sign** with the
band envelope in **43 of 293 samples**, under-reports instability (96 unstable
against 139), and places the stability boundary optimistically high by 0.002 to
0.022 in load, the error growing with availability.

**Consequence for the poster:** the v2B `alpha_Q = 0` and `delta_alpha_Q = 0`
curves may not be quoted as physical boundaries. The interaction coefficient and
the conditional repair result are unaffected. Recorded as N11 in `CLAIMS.md`.

## F8 — H5B was never executed (E28)

The preregistered near-boundary stratum, `abs(alpha_Q) <= 0.05`, holds **0 of 293
samples**, because of F7: the smallest `abs(alpha_Q)` reachable by the chattering
selector is `0.0980`. The Mann-Whitney statistic is undefined and the frozen
protocol returns FAIL.

This is **not** evidence against the port-closure mechanism. Unlike v2A's H5, the
observable is not degenerate: `m_4` ranges over `0.0067` to `0.3570` at condition
number below `9e+03`. The stratum was empty, not the measurement.

A post-hoc recomputation on the band envelope is recorded in `V2B_RESULTS.md`
and is **descriptive only**. It passes nothing. A hypothesis tested on an
observable chosen after seeing the data is not tested. Recorded as N11b.

## F9 — Forward single-branch tracking gives a FALSE NEGATIVE on the flagship (E31)

The v2C order-4 audit re-ran all 16 subsets with the argmax selector anchored on
the base inter-area mode. At the full flagship it returns `-0.734`, i.e. it
declares a portfolio stable that has a right-half-plane mode in the band. It
follows the well-damped 0.66 Hz descendant while the 0.575 Hz one crosses.

E15 did not hit this because it anchored BACKWARDS, on the flagship's own
critical mode, which by construction cannot mis-select at the flagship. That is
why the defect survived until v2B: the original audit could not have exposed it.

The family envelope returns `+0.14467` at the flagship, agreeing with E15, and is
well defined at every subset in either direction. Recorded as N12 in `CLAIMS.md`.

**Consequence:** no result on this system may use a forward argmax single-branch
selector as a stability endpoint. The order-4 claim itself survives, O45.

## F10 — The naive fixed-structure Q/V coordinate is singular at zero gain (F7 Safeguard A)

Scaling the PI voltage-regulator gains to zero leaves one decoupled marginal
integrator per converter: four exact zero eigenvalues on the flagship and a
continuum of equilibria at `g = 0`, and slow real modes collapsing onto the
origin as `g -> 0`. Zero gain cannot sit inside `Theta_reg` in that coordinate.
Replaced by the leaky regulator (`voltage_leak`), regular on `[0, inf)`.
`results/F7/F7_safeguard_A_audit.json`, CLAIMS O64.

## F11 — Corner-label quadtree refinement was abandoned before completion (F7)

The first F7 run refined a cell only when its four corner labels differed. The
exploration had just shown a narrow instability region near 0.68 Hz, which a
corner-only rule can miss entirely. The run was stopped before any map finished
and replaced by the nine-point rule with near-axis monitoring (Safeguard C). No
number from the stopped run exists or was used.

## F12 — A gradient-cosine tangency test depended on axis scaling (F7)

The first boundary classifier called a crossing a tangency when the path
direction made a small angle with the in-plane gradient of `Re lambda`, in
range-normalised coordinates. It labelled the F2 crossing at `g = 0` along `k` a
tangency only because the policy gradient dwarfs the excitation gradient in
those units. Replaced before the maps ran by a scale-free rule: a tangency is a
crossing whose local parabola along the path has its stationary point within
one edge length.

## F13 — Four-connectivity overcounted region components (F7)

Counting components with four-connectivity split thin diagonal strips (a
`kappa = 3` strip at `g ~ 0.36`, `k ~ 1.77`) into single-node components.
Components are counted with eight-connectivity and only those of at least 25
lattice nodes are called substantial.

## F14 — Regulator-leak sensitivity of boundary locations (F7, a limitation, not a bug)

Label agreement with the frozen leak `w = 0.05` is 89 % at `w = 0.02` and 72 % at
`w = 0.2`, concentrated at points adjacent to a boundary. Boundary positions and
area fractions are therefore conditional on `w`; every qualitative F7 result was
re-traced at all three leaks and survives (`results/F7/F7_leak_qualitative.txt`).

## Not attempted

Everything from IEEE-39 onward: the flagship model, the twelve-action census,
the Monte Carlo campaigns, negative controls, the detuning experiment, the
nonlinear time-domain check and the repair comparison. Gates 4, 6, 7 and 8 are
open, and no poster wording may depend on them.
