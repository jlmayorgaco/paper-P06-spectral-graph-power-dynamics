# EXPERIMENTS

One row per experiment: question, hypothesis, variables, controls, result,
interpretation, limitation. Experiments not yet run are listed with their design
so that the design precedes the data.

## E00 — five-state regression (GATE 0) — DONE

- **Question.** Does the framework reproduce the frozen five-state envelope and
  do its identities hold to machine precision?
- **Hypothesis.** They do; any deviation is an implementation error.
- **Independent variable.** the eight subsets of three actions.
- **Dependent.** spectral abscissa, identity residuals, winding split, Moebius
  orders, repair endpoints.
- **Controls.** direct determinant reference for every lemma; random basis
  changes for gauge invariance; four contour radii.
- **Sample size.** deterministic, 8 portfolios.
- **Result.** PASSED. Max abscissa deviation `4.2e-07`; all identities at
  roundoff; `kappa = 3`; provenance COLLECTIVE on admissible regions;
  `k_crit = 1.351826`, `k_repair = 1.259557`.
- **Interpretation.** The calculus is implemented correctly.
- **Limitation.** A constructed model with no algebraic variables and no power
  flow. It cannot corroborate anything physical.

## E01 — IEEE-9 GFL regression (GATE 1) — DONE

- **Question.** On a converter-dominated case with a real AC operating point,
  does the `kappa = 3` structure appear, and which determinant factor carries the
  created mode?
- **Hypothesis.** The collective factor carries it.
- **Independent variable.** the eight subsets of A (SG damping to zero),
  B (`kp_P` 0.20 to 2.00), C (PLL to `wn = 12`, `zeta = 0.15`).
- **Dependent.** spectral abscissa, frequency, damping, winding split, Moebius
  orders, action graph and cycle scores.
- **Controls.** equilibrium residual certified per portfolio; power-flow
  mismatch and voltage range; additivity of the reduction measured; operating
  point drift measured; every identity against a direct determinant; contour
  admissibility; clearance-tied winding step.
- **Sample size.** deterministic, 8 portfolios, one operating point.
- **Result.** PASSED at the frozen tolerance and corroborating the outside
  envelope to `5.3e-4`. `kappa = 3`. Collective provenance on three admissible
  radii. **Mechanism verdict: PAIRWISE_ACCUMULATION**, not genuine third order.
- **Interpretation.** Non-additivity is real and the created mode belongs to the
  collective factor, but on this case the third-order term is small and
  stabilizing; pairwise interaction is what destabilizes.
- **Limitation.** One operating point, one amplitude vector, three actions, no
  uncertainty, no negative control, no time-domain check.

## Designed, not run

| id | question | why it is blocked |
|---|---|---|
| E02 | do two independent Jacobian methods agree on the IEEE-9 model? | complex-step needs analytic subfunctions; central differences are the only estimate today |
| E10 | IEEE-39 baseline on both the frozen ANDES case and a fresh frozen-point model | needs the ANDES environment and the model decision of METHODS.md §10 executed both ways |
| E12 | exhaustive 2^12 portfolio census at a nominal operating point | needs E10 |
| E13 | map of `kappa` over amplitudes; does `kappa` move continuously? | needs E12. Required before any `kappa` claim, because `kappa` is amplitude dependent |
| E14 | negative controls: fast damped PLL, GFM instead of GFL, strong damping, electrically strong connection, weak amplitudes | needs E10. **Required before any poster wording**: without it the method is an "everything is dangerous" detector |
| E15 | does a repair confined to the feedback core beat asset removal and global retuning at equal change norm? | needs a localized mechanism, which E01 has not established |
| E16 | operating-point Monte Carlo, controllers frozen | needs a frozen flagship |
| E17 | controller-parameter Monte Carlo, operating point frozen | as above; run separately from E16 so the two sources can be told apart |
| E18 | inertia and damping retirement sweep | needs E10 |
| E19 | detuning: does `abs(H_p)` scale with the distance from a network mode to a controller pole? | blocked by R7: the natural partition has no simple-pole expansion |
| E20 | held-out frozen validation | last |
