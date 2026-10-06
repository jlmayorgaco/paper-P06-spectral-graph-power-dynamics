# POSTER RESULTS

What the gate phase actually supports, with conservative wording. Nothing below
depends on an experiment that has not run.

## Result 1 — A portfolio can fail while every part of it is safe

**Numbers (IEEE-9, one frozen AC operating point, three stress actions).**

| portfolio | Re(lambda_crit) | zeta | verdict |
|---|---|---|---|
| base | `-0.105890` | `0.0123` | stable |
| A (SG damping to zero) | `-0.027533` | `0.0032` | stable |
| B (`kp_P` 0.20 to 2.00) | `-0.086124` | `0.0101` | stable |
| C (PLL `wn=12`, `zeta=0.15`) | `-0.112056` | `0.0131` | stable |
| A+B | `-0.008014` | `0.0009` | stable |
| A+C | `-0.034409` | `0.0040` | stable |
| B+C | `-0.073377` | `0.0086` | stable |
| **A+B+C** | **`+0.003217`** | `-0.0004` | **unstable, 1.3577 Hz** |

**Safe wording.** "At the frozen action amplitudes and operating point of this
case, every single intervention and every pair leaves the system stable while the
full portfolio does not. Minimum destabilizing order is three."

**Do not write.** "Pairwise screening is always insufficient." One case, one
amplitude vector.

**Figure.** F02 bar chart of `Re(lambda_crit)` per portfolio; the source data is
`results/tables/E01_regression_ieee9_portfolios.csv`.

## Result 2 — The created mode belongs to the collective interaction factor

The portfolio determinant factors exactly as
`det(I+M) = [prod_a det(I+M_aa)] det(I+Q)`. At the created eigenvalue:
`|full| = 1.1e-13`, `|collective| = 2.2e-12`, `|individual| = 0.0489`,
`sigma_min(I+Q) = 2.4e-16`, `rho(Q) = 1.0000`.

Winding on circles of radius 0.005, 0.010 and 0.020 around the created mode:
full `+1`, individual `0`, collective `+1`.

**Safe wording.** "On a region that isolates the created mode from every base and
lower-order mode, the collective interaction factor carries it, and no isolated
action does."

**Must accompany it.** "The statement is conditional on that region: a larger
region that also encloses a single-action mode moves the count to the individual
factor. The framework refuses such a region rather than reporting a label."

**Figure.** F09 winding decomposition with the admissibility band drawn, plus a
deliberately inadmissible radius shown as refused.

## Result 3 — Non-additive does not mean irreducibly three-way

Moebius decomposition of the spectral abscissa:

| | order-1 truncation | order-2 truncation | exact | irreducible 3rd |
|---|---|---|---|---|
| five-state case | `-0.064298` stable | `-0.023223` stable | `+0.046781` unstable | `+0.070005` |
| IEEE-9 case | `-0.013933` stable | `+0.004023` unstable | `+0.003217` unstable | `-0.000806` |

**Safe wording.** "Two cases with the same minimum destabilizing order reach it
by different mechanisms. In the constructed case the irreducible third-order term
is what crosses the axis. In the converter case it is an order of magnitude
smaller and stabilizing, and the crossing is driven by accumulated pairwise
interaction. Reporting `kappa = 3` without this decomposition would conflate the
two."

This is the strongest methodological result of the phase and it should be on the
poster, not hidden. It is what distinguishes the work from a screening heuristic.

**Figure.** F10 interaction magnitude by order, both cases side by side.

## Result 4 — The identities are exact, so the decomposition is not a fit

Schur `1.8e-14`, determinant lemma `6.8e-15`, individual-times-collective
`8.9e-16`, action numerical rank exactly 1 with reconstruction error `5.5e-16`,
reduced-matrix additivity `2.3e-17`, operating-point drift `9.1e-14`.

**Safe wording.** "Every step from the reduced model to the interaction factor is
an exact algebraic identity verified to machine precision; nothing in the
decomposition is fitted or approximated."

## Result 5 — Two numerical requirements without which the method reports nonsense

1. A uniform contour raster returns a clean integer winding of zero for a mode
   at `zeta = 0.012`. The sample step must be tied to the distance from the
   contour to the nearest mode.
2. The individual/collective split is only interpretable on a region that
   excludes every base and lower-order mode.

**Safe wording.** "Two conditions govern whether a provenance number means
anything: the contour step must resolve the clearance, and the region must
isolate the created mode. Both are enforced and reported, not assumed."

## What is NOT ready for the poster

- Any IEEE-39 result. Not built.
- Any cycle or SCC causality claim. The IEEE-9 feedback core is the trivial
  three-action component and the strongest cycle score is `1.6e-3`.
- Any repair claim beyond the five-state case, and no comparison against asset
  removal or global retuning.
- Any robustness, Monte Carlo or negative-control statement.
- Any reuse of TX3 numbers. The frozen TX3 envelope rejected materiality on the
  margin-setting mode of the IEEE-39 benchmark; see `CLAIMS.md`.
