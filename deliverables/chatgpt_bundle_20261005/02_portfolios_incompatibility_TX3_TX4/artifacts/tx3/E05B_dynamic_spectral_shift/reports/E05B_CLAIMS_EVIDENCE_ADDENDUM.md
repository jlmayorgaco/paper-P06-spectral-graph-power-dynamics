# E05B claims--evidence addendum

The accepted E05 claims matrix is retained byte-for-byte because it belongs to the frozen E05 review
package. These new rows are an E05B addendum pending independent review.

| Claim | Exact statement under test | Required evidence | Preregistered tolerance | Experiment | Result | Pass/fail | Permitted paper wording |
|---|---|---|---|---|---|---|---|
| C3b | Full descriptor externality separates exactly into algebraic, mass-scaling, and finite-pole externalities. | Schur-factor audit against frozen full E03 logdets; all pair/triple components at 24 development plus 16 independent holdout points. | Pointwise and integrated residual at most 1e-7; component closure at most 1e-10. | E05B-D1 | Maximum residuals are below 7.4e-13. | PASS | State the exact finite-DAE factor separation; do not infer modal materiality from the full descriptor term. |
| C3c | Finite-pole externalities occur reproducibly on an independent operating-point ensemble. | All 28 pairs and 56 triples at 16 new Sobol points with conservative uncertainty. | At least two pairs and two triples resolved on at least 8/16 points with same-sign fraction at least 0.75. | E05B-D2 | 10 pairs and 5 triples pass. | PASS | Distinguish genuine finite-pole non-additivity from algebraic re-equilibrium interaction. |
| C3d | Connected contour moments localize finite-pole motion for frozen pair and triple candidates. | Frozen OP00 pole families; direct Xi traces; exact residues; independent 128-node contour audit. | At least one pair and one triple valid/resolved on at least 8/16 holdout points with same-sign fraction at least 0.75. | E05B-D3 | A7-A8 and both triples pass; A3-A6 fails its frozen numerical contour rule. | PASS | Claim existential localized pole motion only; no cycle, causality, or material damping-margin claim. |
