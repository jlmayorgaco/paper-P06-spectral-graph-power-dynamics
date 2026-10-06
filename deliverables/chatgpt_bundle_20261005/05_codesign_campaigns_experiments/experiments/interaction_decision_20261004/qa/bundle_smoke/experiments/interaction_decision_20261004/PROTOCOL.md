# Frozen discovery and validation contract — 4 October 2026

Problem: determine whether simultaneous SG/GFL interventions violate a spectral
margin although a sum of exact single-site effects predicts that same margin.
Gap: previous identities and tighter bounds have not changed an IEEE-39 decision.
Insight: remove each site's finite self-interaction exactly before bounding the
remaining cross-site closed walks. This is a specialization of classical
determinant/interaction theory, not a claim to invent decentralized control.
Falsifier: no material decision reversal, invalid contour assumptions, or a
remainder larger than the decision margin defeats the proposed certificate.
Consequence sought: a signed interaction correction with a checkable error budget
and an actionable finite replacement/tuning rule.

## Immutable inputs and semantics
Use the exported fixed-equilibrium IEEE-39 physical model in
experiments/graph_gsp_codesign_20261003/model and its baseline.toml.
The preserved sharing convention, P0 dispatch, controller channels and bounds
apply unchanged. Delays are physical, fixed at 40 ms during each design search.
Small-signal requirement is Re(lambda) <= -0.05 /s, excluding gauge.
Single-site intervention means changing one site's rho/Kp/Ki in the FULL network;
it does not mean disconnecting the physical network or deleting model coupling.
The additive predictor is a counterfactual surrogate, not a competing published
method and not itself a realizable independent system.

## Discovery (explicitly exploratory; save every case and failure)
Family R: common rho = .875 + .005*k, k=0,...,17 (through .960), original gains.
Family P: rho=.875, common Kp factor 1+.02*k, k=0,...,15 (through 1.30), Ki fixed.
Family I: rho=.875, Kp fixed, common Ki factor 1+.02*k, k=0,...,15.
Track all 11 pre-existing physical positive-imaginary roots from the baseline
catalog, with intermediate steps <= .005 rho / .02 relative gain; do not label
this seed catalog a complete DDE spectrum. Compare each joint tracked root to
z_N=z_0+sum_i(z_i-z_0), using the corresponding single-site branches.
Report unresolved tracking and mode collisions, never relabel a jumped root.
Reject a branch comparison if a per-step root jump exceeds 2 /s or roots merge
within 1e-5 /s. Candidate screening is numerical, not interval-certified.

## Selection and closure
Select the earliest increasing-coordinate case per family with nodal predicted
margin <= -1e-4 /s and joint violation >= 1e-4 /s. If none, report negative;
do not alter thresholds/grid retroactively. A bracketed transition may be
refined by bisection; refinement is explicitly adaptive and recorded.
For the strongest selected gap, derive the block-resummed contour correction
and attempt a continuous interval exclusion proof on a contour containing the
relevant cluster. A proof of a root to the right of -0.05 rejects the candidate;
it is not a proof about all allowable gains or global replacement capacity.
Follow-up designs/campaigns require a separately timestamped addendum before
their evaluations; discovery examples are not statistical generalization.

## Comparison and claim gates
Compare first-order sensitivity, exact singleton-additive predictor, ordinary
quadratic expansion, and finite block-resummed correction under the same model.
No global optimum, safe nonlinear replacement, hardware validation, physical
line-cycle law, necessity of communication, or new fundamental law follows.
Nonlinear trajectories are deferred until a meaningful spectral decision exists.
All new evidence is written here. Historical frozen files remain unchanged.
