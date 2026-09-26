# Final rejected wording (IAS poster and TPWRS manuscript)

Every item below is **prohibited**. Each is paired with the reason and, where one
exists, the permitted replacement. The earlier forbidden lists
(`POSTER_FINAL_SAFE_CLAIMS.md`, v2, `G5_NOVELTY_STATEMENT.md`) remain in force;
this list supersedes them where the two differ.

## Stability language (Decision 1)

| rejected | why | permitted |
|---|---|---|
| "the IEEE-39 system is stable" / "Hurwitz" / "globally small-signal stable" | D = 0 and no governor leave an exact neutral common-frequency mode (a Jordan chain with the rotation) | "transversely (relative) stable: alpha(A_perp) < 0" |
| "whole-RHP" for IEEE-39 without qualifier | the rotation and drift eigenvalues sit at 0 | "transverse whole-RHP" |
| "eigenvalues with abs(lambda) <= 1e-3 are ignored" | a physical eigenvalue may be arbitrarily close to 0 | "quotient by the exact center subspace C = span{R_x, w}" |
| "the two zero modes are gauge symmetries" | only the rotation is a gauge; the drift is a missing-restoration property | the separate statements of `TRANSVERSE_STABILITY_QUOTIENT.md` §2 |

The model-scope sentence is mandatory: *"The frozen IEEE-39 model has no primary
frequency restoration and is analysed in relative / transverse coordinates.
Absolute common-frequency restoration is outside this benchmark."*

## Units (Decision J)

| rejected | permitted |
|---|---|
| "4270.7 MW", "4.27 GW of PV", "keeping all 4270.7 MW of PV" | "2096.6 MW of active dispatch on 4270.7 MVA of converter rating" |
| "restoring SG30 costs 1040 MW of PV" | "gives up 436.1 MW of converter active dispatch (1040 MVA machine kept)" |
| "restoring SG37 costs 970.2 MW" | "gives up 321.5 MW (970.2 MVA machine kept)" |
| "X MW of PV capacity" (any X) | active dispatch [MW] or rating [MVA]; no PV nameplate exists |
| "replaced megawatts rank instability" / "MW and inertia are informative predictors" | "aggregate removed active dispatch does not identify the minimum incompatible coalition at the first failing order; rating / inertia may rank risk but do not recover the policy-dependent witness structure" |
| "matched on MW" (E14, E18 controls) | "matched on active dispatch Pg (UC02 rerun)" |
| the frozen `E39_ROC_PR.png` | `results/UC/UC02/UC02_ROC_PR_unit_corrected.png` |

## Novelty

The following are never claimed as new:

- generalized Nyquist, the return difference, the −1 crossing, impedance / port
  criteria;
- Schur complements, "self-energy" terminology;
- hypergraphs, clutters, minimal failure sets (N-k contingency analysis
  already uses minimal defective k-sets on hypergraphs);
- piecewise-constant integer counts;
- Taylor second-order accuracy;
- ordinary Lyapunov stability, sum-of-squares, small gain, Perron–Frobenius
  monotonicity;
- ordinary (persistent) homology, Stanley–Reisner / minimal non-faces;
- the NP-hardness of cardinality-constrained spectral destabilization. It is
  the sparse-PCA clique reduction (Magdon-Ismail 2017), and robust stability is
  NP-hard (Poljak–Rohn 1993; Nemirovskii 1993).

## Mechanism and scope

| rejected | why |
|---|---|
| "Track A is caused by converter control" | E14 N5/N1: not controller-caused, but controller-repairable |
| "cycle holonomy explains the failure" | falsified (E18; the UC02 rerun confirms, ratio 0.72) |
| "a fixed converter retune repairs the failure" as a general claim | bounded authority (E35: 117/240) |
| "the undamped condenser restores composability" | it creates its own unstable swing mode (G1) |
| "bus 30 is the collective enabler" | not supported after adjustment (E38; UC02 rerun p = 0.18) |
| "the incompatibility persists with governors" or "is robust to frequency control" | FC03: with the documented TGOV1N governors the P4 flagship is transversely stable; the kappa = 4 coalition is conditional on the absence of primary frequency control |
| "independently validated by ANDES" for anything beyond the base inter-area mode and the power flow | ANDES does not implement the L0 GFL; the documented IEEEX1 is itself unstable in ANDES (six real modes near +1.03) |
| "EMT", "electromagnetic transient", "hardware", "field" | the TDS is a nonlinear phasor-domain DAE |
| any nonlinear composability claim from a disturbance beyond the model-scope guards | FC05: the finite-disturbance thresholds are scope-censored (see FINAL_TRANSACTION_THEORY_AND_EVIDENCE, §E–G) |
| "the relocated port was validated on the 194 Kundur crossings" | that was method development; the validation is the holdout (28/29, the miss explained) |
| "certificate" for any frequency-grid result | SCREENING only (FC09) |
| "the F7 labels are certified in the g <= 0.005 strip" | only after the direct path; 16 points there remain unresolved |
| a Hessian or second-order claim on IEEE-68 | NL02 accuracy there is only about 1e-3 |
| "the structure is universal" | IEEE-68: the preregistered 4-candidate map is empty |

## Final-campaign additions (after FC03–FC13)

| rejected | why | permitted |
|---|---|---|
| "nonlinear composability is stricter than small-signal composability" / "kappa_NL < kappa_RHP" | FC05/06: inside model scope `H_NL = H_RHP_perp` in 6/6 cases | "within the declared model validity, finite-disturbance composability coincides with transverse spectral composability" |
| "the flagship fails at 130 MW" / any finite `r_S` quoted as a stability margin | 84/96 thresholds are OUTSIDE_MODEL_SCOPE guard events (GFL voltage band, PSS limit) | "first reaches the declared converter voltage envelope at …" |
| "the first-event complex shows nonlinear incompatibility" | it records limiter activation of limits the model omits (rule 24) | "order in which portfolios leave the declared envelope (scope diagnostic)" |
| "684–792 MVA of condenser are needed for stability" | the requirement is to stay inside the guards, not stability | "… to keep a 200 MW disturbance within the declared converter voltage envelope, in this model" |
| "nonlinear planning changes the optimal portfolio" | P2 is unchanged where decidable; out of scope at P4 | none |
| "the incompatibility persists with primary frequency control" (unqualified) | the P4 coalition is stabilized by the documented governors | "policy-dependent incompatibility persists in a smaller region; the P4 four-bus coalition is stabilized" |
| "Hopf bifurcation causes the finite-disturbance failure" / "subcritical everywhere" | l1 was computed at three boundaries; finite thresholds are scope-censored | "the three examined boundaries are subcritical" |
| "we prove NP-hardness of portfolio destabilization" as a contribution | a corollary of the sparse-PCA clique reduction | "a known consequence of …" |
| "topological obstruction" / "new resilience homology" | Stanley–Reisner / filtration of the definition; local spheres are tautological | "the tolerated portfolios form a simplicial complex whose minimal non-faces are the minimal coalitions" |
| "certified" for BC03, SG/Perron, top-sum, Lyapunov or energy bounds | sampled screens (FC09); Lyapunov ratio ≤ 1e-11 (FC08); energy NOT_APPLICABLE | none (report as failed) |
| "the second-order response requires the AC curvature" as new | Taylor consistency is classical; see FC07 for the measured magnitude | the FC07 wording in `FINAL_TRANSACTION_THEORY_AND_EVIDENCE.md` B.6 only |
| "validated by ANDES" for governed-model portfolio results | ANDES reproduces only the governor damping direction on the base case | "consistent in direction with an independent ANDES model of the same governors" |
