# TX4 Blind Portfolio Prediction + Contextual Return

## 1. Executive verdict

The frozen TX4 IEEE-39 model contains a clean finite witness: H4 =
{30,33,35,37} is unstable at P4 while all 15 proper subsets are stable. One
all-SG network kernel plus four local SG-to-GFL substitutions predicted this
V4 verdict before full-order reveal: 16/16 verdicts, zero false-safe cases,
exact H, and exact kappa=4. The contextual return reaches unity at the same
g-only boundary as the full DAE, and the nonlinear phasor-domain TDS preserves
the stable/unstable/stable ordering.

The result is scoped. The V9 census is not an exact scalable predictor:
395/512 verdicts, 116 false-safe cases, and incorrect H recovery, although
kappa=4 is recovered. The reduced closure is slower than the direct
full-order benchmark. The g-boundary run was reduced-only but executed after
the reveal, so it is a reduced reproduction rather than a fully blind
boundary claim. Final decision: CASE B, a good IAS poster with a limited
predictive claim.

## 2. Scientific question

Can one reusable network kernel and local device substitutions predict a
minimal incompatible SG-to-GFL portfolio, and can an exact contextual return
explain the failure? The answer is yes for the preregistered four-bus
flagship at P4, but not yet for the nine-bus transfer census.

## 3. Prior evidence and negative lessons

The frozen TX4 record already established the full-order H4 witness, a
controller-dependent g boundary, nonlinear phasor TDS behavior, and limited
independent ANDES phasor comparisons. This campaign did not modify that
record or any prior negative PowerDynamics/CDW line. V9 and runtime results
are retained as negative evidence against unqualified scalability and
speed claims.

## 4. Frozen TX4 model

The model is the IEEE-39 phasor-domain DAE with the frozen synchronous-machine
and AVR semantics and the ten-state custom GFL replacement. P4 is
(g,k,t,h)=(0.03625,1.425,1.5,1). The critical band is 0.3--1.5 Hz, positive
frequency. Equilibria are resolved for each full-order case. No new
controller family, governor, limit, EMT interface, or network model was
introduced.

## 5. Blind-prediction firewall

V4={30,33,35,37}; V9={30,31,32,33,34,35,36,37,38}. Before reveal, the
predictor imported only model builders, the all-SG baseline, local replacement
models, controller parameters, and port construction. It did not import
FC01, PCV02, FC10, H, kappa, or historical boundary tables. It constructed
one T0/K(s) and one local Delta-Y model per candidate. Blind files were
hashed and committed before reveal. The firewall is documented in
docs/TX4_BLIND_FIREWALL.md.

## 6. Baseline network kernel

The kernel is K(s)=E^T T0(s)^(-1) E. For portfolio S, D_S contains local
device changes, M_S=D_S K_SS, and Q_S=(I+D_S K_d)^(-1)D_S K_o. The audited
port convention gives the exact determinant factorization under its regularity
assumptions. It is not claimed to reproduce the full spectral abscissa
uniformly away from a root.

## 7. Local SG-to-GFL replacement models

Each local model is equilibrium-matched at its candidate bus using the frozen
GFL equations and P4 parameters. The same local library is reused across V4
and V9; no portfolio coefficient is fitted after reveal.

## 8. Exact reduced portfolio closure

The closure search uses the fixed frequency band, deterministic grids, fixed
seeds, root refinement, and local-factor floor. An accepted right-half-plane
closure root yields UNSTABLE_PREDICTED; otherwise the result is
STABLE_PREDICTED or unresolved. Full-order spectra are reveal-stage ground
truth.

## 9. Blind V4 prediction

The predictor returned 15 stable predictions and one unstable prediction,
exactly H4. The full-order reveal returned 15 stable and one unstable cases:
accuracy 16/16, false-safe 0, false-unstable 0, exact H yes, exact kappa 4.

For H4, the reduced root was Re(s)=0.12700633761 s^-1 at
0.62227967095 Hz. Full order gave alpha=+0.12700646714 s^-1 at
0.62227966958 Hz. Root-real error was -1.30e-7 s^-1 and frequency error was
1.37e-9 Hz.

## 10. Full-order reveal

The answer tables were loaded only after the blind commit. The PCV02 core
table omitted BASE; comparison code added BASE from the frozen FC01 reference
without changing predictions. V9 truth contains 327 stable and 185 unstable
portfolios. The reduced predictor classified 395/512 correctly, with 116
false-safe and one false-unstable case. Its H family differs from the
full-order 42-edge family, although minimum blocker cardinality is correctly
4.

## 11. Clean controller-boundary reproduction

The reduced-only g search gave g_return=0.20768138607, frequency
0.706424783 Hz, and collective sigma_min=1.38e-8. The frozen full-order
boundary is g_full=0.20768140519; the difference is 1.91e-8, below the
preregistered 1e-4 tolerance. At g plus or minus 1e-4, collective sigma_min
was about 6.02e-5 and local sigma_min about 0.489.

This search ran after reveal. The solver itself did not import the answer key,
but the output is labeled REDUCED_REPRODUCTION_AFTER_REVEAL. The strict
fully-blind C05 claim is therefore limited.

## 12. Contextual-return identity

For H and R=H minus {i}, define C_H=I+Q_H, G_R=(I+Q_RR)^(-1), and
R_i|R=Q_iR G_R Q_Ri. Classical block algebra gives:

det(I+Q_H) = det(I+Q_RR) det(I-R_i|R).

This is standard Schur/block determinant algebra, not a new theorem. At the
boundary, the nearest return eigenvalue in the archived convention was
-0.9999999649 - 5.1e-9 i; the largest Schur residual was 3.44e-16.

## 13. Numerical return validation

The nine-case audit used BASE, representative singleton/proper cases, and H4
with finite-difference factors 0.5, 1, and 2: 27 rows. NumPy and SciPy
eigensolvers agreed to displayed precision. The maximum equilibrium algebraic
residual in the mode audit was 1.30e-12. A valid E,A descriptor pair was not
available from the frozen DAE, so that check is BLOCKED/NOT_AVAILABLE.

The contextual-return derivative at a=g was checked analytically and by
central finite difference; maximum error was 7.20e-6 in the recorded
realization. This is an implementation check, not new eigenvalue-sensitivity
mathematics.

## 14. Minimality separation

| retained GFL buses | missing bus | alpha s^-1 | frequency Hz | verdict |
|---|---:|---:|---:|---|
| 30+33+35 | 37 | -0.2046876902 | 0.916019 | STABLE |
| 30+33+37 | 35 | -0.1755083817 | 0.918500 | STABLE |
| 30+35+37 | 33 | -0.1678557851 | 0.914639 | STABLE |
| 33+35+37 | 30 | -0.2030092264 | 0.917627 | STABLE |
| 30+33+35+37 | -- | +0.1270064671 | 0.622280 | UNSTABLE |

All 15 proper subsets are stable at P4. H4 is a minimal dynamically unstable
portfolio in this model and policy. It is not a universal SG-to-GFL claim.

## 15. Local versus collective failure

At the contextual boundary, proper-subset separation eta_H was 0.120831 in the
archived realization, while collective sigma_min was 2.24e-8. The reduced
boundary run showed local sigma_min 0.489 and local determinant minimum 0.433.
The accepted interpretation is collective closure failure with regular local
factors, not a single-device local singularity.

## 16. Fifteen proper subsets

results/TX4_15_PROPER_SUBSETS.csv contains BASE, all singletons, all pairs,
all triples, and H4. Every proper row has stable full-order alpha at P4, and
the closure separation columns remain nonzero at the H4 root frequency.

## 17. V9 512-portfolio benchmark

V9 is a scalability stress test, not a success claim. The same kernel and
local library classified all 512 subsets. Metrics are accuracy 0.771484375,
false-safe 116, false-unstable 1, H recovery false, and kappa recovery true.
This limits poster language to flagship prediction rather than general exact
portfolio prediction.

## 18. Runtime

| set | full total s | reduced total s | full/reduced |
|---|---:|---:|---:|
| V4, 16 cases | 2.397 | 4.360 | 0.550 |
| V9, 512 cases | 72.490 | 309.649 | 0.234 |

The reduced approach is slower here. Its value is mechanistic decomposition and
finite flagship prediction, not speed.

## 19. Aggregate MW/MVA counterexample

The fixed V9 rule found a same-cardinality pair with 0.2513% normalized
aggregate mismatch:

| portfolio | cardinality | MW | MVA | alpha s^-1 | status |
|---|---:|---:|---:|---:|---|
| 30+31+32+34+36+38 | 6 | 3659.870 | 6509.2 | +324.249 | UNSTABLE |
| 31+32+34+35+37+38 | 6 | 3652.305 | 6499.9 | -0.1664 | STABLE |

There is no exact MW/MVA match under the stated tolerance. This is a
near-match supporting the limited statement that aggregate penetration and
cardinality do not fully determine compatibility.

## 20. Lower-order screens

Singleton, pair, and triple screens all pass while H4 fails. H4 values are
B3 first-order -0.203249, B4 additive -0.212867, B5 pairwise -0.210389, and
B5b third-order alpha -0.183673, all predicting stable H4. Exact reduced
closure predicts unstable and agrees with full DAE. Third-order alpha
approximation is kept distinct from third-order determinant truncation.

Correct statement: lower-order safe screens do not imply portfolio
composability. It is not that conventional full eigenanalysis fails.

## 21. Controller remediation

Only g was changed. H4 has alpha=+0.1270064671 s^-1 at P4 and
alpha=-0.0173846763 s^-1 at fixed g=0.25. The derivative, direction, and
root-iteration record are in results/TX4_CONTROL_REMEDIATION.csv. This is a
controller-only removal of the identified blocker, not an economic optimum.

## 22. Nonlinear TDS

The frozen G2 disturbance is a +2% active-load pulse at bus 20 for 0.2 s,
using the archived nonlinear phasor-domain BDF solver and common observable.
The three declared runs were proper 30+33+35, H4 original, and H4 at g=0.25.
Outcomes were STABLE, UNSTABLE, and STABLE. All growth signs agreed with
linear alpha signs. Primary exponent errors were below 3.51e-4 s^-1. This
is not EMT.

## 23. ANDES

The existing E31 independent ANDES comparison was reused without retuning. It
supports limited phasor sign comparison for base and selected portfolios, but
its SG/IEEEX1/IEEST/network assembly is not the custom TX4 GFL model. The
parameterized custom-GFL g boundary is NOT_TESTED and cannot be advertised as
independent custom-GFL validation.

## 24. EMT status

No EMT portfolio claim was run. The ParaEMT path remains unresolved and is not
needed for the narrow phasor-domain poster story.

## 25. Literature novelty

The focused 2020--2026 audit found substantial prior work on grid strength,
multi-inverter return ratios, impedance stability, controller stability
regions, and stability-constrained IBR planning. It did not identify in the
checked sources the exact combination of exhaustive finite SG-to-GFL
replacement minimality, pre-reveal reduced prediction, and a contextual
return conditioned on remaining devices. This is a scoped distinction, not a
proof of priority. See docs/TX4_BLIND_PREDICTION_NOVELTY_AUDIT.md.

## 26. Limitations

This is one documented IEEE-39 phasor-domain implementation with a custom
ten-state GFL device and no EMT portfolio validation. The strict blind
g-boundary timing deviation is disclosed. V9 H recovery fails and runtime
does not improve. ANDES does not reproduce the custom GFL boundary. No
probabilistic generalization, global radius, universal minimality, or
economic-optimal intervention claim is allowed.

## 27. Reviewer 1 assessment

Novelty 7/10; correctness 8/10; model adequacy 6/10; evidence 8/10;
usefulness 7/10; IAS impact 8/10; TPWRS readiness 5/10.

The contribution is the finite blind-predicted mechanism package, not a new
Schur theorem or stability-radius theory. V9 and boundary timing prevent
stronger predictive wording.

## 28. Reviewer 2 assessment

Novelty 7/10; correctness 8/10; model adequacy 5/10; evidence 7/10;
usefulness 7/10; IAS impact 8/10; TPWRS readiness 4/10.

The critical mode is an inter-area electromechanical oscillation with
network-mediated converter participation. Residual, conditioning, and
eigensolver checks do not indicate an index or numerical pathology. Richer
inverter physics, limits, and independent custom-GFL validation remain open.

## 29. IAS poster verdict

CASE B: GOOD IAS POSTER, LIMITED PREDICTIVE CLAIM. Recommended title:
SAFE ALONE, UNSAFE TOGETHER — Predicting Minimal Dynamic Incompatibility in
SG-to-GFL Replacement Portfolios.

The poster may show the V4 16/16 prediction, contextual-return mechanism,
g-only repair, and TDS. It must show V9 as a negative scope result and label
the boundary timing deviation.

## 30. Transactions path

The campaign is a focused IAS poster package and a basis for a new six-page
paper only if the limitations remain visible. It is not TPWRS-ready because
V9 generalization fails, the reduced method is slower, and model/independent
validation scope is narrow. All raw tables, scripts, hashes, figures, report,
and handoff are delivered on research/tx4-blind-portfolio-prediction-final
with no push.

## Reproducibility index

Key outputs are:
results/TX4_V4_BLIND_PREDICTIONS.csv
results/TX4_V4_BLIND_VS_FULL.csv
results/TX4_V9_BLIND_PREDICTIONS.csv
results/TX4_V9_BLIND_VS_FULL.csv
results/TX4_G_BOUNDARY_BLIND_PREDICTION.csv
results/TX4_CONTEXTUAL_RETURN_BOUNDARY.csv
results/TX4_MINIMALITY_SEPARATION.csv
results/TX4_LOCAL_VS_COLLECTIVE.csv
results/TX4_RETURN_DERIVATIVE_CHECK.csv
results/TX4_TDS_FINAL.csv
results/TX4_ANDES_FINAL_CHECK.csv
results/TX4_FINAL_CLAIM_MATRIX.csv

The claim hierarchy is in results/TX4_FINAL_CLAIM_MATRIX.csv and the headline
summary will be in results/TX4_BLIND_PREDICTION_HEADLINE.json.
