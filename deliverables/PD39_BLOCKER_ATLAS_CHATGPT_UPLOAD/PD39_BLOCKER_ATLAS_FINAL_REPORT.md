# PD39 Physical Blocker Validation + First Compatibility Atlas

## Executive verdict

This campaign is **CASE C - STOP**.

The fresh-holdout census supplied three low-order C12 portfolios that are
candidate blockers on the default PowerDynamics initialization path:

| candidate | portfolio | default-path alpha [s^-1] | strict Gate A alpha [s^-1] | C12 TDS |
|---|---:|---:|---:|---|
| H1 | 35;36 | +10.1875267 | -0.0979116 | MaxIters |
| H2 | 37;38 | +22.5280341 | +11.8934740 | MaxIters |
| H3 | 30;32;33 | +21.1535104 | -0.1088351 | MaxIters |

The decisive result is not a physical-blocker discovery. It is a physicality
warning: the candidate signs and magnitudes depend strongly on the
equilibrium initialization/tolerance path, and all three C12 integrations
hit the explicit TDS `MaxIters` guard. Their immediate predecessors and
matched controls solve and decay. No selected candidate passed both the
numerical/DAE physicality gate and the nonlinear TDS gate.

The conditional compatibility atlas was therefore not run. No physical
minimal-blocker, compatibility-boundary, local-remediation, contextual-return,
or atlas-theorem claim is permitted from this branch.

## 1. Scope and frozen decision rule

The branch was created from the latest authoritative PD39 result at
`e81d18ab945b` and the preregistration package was committed before new
numerics. The final branch is
`research/pd39-physical-blockers-atlas-v1`; no push was performed.

This narrow campaign allowed only:

1. Gate A numerical/DAE physicality audit.
2. Gate B common nonlinear TDS.
3. Gate C aggregate-matched composition check, conditional on physical
   blocker confirmation.
4. Exact network-closure feasibility/validation audit.

The atlas, compatibility boundaries, weak-node/link ranking, structured
radius campaign, repair optimization, planners, co-design, IEEE-68, EMT, and
new controller search were not run.

The frozen Gate A minimality statistic was

`delta_H = min(alpha(H), -max(alpha(G) for proper G subset H))`,

with a primary separation threshold of `0.02 s^-1`. A positive alpha is a
candidate true-instability signal only after equilibrium and independent
numerical checks are reproducible. A large positive alpha alone is not
physicality evidence.

## 2. Archaeology and selection

The authoritative input was the committed complete census
`results/PD39_255PLUS1_HOLDOUT_CENSUS.csv`, containing 256 portfolios x 24
fresh conditions = 6144 portfolio-condition evaluations. The prior frozen
census summary reported 42 minimal H0 cases and 47 minimal H0.05 cases over
the 24 conditions. Those counts are retained here as prior context; this
narrow campaign did not rerun the complete census.

The deterministic selection rule chose the first two cardinality-2 H0
candidate entries and the first cardinality-3 candidate in the frozen order:

| id | blocker | condition | matched control | blocker/control aggregate distance |
|---|---|---|---|---:|
| H1 | 35;36 | C12 | 32;36 | 0.0013860 |
| H2 | 37;38 | C12 | 36;38 | 0.0148840 |
| H3 | 30;32;33 | C12 | 30;33;35 | 0.0014786 |

The controls were selected before the new numerical results using same
cardinality, then MW, MVA, remaining inertia, and lexicographic tie breaks.

## 3. Gate A: numerical and DAE physicality

The audit covered 17 unique portfolios: the three selected candidates, all
immediate predecessors and cheap proper subsets, and the three matched
controls. It used C12 and nominal conditions, tolerances `1e-8`, `1e-10`, and
`1e-12`, the native PowerDynamics reduction, an independent dense spectrum,
and a descriptor/generalized check where accessible. Central finite-difference
derivative checks used the preregistered fixed relative steps.

Within a fixed solved path, dense and descriptor alpha agreed with the native
alpha to numerical precision and their signs agreed. This is a limited
within-path check, not confirmation that the equilibrium path is the correct
physical branch.

The default census initialization audit reproduced the census candidate signs:

| id | default alpha | default frequency [Hz] | default mechanism label | g_z condition | smallest singular value | residual |
|---|---:|---:|---|---:|---:|---:|
| H1 | +10.1875267 | 0.000000 | electromechanical/control | 1.526e3 | 3.91e-16 | 3.54e-12 |
| H2 | +22.5280341 | 0.549890 | converter-control oscillatory | 1.504e3 | 1.47e-16 | 7.04e-13 |
| H3 | +21.1535104 | 0.000000 | converter-control oscillatory | 1.998e3 | 8.80e-18 | 7.04e-13 |

The strict Gate A path changed the interpretation:

| id | alpha at 1e-8 | alpha at 1e-10 | alpha at 1e-12 | Gate A interpretation |
|---|---:|---:|---:|---|
| H1 | approximately -0.0979 | -0.0979 | fixed-point gate failed at residual 3.54e-12 | stable on the strict solved path; not reproducible as default unstable candidate |
| H2 | approximately +11.8935 | +11.8935 | +22.5280 | positive but extremely path-sensitive in magnitude |
| H3 | approximately -1.18e-8 | -0.1088 | +21.1535 | sign changes with tolerance/path; unresolved |

The surrounding proper subsets and controls remained negative on their solved
paths. At default-path C12, the immediate-predecessor alphas were:

| blocker | immediate predecessor alpha values [s^-1] | default-path delta_H [s^-1] |
|---|---|---:|
| H1 | 36: -0.0971024; 35: -0.0977616 | 0.0971024 |
| H2 | 38: -0.0975787; 37: -0.1378248 | 0.0975787 |
| H3 | 32;33: -0.0978716; 30;33: -0.1081148; 30;32: -0.1051621 | 0.0978716 |

These separations exceed `0.02 s^-1` only on the default initialization
path. They are not physical minimality results because the candidate
equilibria are not path-reproducible.

The audit also exposed very large reduced-Jacobian conditioning, near-zero
singular values, and large state/Jacobian excursions on candidate paths.
PowerDynamics emitted vertex-batch compatibility warnings. The combination
is consistent with a numerical/DAE branch problem being plausible. It does
not identify a unique physical mechanism.

**Gate A: FAIL_UNRESOLVED.**

## 4. Gate B: common nonlinear TDS

The frozen disturbance was used without modification: bus 39, +1% P and Q at
constant power factor from 1.0 to 1.1 s, restoration at 1.1 s, integration to
20 s, common Rodas5P settings and common observables. The runtime guard
`maxiters=100000` was added after an unbounded selected case exceeded the
default iteration budget; this did not change the disturbance, tolerances,
acceptance thresholds, or axes. The return is classified as a TDS failure.

The run contained 26 rows: all selected blockers, immediate predecessors and
matched controls at C12 and nominal.

At C12, all three selected blockers returned `MaxIters`:

| id | linear alpha [s^-1] | TDS status | nonlinear rate | interpretation |
|---|---:|---|---:|---|
| H1 | +10.1875 | failed / MaxIters | unavailable | no nonlinear confirmation |
| H2 | +22.5280 | failed / MaxIters | unavailable | no nonlinear confirmation |
| H3 | +21.1535 | failed / MaxIters | unavailable | no nonlinear confirmation |

The immediate predecessors and controls completed at C12 and decayed. Typical
estimated rates were approximately -0.207 to -0.258 s^-1 for predecessors
and -0.234, -0.236, and -0.258 s^-1 for the three controls. The nominal rows
for the selected blocker labels completed with negative rates because nominal
is not the C12 candidate operating point.

The C12 failed integrations showed very large Jacobian/state diagnostics
before the solver stopped. The correct interpretation is solver divergence or
unresolved trajectory behavior, not evidence of a physical unstable mode.

**Gate B: FAIL.**

## 5. Gate C: composition

The preregistered strong-composition gate was conditional on physical blocker
confirmation. Because Gate A/B did not confirm any physical blocker, no
composition result is promoted and no post hoc threshold is applied.

For transparency, the raw default-path candidate/control comparisons were
written to `results/PD39_BLOCKER_COMPOSITION_TEST.csv`. H2 has a large sign
contrast on the default path, but that contrast is not admissible evidence
for a physical composition effect while the candidate path is unresolved.

**Gate C: NOT APPLICABLE; composition strong = false.**

## 6. Penetration, composition, and mode verdict

The narrow result is **no physical blocker and no compatibility atlas**.

The prior holdout census shows that H0/H0.05 cases are distributed across
multiple conditions rather than defining a universal V8-like identity. The
selected C12 candidates are low-cardinality default-path candidates, not
general physical blockers.

The apparent default-path mode labels are heterogeneous: H1 is labeled
electromechanical/control and H2/H3 are labeled converter-control under the
default path. The tolerance/path changes and TDS noncompletion prevent the
campaign from selecting a physical 7/8 -> 8/8 mechanism. No continuous
physical SG-to-GFL homotopy was run or claimed.

## 7. Exact network-closure feasibility

The installed PowerDynamics model/API does not expose the exact device-port
objects, local admittance partition, and consistently dimensioned `K`, `D`, and
`Q` objects required for the requested network-closure identity. Therefore
the PD contextual return is `BLOCKED`.

The conditional TX4 contextual-return validation was not run after the
physicality STOP. No proxy was constructed, no determinant identity was
claimed, and no `Q -> -1` boundary was manufactured. The exact status is in
`results/TX4_CONTEXTUAL_RETURN_VALIDATION.csv`.

**Closure: BLOCKED.**

## 8. Compatibility atlas and conditional theory

The two-coordinate atlas (`PLL gain scale`, `filter scale`) was not run. No
11 x 11 grid, adaptive boundary, `kappa_0`, `kappa_0.05`, mode-MAC boundary
trace, or contextual-return theorem validation is present in this campaign.

This is intentional. A compatibility atlas built on path-sensitive
equilibria and failed TDS cases would overstate the evidence.

## 9. Figures and retained evidence

The STOP figures are audit figures, not atlas figures:

- `figures/PD39_blocker_atlas/F1_path_dependence_selected_cases.png` and PDF:
  default versus strict-path alpha for candidates and controls.
- `figures/PD39_blocker_atlas/F2_tds_stop_summary.png` and PDF: C12 TDS
  completion/failure pattern.
- `figures/PD39_blocker_atlas/F3_gate_a_tolerance_sweep.png` and PDF: selected
  C12 alpha across the frozen tolerance sweep.

All raw Gate A tables, Gate B traces, deviations, preregistration documents,
and selection records are retained in the repository and in the reproduction
bundle.

## 10. Negative results and forbidden claims

The following claims are not supported:

- that `35;36`, `37;38`, or `30;32;33` is a physically validated minimal H0
  blocker;
- that the default-path positive alpha is a clean electromechanical or
  converter-control instability;
- that a 7/8 -> 8/8 physical mode transition has been established;
- that aggregate composition causes the blocker;
- that a compatibility atlas or blocker boundary exists;
- that a network-closure determinant or contextual-return theorem has been
  validated;
- that any controller, line, topology, or mixed repair works;
- that the result supports generic weak-node, weak-link, hidden-radius,
  holonomy, graph-Laplacian, or broad robust-co-design novelty.

The strongest defensible positive statement is narrower: the stock PD39
workflow produces fresh-holdout low-order instability candidates whose
default-path signs are not reproducible under the frozen numerical audit and
whose common-disturbance trajectories can terminate in solver iteration
failure, while neighboring controls and predecessors decay. This is a
validation warning requiring model/equilibrium-path repair.

## 11. Practical interpretation

For an SG-to-IBR planner, the immediate practical lesson is a quality-control
requirement: do not promote a portfolio to a physical blocker from a single
equilibrium initialization and a single rightmost eigenvalue. Require branch
reproducibility, independent spectrum agreement, conditioning diagnostics,
and a common nonlinear disturbance that actually completes.

The result does not yet justify redesigning the IAS poster around a physical
blocker or compatibility atlas. It is not Transactions-ready as a new
mechanism result.

## 12. Final answers

1. The new campaign did not establish a physical list of blockers. It found
   three default-path C12 candidates: `35;36`, `37;38`, and `30;32;33`.
2. Their physical status is unresolved; the default-path positives are not
   sufficient evidence.
3. A physically validated H0 is not established.
4. A physically validated H0.05 is not established.
5. Default-path immediate-predecessor separations are 0.0971024, 0.0975787,
   and 0.0978716 s^-1 for H1, H2, and H3, respectively.
6. No physical mode change was established.
7. The apparent default modes are path-sensitive and therefore unresolved.
8. Conventional diagnostics were not advanced in this narrow campaign; the
   conditioning/branch issue must be resolved first.
9. `rho_0` was not computed; the atlas was stopped.
10. `rho_0.05` was not computed; the atlas was stopped.
11. No matched-alpha hidden-fragility result was tested in this campaign.
12. No strongest hidden-fragility pair is admissible.
13. No alpha-versus-radius correlation is reported.
14. No dynamic weak node was computed.
15. No dynamic weak link was computed.
16. No static ranking comparison was advanced.
17-20. No 25/50/75/100% designs were run; the scope explicitly prohibited
   repair optimization and co-design.
21. Controller-only repair: not run.
22. Line-only repair: not run.
23. Topology-only repair: not run.
24. Smallest intervention: none found; no intervention search was authorized.
25. Load service for a proposed repair: not applicable.
26. Fresh 24-condition validation of a repair: not applicable.
27. TDS does not support a physical 7/8 -> V8 -> repaired-V8 story because
   the selected C12 blocker integrations fail `MaxIters`.
28. Genuinely new: a narrow, reproducible warning about path-sensitive
   low-order blocker candidates in this stock PD39 workflow.
29. Classical: eigenvalue checks, independent spectra, TDS, conditioning, and
   equilibrium reproducibility are classical validation requirements.
30. Planner usefulness: currently limited to a numerical quality-control
   warning; no planning mechanism is validated.
31. IAS poster redesign: NO.
32. Research status: STOP and repair the equilibrium/model path before any
   blocker or atlas theory is revisited.

## 13. Decision record

| item | result |
|---|---|
| numerical audit | FAIL_UNRESOLVED |
| high-PLL / selected-blocker TDS | FAIL: 3 C12 MaxIters |
| prior H0 counts over 24 holdouts | 42 minimal cases |
| prior H0.05 counts over 24 holdouts | 47 minimal cases |
| penetration/composition/mixed verdict | no physical blocker; composition not promoted |
| mode mechanism | unresolved, path-dependent default labels |
| closure | BLOCKED |
| atlas | NOT RUN_STOP |
| IAS poster redesign | NO |
| Transactions-ready | NO |

**FINAL CASE: C**

