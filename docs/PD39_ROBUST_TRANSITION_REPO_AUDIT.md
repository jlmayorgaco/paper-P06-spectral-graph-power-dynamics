# PD39 Robust SG-to-IBR Transition: Repository Archaeology Audit

Date: 2026-09-13

Status: archaeology completed before any PD39 numerical computation. No
PowerDynamics computation has been run by this audit. No frozen TX4 result was
modified. No push is authorized.

## 1. Scope and current repository state

The requested research direction is a new, independent IEEE-39 study:

> robust SG-to-IBR transition co-design under physically meaningful
> controller, network and operating-point perturbations.

The current repository is primarily a Python/NumPy/SciPy/ANDES research
repository. It contains no committed Julia source, no `Project.toml` for
PowerDynamics, and no `pd39` implementation. The existing IEEE-39 model is an
internal/custom phasor-domain model with selected ANDES reconciliation; it is
not a PowerDynamics model.

The active checkout is:

| item | value |
|---|---|
| branch | `research/series-planning-design` |
| HEAD | `7727df7e00082dd419558b47c92987c74dfd421e` |
| HEAD parent | `70318091b9e04a00ae1deb201bde54328f768182` |
| parent of the preregistration line | `e03a5312` |
| worktree | `C:/Users/walla/Documents/Github/paper-P06-spectral-graph-power-dynamics` |
| push status | no push performed |

The checkout is intentionally not clean because it contains user-owned
research artifacts: 888 staged paths, 2 tracked unstaged paths, and roughly
18,929 untracked files under `reports/`. These files must not be reset, cleaned,
or mixed into the PD39 branch.

Existing worktrees:

```text
C:/Users/walla/Documents/Github/paper-P06-spectral-graph-power-dynamics
  research/series-planning-design @ 7727df7e

C:/Users/walla/Documents/Github/paper-P06-cdw68
  research/cdw68-em-reinforcement-confirmation @ 6ccf7d2b
```

## 2. Branch, tag and commit map

| line | ref / commit | scientific role | status for PD39 |
|---|---|---|---|
| main | `main` @ `c7998eac` | repository integration branch | historical context only |
| TX4 | tag `TX4_FINAL_MANUSCRIPT_FREEZE` @ `69f200df` on `research/paremt-emt-validation` | frozen policy-dependent minimal incompatibility and network-closure manuscript | read-only evidence; do not rewrite |
| CDW original | `research/contextual-dynamic-weakness` @ `e73dd355` | original contextual dynamic weakness campaign and final report | read-only predecessor |
| CDW hardening | `research/contextual-dynamic-weakness-hardening` @ `e03a5312` | hardening, cross-model checks, corridors and final reviewer record | clean shared-history candidate |
| CDW68 preregistration | `6f188531` | IEEE-68 confirmatory preregistration before converter portfolios | negative replication context |
| CDW68 results | `research/cdw-ieee68-replication` and worktree alias `research/cdw68-em-reinforcement-confirmation` @ `6ccf7d2b` | IEEE-68 census, two converter models, ranking and replication decision | read-only cross-benchmark context |
| planning preregistration | `70318091` | PD1 adaptive retune and PD2 census-scale design preregistration | predecessor to current planning results |
| planning results | `research/series-planning-design` @ `7727df7e` | PD1 completed; PD2 completed but inconclusive | relevant planning context; not the PowerDynamics benchmark |

Local tags found:

```text
TX4_FINAL_MANUSCRIPT_FREEZE
IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE
IAS2026_PRE_FINAL_VALIDATION
IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE
IAS2026_TRACKA_F8_F12_POST_F7_FREEZE
```

## 3. What the previous research established

### TX3

The frozen TX3 envelope established finite re-equilibrated externalities and
the connected calculus on the canonical matched IEEE-39 / three-GFL benchmark.
The margin-setting mode did not show a material effect in the relevant TX3
tests: C3a and C3b were rejected and C3c remained unresolved. Cycle/SCC
localization was not assessed. TX3 numbers cannot be reused as evidence for a
PD39 weak-node, robust-radius or repair claim.

### TX4

TX4 froze the following model-conditional results:

- transverse stability labels and a policy-dependent minimal incompatibility
  hypergraph;
- exact network-closure / principal-minor anatomy at the observed boundaries;
- actionable boundary motion and controller/network retuning;
- a nominal four-unit witness `{30,33,35,37}` under the frozen model.

The final author decision explicitly limits the witness: `{30,33,35,37}` is a
nominal frozen-model witness, not a robust weak-bus set. Likewise, `kappa = 4`
is a minimum failing-set cardinality / closure order, not evidence of an
irreducible four-device interaction. The TX4 manuscript and evidence remain
frozen at `TX4_FINAL_MANUSCRIPT_FREEZE`.

Important TX4 model limitations are part of the evidence, not implementation
details to silently remove: the harmonized IEEE-39 configuration substitutes a
documented first-order AVR and uses a custom ten-state GFL; the distributed
IEEEX1/IEEEST configuration is not an equivalence point for the journal claims;
the custom GFL is not independently reproduced by stock ANDES.

### CDW original and hardening

The original CDW campaign tested the distinction between static, dynamic,
contextual and intervention weakness. Its final report records GOLD-A and
GOLD-B passes, while GOLD-C and GOLD-D failed informatively and GOLD-E/GOLD-F
were blocked by missing material. The strongest result was total,
re-equilibrated dynamic branch sensitivity versus finite interventions, not a
universal weakness index. No weakness index was constructed.

The hardening branch strengthened the IEEE-39 evidence. The key result was
contextual same-mode reversal on 17/19 base-stable holdout policies and a
large advantage of total branch sensitivity over static baselines. This remains
benchmark/model-conditional evidence.

### IEEE-68 replication

The preregistered IEEE-68 campaign ended in CASE B:

- the nested reversal failed to replicate: 0/16 policies for both converter
  models;
- the portfolio-conditioned total branch ranking passed, with median Spearman
  `0.996` for Model A and `0.981` for Model B;
- the preregistered rightmost mode was a slow control pole with small branch
  effects, so the material electromechanical ranking was post hoc;
- the preregistered TDS holdout had no applicable reversal case.

Therefore, contextual reversal is not a general law. PD39 must preserve this
negative replication and must not phrase IEEE-39 observations as universal.

### Existing planning/design line

The current planning branch contains a preregistered Python/ANDES campaign:

- PD1 adaptively retuned five converter physical log-coordinates at each of
  240 unstable E35 operating points. It reproduced the flagship value at
  240/240, restored transverse stability at 240/240, and reached the stricter
  `-0.02 s^-1` margin at 234/240. The severe tercile was 80/80 stable.
- PD2 re-ran census-scale plan-level versus single-boundary design. All six
  primary T2/T4 tasks stopped without a safe census lattice; the campaign
  verdict is `INCONCLUSIVE`, not a failure of the underlying research question.

These results are useful planning context, but they are not PowerDynamics
results and cannot be presented as the PD39 benchmark.

## 4. Reusable code and data

The main reusable prior-analysis package is under
`reports/poster/ias2026/research/src/ibr_cycles/`:

- `dynamics/`: descriptor DAE, equilibrium, linearization, eigenanalysis and
  modal tracking;
- `certification/`: transverse quotient, symmetry, physical labels and gates;
- `models/`: toy, IEEE-9, IEEE-39 and IEEE-68 model implementations;
- `reduction/` and `actions/`: Schur/self-energy, low-rank action and Green
  operator calculations;
- `nonlinear/`: prior phasor-domain network/TDS components;
- `io/`: manifests and provenance utilities.

Relevant prior experiments include:

- `experiments/E34_mitigation_frontier.py` — controller/network mitigation
  frontier;
- `experiments/E35_mc_operating.py` — operating-point uncertainty envelope;
- `experiments/final_closure/FC01_transverse_quotient.py` — transverse
  re-audit;
- `experiments/planning_design/` — PD1/PD2 planning extensions;
- the CDW hardening and IEEE-68 code under their respective branch-specific
  research directories.

Relevant frozen inputs include:

- `configs/ias2026/ieee39_network.json`;
- `configs/ieee39_harmonized_dynamic_model.yaml`;
- the TX4/CDW preregistrations and claim matrices;
- the stored IEEE-39 operating-point, portfolio and validation artifacts.

These are archaeology and cross-check resources for PD39. They are not a
license to translate the old model into PowerDynamics and call it an
independent benchmark.

## 5. Frozen material and non-reusable claims

The following must remain read-only for PD39:

1. TX4 manuscript, evidence ledger, tag and frozen result tables.
2. TX3 final claim freeze and negative results.
3. CDW and CDW68 final reports as historical evidence, including their
   negative replication and blocked-gate records.
4. Existing result files used to support already-frozen claims.
5. Any old claim that a nominal witness is a robust weak-bus set, that `kappa`
   is an irreducible interaction order, or that contextual reversal is
   universal.

PD39 must not tune a PowerDynamics controller after seeing an old IEEE-39
   outcome. All equilibrium modifications must be re-solved, load service must
   be preserved, and every confirmatory threshold must be frozen before its
   results are computed.

## 6. Contradictory, obsolete or dangerous paths

- The distributed IEEE-39 IEEEX1/IEEEST data and the harmonized substitute are
  different dynamic models. They must not be conflated.
- The stock/custom ANDES converter limitations mean that prior ANDES parity is
  not PowerDynamics parity.
- Old TX3/TX4/CDW narratives contain superseded wording and model-specific
  interpretations. Claim matrices and final reports take precedence over old
  manuscript prose.
- The existing `part5_planning_design` paper and PD1/PD2 results are Python /
  ANDES planning evidence, not evidence that PowerDynamics IEEE-39 has been
  reproduced.
- Generated PDFs, logs, browser/resource captures, scratch files and large
  untracked `reports/` trees are not source inputs for PD39 unless explicitly
  identified in a manifest.

## 7. PD39 readiness and required next boundary

PD39 is not numerically started. Before Phase 1, the dedicated branch must
freeze at least:

- PowerDynamics and Julia versions and a reproducible project environment;
- stock IEEE-39 network provenance and the exact equilibrium convention;
- candidate SG set and MW-based replacement semantics;
- the exact GFL model, parameters, controller defaults and limits;
- load delivery, dispatch, reactive-power and feasibility rules;
- transverse stability definition, mode handling and admissibility rules;
- perturbation classes and normalized envelopes for controller, network and
  operating-point actions;
- structured-radius search, bracketing, multistart and audit rules;
- hidden-fragility thresholds, co-design targets and TDS selection rules;
- deterministic manifests and the three preregistration documents.

The first PD39 implementation should therefore be a model/equilibrium audit
and qualification gate, not a portfolio sweep. No PD39 final question can yet
be answered from the repository, including candidate eligibility, feasible
portfolio counts, structured-radius distributions, weak nodes/links or TDS
confirmation.

## 8. Dedicated-worktree decision

The cleanest practical starting point is the committed current HEAD
`7727df7e` (which descends from `e03a5312` and contains the latest planning
context), with this audit committed separately before the dedicated branch is
created. The dedicated branch name required by the research brief is:

```text
research/pd39-robust-ibr-transition-codesign
```

The branch must be created in a separate clean worktree. The dirty checkout
must remain untouched. All subsequent PD39 code, preregistration, model audit,
results and reports belong to the dedicated branch. No push is permitted.

