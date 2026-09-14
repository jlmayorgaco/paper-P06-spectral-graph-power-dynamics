# ChatGPT handoff - PD39 physical blocker validation

## One-line verdict

**CASE C - STOP.** Three C12 low-order default-path candidates were found,
but none is a physically confirmed blocker. Gate A is unresolved because the
candidate signs are path/tolerance-sensitive; Gate B fails because all three
C12 candidate TDS integrations return `MaxIters`. The compatibility atlas was
not run.

## Repository and provenance

- Final branch: `research/pd39-physical-blockers-atlas-v1`
- Parent authoritative result: `e81d18ab945b` on
  `research/pd39-255plus1-mechanism-validation`
- Preregistration commit: `1b9ad5fe`
- Gate A script commit: `02381549`
- Default-path audit commit: `d24f0ad3`
- Gate B/TDS guard commit: `c95bab91`
- No push was performed.
- Runtime: Julia 1.11.9, PowerDynamics 5.0.0.

## Frozen selection

The complete prior census was 256 portfolios x 24 holdout conditions = 6144
evaluations. The frozen prior summary had 42 minimal H0 cases and 47 minimal
H0.05 cases. This campaign did not rerun that census.

Selected in frozen order:

| id | candidate | condition | matched control | distance |
|---|---|---|---|---:|
| H1 | 35;36 | C12 | 32;36 | 0.0013860 |
| H2 | 37;38 | C12 | 36;38 | 0.0148840 |
| H3 | 30;32;33 | C12 | 30;33;35 | 0.0014786 |

## Gate A: what was actually found

The audit covered 17 portfolios at C12 and nominal with tolerances 1e-8,
1e-10, 1e-12, native reduction, independent dense spectrum, descriptor
check where accessible, central-FD diagnostics, conditioning, singular
values, and participation.

Default census path at C12:

| id | alpha | frequency | label |
|---|---:|---:|---|
| H1 | +10.1875267 | 0.000000 Hz | electromechanical/control |
| H2 | +22.5280341 | 0.549890 Hz | converter-control oscillatory |
| H3 | +21.1535104 | 0.000000 Hz | converter-control oscillatory |

Strict Gate A path at 1e-10:

- H1: -0.0979116 s^-1, 0.067999 Hz.
- H2: +11.8934740 s^-1, 0 Hz.
- H3: -0.1088351 s^-1, 0.074147 Hz.

Within each solved path, dense and descriptor spectra agree. The issue is that
the equilibrium/path itself is not reproducible across the prescribed audit
paths. Large conditioning and near-zero singular values reinforce the
unresolved DAE/branch concern.

Default-path immediate-predecessor separation is 0.0971024, 0.0975787, and
0.0978716 s^-1 for H1-H3, but it is not a physical minimality result while
the candidate branch is unresolved.

**Gate A: FAIL_UNRESOLVED.**

## Gate B: common TDS

Disturbance: bus 39, +1% P/Q constant power factor from 1.0 to 1.1 s,
restore, integrate to 20 s, common Rodas5P settings. There were 26 rows.

All three selected C12 candidates returned `MaxIters` after the explicit
100000-iteration guard. The guard was recorded as a deviation only to prevent
unbounded runtime; no physics or threshold was changed. Immediate
predecessors and matched controls completed and decayed, with representative
rates between about -0.207 and -0.258 s^-1. A failed TDS is not a physical
instability confirmation.

**Gate B: FAIL.**

## Gate C and atlas

Gate C was conditional on physical blocker confirmation, so it was not entered.
The transparent raw default-path comparisons are in
`results/PD39_BLOCKER_COMPOSITION_TEST.csv`; `composition_strong=false`.

The atlas, 11 x 11 PLL/filter grid, adaptive boundaries, kappa counts,
contextual-return theorem validation, local remediation, and all excluded
campaigns were not run.

Exact PD contextual return is `BLOCKED` because the installed API does not
expose exact K/D/Q port objects and dimensions. TX4 was `NOT_RUN_STOP`; no
proxy was used.

## What ChatGPT should conclude

Allowed:

- The stock workflow can generate low-order C12 candidate blockers that are
  highly sensitive to equilibrium path/tolerance.
- The selected candidates cause TDS solver failure while neighbors decay.
- The correct next task is to repair and independently validate the
  equilibrium/model branch.

Not allowed:

- Calling any candidate a physical minimal H0/H0.05 blocker.
- Naming a physical 7/8 -> 8/8 mode mechanism.
- Claiming compatibility-atlas boundaries, hidden-radius evidence,
  contextual-return closure, or any intervention/repair result.
- Reviving Shapley, cycle/holonomy, graph-Laplacian/Fiedler, generic weak-node,
  weak-link, or broad co-design novelty from this STOP result.

## Final machine-readable status

- numerical audit: `FAIL_UNRESOLVED`
- selected-blocker TDS: `FAIL`
- prior H0 counts over 24 holdouts: `42`
- prior H0.05 counts over 24 holdouts: `47`
- physical blockers confirmed: `0`
- composition strong: `false`
- atlas: `NOT_RUN_STOP`
- closure: `BLOCKED`
- IAS poster redesign: `NO`
- Transactions-ready: `NO`
- final case: `C`

## Key files

- `docs/PD39_BLOCKER_ATLAS_FINAL_REPORT.pdf`
- `docs/PD39_BLOCKER_ATLAS_FINAL_REPORT.md`
- `results/PD39_BLOCKER_ATLAS_HEADLINE.json`
- `results/PD39_BLOCKER_ATLAS_MASTER_CLAIMS.csv`
- `results/PD39_BLOCKER_ATLAS_MASTER_CASES.csv`
- `results/PD39_BLOCKER_NUMERICAL_AUDIT.csv`
- `results/PD39_BLOCKER_DEFAULT_PATH_AUDIT.csv`
- `results/PD39_BLOCKER_TDS_VALIDATION.csv`
- `deliverables/PD39_BLOCKER_ATLAS_CHATGPT_UPLOAD.zip`

