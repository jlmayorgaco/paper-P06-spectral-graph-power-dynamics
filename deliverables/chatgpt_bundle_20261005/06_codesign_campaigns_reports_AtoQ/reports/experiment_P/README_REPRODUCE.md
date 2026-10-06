# Reproducing ExpP

Run from the repository root in PowerShell, on the existing Julia 1.11.9 environment. The project and manifest are frozen to ExpN hashes; do not run `Pkg.update`, regenerate the manifest, or change the ExpN candidate. All commands use the installed environment and checked-in source.

## Inputs and provenance

- Branch: `research/expN-pd-exact-zstar`; parent commit: `d0fecb3264aeb855ab6700fc6ef66120d4b6c33c`.
- `Project.toml` SHA-256: `e170a2f9fca57aa3f8f3f053b318f0d5e3d434679dbdd4a4b8eef5481cfc0e4f`.
- `Manifest.toml` SHA-256: `92bc7849ea4445a4bedd18f732611ec271147e238c349aa5eeefc6d02023e911`.
- Julia 1.11.9, PowerDynamics 5.0.0, NetworkDynamics 1.3.0.
- Model SHA and ExpN candidate SHA are recorded in [P0_BASELINE.json](P0/P0_BASELINE.json) and [MODEL_FREEZE.json](../experiment_N/MODEL_FREEZE.json).

## Ordered run

The single stage runner is `experiments/bnd_expP/run_expP.jl`. On a pristine ExpP output directory, stages P0–P4 and P6 can be invoked in order:

```powershell
julia --project=. --startup-file=no experiments/bnd_expP/run_expP.jl --stage P0
julia --project=. --startup-file=no experiments/bnd_expP/run_expP.jl --stage P1
julia --project=. --startup-file=no experiments/bnd_expP/run_expP.jl --stage P2
julia --project=. --startup-file=no experiments/bnd_expP/run_expP.jl --stage P3
julia --project=. --startup-file=no experiments/bnd_expP/run_expP.jl --stage P4
```

`--through P2`, `--through P5`, `--resume` are implemented. The frozen P5 candidate already exists in this checkout; its pre-freeze action protects against overwriting it. For this completed run, rerun post-freeze validation using the dedicated scripts below, then P6. Do not use `--stage P5` to overwrite frozen results.

```powershell
julia --project=. --startup-file=no experiments/bnd_expP/validate_expP_tds.jl
julia --project=. --startup-file=no experiments/bnd_expP/compare_frozen_linear_tds.jl
julia --project=. --startup-file=no experiments/bnd_expP/validate_expP_declared_event.jl
julia --project=. --startup-file=no experiments/bnd_expP/run_expP.jl --stage P6
```

The specific P2 new-design guard calculation is reproducible without touching the frozen candidate:

```powershell
julia --project=. --startup-file=no experiments/bnd_expP/audit_expP_numerical_guard.jl
```

## Tests

```powershell
julia --project=. --startup-file=no test/bnd_expP/runtests.jl
```

Tests recompute ExpN freeze checks, PLL/SG determinant identities, conditional full-spectrum retention roots, and assert the recorded event metrics, guard root, frozen candidate status, rejected CARE attempt, conditional interval-bound margin/caveat, and open globality ledger. The final suite passed 56 checks across four test sets; see [TEST_RESULTS.json](TEST_RESULTS.json). It emits a benign Julia namespace warning because `ExpP.ROOT` conflicts with an existing identifier in `Main`. The suite does not include fault-injection mutation tests for every corrupted-code scenario listed in the specification.

## Outputs

- `P0`–`P6`: stage summaries, machine-readable results, CSV tables, and the frozen candidate/hash.
- [FUNCTION_INVENTORY.csv](P0/FUNCTION_INVENTORY.csv): reused ExpN APIs and PD-only validation APIs.
- [CLAIM_LEDGER_EXP_P.csv](CLAIM_LEDGER_EXP_P.csv): claim, evidence type, and limits.
- [EQUATIONS_SOURCE_MAP.md](EQUATIONS_SOURCE_MAP.md): source functions, hashes, and sign/unit conventions.
- [PROVENANCE.toml](PROVENANCE.toml): branch/parent, dependency and model hashes, and the supplied specification hashes.
- [STAGE_STATUS.json](STAGE_STATUS.json): most recently completed stage plus the independent P5 subgate outcomes.
- `FIG_P01_conditional_retention_roots.png`, `FIG_P02_robust_pointwise_bounds.png`, and `FIG_P03_declared_step_tds.png` are regenerated from saved CSVs with `python experiments/bnd_expP/render_expP_figures.py`.

Stage elapsed values are recorded in each result. Julia compilation time and peak memory were not measured separately. The P3 wrapper did not instrument direct sensitivity or QP subcalls; its counters therefore report wrapper-level work only.
