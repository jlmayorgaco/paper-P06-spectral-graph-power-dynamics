# Reproduction guide

Run from the repository root with the checked-in Julia environment (Julia 1.11.9 was used). The order below matches the executed campaign: D0 first, freeze all delay assignments, then inspect positive-delay roots.

```powershell
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/baseline_reproduction.jl
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/physical_baseline_reproduction.jl
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/physical_holdout_reproduction.jl
python experiments/delay_dressed_replacement_frontier_20261002/audit_d00.py
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/run_d01_tau0_identity.jl
python experiments/delay_dressed_replacement_frontier_20261002/audit_d01_tau0.py
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/freeze_delay_patterns.jl
python experiments/delay_dressed_replacement_frontier_20261002/freeze_delay_patterns.py
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/run_d01_delay_roots.jl
python experiments/delay_dressed_replacement_frontier_20261002/audit_d01_delay.py
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/run_d02_local_sensitivities.jl
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/run_d03_taylor_validation.jl
julia --project=. experiments/delay_dressed_replacement_frontier_20261002/run_d04_graph_operator.jl
python experiments/delay_dressed_replacement_frontier_20261002/finalize_tables.py
python experiments/delay_dressed_replacement_frontier_20261002/generate_figures.py
python experiments/delay_dressed_replacement_frontier_20261002/package_experiment_zip.py
```

The D1 contour counter may exceed practical runtime for its conservative norm-derived contour. That outcome is recorded as `BLOCKED_EXACT_DDE_SPECTRUM`; do not replace it with Padé or treat a few tracked Newton roots as complete. D4 exits successfully after recording the singular Schur gate and keeps pre-gauge values only under `*_PRE_GAUGE_DIAGNOSTIC` names.

## What reproduces and what does not

- D0 reproduction is complete and reported in `baseline_reproduction/D00_GATE_REPORT.md`.
- The τ=0 characteristic identity reproduces the full finite ODE spectrum for the stored baseline and joint designs.
- Three positive-delay patterns have locally continued dominant ODE pairs with small nonlinear-eigenvalue residuals. The rightmost DDE root count is absent.
- The 104 heterogeneous delay vectors were frozen from one multiset and hashed before any positive-delay root calculation. Four graph-selected patterns are heuristic selections among 100,000 seeded candidate permutations, not certified extrema; 100 seeded random assignments are also recorded.
- Lossless graph-signal descriptors, including the commutator identity check, are recomputed from the frozen delay vectors and graph matrices. There is no validated Schur damping operator to compare to them.
- The lossless branch Taylor test has 60 deterministic direction/amplitude cases. It is neither an AC-loss Taylor validation nor a nonlinear stability proof.
- Delayed co-design, nonlinear DDE events, shadow prices, and global/box upper bounds were not run because required spectral and Schur gates are blocked.
- `package_experiment_zip.py` refreshes the SHA-256 inventory, creates the complete local ZIP bundle, and verifies its archive integrity. The ZIP is excluded from the inventory to avoid self-reference.

Inputs, scripts, status-only tables and all new outputs live under this experiment directory. `BASELINE_MANIFEST.json` records the original Git SHA, dirty state, model-source hashes and package versions. `EXPERIMENT_ARTIFACT_HASHES.json` records the experiment files. Do not edit or regenerate the historical directories referenced by those manifests.
