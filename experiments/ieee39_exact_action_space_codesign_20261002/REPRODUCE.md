# Reproduce Experiment Q

Run from the repository root with Julia 1.11 and the checked-in project environment.

```powershell
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q0_validate_five_events.jl experiments/ieee39_exact_action_space_codesign_20261002/seed_uniform_875.toml
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q12_validate_low_rank.jl
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q3_count_dde_roots.jl
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q3_sample_nyquist_diagnostic.jl
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q4_validate_pi_boundary.jl
julia --startup-file=no --project=. experiments/ieee39_exact_action_space_codesign_20261002/q56_validate_action_space.jl
```

`q3_count_dde_roots.jl` performs the Float64 contour diagnostic. Positive-delay output must be treated as indeterminate unless independent root validation is completed; it must not be promoted to a safe/unsafe certified count. `TABLE_Q03_CONTOUR_ATTEMPT.csv` preserves the original count attempt, while `TABLE_Q03_DDE_ROOT_COUNTS.csv` reports only the accepted classification and stores positive-delay winding indices separately as unverified diagnostics.

The plotted figures are generated from the saved Q3 diagnostic locus and Q4/Q5 CSVs by `make_q_figures.py` after the Julia runs. Matplotlib, pandas, and NumPy are required. All data, scripts, and reports generated for this experiment belong in this directory. Prior experiment directories are read-only inputs.

The current experiment reaches Q0 and the structural Q1/Q2/Q4/Q5/Q6 tests. The Q3 gate remains indeterminate, so the conditional later stages Q7 and Q8, delayed event validation, and the delay frontier have not been executed.
