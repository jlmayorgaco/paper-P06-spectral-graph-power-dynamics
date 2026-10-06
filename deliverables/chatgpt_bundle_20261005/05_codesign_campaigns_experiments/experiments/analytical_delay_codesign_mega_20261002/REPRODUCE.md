# Reproduction sequence

Run from the repository root with its existing Julia project and Python environment. Scripts write only in this directory. Do not run the M1 candidate generator after changing the frozen reference designs; the candidate TOMLs and hashes are already saved.

```powershell
python experiments/analytical_delay_codesign_mega_20261002/write_manifest.py
julia --project=. experiments/analytical_delay_codesign_mega_20261002/q0_validate_five_events.jl experiments/analytical_delay_codesign_mega_20261002/seed_uniform_875.toml
julia --project=. experiments/analytical_delay_codesign_mega_20261002/q12_validate_low_rank.jl
julia --project=. experiments/analytical_delay_codesign_mega_20261002/q4_validate_pi_boundary.jl
julia --project=. experiments/analytical_delay_codesign_mega_20261002/q56_validate_action_space.jl
julia --project=. experiments/analytical_delay_codesign_mega_20261002/m0_recheck_zero_schur.jl
python experiments/analytical_delay_codesign_mega_20261002/m0_aggregate.py
```

The M1 inputs were frozen in `M1_FROZEN_CANDIDATES.json` before evaluation; the sole midpoint is frozen in `M1_BISECTION_FREEZE.json`. Run each existing candidate with `m1_validate_candidate.jl`, then `m1_aggregate.py`. Each run takes roughly minutes because it completes all five nonlinear 61-second events.

```powershell
julia --project=. experiments/analytical_delay_codesign_mega_20261002/m1_validate_candidate.jl experiments/analytical_delay_codesign_mega_20261002/m1_candidates/eta_0375_bisect.toml
python experiments/analytical_delay_codesign_mega_20261002/m1_aggregate.py
julia --project=. experiments/analytical_delay_codesign_mega_20261002/m2_random_action_space.jl
```

For M3, run `m3_a_trace_integral.jl` with delay in milliseconds and optional second argument design TOML. The executed runs used 20, 30 and 40 ms on the seed and on `M1_ZERO_DELAY_DESIGN.toml`. The 40 ms roots were discovered and refined with `m3_c_refine_roots.jl 40` and `m3_c_refine_roots.jl 40 experiments/analytical_delay_codesign_mega_20261002/M1_ZERO_DELAY_DESIGN.toml`. The contour is large (~10,414 rad/s), so these runs are computationally substantial. `m3_a_trace_integral.jl` uses Float64 adaptive quadrature; results are numerical, not interval certified.

```powershell
python experiments/analytical_delay_codesign_mega_20261002/finalize_outputs.py
python experiments/analytical_delay_codesign_mega_20261002/hash_outputs.py
```

`finalize_outputs.py` recreates derived tables and figures from primary executed outputs. It leaves status-bearing BLOCKED tables for gated stages; those rows are not fictitious experiment results. The old phase-only winding indices 64/62 must not be imported as verified DDE root counts.
