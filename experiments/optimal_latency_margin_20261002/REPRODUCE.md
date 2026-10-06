# Reproduce the targeted calculations

Run at the repository root with Julia 1.11.9, the checked-in `Project.toml`/`Manifest.toml`, and Python 3.13.9 with Matplotlib 3.10.6 and NumPy. The scripts import read-only source models and frozen designs from adjacent experiment directories. They write only into this directory. Some complete contour passes take several minutes.

The exact script sequence used to generate the results is:

```powershell
python experiments/optimal_latency_margin_20261002/freeze_protocol.py
julia --project=. experiments/optimal_latency_margin_20261002/run_t01_counts.jl
julia --project=. experiments/optimal_latency_margin_20261002/run_slow_mode_candidates.jl
julia --project=. experiments/optimal_latency_margin_20261002/run_rightmost_exclusion.jl
julia --project=. experiments/optimal_latency_margin_20261002/run_fast_root_locus.jl
julia --project=. experiments/optimal_latency_margin_20261002/run_t02_bisect.jl
julia --project=. experiments/optimal_latency_margin_20261002/refine_t02_crossings.jl
julia --project=. experiments/optimal_latency_margin_20261002/verify_t02_boundary.jl
julia --project=. experiments/optimal_latency_margin_20261002/run_t03_mode_sensitivities.jl
julia --project=. experiments/optimal_latency_margin_20261002/validate_t03_tau.jl
julia --project=. experiments/optimal_latency_margin_20261002/compare_mode_families.jl
python experiments/optimal_latency_margin_20261002/freeze_fixed_gain_designs.py
julia --project=. experiments/optimal_latency_margin_20261002/run_t04_fixed_gain_frontier.jl
julia --project=. experiments/optimal_latency_margin_20261002/refine_t04_crossings.jl
julia --project=. experiments/optimal_latency_margin_20261002/verify_t04_endpoint_crossings.jl
python experiments/optimal_latency_margin_20261002/augment_t01_candidates.py
python experiments/optimal_latency_margin_20261002/make_poster_outputs.py
```

`verify_t04_endpoint_crossings.jl` currently targets ±0.01 ms and writes `T04_ENDPOINT_BOUNDARY_VERIFICATION_0P01MS.csv`. The earlier ±0.001 ms calls are preserved separately in `T04_ENDPOINT_BOUNDARY_VERIFICATION.csv`; each was INDETERMINATE from unresolved contour quadrature. The rerun must never reinterpret them as safe.

The frozen input hashes are in `FROZEN_PROTOCOL.json`; derived file hashes are in `RESULT_HASHES.json`. Reexecuting `freeze_protocol.py` updates its timestamp and should be done only when intentionally starting a fresh run, because the original preregistration is part of this evidence. To inspect the existing results, rerun only `make_poster_outputs.py` and the verification code; do not overwrite the frozen protocol.
