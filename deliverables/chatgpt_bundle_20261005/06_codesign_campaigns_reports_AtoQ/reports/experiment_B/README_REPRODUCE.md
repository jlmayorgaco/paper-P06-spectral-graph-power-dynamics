# Reproducing Experiment B-pre

From the repository root, run:

```powershell
julia --project=. experiments/bnd_expB/run_experiment_B.jl
julia --project=. test/bnd_expB/runtests.jl
```

The first command regenerates all deterministic synthetic cases, ten required CSV tables, auxiliary spectral tables, `RESULTS_EXP_B.json`, the report, and figures in PNG, PDF, and SVG formats. It records the Julia version, BLAS configuration, Git HEAD, and seeds. Figure rendering uses the configured `python` executable (or `PYTHON`) with Matplotlib installed.

The second command runs the focused unit and identity checks. It does not invoke PowerDynamics or use any Experiment-A output.

The optional `S5` case is included only when its deterministic positive-diagonal/indefinite-Hermitian condition is found. The adapter file format is versioned TOML; see `EXPA_INTERFACE_SPEC.md`.
