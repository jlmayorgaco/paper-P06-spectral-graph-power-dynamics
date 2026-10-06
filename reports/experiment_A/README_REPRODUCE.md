# Experiment A reproduction

Run from the repository root in the existing Julia environment:

```powershell
julia --project=. experiments/bnd_expA/run_experiment_A.jl
```

The runner uses the PowerDynamics.jl 5.0.0 and Julia 1.11.9 versions pinned by
the repository `Manifest.toml`. It does not update or install packages. Julia
writes the CSV tables, matrices, JSON results, and report. The runner then
creates `TABLES_EXP_A.md` from the CSV tables with `render_markdown_tables.py`
and the figures with `plot_experiment_A.py` using the Python `matplotlib`
available on `PATH`.

The lightweight algebra/indexing checks requested for this experiment run with:

```powershell
julia --project=. experiments/bnd_expA/test/runtests.jl
```

No random sampling is used. Outputs are written under `reports/experiment_A/`.
The raw descriptor matrices and eigenvectors are retained under
`reports/experiment_A/matrices/`; the machine-readable tables are under
`reports/experiment_A/tables/`.
The complete named equilibrium vectors are in `reports/experiment_A/matrices/`.

The raw repository working tree already contains unrelated local changes.
This experiment writes only to its new `experiments/bnd_expA/` and
`reports/experiment_A/` paths, and does not commit or push.
