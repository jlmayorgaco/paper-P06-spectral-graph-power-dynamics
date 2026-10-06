# Experiment C reproduction

From the repository root, run:

```powershell
julia --project=. experiments/bnd_expC/run_experiment_C.jl
```

The run reads the Experiment-A bus-33 `A_reduced` matrix/state partition and Experiment-A/B machine-readable status artifacts without modifying them. It rebuilds the all-SG baseline and bus 30/33/35/37 cases with the pinned PowerDynamics 5.0.0 model and frozen nominal SimpleGFLDC template on the first run. It retains those exact reduced matrices/state maps in `reports/experiment_C/matrices/`, keyed by hashes of the project, manifest, model adapter, and ExpA inputs. Later runs reuse only a matching cache; pass `--rebuild-models` to force a fresh PowerDynamics extraction. No tuning, package update, or network push occurs.

The runner regenerates CSV tables, JSON results, Markdown report and tables, matrices, and figures under `reports/experiment_C/`. Julia and package versions, BLAS configuration, Git HEAD, source artifact hashes, timestamp, and the bus-33 reconstruction residual are recorded.

Focused unit checks are in `test/bnd_expC/runtests.jl`; ExpB checks remain in `test/bnd_expB/runtests.jl`.

