# Reproduce Experiment F2

From the repository root, run:

```powershell
julia --project=. --startup-file=no experiments/bnd_expF2/run_experiment_F2.jl
python experiments/bnd_expF2/finalize_results_F2.py
python experiments/bnd_expF2/render_figures_F2.py
```

The analytic assembly includes the audited nine-state SimpleGFLDC equations
and the frozen IEEE-39 branch/load data. It imports the independent analytic
ExpE model, never PowerDynamics. The F2 tolerance and basis were preregistered
in `PREREGISTRATION.md`; a failed modal gate intentionally stops F3/F4 spectral
gain-law work.
