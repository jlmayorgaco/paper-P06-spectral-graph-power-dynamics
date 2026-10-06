# Reproduce Experiment F1

From the repository root, run:

```powershell
julia --project=. --startup-file=no experiments/bnd_expF1/run_experiment_F1.jl
python experiments/bnd_expF1/render_figures_F1.py
```

This is a deterministic synthetic second-order reference. It uses no GFL
model and makes no claim about globally optimal noncommuting damping.
