# Reproduce Experiment F0

From the repository root, run:

```powershell
julia --project=. --startup-file=no experiments/bnd_expF0/run_experiment_F0.jl
python experiments/bnd_expF0/run_mode_overlap_F0.py
python experiments/bnd_expF0/finalize_results_F0.py
python experiments/bnd_expF0/render_figures_F0.py
```

Inputs are the frozen IEEE-39 CSV copies in `reports/experiment_D/inputs/`.
The analytic code assembles the passive branch Ybus directly and does not
import PowerDynamics. The reports include input hashes, the Kron identity
residual, graph spectrum, candidate classifications, and the ExpC mode
projection diagnostic when the frozen q-state map is compatible.
