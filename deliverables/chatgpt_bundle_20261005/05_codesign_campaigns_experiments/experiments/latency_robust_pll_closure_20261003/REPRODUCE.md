# Reproducing the latency-robust PLL closure

The bundle preserves the repository layout required by the Julia scripts. It
includes the closure experiment, its frozen parent, model-source snapshots,
the required IEEE-39 report inputs, and the repository's `Project.toml` and
`Manifest.toml`. Julia dependencies must already be available or instantiated
from the manifest.

From the extracted repository root:

```powershell
python experiments/latency_robust_pll_closure_20261003/audit_parent.py
julia --project=. experiments/latency_robust_pll_closure_20261003/revalidate_gain_rank.jl
julia --project=. experiments/latency_robust_pll_closure_20261003/map_closure_families.jl
python experiments/latency_robust_pll_closure_20261003/finalize_closure.py
python experiments/latency_robust_pll_closure_20261003/audit_predictor_correction.py
python experiments/latency_robust_pll_closure_20261003/write_final_report.py
```

For a specific stored design ID, the exact characteristic and event reruns are:

```powershell
julia --project=. experiments/latency_robust_pll_closure_20261003/evaluate_design.jl experiments/latency_robust_pll_closure_20261003/designs/<ID>.toml
julia --project=. experiments/latency_robust_pll_closure_20261003/discover_boundary_roots.jl experiments/latency_robust_pll_closure_20261003/designs/<ID>.toml
python experiments/latency_robust_pll_closure_20261003/augment_root_catalogue.py <ID>
python experiments/latency_robust_pll_closure_20261003/correct_spectral_from_coverage.py <ID>
julia --project=. experiments/latency_robust_pll_closure_20261003/run_event_local.jl experiments/latency_robust_pll_closure_20261003/designs/<ID>.toml
julia --project=. experiments/latency_robust_pll_closure_20261003/export_linear_dde.jl experiments/latency_robust_pll_closure_20261003/designs/<ID>.toml
python experiments/latency_robust_pll_closure_20261003/linear_dde_method_steps.py <ID>
```

The last two commands run a **linear** positive-delay DDE realization. They
do not reproduce nonlinear delayed events. All five physical event scripts
run at zero added PLL measurement delay. The model's intrinsic PLL low-pass
filter is retained in both calculations.

The optimization predictor uses measured actuator sensitivities from stored
one-event probes and an exact-root local gradient. Its outputs are proposal
files; the subsequent characteristic and physical-event checks determine
acceptance. Existing frozen parent results are never rewritten.

`TABLE_F5_MULTISTART.csv`, `TABLE_F6_GENERIC_OPTIMIZER_COMPARISON.csv`, and
`TABLE_F7_PARETO_FRONTIER.csv` explicitly mark work that was not executed.
Their presence is not evidence of a completed comparison.
