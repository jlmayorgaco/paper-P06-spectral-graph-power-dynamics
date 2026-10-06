# Reproduce the current ExpE analytical work

Run from the repository root with the pinned `Project.toml`/`Manifest.toml`. The frozen domain must be created before any stability evaluation:

```powershell
julia --project=. experiments/bnd_expE/freeze_domain.jl
julia --project=. experiments/bnd_expE/audit_collective_model.jl
julia --project=. experiments/bnd_expE/evaluate_full_replacement.jl
julia --project=. experiments/bnd_expE/audit_closure.jl
julia --project=. experiments/bnd_expE/continuation_probe.jl
julia --project=. experiments/bnd_expE/one_anchor_branches.jl
$env:EXP_E_ANCHOR_BUS='38'; julia --project=. experiments/bnd_expE/run_analytic_continuation.jl
python experiments/bnd_expE/summarize_provisional.py
python experiments/bnd_expE/audit_electrical_balance.py
julia --project=. test/bnd_expE/runtests.jl
python experiments/bnd_expE/print_status.py
```

The continuation script supports `EXP_E_RESUME=1` to resume its saved branch. The current analytical result is provisional and has not been frozen for independent detailed-model validation.
