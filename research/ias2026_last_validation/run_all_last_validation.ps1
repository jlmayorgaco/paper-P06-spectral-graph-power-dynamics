$ErrorActionPreference = 'Stop'
$campaign = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $campaign
try {
    python code/python/run_true_same_model_oracle.py --campaign-root .
    julia --project=env/julia code/julia/run_true_same_model_ieee39.jl --campaign-root .
    python code/python/reconcile_true_same_model.py --campaign-root .
    python code/python/run_true_same_model_gate.py --campaign-root .
    julia --project=env/julia code/julia/run_true_same_model_mechanism.jl --campaign-root .
    julia --project=env/julia code/julia/run_true_same_model_tds.jl --campaign-root .
    julia --project=env/julia code/julia/run_p5_simplegfldc_ieee39.jl --campaign-root .
    julia --project=env/julia code/julia/run_second_model_discovery.jl --campaign-root .
    python code/python/freeze_mixed_holdout.py --campaign-root .
    python code/python/package_last_validation.py --campaign-root .
}
finally {
    Pop-Location
}
