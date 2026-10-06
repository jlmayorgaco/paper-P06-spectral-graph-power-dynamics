# Reproduction

For an independent **best-design recheck**, run from the repository root with the Julia project and Python packages recorded in `PROVENANCE.md`:

```powershell
julia --project=. experiments/latency_robust_pll_codesign_20261003/evaluate_design.jl experiments/latency_robust_pll_codesign_20261003/BEST_FOUND_DESIGN.toml
julia --project=. experiments/latency_robust_pll_codesign_20261003/run_event_local.jl experiments/latency_robust_pll_codesign_20261003/BEST_FOUND_DESIGN.toml
```

For the entire campaign, use the archived stage scripts in order: L0 reproduction and both event validations; L1 frozen probes and action validation; L2 gradients, finite differences and family scan; L3 initial predictor/corrector and each design's `evaluate_design.jl`, root-coverage discovery, and five-event correction; L4 family map; then `finalize_results.py` and `write_report.py`. `L3_OPTIMIZATION_TRACE.csv` and each `designs/*.toml` give the executed candidate order. Some freeze scripts refuse to overwrite existing artifacts: rerun the full pipeline in a **separate clean copy**, retaining `EXPERIMENT_MANIFEST.json` and the frozen probe tables. Do not regenerate preregistration after seeing results.

The original repository structure and prior experiment modules are required. `source_snapshot/` and `input_snapshot/` preserve source and data inputs captured after the runs for comparison or recovery in a separate checkout. A five-event validation takes several minutes. Rechecking with `BEST_FOUND_DESIGN.toml` writes a separate `BEST_FOUND_DESIGN` output directory under this experiment.
