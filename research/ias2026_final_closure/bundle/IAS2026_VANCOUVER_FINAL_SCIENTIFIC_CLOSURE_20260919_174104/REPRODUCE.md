# Fresh extraction reproduction

This package excludes FINAL_PAPER.pdf and FINAL_POSTER_DRAFT.pdf. Run from the extracted bundle root.

```powershell
python code/python/run_p1_canonical.py --campaign-root .
python code/python/run_p1_canonical_transfer.py --campaign-root .
python code/python/compare_p1_transfer.py --campaign-root .
python code/python/run_true_same_model_gate.py --campaign-root .
julia --project=env/julia code/julia/run_p1_gfl11.jl --campaign-root .
python code/python/build_executive_report_pdf.py --campaign-root .
```

The true same-model gate is intentionally STOPPED_BY_GATE until the exact custom synchronous-machine/AVR/PSS Julia port is completed. Alternative-model results remain model-conditioned negative evidence.
