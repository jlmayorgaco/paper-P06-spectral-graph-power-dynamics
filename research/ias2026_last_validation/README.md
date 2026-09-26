# IAS 2026 Vancouver final scientific closure

This is the execution root for the final-closure campaign requested by the
IAS2026 master prompt. It is intentionally separate from the earlier
`research/ias2026_bulletproof/` corrective campaign. The earlier campaign is
preserved as provenance; this root is the adjudication layer.

## Decision state

The package is scientifically useful but does not clear the award-level claim
gates. The canonical GFL11 device/port parity is numerically verified, the
official PowerDynamics and SimpleGFLDC campaigns are valid alternative-model
negative results, and the official-model TDS solves are reproducible. The true
same-model frozen IEEE-39 cross-code gate is `STOPPED_BY_GATE` because the exact
Julia port of the frozen custom SynchronousMachine/AVR/PSS model was not
completed. Consequently the collective mechanism, common physical robustness,
new independent holdout, and scaling gates are not promoted.

No final paper or poster rewrite is authorized by this package. The only PDF
created here is the executive audit report, and it is not a submission paper.

## Layout

- `raw/`: immutable copied evidence, separated into canonical, same-model,
  second-model, TDS, holdout, robustness, and scaling scopes.
- `code/`: reproducibility scripts for the campaign root.
- `reports/`: machine-readable dashboard and gate reports.
- `docs/`: claim ledger, deviations, model reconciliation, negative results,
  preregistration, and reviewer attack matrix.
- `derived/tables/`: derived tables copied from the corrective campaign.
- `figures/`: evidence figures only; no paper/poster PDFs.
- `bundle/`: generated ZIP packages.

## Reproduction

Run from this directory:

```powershell
python code/python/run_p1_canonical.py --campaign-root .
python code/python/run_p1_canonical_transfer.py --campaign-root .
python code/python/compare_p1_transfer.py --campaign-root .
python code/python/run_true_same_model_gate.py --campaign-root .
julia --project=env/julia code/julia/run_p1_gfl11.jl --campaign-root .
python code/python/build_executive_report_pdf.py --campaign-root .
python code/python/package_final_closure.py --campaign-root .
```

The status vocabulary is restricted to the values listed in
`docs/STATUS_VOCABULARY.md`. Every stopped, negative, single-class, or
not-tested result is retained rather than converted into a positive claim.
