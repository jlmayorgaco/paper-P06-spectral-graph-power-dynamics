# IAS 2026 Bulletproof Closure Campaign

This directory is the new, post-freeze execution layer for the IAS 2026 closure campaign. It is intentionally separate from the historical TX4/CDW and poster result directories.

## Provenance

- Canonical base tag: `IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE`
- Canonical base commit: `b9f274e241a793bd2694d488f6e53c9aca6e6ac5`
- Parent TX4 freeze: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Campaign branch: `research/ias2026-bulletproof-closure-v1`
- Historical evidence is read-only input. New outputs belong here.

## Operating rule

This campaign audits and falsifies claims. A missing package, failed gate, unsupported uncertainty structure, or negative holdout is recorded as evidence and stops only dependent phases. The phrase “full-order reference within the frozen model” is used instead of “ground truth”.

## Status

Corrective-audit status: corrected Gate 0 reproduces the nominal frozen
IEEE-39 P4 result and separates the global V9 census (402/512) from the
targeted 0.3–1.5 Hz census (511/512); property tests pass; and the official
PowerDynamics Gate A passes. Same-model GFL parity, robustness, a second model,
new blind holdout, and new Julia TDS remain incomplete. No final paper/poster
bundle is authorized yet. The dashboard, claim ledger, and audit report are in
`reports/` and `docs/`.

To reproduce the core gates, run `python research/ias2026_bulletproof/src/python/run_all_python.py --phase gate0` and `python research/ias2026_bulletproof/src/python/run_all_python.py --phase theory-tests`. Run Gate A with `julia --project=research/ias2026_bulletproof/env/julia research/ias2026_bulletproof/src/julia/run_pd39_gate.jl`.

The current status dashboard is written to `reports/GATE_DASHBOARD.md` after each orchestrated phase.

Execution is isolated from the primary dirty checkout.
