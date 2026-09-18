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

Executed closure status: Gate 0 reproduces the nominal frozen IEEE-39 P4
result, theorem tests pass, the official PowerDynamics IEEE-39 equilibrium gate
passes, and the stronger V9 transfer expectation is refuted. The dashboard,
claim ledger, and final audit are in `reports/` and `docs/`.

To reproduce the core gates, run `python research/ias2026_bulletproof/src/python/run_all_python.py --phase gate0` and `python research/ias2026_bulletproof/src/python/run_all_python.py --phase theory-tests`. Julia phases use the isolated project in `research/ias2026_bulletproof/env/julia`.

The current status dashboard is written to `reports/GATE_DASHBOARD.md` after each orchestrated phase.

Execution is isolated from the primary dirty checkout.
