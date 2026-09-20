# Repository Inventory and Freeze Record

This inventory records the execution checkout, the historical evidence freeze,
and the environment used for the closure campaign.

## Execution checkout

- Repository: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics`
- Campaign worktree: `C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics-ias2026-bulletproof`
- Branch: `research/ias2026-bulletproof-closure-v1`
- Starting HEAD: `b9f274e241a793bd2694d488f6e53c9aca6e6ac5`
- Starting tag: `IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE`
- User checkout: preserved; its pre-existing dirty state was not reset, cleaned,
  rebased, or overwritten.
- Campaign policy: no post-hoc parameter retuning; all negative and stopped
  gates are retained.

## Relevant historical anchors

- `IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE` -> `b9f274e...`
- `IAS2026_PRE_FINAL_VALIDATION` -> `5b536729...`
- `TX4_FINAL_MANUSCRIPT_FREEZE` -> `f003ee2b...`
- TX4 parent evidence commit: `69f200dfe1dfb74cc4f678ad25a0c3b1751d62a6`
- Historical blind worktree: `C:\w\tx4blind`

The campaign reads the frozen final-closure tables and the historical blind V9
reveal. Those source paths are recorded in each derived table and in
`derived/manifests/gate0_manifest.json`.

## Runtime

- Windows 11 Pro, build `10.0.26200`
- Python `3.13.14`
- Campaign runtime imports: NumPy `2.3.5`, SciPy `1.16.3`, pandas `2.3.3`,
  matplotlib `3.10.8`, pytest `8.4.2`
- Frozen-record runtime: Python `3.13.14`, NumPy `2.5.2`, SciPy `1.18.1`,
  pandas `3.0.5`, ANDES `2.0.0`
- Julia `1.11.9`
- Isolated Julia project: `env/julia/Project.toml`
- PowerDynamics `5.0.0`
- Hardware observed: Intel Core Ultra 9 185H; approximately 31.5 GB visible
  RAM; 22 logical processors.

The Python package mismatch between the frozen record and the present runtime
is documented as a deviation. Gate 0 compares exported numerical artifacts;
it does not silently relabel this runtime as the original environment.

## Execution provenance

- Gate 0: `reports/GATE0_REPORT.md`
- Theory tests: `logs/theory_tests.log`
- PowerDynamics gate: `raw/powerdynamics/pd39_equilibrium_gate.md`
- Dashboard: `reports/GATE_DASHBOARD.md` and `.csv`
- Final claim ledger: `docs/CLAIM_LEDGER.md` and
  `derived/tables/TABLE_CLAIM_LEDGER.csv`
