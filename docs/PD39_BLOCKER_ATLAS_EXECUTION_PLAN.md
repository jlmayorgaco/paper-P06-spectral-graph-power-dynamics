# PD39 Blocker Atlas v1 — execution plan

1. Verify branch, parent HEAD, runtime and authoritative input hashes; commit
   the current-state note and this preregistration package.
2. Parse the complete holdout census and write the deterministic three-blocker
   selection table before any new solve.
3. Run Gate A on C12 and nominal for selected blockers, all immediate
   predecessors, all cheap proper subsets, and matched controls. Preserve all
   tolerance, descriptor, independent-Jacobian and conditioning rows.
4. Revalidate minimality separation and classify each mode using the frozen
   physicality rules.
5. Run Gate B with the common 1% P/Q, 100 ms bus-39 disturbance. Store table
   and raw traces.
6. Apply Gate C exactly; issue the fixed GO/PARTIAL GO/STOP decision.
7. If GO only, run the two-coordinate compatibility atlas and boundary
   continuation, then the exact-return audit/theory, diagnostics, and one-shot
   remediation proof of concept. If not GO, do not run these numerics.
8. Generate figures, final report, handoff, machine-readable claims and the
   compact upload and reproducibility bundles.
9. Verify required files, commit, check clean branch, and print the final
   CASE A/B/C line. Never push.

Expected primary output files include:

- `results/PD39_SELECTED_BLOCKERS.csv`
- `results/PD39_BLOCKER_NUMERICAL_AUDIT.csv`
- `results/PD39_BLOCKER_TDS_VALIDATION.csv`
- `results/PD39_COMPATIBILITY_ATLAS.parquet` when and only when GO
- `results/PD39_BLOCKER_BOUNDARIES.csv` when and only when GO
- `results/TX4_CONTEXTUAL_RETURN_VALIDATION.csv`
- `results/PD39_BLOCKER_ATLAS_HEADLINE.json`
- `results/PD39_BLOCKER_ATLAS_MASTER_CLAIMS.csv`
- `results/PD39_BLOCKER_ATLAS_MASTER_CASES.csv`
