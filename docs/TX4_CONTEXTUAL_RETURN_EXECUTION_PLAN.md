# TX4 contextual-return execution plan

This plan is subordinate to `TX4_CONTEXTUAL_RETURN_PREREG.md` and is
committed before new numerics.

## Sequence

1. Verify the parent/tag/worktree and perform read-only archaeology against
   FC01, PCV04/FC13, E31, G2, and the exact port implementation.
2. Run the numerical truth audit on H4, its four one-device removals, and the
   four distinct fixed controls BASE, singleton 30, singleton 33, and pair
   30+33.
   Use tolerance factors 0.5, 1, and 2 and independent Jacobian/eigensolver
   paths; run a descriptor check only if a valid pair is exposed.
3. Recompute the full H4 g-only boundary and frozen sweep. Save alpha,
   critical mode, local factors, collective factor, contextual returns,
   determinant residuals, and equilibrium diagnostics.
4. Evaluate all 15 proper subsets at the matched root; construct the
   minimality and local-versus-collective tables.
5. Compute the contextual-return derivative at the root versus the fixed
   finite difference.
6. Run the fixed G2 D2 nonlinear phasor TDS on H4 at P4 and at g_after=0.25.
   Also run the four P4 triples only if runtime is bounded and the same rules
   can be retained.
7. Reuse/audit E31 ANDES outputs; do not rebuild a model.
8. Build the compact conventional-screen comparison, perform the focused
   2020--2026 primary-source literature audit, and write the claim hierarchy.
9. Generate the required figures, theorem note, report, poster content,
   reviewer attacks, handoff, and upload ZIP.
10. Run manifest/PDF/render/output-count checks, commit all campaign artifacts
    on `research/tx4-contextual-return-final`, and do not push.

## Required new outputs

- `results/TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv`
- `results/TX4_CONTEXTUAL_RETURN_SWEEP.csv`
- `results/TX4_PROPER_SUBSET_CLOSURE.csv`
- `results/TX4_CONTROL_BOUNDARY_REMEDIATION.csv`
- `results/TX4_FINAL_CLAIM_MATRIX.csv`
- figures in `figures/ias_final/`

## Explicit exclusions

Do not run weak-node/link ranking, structured-radius directions, repair
optimization, planning comparisons, IEEE-68, EMT, new controller searches,
PD39/PowerDynamics campaigns, or unrelated branches. Do not change grids,
thresholds, or claims after inspecting outcomes.
