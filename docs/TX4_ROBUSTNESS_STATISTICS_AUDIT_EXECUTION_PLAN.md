# TX4 Robustness Statistics Audit — Execution Plan

1. Freeze and document the source campaign state.
2. Implement blocker-set logic and run six synthetic unit cases, including
   empty U0, singleton H4, a proper singleton, a composite blocker, mixed
   minimal blockers, and a noncomposable failure case.
3. Audit legacy H4 flags, NONCOMPOSABLE handling, eta semantics, and row-level
   evaluator provenance without rerunning the campaign.
4. Compute exact all-16 blocker truth for the existing exact calibration rows;
   validate the reconstructed ExtraTrees tier with grouped held-out metrics.
5. Apply the strict surrogate gate. On failure, run checkpointed exact QMC and
   MC all-16 recomputation using the fixed coordinates.
6. Build corrected primary statistics, delta_H4 and eta diagnostics, exact
   provenance, invariant checks, figures, report, handoff manifest, and
   reproducibility/upload archives.
7. Stop only when every invariant is PASS; if any final invariant fails, do
   not claim completion.

