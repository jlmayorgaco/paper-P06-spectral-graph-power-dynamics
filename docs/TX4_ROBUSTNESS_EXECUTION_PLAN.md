# TX4 Robustness Execution Plan

1. Freeze this state, preregistration, claim matrix, and deviation log in one
   commit before retaining new campaign outputs.
2. Run the exact reduced-DAE evaluator over all 16 target-family portfolios
   for deterministic sweeps, phase maps, QMC, MC, Morris, and Saltelli
   conditions. Use deterministic condition IDs and resume-safe partition files.
3. Derive confidence intervals, sensitivity, logistic, g-star, return,
   mode-family, and TDS summaries from the master condition table only.
4. Run the fixed 32-condition Julia cross-code spot-check and record the
   exact script hash, Julia version, and paired verdicts.
5. Generate the claim matrix, headline JSON, statistical appendix, final
   report, poster summary, handoff, and both ZIP bundles. Validate ZIP contents,
   CSV row counts, Parquet readability, and PDF page counts.
6. Commit all campaign artifacts on this branch, verify the frozen sibling
   worktree remains unchanged, and stop without pushing.
