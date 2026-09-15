# TX4 blind-prediction execution plan

1. Freeze this worktree from `69f200df`; complete read-only archaeology.
2. Commit current state, preregistration, claim matrix, deviations, and
   firewall before any prediction/evaluation.
3. Build only the baseline `T0/K` and local Delta-Y library, then enumerate and
   commit immutable V4 blind predictions. Hash the prediction files.
4. Reveal the frozen V4 full-order answer key and compare; do not tune settings.
5. Using the same frozen library, enumerate and commit all 512 V9 predictions;
   then reveal and compare to the frozen V9 answer key.
6. Measure full versus reduced runtime using pinned BLAS threads.
7. Recompute the g-only reduced boundary blindly, then reveal the full-order
   boundary and report both errors.
8. Validate contextual returns, minimality separation, local/collective
   regularity, proper subsets, aggregate-matched pairs, and screen comparison.
9. Run fixed g-only remediation and the frozen nonlinear phasor TDS.
10. Reuse the E31 ANDES validation, perform the focused literature audit, and
    write the claim hierarchy, figures, reviewers, report, poster content, and
    handoff.
11. Validate PDFs, manifests, hashes, upload bundles, and commit all artifacts
    on this branch. Never push or start another campaign.

Required top-level result files include the V4/V9 blind and reveal tables,
g-boundary and contextual-return tables, minimality/locality/screen/runtime
tables, final claim matrix, headline JSON, reports, figures, and upload ZIPs.
