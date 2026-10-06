# Latency-robust PLL co-design at fixed IEEE-39 replacement

This self-contained **experiment output** reuses the frozen repository model; see `PROVENANCE.md` and `EXPERIMENT_MANIFEST.json`. At fixed 88.455141% GFL replacement, it tunes ten independent PLL proportional and integral gains against a common, exogenous delay on each PLL measurement/error path. No files from prior experiments were edited.

The best fully validated design **found in the executed search** is `N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp`. Its numerical local uniform-delay threshold is **40.844904 ms**, compared with **39.383359 ms** for nominal N and **37.387188 ms** for zero-delay-tuned Z. It passes the full zero-delay physical spectrum and all five frozen nonlinear 61-second events. The corresponding `rho`, `Kp`, and `Ki` vectors are in `BEST_FOUND_DESIGN.toml`.

Read `STATUS.md` for claim gates, `THEORY_LATENCY_ROBUST_CODESIGN.md` for the exact rank-10 action and simple-root sensitivity, `L3_OPTIMIZATION_TRACE.csv` for accepted and rejected steps, `L4_MODAL_FAMILY_MAP.csv` for modal identity, and `POSTER_CLAIM_LEDGER.md` for the allowed claims. `figures/` contains diagnostic plots generated from the CSVs.

The word “exact” describes the exponential DDE characteristic and algebraic gain action. Numerical roots, contour counts, and optimization are Float64 calculations. The threshold is a **numerically observed local first crossing**, without an interval certificate of all characteristic roots or a global optimality proof. The nonlinear events use **zero delay**; the experiment does not establish nonlinear safety under positive delay.
