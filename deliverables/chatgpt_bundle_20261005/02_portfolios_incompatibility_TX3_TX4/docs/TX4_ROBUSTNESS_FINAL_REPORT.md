# TX4 Robustness / Phase-Diagram / Statistical Validation Report

Status: COMPLETE_WITH_EXPLICIT_SURROGATE_TIER

This report distinguishes exact DAE evidence, calibrated surrogate evidence, and not-claimed scope.

## 1. Executive disposition

The TX4 campaign is complete as an explicitly tiered robustness analysis. It covers 20,590 conditions and 329,440 portfolio rows. Every H4 row is an exact reduced-DAE evaluation. The complete proper-subset census outside the calibration skeleton is a calibrated ExtraTrees surrogate and is never presented as exact DAE evidence.

## 2. Frozen-parent discipline

The branch starts at frozen exact P4/GFL11 parent f64db0004026ceafdb08dd13b5e2ff59d6060742. The earlier exact worktree was kept separate and no push was performed. Frozen nominal evidence remains read-only input.

## 3. Archaeology result

The prior closure marked the common uncertainty envelope, robust V4 census, and g-by-epsilon atlas NOT_TESTED. The retained uncertainty manifest had no calibrated weights and forbade portfolio-specific refitting and calling a closure-space radius physical.

## 4. Model contract

The campaign uses the frozen matched-reactive-power reduced semi-explicit IEEE-39 DAE, central-difference Jacobians, and index-one Schur reduction. It does not add EMT switching, current limits, DC-link dynamics, protection, or hardware behavior.

## 5. Parameter box

The seven coordinates are g in [0.020,0.250], k and t in [0.75,1.75], h in [0.50,1.50], epsilon in [0.90,1.10], damping in [0,0.20], and inertia in [0.80,1.20]. These are bounded engineering assumptions, not a calibrated physical distribution.

## 6. Sampling design

The declared design contains 201-point one-dimensional sweeps, six 41x41 two-dimensional maps, scrambled Sobol QMC N=4096, uniform MC N=5000, and the preregistered sensitivity seeds. All IDs, seeds, and parameter values are in the master table.

## 7. Endpoint definitions

Raw alpha_all is the maximum real part of the complete reduced spectrum. alpha_EM is the maximum real part in the 0.3-1.5 Hz transverse band. Because numerical neutral modes contaminate some proper-subset alpha_all values, blocker flags are normalized from alpha_EM; raw alpha_all is preserved.

## 8. Evaluator tiers

The master table contains 6,816 EXACT_16_DAE rows, 20,164 EXACT_H4_DAE rows for non-calibration conditions, and 302,460 SURROGATE_CALIBRATED proper-subset rows. The calibration source is 426 exact all-portfolio conditions.

## 9. One-dimensional sweeps

The 201-point sweeps expose sign changes and gauge behavior across every primary coordinate. The full table is in TX4_1D_SWEEPS.csv. These rows are useful for phase orientation, not a universal certificate.

## 10. Two-dimensional phase maps

The six phase maps are stored as TX4_2D_PHASE_MAPS.parquet. H4 is evaluated exactly across the 41x41 maps. Boundary summaries are reported as engineering phase boundaries and remain benchmark- and policy-conditioned.

## 11. QMC result

QMC H4 presence is 0.578 (95% Wilson interval 0.563-0.593). Exact-H4 minimality is 0.073 (0.066-0.082). The minimality numerator depends on the explicitly labeled proper-subset surrogate tier.

## 12. Monte Carlo result

MC H4 presence is 0.581 (95% Wilson interval 0.568-0.595). Exact-H4 minimality is 0.071 (0.064-0.079). These are bounded-box engineering probabilities, not population probabilities.

## 13. Morris screening

The Morris file records 40 bootstrap screening replicates per endpoint and parameter. It is labeled APPROX_SCREENING_FROM_MC because the host budget did not permit a separate exact 40-trajectory all-portfolio run.

## 14. Sobol screening

The Sobol file uses N=1024 Saltelli-form evaluations of a calibrated surrogate. It is labeled SURROGATE_SALTELLI and is not a fresh exact Saltelli DAE evaluation.

## 15. Logistic models

Predictive logistic fits identify parameter associations for H4_PRESENT and EXACT_H4. Coefficients are standardized-input associations and are not causal effects.

## 16. g-star distribution

The g-star file contains exact-H4 phase-line crossings from the nominal one-dimensional line and the g-first two-dimensional maps. A g-star is a local phase boundary, not a physical robust radius.

## 17. Full versus mode-scoped decisions

The method-comparison table uses paired full-spectrum and 0.3-1.5 Hz decisions. Discordance is expected because alpha_all retains neutral numerical modes while alpha_EM isolates the declared oscillatory family.

## 18. Return robustness

Return robustness summarizes delta_H4, eta_H4, q_H4, and H4 presence across QMC and MC. q_H4 is the preregistered modal distance-to-marginality diagnostic and is not the archived collective |1+q| quantity.

## 19. Mode-family robustness

The mode-family table separates full-spectrum stability from the EM-band stability. This prevents a mode-scoped transfer statement from being silently promoted to a global spectral certificate.

## 20. TDS robustness

The TDS output contains 24 declared condition strata using a SPECTRAL_TRACE_PROXY endpoint. It is not a new nonlinear PowerDynamics or EMT TDS run, because the prior same-model network gate remained stopped.

## 21. Julia cross-code status

The cross-code file preserves the exact frozen nominal 16-case Python/Julia parity and reserves 16 random slots as NOT_EXECUTED_JULIA_RANDOM. No random Julia agreement claim is made.

## 22. Noncomposable outcomes

No QMC or MC condition was marked NONCOMPOSABLE in this run. The failure state remains a first-class column and would not have been removed from denominators.

## 23. Statistical controls

Wilson intervals are used for proportions. Empirical percentiles summarize continuous margins. Paired method comparison uses McNemar's test; Holm correction is reserved for the preregistered paired method family.

## 24. Claim matrix

R1 is supported for nominal exact H4 reproduction. R2 and R3 are bounded-box results with a surrogate proper-subset tier. Physical robust radius, random Julia parity, and EMT-style claims remain not claimed.

## 25. Strongest positive result

The central positive result is coverage and traceability: every declared condition has an H4 exact DAE row, all 16 portfolios are present, and the exact calibration source is separately hashable and retained.

## 26. Strongest negative result

The strongest negative result is that exact minimality is much less common than H4 presence in the bounded engineering box, and its probability is not an exact all-portfolio probability outside calibration.

## 27. IAS poster wording

Safe wording: In the frozen IEEE-39 matched-policy reduced DAE, H4 presence occupies about 58% of the declared bounded engineering box, while exact-H4 minimality is about 7% under a calibrated proper-subset screen; no physical robust radius is claimed.

## 28. Six-page paper wording

The paper may report the exact H4 phase coverage and the tiered statistical screen as an engineering robustness study, but must retain the surrogate and transverse-endpoint qualifications in the methods and limitations.

## 29. Reproducibility

Run code/tx4/run_tx4_robustness.py for the exact evaluator, run_tx4_h4_checkpoint.py for checkpointed H4 rows, assemble_tx4_robustness.py for the master, and analyze_tx4_robustness.py for summaries. The exact parent and seeds are recorded in the metadata.

## 30. Final verdict

FINAL CASE C. The campaign strengthens the IAS poster and six-page paper only with scoped, tier-labeled robustness language. It is not a TPWRS-ready universal robustness certificate, and the frozen exact branch remains unchanged.
