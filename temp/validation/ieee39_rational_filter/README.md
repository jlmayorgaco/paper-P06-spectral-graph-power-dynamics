# IEEE 39 Validation Protocol for Rational Graph-Filter Damping Estimator

This folder is the validation harness for the long paper version. It is not a
completed validation result. It defines what must be run once ANDES and a dynamic
IEEE 39-bus case are available.

## Current Environment Status

Checked on 2026-06-08:

- ANDES 2.0.0 is available in the active Python environment used by this
  workspace.
- ANDES includes `ieee39_full.xlsx`, `ieee39.xlsx`, and `ieee39.raw`.
- The packaged `ieee39_full.xlsx` case has 10 GENROU machines and no active
  REGCA1/REGCP1/REGF1/REGF2/PLL1 devices.
- A case generator now creates derived IEEE 39 cases with active
  `REGCP1`/`PLL1` grid-following and `REGF1` grid-forming devices.

Therefore no final full-IBR IEEE 39 estimator claim can be made yet. A
synchronous ANDES baseline, an IEEE39-derived surrogate experiment, and first
full-ANDES IBR replacement pilot cases have now been run.

## Final Claim-Readiness Audit

Script:

```text
validation/ieee39_rational_filter/final_claim_audit.py
```

Command:

```text
python validation/ieee39_rational_filter/final_claim_audit.py
```

Outputs:

```text
outputs/final_claim_audit/claim_audit.csv
outputs/final_claim_audit/claim_audit.json
outputs/final_claim_audit/final_readiness_report.md
```

Interpretation:

- The paper is ready as a rigorous methodological/master-thesis package.
- The second-order correction, rational graph-filter interpretation, reduced
  PLL check, IEEE39-derived surrogate, and ANDES case-construction workflow are
  documented.
- A final IEEE Transactions empirical claim is still blocked until calibrated
  ANDES IBR cases are compared against full eigensolve baselines.

## Phase E0 Inertia-Control Cross-Term Gate

Script:

```text
validation/ieee39_rational_filter/phase_cross_inertia_control.py
```

Representative commands:

```text
python validation/ieee39_rational_filter/phase_cross_inertia_control.py --mode network --eps 0.01
python validation/ieee39_rational_filter/phase_cross_inertia_control.py --mode network --eps 0.01 --sweep-singletons --out outputs/phase_cross_inertia_control_sweep_network
```

Outputs:

```text
outputs/phase_cross_inertia_control/phase_cross_inertia_control_status.json
outputs/phase_cross_inertia_control/phase_cross_inertia_control_report.md
outputs/phase_cross_inertia_control_sweep_network/sweep_summary.csv
outputs/phase_cross_inertia_control_sweep_network/phase_cross_inertia_control_status.json
```

Current result:

- Mix60/no-PSS network-family all-active test: `NOT_MEASURABLE`.
- Singleton GENROU.M x PLL1.Kp/Ki sweep: `0/13` measurable cross cases.
- The strongest cross fraction was the all-active case, `2.23882e-04` of the
  true joint pole shift.

Interpretation:

- The inertia-control cross term is mathematically present in the local
  expansion, but it is not an empirical headline in this initial ANDES gate.
- Keep it as a theoretical local interaction unless a physically justified
  virtual-inertia/control case makes it measurable.
- Continue to Experiment A for the estimator advantage gate.

## Completed Preliminary Run

Script:

```text
validation/ieee39_rational_filter/run_ieee39_surrogate_experiments.py
```

Command:

```text
python validation/ieee39_rational_filter/run_ieee39_surrogate_experiments.py --trials 60 --out outputs/ieee39_rational_filter
```

Outputs:

```text
outputs/ieee39_rational_filter/andes_baseline_summary.json
outputs/ieee39_rational_filter/surrogate_summary.json
outputs/ieee39_rational_filter/estimator_table.csv
outputs/ieee39_rational_filter/weak_node_table.csv
outputs/ieee39_rational_filter/weak_link_table.csv
outputs/ieee39_rational_filter/network_metadata.json
outputs/ieee39_rational_filter/fig_estimator_errors_by_rho.pdf
outputs/ieee39_rational_filter/fig_line_frequency_vs_damping.pdf
```

Key preliminary results:

- ANDES baseline: 39 buses, 46 lines, 10 GENROU, 10 IEEEST/PSS, no active IBR
  devices.
- Full ANDES baseline least-damped stable oscillatory pole:
  `-1.3460 +/- j8.6115`, damping ratio `0.1544`, frequency `1.3706 Hz`.
- IEEE39-derived Kron-reduced surrogate, 300 scenarios:
  - diagonal median/p95 relative error: `0.468% / 24.15%`,
  - second-order median/p95: `0.032% / 9.91%`,
  - reduced QEP r=6 median/p95: `0.037% / 5.44%`,
  - adaptive median/p95: `0.031% / 7.43%`.
- Representative 60% surrogate scenario:
  - best weak-node damping intervention: bus 30, `Delta zeta=0.0112`,
  - best line reinforcement: line 25--37, `Delta zeta=0.00631`,
  - 35 of 46 line reinforcements satisfy the surrogate Braess-like diagnostic
    `Delta nu_c > 0` and `Delta zeta_min < 0`.

These are IEEE39-derived surrogate results, not full ANDES IBR results.

## Packaged IEEE39 Positive-Mode Diagnostic

Script:

```text
validation/ieee39_rational_filter/diagnose_ieee39_positive_modes.py
```

Command:

```text
python validation/ieee39_rational_filter/diagnose_ieee39_positive_modes.py --out outputs/ieee39_mode_diagnostics
```

Outputs:

```text
outputs/ieee39_mode_diagnostics/baseline_positive_mode_diagnostics.csv
outputs/ieee39_mode_diagnostics/positive_mode_participation.csv
```

Result:

| Variant | Removed model sheets | Positive modes | Critical stable oscillatory mode |
| --- | --- | ---: | --- |
| base | none | 9 | zeta 0.1544, 1.3706 Hz |
| no_tgov | `TGOV1N` | 4 | zeta 0.0568, 0.6224 Hz |
| no_pss | `IEEEST` | 0 | zeta 0.1545, 1.3705 Hz |
| no_exciter_pss | `IEEEX1`, `IEEEST` | 0 | zeta 0.1538, 1.3725 Hz |
| genrou_only | `TGOV1N`, `IEEEX1`, `IEEEST`, `ACEc`, `Toggler` | 0 | zeta 0.0606, 0.6285 Hz |

Interpretation:

- The positive-real modes in the packaged `ieee39_full.xlsx` case are
  non-oscillatory.
- Participation factors concentrate in `TGOV1N`, `IEEEX1`/`IEEEST`, and
  internal `GENROU` transient-voltage states, not in angle-frequency
  electromechanical pairs.
- Removing only `IEEEST` removes all positive-real modes while leaving the
  1.37 Hz least-damped stable oscillatory pair essentially unchanged.
- Therefore the packaged case is useful for reproducible extraction tests, but
  the positive modes should be treated as a model-quality warning rather than
  as a physical instability claim.

## Full ANDES IBR Replacement Pilot

Script:

```text
validation/ieee39_rational_filter/create_ieee39_ibr_cases.py
```

Command:

```text
python validation/ieee39_rational_filter/create_ieee39_ibr_cases.py --out-cases validation/ieee39_rational_filter/cases --out-results outputs/ieee39_ibr_cases
```

Generated cases:

```text
validation/ieee39_rational_filter/cases/ieee39_ibr_gfl20.xlsx
validation/ieee39_rational_filter/cases/ieee39_ibr_gfm20.xlsx
validation/ieee39_rational_filter/cases/ieee39_ibr_mix20.xlsx
validation/ieee39_rational_filter/cases/ieee39_ibr_mix40.xlsx
validation/ieee39_rational_filter/cases/ieee39_ibr_mix60.xlsx
```

Outputs:

```text
outputs/ieee39_ibr_cases/case_generation_summary.json
outputs/ieee39_ibr_cases/case_generation_summary.csv
```

Construction rule:

- retain `PV` and `Slack` devices to preserve the original power-flow dispatch,
- keep slack generator 10 synchronous,
- remove `GENROU`, `TGOV1N`, `IEEEX1`, and `IEEEST` rows attached to converted
  machines,
- add GFL dynamic replacements as `REGCP1 + REECA1 + REPCA1 + PLL1`,
- add GFM dynamic replacements as `REGF1`.

Pilot results:

| Case | GFL gens | GFM gens | PFlow/EIG | full DAE stable? | stable oscillatory zeta |
| --- | --- | --- | --- | --- | --- |
| `ieee39_ibr_gfl20` | 9, 8 | none | yes | no, 1 positive real mode | 0.1148 |
| `ieee39_ibr_gfm20` | none | 9, 8 | yes | no, 6 positive modes, 2 oscillatory | 0.1438 |
| `ieee39_ibr_mix20` | 9 | 8 | yes | no, 3 positive real modes | 0.1370 |
| `ieee39_ibr_mix40` | 9, 7 | 8, 6 | yes | no, 2 positive real modes | 0.1266 |
| `ieee39_ibr_mix60` | 9, 7, 5 | 8, 6, 4 | yes | yes under `1e-7` positive-real tolerance | 0.1236 |

Interpretation:

- These cases prove that IEEE39 can be modified inside ANDES to replace selected
  SG dynamics with GFL/GFM inverter dynamics.
- They do not yet prove the paper's estimator. Four cases are unstable under
  uncalibrated generic controller parameters, and the stable mixed 60% case
  still needs a physically matched surrogate extraction before estimator error
  can be reported.
- The next technical step is controller calibration plus extraction of the
  reduced `(M,D,L)`/control-coupling surrogate from the full ANDES state matrix
  `EIG.As`.

## Phase-0 Bridge Artifacts and IEEE39 IBR Diagram

Script:

```text
validation/ieee39_rational_filter/phase0_bridge_and_diagram.py
```

Command:

```text
python validation/ieee39_rational_filter/phase0_bridge_and_diagram.py
```

Outputs:

```text
outputs/ieee39_phase0_bridge/phase0_bridge_status.json
outputs/ieee39_phase0_bridge/phase0_gate_table.csv
outputs/ieee39_phase0_bridge/phase0_bridge_report.md
outputs/ieee39_phase0_bridge/fig_ieee39_mix60_ibr_map.pdf
outputs/ieee39_phase0_bridge/fig_ieee39_mix60_ibr_map.png
paper_ieee_transactions/figures/fig5_ieee39_ibr_map.pdf
paper_ieee_transactions/figures/fig5_ieee39_ibr_map.png
```

The Mix60 diagram maps the actual replacement used in the stable pilot:

- GFL: generators 5, 7, and 9 at buses 34, 36, and 38,
- GFM: generators 4, 6, and 8 at buses 33, 35, and 37,
- SG retained: generators 1, 2, 3, and 10 at buses 30, 31, 32, and 39.

The phase-0 gate table is intentionally conservative. It marks baseline
ANDES simulation, positive-mode diagnosis, and case construction as completed;
generic IBR stability as partial; and both the physical ANDES-to-`(L,M,D)`
bridge and final estimator validation against full ANDES IBR poles as blocked.

## Strict Phase-0 Physical Bridge Audit

Script:

```text
validation/ieee39_rational_filter/phase0_physical_bridge_audit.py
```

Command:

```text
python validation/ieee39_rational_filter/phase0_physical_bridge_audit.py
```

Outputs:

```text
outputs/phase0_physical_bridge/phase0_physical_bridge_status.json
outputs/phase0_physical_bridge/phase0_physical_bridge_report.md
outputs/phase0_physical_bridge/base_mode_match.csv
outputs/phase0_physical_bridge/no_pss_mode_match.csv
outputs/phase0_physical_bridge/base_matrices.json
outputs/phase0_physical_bridge/no_pss_matrices.json
```

Acceptance rule:

- critical mode frequency error `<= 5%`,
- first-six mean frequency error `<= 5%`,
- positive damping sign reproduced by the physically extracted scalar `D`.

Result:

- `BLOCKED`.
- The physical `L` and `M` extraction places several modes in the
  electromechanical band, but the extracted scalar machine damping is zero
  (`GENROU.D = 0`) in the packaged case.
- The full ANDES damping is supplied by controller/exciter/PSS states, so a
  scalar nodal `D` cannot be claimed without a validated controller-state
  reduction.

Interpretation:

This is a useful negative audit. It prevents the surrogate estimator from being
misreported as validated against full ANDES poles. The dependent estimator,
reversal, and planning phases remain blocked as full-ANDES evidence until this
bridge is replaced by a validated reduced model.

## Phase-1 Controller Calibration Audit

Script:

```text
validation/ieee39_rational_filter/phase1_controller_calibration_audit.py
```

Command:

```text
python validation/ieee39_rational_filter/phase1_controller_calibration_audit.py
```

Outputs:

```text
outputs/ieee39_phase1_calibration/phase1_calibration_sweep.csv
outputs/ieee39_phase1_calibration/phase1_calibration_sweep.json
outputs/ieee39_phase1_calibration/phase1_mode_participation.csv
outputs/ieee39_phase1_calibration/phase1_best_by_case.csv
outputs/ieee39_phase1_calibration/phase1_profile_summary.csv
outputs/ieee39_phase1_calibration/phase1_calibration_report.md
validation/ieee39_rational_filter/cases/phase1_calibrated/
```

Profiles tested:

- baseline generated parameters,
- slower PLL,
- softened GFL voltage/reactive controls,
- softened GFM controls,
- combined soft profiles,
- diagnostic removal of retained `IEEEST`/PSS,
- diagnostic removal of retained `IEEEX1` plus `IEEEST`.

Key result:

- Controller-only profiles do not stabilize the low-penetration cases.
- `ieee39_ibr_mix60` remains stable under the baseline and `soft_gfl`.
- Diagnostic removal of retained synchronous-machine controls stabilizes
  `GFL20`, `Mix20`, and `Mix40`, but `GFM20` remains unstable.
- Participation-factor audit confirms that the near-13.7 Hz least-damped stable
  modes are high-frequency control or auxiliary modes, not low-frequency
  inter-area electromechanical modes.

Interpretation:

The blocker is not just IBR gain tuning. A final benchmark needs a defensible
calibrated dynamic library for retained synchronous-machine controls and IBR
controls, followed by a validated reduced surrogate before estimator errors are
claimed against full ANDES IBR poles.

## IEEE39 Evidence Map

Script:

```text
validation/ieee39_rational_filter/ieee39_result_map.py
```

Command:

```text
python validation/ieee39_rational_filter/ieee39_result_map.py
```

Outputs:

```text
outputs/ieee39_result_map/fig_ieee39_result_map.pdf
outputs/ieee39_result_map/fig_ieee39_result_map.png
outputs/ieee39_result_map/ieee39_result_map_summary.json
paper_ieee_transactions/figures/fig6_ieee39_result_map.pdf
paper_ieee_transactions/figures/fig6_ieee39_result_map.png
```

The evidence map overlays the IEEE 39 topology with:

- Mix60 SG/GFL/GFM replacement classes,
- top weak-node damping-placement buses from the surrogate,
- the best damping-aware line reinforcement,
- the strongest Braess-like damping reversals,
- the current Phase-1 ANDES audit status.

Interpretation:

This is a traceability figure for the current draft. It combines validated
surrogate results and ANDES case-construction/calibration-audit artifacts, but
it is not final estimator validation against calibrated full-ANDES IBR poles.

## Required Inputs

Place the case files here or pass paths via CLI:

```text
cases/ieee39.raw
cases/ieee39.dyr
cases/ieee39_ibr_config.json
```

The dynamic case should document:

- synchronous generator models, preferably GENROU or the closest ANDES-supported
  equivalent,
- excitation/governor/PSS models if used,
- GFL IBR models with PLL, e.g. REGCA1/REGCP1/REECA1/REPCA1 where supported,
- GFM models, e.g. REGF1/REGF2/REGF3 where supported,
- penetration schedule,
- controller gains and tuning rules.

## Experiments

### A. Estimator Accuracy

For each penetration level:

1. run power flow,
2. run full small-signal eigenanalysis,
3. extract the critical oscillatory pole and true damping margin,
4. build the reduced `(M,D,L)` surrogate when physically meaningful,
5. compute:
   - homogeneous damping formula,
   - Fiedler-only formula,
   - diagonal all-mode formula,
   - second-order correction,
   - reduced QEP,
   - adaptive selected estimate.

### B. Network-vs-Control Split

Compare:

```text
zeta_net^(0)
Delta zeta_ctrl^(2)
zeta_full
```

against IBR penetration and PLL/GFM control settings.

### C. Trigger Precision/Recall

A trigger is "safe" if the selected cheap estimate is within a predefined
relative error tolerance against the full eigensolve.

Report:

- true safe,
- false safe,
- true reject,
- false reject.

False-safe cases are the serious failures.

### D. Weak Node Ranking

For each candidate bus and intervention type:

1. compute analytic or finite-difference sensitivity,
2. apply the intervention,
3. run full eigensolve,
4. compute ranking regret.

Controls:

- random bus,
- inertia-only ranking,
- Fiedler participation,
- DRI or published damping-placement index where implementable.

### E. Weak Link and Braess-Like Reversal

For each line reinforcement:

1. compute frequency sensitivity,
2. compute damping-margin sensitivity,
3. apply finite reinforcement,
4. run full eigensolve,
5. flag:

```text
Delta nu_c > 0 and Delta zeta_min < 0.
```

Controls:

- frequency-only line ranking,
- finite-difference full eigensolve,
- published QEP sensitivity formulas.

## Outputs

Expected output directory:

```text
outputs/ieee39_<timestamp>/
```

Expected files:

```text
case_metadata.json
estimator_table.csv
decomposition_table.csv
trigger_confusion.csv
weak_node_table.csv
weak_link_table.csv
braess_audit_table.csv
runtime_table.csv
fig_damping_region.pdf
fig_estimator_errors.pdf
fig_network_vs_control.pdf
fig_weak_nodes.pdf
fig_weak_links.pdf
fig_frequency_vs_damping.pdf
fig_trigger_confusion.pdf
```

## Success Criteria

The paper can make a strong benchmark claim only if:

1. the surrogate tracks the relevant full ANDES poles over a documented regime,
2. the second-order correction improves the diagonal screen in a nontrivial
   subset of IEEE 39 cases,
3. the trigger has no or very few false-safe cases,
4. weak-node/link rankings have lower regret than registered controls,
5. any Braess-like reversal is audited for pole switching and mechanism.

If these criteria fail, report the failure and reframe the paper as a
reduced-model method or theoretical diagnostic.
