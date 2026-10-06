ALL-PLL ANALYTICAL GAIN MAP — ONE PREREGISTERED VALIDATION

Start with REPORT_ES.txt and SUMMARY.json. Protocol and prediction hashes
precede independent model validation. No replacement optimization was run.

Reproduction from the repository root (requires its frozen source inputs):
  python -B experiments/all_pll_gain_map_validation_20261003/predict.py
  julia --startup-file=no --project=experiments/physical_collective_damping_20261003/source_snapshot experiments/all_pll_gain_map_validation_20261003/validate.jl
  python -B experiments/all_pll_gain_map_validation_20261003/check_independent_roots.py
  python -B experiments/all_pll_gain_map_validation_20261003/finalize.py

These commands regenerate outputs. Preserve this executed experiment before
rerunning it; use a separate repository copy for independent reproduction.
No Julia package installation is needed with the recorded environment.
PACKAGE_VERSIONS.csv and FINAL_MANIFEST.json record runtime and source hashes.
Re-aggregation alone uses finalize.py and does not run the model.

TABLE_01: reconstruction for every nonreal pole in the frozen catalog.
TABLE_02: nine vector gain derivatives versus centered finite differences.
TABLE_03: baseline, first-order prediction and exact calculated nodal gains.
TABLE_04: parametric-model modal predictions and checks. Its PLL_pattern_MAC
          compares every mode to the SINGLE TARGET pattern. It is meaningful
          as a preservation gate only on the target row; it is NOT the MAC
          against each mode's own baseline. This frozen table is preserved.
TABLE_05: two prescribed poles AND patterns, restricted assignment residuals.
TABLE_06: independent nonlinear-equation equilibrium and AD Jacobian parity.
TABLE_07: numerical all-root count to the right of the declared spectral margin.
TABLE_08: five nonlinear delayed event results for the analytical candidate.
TABLE_09: independent full-model roots, each with its OWN baseline-pattern MAC.
PREDICTED_ROOTS and TARGET_PATTERN: sealed before independent validation.
designs/: complete rho, Kp, Ki and fixed delay.
independent/: Jacobians exported directly from nonlinear equations.
nonlinear/: all five candidate trajectories.

FIG_01: analytical reconstruction, other-mode predictions and nodal retuning.
FIG_02: historical baseline versus the new analytical candidate. Replacement
        and gains both change; this is NOT a same-rho causal tuning comparison.

The open THEORY.tex may receive the generated validation section through
update_theory.py, preserving its pre-validation bytes under source_snapshot/.
Previous theory provenance records and frozen experiment results are unchanged.
The native LaTeX compiler status is recorded separately.

Boundaries: numerical consistency and one declared-event feasible candidate;
no global optimality, no maximum replacement, no heterogeneous-delay result,
no proof that neighboring feedback is necessary, no certified nonlinear safety.
