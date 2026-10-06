# Experiment M — final summary

**EXP_M_STATUS: FAIL / MODEL_NOT_RECONCILED.** The unchanged analytical model used in Experiments A–K is not identical to the installed PowerDynamics model at mixed SG/GFL operating points. It must not be used for another `Z=(ρ,Kp,Ki)` optimization. No optimization, retuning, commit, or push was performed in this audit.

## What the measurements establish

- The installed IEEE-39 example and constructed devices use **100 MVA, 60 Hz**. The proposed 50/60 Hz explanation is false. Pure GFL replacements do carry a wrong **1.0 kV** busbar base where the input specifies **16.5 kV**; setting it to 16.5 kV changes neither the ExpK equilibrium residual nor any finite pole in a direct ablation ([bases](tables/TABLE_M01_component_bases.csv), [ablation](tables/TABLE_M13_voltage_base_ablation.csv)).
- Buses 31 and 39 retain their initialized ZIP injections. Direct device P/Q balances close below `9e−11 MW`. The mixed candidate does **not** retain the analytical per-device P/Q allocation: the initialized ExpK bus-39 SG produces `78.405 MW`, versus `2.710 MW` in the frozen analytical sharing. ExpK bus-38 SG produces `−0.230 MW`, versus `+0.598 MW` ([device ledger](tables/TABLE_M04_direct_device_PQ_balance.csv)). Fractional SG rating and GFL current scaling alone do not enforce the target dispatch.
- Complete reduced-matrix relative errors are `2.92e−15` (all-SG), `2.48e−3` (ExpG), and `2.51e−4` (ExpK). The initialized PD ExpK rightmost pole is `+0.005106001794 s⁻¹`; the old analytical model reports `−0.050001000474 s⁻¹` ([matrix identity](tables/TABLE_M11_matrix_identity.csv), [poles](tables/TABLE_M16_pole_matching_summary.csv)). ExpG has a second near-zero pole whose physical/gauge classification remains unresolved.
- A separate positive-feedback Schur closure built from the **initialized PD open-loop ports** reproduces direct PD matrices within `7.7e−18` relative A error and complete finite poles within `1.65e−9 s⁻¹` for ExpK ([corrected identity](tables/TABLE_M21_corrected_closure_identity.csv)). This is a fixed-point reconstruction, not an optimizer-ready parameterized operating-point map.
- The corrected closure and direct PD linearization have relative trajectory differences at most `2.01e−12` for the tested bus-16 load and SG-39 state perturbations. Nonlinear PD pulse errors decrease with amplitude; the upper two amplitudes show approximately quadratic absolute error. The smallest amplitude reaches a numerical floor, making the three-point exponent `1.56` inconclusive as an asymptotic law ([linear](tables/TABLE_M19_linear_trajectory_identity.csv), [nonlinear](tables/TABLE_M20_nonlinear_linear_scaling.csv)).

## Decision

**MODEL_IDENTITY_CERTIFIED: NO. MODEL_SAFE_FOR_Z_OPTIMIZATION: NO.** The next model reconstruction must enforce initialized per-device SG/GFL P/Q shares, reject machine states outside declared bounds, set pure-GFL voltage bases from the frozen bus input, and explicitly classify every near-zero mode. Only then should the support and `ρ` identity gates be repeated.

The detailed evidence, methods, provenance, figures, and limitations are in [REPORT_EXP_M.md](REPORT_EXP_M.md), [DERIVATION_EXP_M.md](DERIVATION_EXP_M.md), and [TERMINAL_SUMMARY_EXP_M.txt](TERMINAL_SUMMARY_EXP_M.txt).
