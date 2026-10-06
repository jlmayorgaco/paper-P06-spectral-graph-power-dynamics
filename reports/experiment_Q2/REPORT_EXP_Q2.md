# Experiment Q2 — Secure SG→GFL co-design gate report

## Result

**Primary status: `FAIL_FREQUENCY_METRIC`.** The baseline and the local analytical frequency identity pass, but the project does not contain a frozen bandwidth or event rise-time definition for applying the 0.5 Hz/s RoCoF limit to a bus-voltage frequency. The stipulated bus output `f_grid=(1/2π)d arg(V)/dt` has a nonzero algebraic phase jump under the ideal load step. Its unfiltered continuous RoCoF therefore has no finite peak. Finite sampled and filtered values vary with estimator bandwidth. Selecting one of them would silently redefine a hard constraint, so no secure KKT problem or candidate was frozen.

This is an input-definition failure, not evidence that the SG/GFL design is infeasible under every physically measured RoCoF definition. The previously frozen project limit was used with retained-SG COI; the audit found no technology-independent grid-frequency estimator or filter in ExpG, ExpP, or ExpQ.

## Gates completed

| Stage | Status | Evidence |
|---|---|---|
| Q0 baseline | PASS, 14/14 | [baseline table](BASELINE/TABLE_Q2_BASELINE_REGRESSION.csv) |
| F0 frequency output | PASS_CONTINUOUS; RoCoF METRIC_DEPENDENT | [continuity](F0/TABLE_Q2_F0_frequency_metric_continuity.csv), [edge errors](F0/TABLE_Q2_F0_continuity_errors.csv), [bandwidth sensitivity](F0/TABLE_Q2_F0_bandwidth_sensitivity_allSG.csv) |
| F1 linear/PD identity | PASS_LOCAL_LINEAR_IDENTITY | [16 independent PD runs](F1/TABLE_Q2_F1_grid_frequency_identity.csv) |
| F2 modal residues | PASS reconstruction; event exposes unfiltered impulse | [residues](F2/TABLE_Q2_F2_modal_frequency_contributions.csv), [ExpN bandwidth sensitivity](F2/TABLE_Q2_F2_ExpN_bandwidth_sensitivity.csv) |
| F3 PLL DC support | PLL_DC_SUPPORT_INDEPENDENT | [30-point centered differences](F3/TABLE_Q2_F3_PLL_DC_support.csv), [dynamic derivatives](F3/TABLE_Q2_F3_dynamic_gain_derivatives.csv) |
| F4 affine stiffness | NONAFFINE | [200 mixed/support samples](F4/TABLE_Q2_F4_steady_frequency_samples.csv), [holdout](F4/TABLE_Q2_F4_affine_holdout.csv) |

## Measurements and evidence

The primary passive output used for model identity and continuity is the derivative of the unwrapped voltage phase at generator buses 30–39. It exists in all-SG, mixed, and all-GFL architectures, adds no current, and uses the same channels in every architecture. For filtered diagnostics, the fixed comparison window is a cubic 11-sample Savitzky–Golay derivative at 100 Hz (0.10 s support, derivative −3 dB response 11.7069 Hz). This window is explicitly a diagnostic because no frozen project bandwidth could be recovered.

At `rho=1e-8`, the largest absolute difference from the `rho=0` grid-frequency metrics over the tested buses/metrics was `5.61856e-11`; the condition-aware tolerance was `2.99701e-5`. Thus the passive grid-frequency metric is continuous across the tested SG-removal boundary. The internal PLL-frequency diagnostic is not used as system frequency.

F1 compared the analytical output with independent PowerDynamics traces at all-SG, corrected ExpG, ExpN, and five deterministic random mixed designs. It used 1 MW and 5 MW steps and the same bus-phase estimator at buses 30–39. Maximum normalized linear/PD error was `0.1005%`; maximum normalized 1-to-5 MW scaling error was `0.0804%`; maximum P/Q trim error was `5.73e-13 pu`. The raw channel-by-time series for all 16 PD runs and analytic traces are in `F1/TABLE_Q2_F1_grid_frequency_sensors.csv` (64,160 rows). This validates local identity only; it is not a large-signal validation of a new candidate.

For the frozen ExpN point, the analytical grid-frequency peak for 100 MW is `21.25465 Hz`, steady offset `23.10445 Hz`, and rightmost pole `-0.050000001105 s^-1`. Its load-step algebraic phase jump is `0.0012662112 rad` at bus 38. The unfiltered continuous RoCoF therefore includes a distributional impulse. Finite estimates vary: sampled 100 Hz difference `2.01815 Hz/s`; 0.04, 0.10, 0.20, and 0.40 s Savitzky–Golay windows yield `1.39992`, `0.80353`, `0.78338`, and `0.74261 Hz/s`. All tested finite-window estimates exceed the project limit at this baseline; those values do not select or validate a project measurement bandwidth.

The complete finite-mode residue sum reconstructs the F2 frequency response at six checked times within `1.53e-11 Hz`; its zero-time smooth RoCoF reconstruction error is `7.08e-10 Hz/s`. The spectral-abscissa mode dominates the frequency peak and steady offset, but rapid local modes dominate the finite smooth RoCoF. This further rules out using the rightmost pole as a RoCoF proxy.

Centered gain differences over 30 withheld designs gave maximum `|dF_inf/dKp|=7.66e-14` and `|dF_inf/dKi|=1.25e-14 Hz/MW per gain unit`. Dynamic derivatives at active GFL bus 38 in ExpG were nonzero for the filtered RoCoF but near zero for F∞, consistent with PLL shaping rather than sustained active-power support. The pointwise beta derivative is evaluated at the frozen ExpG peak frequency; it is not a recomputed full H-infinity optimum.

The 200-design F4 study had 186 finite steady responses. A 100-design fit / 86-design holdout rejected affine stiffness: maximum holdout relative error `1.78024` and maximum absolute error `3298.27 MW/Hz`. No affine formula was used as a global lower bound.

## Optimization and validation disposition

No d=0 robust KKT, 0–100 MW continuation, support search, or candidate freeze was run. The hard `R_peak_grid <= 0.5 Hz/s` constraint has no single finite value for the ideal step under the unfiltered measurement; every filtered sensitivity value corresponds to a different constraint. An optimizer could return a bandwidth-conditional candidate, but it would not be the unique `Z*_secure` requested under the frozen inputs. The fail-closed response is to stop before tuning rather than use a post hoc filter. Accordingly, no post-freeze PowerDynamics validation, 100 MW grid-frequency trajectory validation, or GSP work is claimed.

The ExpN baseline's pre-existing independent PowerDynamics spectrum remains reproduced: analytic alpha `-0.050000001105`, PD alpha `-0.050000000139 s^-1`. This is not a new secure candidate. The two generated figures are [metric continuity](figures/FIG_Q2_01_frequency_metric_continuity.png) and [stable versus frequency-insecure ExpN](figures/FIG_Q2_02_ExpN_stable_not_frequency_secure.png). No frontier/authority/BND optimization figures are generated because there is no frozen secure optimum.

## Reproducibility

- Branch: `research/expQ2-secure-kkt`; no commit or push.
- Model SHA: `e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a`.
- Julia / PowerDynamics / NetworkDynamics: 1.11.9 / 5.0.0 / 1.3.0.
- Frozen inputs, candidates, and dependencies were not modified.
- Scripts are in `experiments/bnd_expQ2/`; source and tests are in `src/bnd_expQ2/` and `test/bnd_expQ2/`.
- Unit tests: 7/7 passed. The F0 numerical architecture tests passed all 1,000 continuity rows.
- F0 runtime: 414.437 s for 60 architecture points. F1: 8 configurations, 16 PD simulations. F3: 30 centered DC cases (60 finite-difference perturbation models), plus six dynamic perturbation evaluations. F4: 200 fixed-support model evaluations. F1/F3/F4 wall time was not instrumented separately.
- No dependency or lockfile changes.
