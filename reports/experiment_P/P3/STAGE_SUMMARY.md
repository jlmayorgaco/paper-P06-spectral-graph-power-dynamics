# P3 — Joint continuous co-design on the incumbent support

Status: **PASS_LOCAL_CERTIFIED**. Fixed support is SG bus 38, with one retention variable and all 20 local controller gains active.

- Equation verified: joint first-order SQP subproblem for retained-MW objective under every spectral pole within `1e-6 s⁻¹` of the rightmost pole; after each trial, the complete gauge-quotiented finite spectrum is recomputed. No one-pole correction or gain grid is used.
- Incumbent/candidate: `1.124344477653483` / `1.124344477653483` retained MW; α=`-0.0500000011050941` s⁻¹. The ExpN incumbent was preserved because no strict improvement exceeding `1e-8 MW` was resolved by this SQP run.
- KKT: primal `0.0`, stationarity `0.0`, complementarity `8.195473598010088e-10`; LICQ `true`, active rank `1`, SOSC `VACUOUS_CRITICAL_CONE_STRICT_GAIN_BOUNDS`. Critical eigenvector conditioning `913.3877667018162`.
- Guard qualification: this solve and the frozen candidate use the historical `1e-9 s⁻¹` guard, not the proposed new-candidate `1e-6 s⁻¹` guard. The separate rank-two audit finds a conditional guarded root at `1.1243522749761832 MW`, but that point was not KKT-corrected or PD-validated. Do not treat the P3 KKT certificate as applying to that point or to the larger guard.
- GSP: actual feature rank `5` of `5`; 10 independent local PLLs remain, no communications; fit residual `0.0`, replacement loss `0.0` MW.
- Evaluation count: two spectrum evaluations were counted by the stage wrapper; direct mode-sensitivity calls and QP linear algebra are not separately instrumented. Accepted predictor-corrector steps `0`; elapsed `34.598 s`.
- Certificate scope: `LOCAL_KKT_CERTIFIED` for the fixed support only. Other architecture regions and disconnected branches remain open; no globality claim. This allows P4 robustness/transient audit and P5 freeze/independent validation.
