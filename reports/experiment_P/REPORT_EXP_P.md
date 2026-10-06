# Experiment P — Algebraic frontiers and incremental ExpN co-design

## Result

**Primary status: `FAIL_TDS`.** The ExpN PD-exact model is reused without rebuilding it from E/H/K. P0–P3 verify the baseline, fixed-support algebraic identities, conditional retention roots, and a nominal local KKT point. Independent PowerDynamics confirms the frozen nominal candidate's equilibrium and full spectrum, but the candidate fails the inherited sustained 100 MW step requirement. Robust co-design and globality remain open.

Work remained on branch `research/expN-pd-exact-zstar`; parent commit: `d0fecb3264aeb855ab6700fc6ef66120d4b6c33c`. No commit or push was made. Frozen ExpN model SHA: `e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a`. ExpN candidate SHA: `3915a8f57da0f552779b11464572f14a86f0b648b532be77124f668eeff9f263`. Project and Manifest hashes remain `e170a2f9fca57aa3f8f3f053b318f0d5e3d434679dbdd4a4b8eef5481cfc0e4f` and `92bc7849ea4445a4bedd18f732611ec271147e238c349aa5eeefc6d02023e911`.

## Frozen incumbent and candidate

The original component dispatch weights sum to `5402.761089978847 MW`. The pre-existing ExpN candidate remains support `{38}`, `rho_38=0.9986453680992127`, all `Kp=7.853981633974483`, all `Ki=986.9604401089358`, retained SG `1.124344477653483 MW`, GFL dispatch `5401.636745501193 MW` (`99.97918944667498%`), and analytical alpha `−0.0500000011050941 s⁻¹`.

ExpP's unchanged nominal candidate freeze is [Z_P_NOMINAL_FINAL.toml](P5/Z_P_NOMINAL_FINAL.toml), SHA-256 `dbe727713acccec9b0802619f4b2e10763a0a88b9799ba24a74411ac88ff7046`. It equals the incumbent; it is a nominal fixed-support candidate, not a robust/transient-feasible solution. The candidate was not edited after freeze.

## Stage record

| Stage | Result | What was verified | Boundary of the result |
|---|---|---|---|
| P0 | PASS | ExpN reproduction, component P/Q split, one withheld mixed case, full physical poles and reduced A | 0 design-time PD calls; two independent PD validations. Does not add an optimum. |
| P1 | PASS | PLL single-device rank-one update, affine same-device gain pair, SG rank-two retention, two-SG rank-four identity, real PI boundary | Conditional fixed-architecture identities. Near-mode conditioning reaches `6.01e11`; measured determinant residuals are recorded. |
| P2 | PASS | Algebraic single-SG roots on all ten generator buses; complete physical spectrum checked | Roots are conditional at frozen gains. The bus-38 root recovers ExpN within `7.21e-9 MW`; it is not a joint optimum. |
| P3 | PASS_LOCAL_CERTIFIED | One joint local SQP/KKT solve over bus-38 retention and 20 gains, with complete-spectrum rechecks | Only support `{38}` was jointly corrected; zero accepted descent steps; historical `1e-9 s⁻¹` guard. The suggested `1e-6` guard root is diagnostic and lacks a joint KKT/PD validation. |
| P4 | INCOMPLETE_ROBUST_OPTIMIZATION | Frozen ExpG full-block uncertainty definition reused; nominal incumbent rejected by a pointwise witness, and a fixed-gain conditional point passes a Float64 Lipschitz interval lower-bound test | The interval computation used no outward rounding and has only `3.21e-11` beta margin. It is provisional numerical evidence, not a formal validated-arithmetic result. No joint robust KKT or robust optimum. Sustained transient was not evaluated in this stage. |
| P5 | FAIL_TDS | Frozen nominal candidate passes independent PD equilibrium/full-spectrum identity and six small-pulse/linear comparison checks | The declared 100 MW sustained bus-16 step fails both inherited limits and does not settle within 60 s. The candidate is rejected for the combined requirement. |
| P6 | GLOBAL_GAP_OPEN | 1,024 SG-support masks enumerated in a ledger; GSP gain representations evaluated at the incumbent | Only one support has a joint KKT point; 1,023 remain open. The `0 MW` lower bound is trivial. No global or robust claim. |

Phase results and stage-specific equations, residuals, time, and outputs are in the `P0`–`P6` folders. `STAGE_STATUS.json` records `last_completed_stage=P6` and the separate P5 nominal-PD, small-pulse, and declared-event subgates.

Figures render only saved observations: [conditional roots vs full-spectrum alpha](FIG_P01_conditional_retention_roots.png), [pointwise robustness bounds](FIG_P02_robust_pointwise_bounds.png), and [the sustained-step PD trajectory](FIG_P03_declared_step_tds.png). No interpolated robust frontier or uncomputed continuation curve is drawn.

## Algebraic results

The PD-exact raw component port convention is current into a device, so `Yinj=−Yraw`; the network closure uses `T=Ystatic+ΣΠYrawΠᵀ`. On a fixed architecture, each GFL PLL update has numerical rank one (`1.01e-17` maximum rank residual), the same-device two-gain determinant is affine to measured residual `2.13e-8`, and cross-device terms are nonzero. The near-mode reference pencil condition is `6.01e11`, so the small determinant residual is conditioned evidence, not exact floating-point arithmetic.

Single-SG retention has the conditional determinant `det(I₂+εN_i)=1+tr(N_i)ε+det(N_i)ε²`; the tested single/two-SG operator residuals are `6.94e-18` / `6.10e-18`. The two-port determinant includes a nonzero cross term; its maximum determinant residual is `7.99e-10`. At zero frequency the PI boundary is treated as a real line, not by inverting a singular complex-to-real 2×2 system.

At the historical guard `1e-9 s⁻¹`, the bus-38 root gives retained SG `1.1243444848621267 MW`, `7.20864e-9 MW` above ExpN. The proposed new-design guard derived from the observed P0 alpha error is `max(1e-6,10×9.66594e-10)=1e-6 s⁻¹`; its conditional bus-38 root is `1.1243522749761832 MW`. It is not a redesigned candidate and is not PD-validated.

## Nominal KKT and globality

P3 uses the dispatch-weighted retained-SG objective, not installed capacity. For fixed support `{38}`, its recorded residuals are primal `0`, stationarity `0`, complementarity `8.19547e-10`, LICQ true with active rank 1, and a strict-gain-bound/vacuous-critical-cone SOSC check. The active pole condition is `913.39`. The incumbent remains the candidate because the local SQP found no resolvable strict decrease. These checks apply to the historical guard and fixed support only.

The simple-pole sensitivity implementation is reused from ExpN, whose `TABLE_N08_derivative_validation.csv` contains 90 centered-difference comparisons over three withheld cases. Maximum absolute discrepancies were `2.75e-7` for rho, `3.21e-8` for Kp, and `7.28e-10` for Ki; relative errors grow for derivatives near zero. ExpP did not add a separate P3 finite-difference sweep at the active KKT point, so this is supporting inherited validation, not a new ExpP gradient gate.

P6's nominal spectral-only valid interval is `0 ≤ J* ≤ 1.124344477653483 MW`; its gap is `1.124344477653483 MW`, far above the study tolerance `0.01 MW`. For robust-plus-transient design, the nominal candidate is not feasible, so this interval provides no upper bound and no combined-problem gap.

The GSP basis has measured rank `5/5` on the ten generator feature rows. Two common gains, feature prefixes with 2–5 features per channel (4–10 coefficients), and the 20-gain representation all reproduce the incumbent's uniform gain vector exactly, with zero fit error and zero replacement loss at that point. This is an incumbent representation comparison, not a GSP optimization or a controlled self-energy truncation certificate. Ten local PLLs are preserved; no remote signals or communications were added.

## Robustness and transient result

The reused ExpG full-block requirement is normalized `beta_req=1.6991206999182038e-6`, not a physical percentage. At the ExpN candidate, the frequency-zero pointwise witness bounds its radius above by `3.746493388818478e-11`, proving that candidate violates the requirement. For a fixed-gain one-coordinate continuation point, ε₃₈=`0.0013695147207744968`, retained SG is `1.136697218242845 MW`, and complete-spectrum α is `−0.05169088501717586 s⁻¹` (101 physical poles). The sampled peak gives only an upper bound on its radius, `1.7329480509001635e-6`; the separate 34,333-node Lipschitz interval calculation reports a lower bound `1.6991528445276632e-6`, exceeding beta_req by `3.214460945939892e-11`. The interval routine uses Float64 arithmetic without outward rounding, so this remains provisional numerical evidence rather than a formal H∞ certificate. The attempt to obtain a separate CARE upper bound on the norm is rejected: relative residual `4.23175e-4`, maximum residual eigenvalue `+0.0661052`. The conditional point is neither a joint robust KKT solution nor a frozen/PowerDynamics-validated candidate.

The inherited event is a 100 MW sustained active-power step at bus 16, beginning at 1 s, with the COI speed measured from retained SG bus 38. In the independent PD trajectory the peak absolute frequency deviation is `39.3334 Hz` and peak sampled RoCoF is `1.07221 Hz/s`, against project limits `0.5 Hz` and `0.5 Hz/s`; no 2% settling is observed within 60 s. The small 0.1 s pulses at buses 8, 16, and 29 scale well and closely match the linear response, but are a different event and do not rescue the sustained-step failure. No hard current-limiter safety claim is made.

## Reproduction and limitations

Use the commands in [README_REPRODUCE.md](README_REPRODUCE.md). Tests recompute the algebraic identities and roots and assert the PD/TDS evidence and candidate hashes. The P5 candidate is immutable; use the checked-in post-freeze validation scripts to reproduce its independent checks without overwriting its contents.

Not completed: joint robust KKT; robust/transient co-design and validation of a robust candidate; transient constraints in the optimizer; a new candidate meeting the `1e-6 s⁻¹` numerical guard; branch-complete support exploration; a nontrivial global lower bound; GSP optimization and controlled self-energy remainder/pivot comparison; generic-solver falsification; exhaustive PI region topology; peak/tail certificate for the sustained nonlinear step; mutation/fault-injection tests for each deliberately corrupted implementation scenario. No result from ExpE/G/H is substituted for missing same-model evidence.
