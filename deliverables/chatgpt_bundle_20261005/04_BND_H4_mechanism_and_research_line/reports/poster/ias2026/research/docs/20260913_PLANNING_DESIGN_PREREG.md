# Preregistration: planning and design campaign (PD1, PD2)

Date: 2026-09-13. Branch `research/series-planning-design` (from `e03a5312`).
Frozen rules: `configs/ias2026/planning_design_prereg_v1.yaml`. Code:
`experiments/planning_design/`. This document, the config and the code are
committed before any campaign computation.

The campaign supplies the missing design evidence for Part V of the paper series
(`reports/papers/series_replacement_portfolios/`). It answers two questions left
open by the frozen record:

- **Gap 3** (`docs/TPWRS_REMAINING_GAPS.md`): E35 found that the frozen RC
  converter retune, optimized once at the nominal point, restores stability at
  117 of 240 unstable held-out operating points (80/80 mild, 37/80 medium,
  0/80 severe), while a 25 % condenser works at all 240. Nobody asked whether an
  *adaptive* retune, re-optimized per operating point, succeeds everywhere.
- **GOLD-D at census scale** (CDW E9): at the four-unit scale single-boundary
  tuning always sufficed; at the nine-unit census (512 subsets) neither method
  converged within 6 iterations, so the result was inconclusive.

## PD1: adaptive per-point retune

**Samples.** The 240 E35 samples with status ACCEPTED and an unstable four-unit
portfolio. Each operating point is regenerated with the E35 rule
(`default_rng(20260913 + 100003·sample)`, the stored load, reactive load and
availability, load scatter 0.05, dispatch jitter 0.10).

**Reproduction gate.** The flagship inter-area abscissa `alpha_IA` must match
the frozen E35 value to 1e-6. If more than 5 % of samples fail, the campaign
stops and reports a defect. The pre-run check on one sample reproduced it
exactly (0.19660086465931648).

**Optimizer.** E34 M1 verbatim: minimum `||v||²` over the five physical
log-coordinates (PLL natural frequency and damping, active and reactive outer
bandwidths, measurement-filter bandwidth), bounds ±ln 8, SLSQP from −0.15,
`maxiter` 90, `ftol` 1e-7. The one change is the constrained quantity: the
transverse spectral abscissa (exact quotient by the centre subspace) replaces
the E34 disc rule `|λ| > 1e-3`. The target is −0.02, the first E34 margin.

**Repairs labelled at each point.** Frozen RC retune; 25 % condenser; adaptive
retune.

**Labels.**
- **S0 (primary):** the transverse whole-right-half-plane status of the FC01
  direct path (`physical_report(case).status_rest`) is STABLE. Unresolved and
  infeasible outcomes count as failures.
- **S1:** S0 and `alpha_perp ≤ −0.02`.
- **SB (comparability only):** no eigenvalue in the right half plane within the
  0.3–1.5 Hz inter-area band. This is the E35 criterion behind 117/240.

**Severity.** The frozen E35 terciles of `alpha_IA` over the 240: mild ≤ 0.12215459,
medium ≤ 0.24183127, severe above.

**Hypotheses and rules (fixed now).**

| id | statement | rule |
|---|---|---|
| H-PD1a | the adaptive retune restores transverse stability across the envelope | S0 rate ≥ 0.95 SUPPORTED; [0.80, 0.95) PARTIAL; < 0.80 NOT SUPPORTED |
| H-PD1b | it does so in the severe tercile | severe S0 rate ≥ 0.80 SUPPORTED |
| H-PD1c | it beats the frozen retune, paired | exact one-sided McNemar p ≤ 0.01, with more adaptive-only than RC-only successes |

**Descriptive outputs.** Per-tercile S0/S1/SB rates for RC, condenser and
adaptive; Spearman correlation between the required change norm and severity
among the adaptive successes; fraction of solutions with an active bound; SLSQP
success fraction.

**Scope notes.**
- E35 predates the MW/MVA unit correction. The solvers never used the
  mislabelled quantity, so the labels are unaffected.
- E35 predates the transverse re-audit. Every PD1 label is recomputed on the
  transverse quotient (S0); the band criterion is reported only to compare with
  117/240.
- E35 is the census-policy Track-A model (fixed reactive power, `k = t = 1`),
  not the TX4 policy planes.

## PD2: census-scale plan-level design (E9 re-run)

**Unchanged from E9** (CDW prereg `b8ae3082`): targets and policies (T2 = V9 at
D01 and T4 = V9 at D03 are primary; T1 = H4 at D01 and T3 = H4 at D11 are
controls), the parameter families (control: Q/V gain of every replaced unit and
AVR-gain scale of every surviving machine; topology: admittance scale of the 46
branches; joint), bounds and ranges, `EPS = 0.02`, the SPR total gradients, the
linearized QP step with its Gordan least-violation fallback, the stopping rules
(single: `alpha_T ≤ −EPS/2`; plan: `Φ ≤ −EPS/2`) and the GOLD-D criterion.

**Changed, and why.** At census scale 170–179 of the 512 subsets are unstable
at the nominal policy, many through fast converter-control modes with abscissae
of 10²–10⁴ s⁻¹. E9 accepted every QP step even when the merit increased (for
example Φ rose from 1983 to 25 701 in one task) and capped the active set at 25.
- **Budget:** 30 iterations, at most 40 lattice evaluations and 6 h per task
  (was 6 iterations).
- **Globalization:** a step is accepted only if the merit decreases (single:
  `alpha_T`; plan: `Φ`). Otherwise the trust radius is halved (from 0.25 × range
  down to 0.25/16 × range) and the QP is re-solved with the same gradients.
  After an accepted step the radius doubles, up to 0.25.
- **Plan active set:** every subset with `alpha ≥ −EPS`, ranked by `alpha`, cap
  60, keeping previous constraints while they remain within `−EPS`. Non-finite
  abscissae are skipped.
- **Start:** the nominal policy, not the E9 final designs.

**Hypothesis H-PD2 (GOLD-D at census scale, T2 or T4, any family).**
- SUPPORTED: single-boundary tuning leaves an unsafe subset (`alpha_T` STABLE,
  `Φ ≥ 0`) while plan-level design is safe (`Φ < 0`, no unresolved subset, no
  unstable subset).
- NOT SUPPORTED: a census plan-level safe design exists, but single-boundary
  tuning did not leave an unsafe subset.
- INCONCLUSIVE: no census plan-level design reached a safe lattice within the
  budget.

The controls T1/T3 are reported against E9 (single-boundary tuning sufficed) and
are not gated.

## Compute budget (measured before freezing)

| operation | time |
|---|---|
| one E34-M1 constraint evaluation (PD1) | 0.11 s |
| one census portfolio evaluation (PD2) | 0.23 s (512 subsets ≈ 117 s) |
| one joint-family gradient (55 parameters) | 7.3 s |

- **PD1:** at most about 550 evaluations per sample, so about 15 min on 16 workers.
- **PD2:** a census plan task with 60 active gradients runs at up to about 7 min
  per iteration, plus lattice evaluations. The worst case is about 5 h per task,
  with 24 tasks on 18 workers, which is below the 10-h ceiling.

## Pre-run checks (not results)

- **PD0 timing:** on one campaign sample it reproduced the frozen E35 flagship
  label and evaluated the M1 constraint once at the optimizer's starting point.
  No optimization was run.
- **Code smoke tests,** on inputs outside the campaign only:
  - a stable E35 sample, which is not among the 240;
  - the four-unit lattice at discovery policy D04, which is not a PD2 target,
    for 2 iterations.
- **One execution fix before freezing:** the per-task checkpoint directory was
  created lazily; it is now created at task start.

## What will not be done

- No threshold, target, tercile or rule changes after any campaign result is seen.
- Execution fixes are logged in `docs/20260913_PLANNING_DESIGN_DEVIATIONS.md`.
- No push.
