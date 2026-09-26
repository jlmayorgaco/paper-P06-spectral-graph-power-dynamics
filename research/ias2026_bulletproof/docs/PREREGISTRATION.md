# Preregistration — IAS 2026 Bulletproof Closure

**Freeze time:** 2026-09-18 local execution start

## Scope

The primary benchmark is the frozen IEEE-39 SG→GFL model and its declared transverse formulation. The flagship portfolio is `{30,33,35,37}` at P4 (`g=0.03625`, `k=1.425`, `t=1.5`, `h=1`). Existing post-freeze evidence is retrospective; any new holdout must be frozen before reveal.

## Labels

Every result receives exactly one label from the plan: `PROVED`, `NUMERICALLY_VERIFIED`, `IEEE39_VALIDATED`, `POWERDYNAMICS_VALIDATED`, `SECOND_MODEL_VALIDATED`, `NONLINEAR_TDS_VALIDATED`, `SYNTHETIC_PILOT`, `CONSTRUCTED_COUNTEREXAMPLE`, `RETROSPECTIVE`, `BLIND_HOLDOUT`, `NOT_TESTED`, `REFUTED`, `UNRESOLVED`, or `STOPPED_BY_GATE`.

## Frozen thresholds

- P4 H4 alpha absolute error: `<= 1e-6 s^-1`.
- P4 H4 frequency absolute error: `<= 1e-6 Hz`.
- All 15 proper-subset stability signs must match.
- P4 boundary `g*` absolute error: `<= 1e-5`.
- Theory identity median residual `< 1e-13`; maximum `< 1e-10` except deliberately ill-conditioned cases, which report condition and backward error.
- TDS verdict equality uses the frozen declaration and reports the declared tolerance.
- No downstream claim is promoted when its prerequisite gate fails.

## Holdout lock

The new holdout prediction file is created before any reveal. Existing branch doublings `0,1,13,43`, the frozen F7 policy planes, and all archived final results are retrospective and cannot be called blind.

## Deviations

Any change to thresholds, selected policy, candidate action, model, or metric after a prediction freeze must be added to `docs/DEVIATIONS.md` before rerunning.
