# Poster-safe claims, v2 (unit-corrected) — supersedes POSTER_FINAL_SAFE_CLAIMS.md

v1 is frozen (FREEZE_SHA256SUMS.txt) and left untouched. v2 changes only what
the MW/MVA correction requires (`docs/PHASE_I_UNIT_CORRECTION.md`). Every other
statement of v1 is carried over verbatim in meaning. The poster sources
(`sections/*.tex`) contain none of the corrected statements. The poster has
**not** been edited.

## Unit rule (new, binding)

- **MW** is used only for active power:
  - dispatch `Pg`, measured;
  - the documented synchronous limit `Pmax`.
- **MVA** is used only for apparent-power ratings `Sn`: machines, converters,
  condensers.
- **Mvar** is used for reactive power.
- No PV nameplate exists in these benchmarks. Never write "X MW of PV
  capacity".
- Flagship, in one sentence: *the four replacements displace 2096.6 MW of
  active dispatch from synchronous machines rated 4270.7 MVA; the converters
  are rated 4270.7 MVA and carry that dispatch.*

## Safe to state

**S1–S7, S9 — unchanged from v1.** They contain no MW or MVA quantity.

**S8 — mitigation, stated as a frontier (unchanged; units verified).**
> At this operating point a converter retune reaches the full base-case margin
> with no converter active dispatch forgone and no synchronous capacity added,
> and a 166–270 MVA synchronous condenser at a single bus does the same without
> touching the converters.

E34 re-scored on measured Pg (UC03): the Pareto front membership is identical.

**S10 — the conventional-baseline comparison (CORRECTED).**
> A model keeping interactions only to third order does not merely rank these
> portfolios badly — under the declared threshold it never predicts instability
> at all, declaring every portfolio stable including the 139 that are not.
> Removed inertia (AUC 0.78–0.87) and the retired machine rating (0.77–0.86,
> MVA) rank instability reasonably well; the removed active dispatch does not at
> the minimum failing order (0.57, permutation p = 0.46) and only moderately at
> larger sizes (0.64–0.73); the short-circuit family is weak, 0.59–0.72.

E39 reproduced exactly (1e-16), plus UC02 (Pg, Sn, Pmax, Q, fresh permutation
seed). Figure: `results/UC/UC02/UC02_ROC_PR_unit_corrected.png` replaces
`E39_ROC_PR.png`, whose curve labelled `replaced_mw` is the MVA rating.

## Must be stated alongside (unchanged): Q1–Q4 of v1

## Forbidden wording (v1 list, plus)

- everything forbidden in v1, except the E39 item, which becomes:
  - "conventional baselines are anti-predictive" — refuted for inertia and
    machine rating; only the short-circuit family and removed dispatch are
    weak;
- "replaced megawatts rank instability" or "megawatts are informative rankers"
  — the old column was MVA; true megawatts of dispatch are weak (UC02);
- "4270.7 MW", "keeping all 4270.7 MW of PV", "costs 970–1040 MW of PV",
  "matched on MW" for the E14/E18 controls (they were matched on MVA; the
  Pg-matched reruns are in UC02);
- any "PV capacity" in MW.

## Figures cleared for use

As in v1, except `E39_ROC_PR.png` → `results/UC/UC02/UC02_ROC_PR_unit_corrected.png`.
