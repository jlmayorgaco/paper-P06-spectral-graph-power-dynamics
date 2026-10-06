# ParaEMT EMT-validation campaign (IAS2026 / TX4) — V2 report

**Status: V2 STOPPED at gate G3a.** Branch `research/paremt-emt-validation`, not
pushed.

**V1 (unchanged, permanent):** G3 FAIL, G4 FAIL, EMT04–EMT18 BLOCKED.

**V2** (preregistered at 0ab8efd4, before any V2 computation):

| gate | result |
|---|---|
| SW | PASS |
| **G3a** | **FAIL** — the preregistered ringdown sub-check failed on the blind case. Every state trajectory agreed to ≤ 2.1e-4 of its excursion. |
| G3b, G4a, G4b | **NOT RUN** (stop rule). The GFL blind holdouts remain blind. |
| EMT04–EMT18 | BLOCKED |

**New finding.** The v1-frozen modal estimator A (multi-channel matrix pencil)
is numerically unfit for its purpose. On synthetic IEEE-39-type signals it
crashes in 3 of 4 cases and is unresolved in the fourth. Any continuation needs a
pre-validated estimator before any portfolio run.

## 1. Required wording

> The first preregistered EMT realization failed its device-interface gates.
> A second preregistered realization, designed after diagnosing those failures
> and evaluated on new holdouts, stopped at its synchronous-machine
> transcription gate. Every machine state trajectory agreed with the frozen
> phasor model to 2e-4 of its excursion. The preregistered ringdown sub-check,
> which reused the frozen v1 matrix-pencil estimator, returned a spurious mode
> for both the EMT run and its phasor reference. The grid-following-converter
> gates were not reached, and no portfolio result was produced.

The V1 FAIL is not replaced, and there is no V2 PASS.

## 2. Side by side

| element | V1 outcome | post hoc V1 diagnosis | V2 prospective design | V2 blind-holdout result | final EMT result |
|---|---|---|---|---|---|
| network / operating point | G1 PASS (Ybus 1.7e-16; V 1.5e-6) | — | unchanged | — | **reproduced** (V19 / I01) |
| SG transcription | not separated (G3) | exact to 1.9e-4 with an algebraic line (D1) | G3a: algebraic network; bus 30 regression + bus 36 blind; v1 rules | trajectories PASS (bus 36: ≤ 2.1e-4); **ringdown sub-check FAIL** (spurious estimator mode in both EMT and reference; Δα 0.016 > 0.005) | **G3a FAIL** (preregistered rule) |
| SG electromagnetic embedding | G3 FAIL (efd 2.26 % > 2 %) | the line's (X/ω0)dI/dt (D3: 1.4e-3 against a dynamic-line reference) | G3b: against a dynamic-line reference | NOT RUN | not established |
| GFL transcription | not separated (G4) | exact to 0.78 % with an algebraic line (D1) | G4a: v1 harness; bus 35 and 37 blind | NOT RUN (holdouts blind) | not established |
| GFL EMT interface | G4 FAIL (ideal current source; Nyquist growth) | current source into a lossless trapezoidal inductor (D4) | G4b: voltage source behind Rf–Lf; 6 checks | NOT RUN (holdouts blind); software gate SW PASS | not established |
| modal estimator (frozen v1 §5) | never exercised on IEEE-39 | — | reused unchanged | unresolved on every unit ringdown; crashes on 3/4 synthetic multi-channel signals (DV2-2) | **unfit as frozen** |
| EMT04–EMT18 portfolio | BLOCKED | — | resume unchanged if all gates pass | — | **BLOCKED** |

## 3. V2-0 ledger reconciliation (commit cf8a0b38)

- **I04 / V21.** The final ledger records the custom-GFL results as VALIDATED:
  reproduced in ANDES with the same equations (post-cumulant Phase 8, GATE 5
  PASS). The v1 EMT report had wrongly said "NOT YET TESTED". Independent phasor
  reproduction and EMT reproduction are now kept apart. EMT reproduction was not
  achieved, in v1 or in V2.
- **T08.** The Kundur holdout scores 28/29 under the preregistered rule. The miss
  is a real-pair coalescence, so the score is 28/28 for origin crossings.
- **Claim matrix.** It was regenerated from the final canonical ledger: 43 rows,
  none with EMT support.

## 4. V2 design

**Documents** (commit 0ab8efd4):
- `docs/20260911_PAREMT_EMT_PREREG_V2.md`;
- `docs/20260911_PAREMT_GFL_VOLTAGE_INTERFACE_DERIVATION.md`.

**Gates.** SW → G3a → G3b → G4a → G4b. Any failure stops the sequence.

**GFL interface.**
- An average-value voltage source behind the physical Rf–Lf filter.
- The measured filter currents are i_d and i_q; the controller integrates 9
  states.
- Command: `e_cmd = e_TX4 + j(xf/ω0)θ′i`, evaluated at the predicted states
  and at linearly extrapolated phasors.
- No numerical damping.

**Blind holdouts**, frozen before any run:
- SG: bus 36 at its canonical dispatch;
- GFL: bus 35 at S = 2 + j0.5;
- GFL: bus 37 at its canonical dispatch.

**Code.** Commit f8b8fb15 holds:
- `overlay/tx4_emt_v2.py`;
- `v2/EMTV2_{refs,sw_tests,gates,G3a_diagnostics}.py`;
- `v2/emtv2_rules.py`.

The v1 kernel is untouched.

## 5. Gate SW: PASS (`results/EMTV2/SW/sw_tests.json`)

| test | criterion | result |
|---|---|---|
| T1 identity: (P) at e_cmd equals f_i,TX4 (2000 random vectors, θ′ ∈ [−5, 5] rad/s) | ≤ 1e-12 relative | 9.1e-17 |
| T1 e_cmd = e_TX4 + j(xf/ω0)θ′i | ≤ 1e-12 | 2.2e-16 |
| T2 abc-level physical filter vs dq ODE (0.2 s, smooth drive) | ≤ 1e-7 | 2.2e-12, against an excursion of 6.3e-3 |
| T3 companion (Req, icf, Gv1) vs ParaEMT `numba_InitNet` | ≤ 1e-12 | 1.8e-16 |
| T4 static bases and signs at t = 0 (3 cases) | ≤ 1e-12 | ≤ 3.0e-16 |

The interface is implemented as derived.

## 6. Gate G3a: FAIL (`results/EMTV2/G3a/`)

**Run.** Algebraic network, the v1 SG kernel, 50 µs, Pm ×1.02 on [1.0, 1.2) s,
10 s. The reference is the quasi-static phasor model (canonical classes, Radau).

| quantity | rule | R-SG30 (regression) | B-SG36 (blind) |
|---|---|---|---|
| finite | — | yes | yes |
| equilibrium ΔP, ΔQ at 0.9 s | ≤ 1e-4 | 0, 0 | 0, 0 |
| worst state error / excursion (δ, ω, e′q, e′d, efd, pss_w, pss_l) | ≤ 0.02 | 1.9e-4 (ω) | 2.1e-4 (ω) |
| frozen `matrix_pencil` on ω − 1: α_EMT, α_ref | \|Δα\| ≤ 0.005 | −0.118368 / −0.118367 | **+6.3134 / +6.2976 (Δα 0.0158)** |
| f_EMT, f_ref | \|Δf\| ≤ 0.005 | 1.177112 / 1.177112 | 1.1051 / 1.1097 (Δf 0.0046) |
| case verdict | all rules | PASS | **FAIL** |

**Gate G3a: FAIL.** The stop rule applies: G3b, G4a and G4b were not run, and
EMT04–EMT18 stay BLOCKED. No threshold was changed and nothing was retuned.

## 7. Diagnosis (post hoc; `results/EMTV2/G3a/G3a_diagnostics.json`)

**DV2-1 — the physical ringdown agrees. The failure is in the estimator.**

The frozen `matrix_pencil` reports **resolved = False** in all four evaluations:
both cases, EMT and reference alike. Its reconstruction residual is 0.9987–0.9989
and its model order is 10 or 16. On B-SG36 it selects a spurious growing mode,
α ≈ +6.3 s⁻¹, for the reference and for the EMT run alike. That mode is not
physical: the trajectories decay.

An independent fixed-order damped-sinusoid fit of ω − 1 over [1.2, 10] s gives:

| case | α_EMT | α_ref | f_EMT = f_ref | fit residual |
|---|---|---|---|---|
| R-SG30 | −0.1182316 | −0.1182316 | 1.177050 Hz | 0.0024 |
| B-SG36 | −0.2439564 | −0.2439563 | 1.250783 Hz | 0.018 |

EMT and the reference agree on the physical mode to about 1e-7 s⁻¹. The
sub-check failed because the frozen estimator is unresolved, which the
preregistered G3a rule did not require it to be. It did not fail because of a
model discrepancy.

- **v1.** The v1 EMT02 ringdown "agreement" came from the same unresolved
  estimator; its value there happened to be the physical one.
- **V2.** The V2 prereg added "where resolved" logic only for the G4b estimator
  R. G3a was preregistered with the unconditional v1 rule, and it is judged by
  that rule.

**DV2-2 — the frozen estimator on synthetic multi-channel data.**

The input imitates the EMT05–EMT18 use: 20 channels, 30 s, a target mode plus a
damped 0.95-Hz mode, and the pulse at 1 s. `classify` is the frozen v1 function.

| true α, f | estimator A (matrix pencil) | estimator B (Hilbert) | classify |
|---|---|---|---|
| +0.127, 0.622 Hz (P4 H4-like) | **crash: LinAlgError** (Vandermonde overflow) | resolved, +0.1266, 0.6216 Hz | **ERROR** |
| −0.144, 0.637 Hz | unresolved (−0.1440, 0.6370 Hz) | resolved, −0.1468 | UNRESOLVED |
| −0.0174, 0.71 Hz (G_S-like) | **crash** | unresolved, −0.0190 | **ERROR** |
| +0.0036, 0.705 Hz | **crash** | unresolved, +0.0021 | **ERROR** |

**Mechanism.** The model order comes from the singular-value threshold 1e-4, and
it admits spurious modes with |z| ≫ 1. `z**k` then overflows in the
least-squares residue fit, which gives inf, then NaN, then an SVD failure. When
it does not crash, the spurious columns drive the residual to about 1, so the
case is never "resolved".

**Consequence.** Estimator A was frozen in v1 §5 and reused unchanged in V2, and
its α_A is the reported α_EMT. It could not have produced a STABLE/UNSTABLE
verdict for EMT05–EMT18 even if every device gate had passed.

This is a defect of the frozen v1 design. It was found only now, because v1 never
reached the portfolio experiments. Under the rules it is not repaired inside V2.

## 8. What was not run, and why

**G3b, G4a and G4b were not run** (V2 prereg §5.1, stop rule).
- The GFL blind holdouts, B1-GFL35 (S = 2 + j0.5) and B2-GFL37 (canonical
  dispatch), have never been simulated in EMT.
- Their phasor references (quasi-static and dynamic-line) were computed before
  the gates, as frozen predictions. Only these exist:
  - `results/EMTV2/refs/`;
  - the 50-µs grids in `external/paremt_runs/v2_refs/` (hashed).
- Both GFL holdouts therefore stay prospective for any continuation.
- The voltage-source GFL interface was verified only at the software level
  (gate SW). It has **not** been run in any network.

**EMT04–EMT18** remain BLOCKED, and no IEEE-39 EMT portfolio run exists.

## 9. Claim ledger

There is no change: `results/20260911_EMT_CLAIM_MATRIX.csv` (V2-0) stands.
- No claim receives EMT support, and none is EMT-refuted.
- V19 (network and operating point) remains the only EMT QUANTITATIVELY
  REPRODUCED claim.
- V20 (SG) and V21 (custom GFL) keep their final-ledger status. That status
  rests on the independent ANDES phasor reproductions, which the EMT outcome
  does not touch.

## 10. For the author: what a V3 would need

This is a decision memo, not executed.

1. **Replace or repair estimator A before any use.** Freeze a numerically safe
   pencil and pre-validate it on known signals, then commit the validation
   before any TX4 EMT result. The validation set:
   - synthetic multi-channel signals, including the four in DV2-2;
   - the frozen phasor TDS of P4 H4 and G_S, where it must reproduce the frozen
     eigenvalues.

   A safe pencil would, for example:
   - evaluate residues with |z| normalized (or with a Vandermonde built from
     time-reversed or scaled powers);
   - use a model order fixed by an information criterion;
   - drop modes with α > 5 s⁻¹ before the energy test.

   Estimator B alone is an alternative, but it is weaker near the boundary
   (unresolved at |α| ≤ 0.02 in DV2-2).
2. **SG ringdown sub-check.** Condition it on the estimator being resolved, or
   use a preregistered fixed-order fit (DV2-1 type). Keep the 2 % trajectory rule
   unchanged.
3. **Blind SG holdout.** B-SG36 has now been seen (G3a), so a V3 needs a new one
   (for example bus 38 at its canonical dispatch). The GFL holdouts B1-GFL35 and
   B2-GFL37 can be kept, because they are still blind.
4. **Everything else unchanged:**
   - the V2 interface, derivation and G4b checks;
   - EMT04–EMT18, their frozen predictions and their success rules.
5. **Stop rule unchanged.** If the voltage interface fails G4b on the blind
   holdouts, STOP. There is no other interface in the same campaign.

## 11. Reproducibility

**Rerun from commit f8b8fb15** (paths relative to `reports/poster/ias2026/research`):

```
.venv/xtool-paremt   experiments/paremt_emt/EMT_setup_tx4.py
.venv/tx3-analysis   experiments/paremt_emt/v2/EMTV2_refs.py
.venv/xtool-paremt   experiments/paremt_emt/v2/EMTV2_sw_tests.py
.venv/xtool-paremt   experiments/paremt_emt/v2/EMTV2_gates.py
.venv/xtool-paremt   experiments/paremt_emt/v2/EMTV2_G3a_diagnostics.py
```

**Byte identity.** The chain was executed twice: before the code commit, and
again from the committed HEAD. The gate, SW and diagnostics JSON/CSV outputs are
identical, apart from the provenance fields (git HEAD, hashes). The reference
arrays are identical. Details are in §12.

**Provenance.**
- `results/EMTV2/EMTV2_gates.json` and `SW/sw_tests.json`: git HEAD, ParaEMT
  upstream SHA d79d735a, working-copy diff SHA, and environment hash.
- `refs/refs_manifest.json`: the sha256 of every reference file.

## 12. Determinism record

Source: `results/EMTV2/EMTV2_rerun_identity.json`. The second execution ran
from the committed HEAD f8b8fb15 and was compared with the first.

| outputs | result |
|---|---|
| G3a checks and trajectories (CSV) | byte-identical |
| EMTV2_gates.json, sw_tests.json, G3a_diagnostics.json | identical, provenance fields excluded |
| all 10 committed reference NPZ files | identical arrays |
| refs_manifest.json | identical; file hashes excluded because the NPZ zip metadata differ |

The committed provenance records git HEAD f8b8fb15.

---

## What EMT did and did not validate (after V1 and V2)

The table separates **theorem validity** from **model-fidelity corroboration**.
The statuses are those of the final canonical ledger (V-ids; old ids in
brackets).

| Theory / claim | Mathematical status | Phasor evidence | Independent phasor evidence | EMT evidence (V1 + V2) | Final status |
|---|---|---|---|---|---|
| V01, V02, V08, V22, V26, T04, T06, T08 (theorems) | proved | yes | where applicable | not applicable | PROVED |
| V03 P4 trap, κ = 4 [B02] | — | PCV02; G2 TDS | ANDES Phase 8, 32/32 | BLOCKED | VALIDATED (model-specific witness) |
| V06, V07 policy dependence [B01] | structure proved | clean g-only counterfactual | ANDES Phase 8, P4 vs G_S | BLOCKED | VALIDATED |
| V10, V11 port derivative, line ranking [T10] | proved (V10) | PCV04 | ANDES R3 11/12, ρ 1.00 | BLOCKED | PROVED + VALIDATED; VALIDATED (local ranking) |
| V13, V14 remediation, governors | — | PCV04; FC03; G2 | partial ANDES | BLOCKED | VALIDATED (model-specific); BENCHMARK-SPECIFIC |
| V15–V18 robustness | — | PCV05 | — | BLOCKED / not targeted | as in the ledger |
| V30, B07, B10 nonlinear scope | normal form | G2 TDS 32/32 | — | BLOCKED | NUMERICALLY VALIDATED |
| V19 network and operating point [I01] | — | reference | ANDES; pandapower after the tap fix | **EMT QUANTITATIVELY REPRODUCED** (V1 G1) | CONDITIONAL (pandapower translation) |
| V20 SG dynamics [I02, I03] | — | internal | ANDES R2/R3, E1 | V1 G3 FAIL; V2 G3a FAIL (trajectories 2e-4; estimator sub-check) → EMT UNRESOLVED | VALIDATED (by ANDES) |
| V21 custom GFL in a second simulator [I04] | — | project DAE | **ANDES Phase 8, GATE 5 PASS** | V1 G4 FAIL; V2 G4a/G4b NOT RUN (blind) → EMT UNRESOLVED | VALIDATED as a reproduction of computation; EMT reproduction not achieved |
| cumulant and closure diagnostics, V05, V12, B04, B06, negatives | as in the ledger | as in the ledger | — | not applicable / BLOCKED | unchanged |
