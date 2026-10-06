# EMT Scientific Qualification Campaign (EMT-SQ1) — report

**Status: SQ1 STOPPED at SQ1-2 (blind synthetic qualification). By the
preregistered rule, the ParaEMT EMT line is STOPPED PERMANENTLY.** Branch
`research/paremt-emt-validation`, not pushed.

**SQ1-2**, run on the frozen V3 estimator (seed 20260913, 210 domain cases):

| gate | result |
|---|---|
| GA — wrong-sign verdicts | **0** (pass) |
| GB — resolved, among cases with \|α\| ≥ 0.01 | **78.4 %** (fail; ≥ 95 % required) |
| GC — p95 α error, resolved cases | **0.0066 s⁻¹** (fail; ≤ 0.005 required) |
| GD — p95 f error, resolved cases | 0.0017 Hz (pass) |
| GE — crashes | none (pass) |

**Not run:**
- the phasor-TDS holdout T1–T9;
- the device gates (the blind SG38, GFL35 and GFL37 were never simulated);
- EMT04–EMT10.

## 1. History (binding; nothing rewritten)

| campaign | preregistration | permanent outcome |
|---|---|---|
| V1 | c2947bd8 (+ A1) | G1 PASS; **G3 FAIL, G4 FAIL** (interface / device gates) |
| V2 | 0ab8efd4 | SW PASS; **G3a FAIL** (caused by the frozen v1 estimator) |
| V3 | 6d339006 (+ AM1) | **E0 FAIL** (criterion A5: single-channel ringdowns, 14/15) |
| SQ1 | 0555f587 | **SQ1-2 FAIL** (GB 78.4 %, GC 0.0066) → ParaEMT line stopped permanently |

**Statement:**

> Four preregistered campaigns attempted an EMT corroboration of the IEEE-39
> results.
>
> - The first three stopped at infrastructure or measurement gates.
> - The final, task-specific qualification campaign froze the estimator and the
>   model interfaces in advance. It tested the estimator on a fresh blind
>   synthetic set derived from the frozen phasor models of the actual
>   scientific domain.
> - The estimator made no wrong-sign call. It did not reach the preregistered
>   resolution and accuracy targets, and the failures concentrate where the
>   target mode is a minority mode beside a dominant, similarly damped
>   neighbouring mode family.
> - By the preregistered rule, the ParaEMT line stopped before any blind device
>   or portfolio holdout was simulated.
>
> No EMT evidence for or against the portfolio, policy or design claims exists.

None of the four failures is described as invalid.

## 2. SQ1-0 — frozen instrument

**Estimator.** `experiments/paremt_emt/v3/emtv3_estimator.py`, commit f8e0d44f,
sha256 `849472d8…3b0817f`. It was unchanged throughout SQ1; every SQ1 script
asserts the hash.

**Frozen interfaces:**
- the v1 SG kernel;
- the V2 GFL voltage interface;
- no damping.

## 3. SQ1-1 — scientific measurement domain (frozen predictions only)

**Target modes.** Taken from the frozen `band_re` and `band_hz` of all 126
EMT04–EMT16 cases.
- **Hull:** α ∈ [−0.229, +0.199] s⁻¹ and f ∈ [0.557, 1.047] Hz.
- **Domain:** α ∈ **[−0.23, +0.20]** s⁻¹ and f ∈ **[0.55, 1.05]** Hz.

**Structure.** Taken from a **linear modal library** of 26 frozen non-holdout
phasor models (`results/EMTSQ1/domain/`):
- the channel residues of every transverse mode for the preregistered pulse,
  through the v1 §5 channel map;
- validated against the development TDS trace P4 BASE: 0.30 % median and
  0.83 % max relative error;
- 28–32 channels;
- per-channel target energy share: median 0.94, and **30 % of channels below
  0.1**;
- 16 nuisance modes of at least 1 % share;
- **in 9–10 of the 26 models,** the rightmost in-band mode, the preregistered
  target (the 0.91-Hz family), is a **minority mode**. A 0.64-Hz family 0.26 Hz
  away decays only 0.003–0.07 s⁻¹ faster and carries 85–92 % of the median
  channel's energy. These are P4 subsets that contain 30, 35 or 37 without the
  full H4.

**V3 S6** remains a stress test that documents a limitation of the estimator.
It is not relabelled.

## 4. SQ1-2 — blind synthetic qualification (seed 20260913)

**Generator.** Each case takes a random template from the library, keeping its
residues (perturbed ×U(0.8, 1.25) and e^{jU(−0.2, 0.2)}) and all of its
nuisance modes. The target is redrawn from the domain. Affine trends and noise
(σ = ν·rms, ν ~ log-U[1e-4, 3e-3]) are added, and f_pred = f + U(±0.03).
Records are 90 s, and the frozen tracker, sentinel and adaptive rule are used.

| stratum | n | resolved | correct conclusive | wrong sign | p95 \|Δα\| | p95 \|Δf\| |
|---|---|---|---|---|---|---|
| SD1 domain | 130 | 82.3 % | 80.8 % | 0 | 0.0080 | 0.0022 |
| SD2 near boundary (\|α\| ≤ 0.01) | 40 | 90.0 % | 60.0 % (the rest UNRESOLVED, allowed) | 0 | 0.0014 | 0.0003 |
| SD3 low observability | 40 | 67.5 % | 67.5 % | 0 | 0.0052 | 0.0014 |
| **ST1 stress** (V3-S6-like, reported only) | 20 | 100 % | 100 % | 0 | 0.0044 | 0.0010 |
| **ST2 stress** (α ∈ [−0.6, −0.3], reported only) | 15 | 0 % | 0 % | 0 | — | — |

| gate (SD1–SD3) | rule | result | verdict |
|---|---|---|---|
| GA | 0 wrong-sign verdicts, \|α\| ≥ 0.005 | 0 | pass |
| **GB** | **≥ 95 % resolved, \|α\| ≥ 0.01** | **131/167 = 78.4 %** | **FAIL** |
| **GC** | **p95 \|Δα\| ≤ 0.005** | **0.0066** | **FAIL** |
| GD | p95 \|Δf\| ≤ 0.005 | 0.0017 | pass |
| GE | no crash | none | pass |

- CI coverage of α_true among resolved cases: 94.7 %. Sixteen resolved cases
  have |Δα| > 0.005.
- Determinism rerun: identical code and seed; the CSV is **byte-identical** (sha256 `30f71465…4ab62f`, 245 rows) and the JSON is identical except `wall_s` (`results/EMTSQ1/SQ1-2/sq1_rerun_identity.json`).

**SQ1-2: FAIL. STOP the ParaEMT line permanently** (prereg SQ1 §3).

## 5. Post hoc diagnosis (not a gate; `sq1_posthoc_breakdown.json`)

**By template.** Among the domain cases with |α| ≥ 0.01:

| template class | n | resolved | p95 \|Δα\| (resolved) | p95 \|Δf\| (resolved) |
|---|---|---|---|---|
| **dominant-nuisance** templates (9 P4 subsets with a 0.64-Hz family carrying > 50 % median share) | 54 | **40.7 %** | **0.059** | **0.127** |
| all other templates | 113 | 96.5 % | 0.0056 | 0.0015 |

**Unresolved reasons** (36 cases): 13 were at the bound together with a wide CI,
11 had a wide CI only, and 12 involved the residual rule R3.

**Reading.** The frozen estimator is a single-mode targeted tracker (V3
specification, as the V3 preregistration required).
- **Isolated or dominant target:** resolution ≈ 96.5 %, with no wrong sign. The
  accuracy is marginal: p95 0.0056 s⁻¹ against the 0.005 target.
- **Minority target beside a dominant, similarly damped mode 0.26 Hz away:**
  - it cannot isolate the target within the effective record length;
  - the jackknife guard withholds most of these cases (UNRESOLVED);
  - the remainder are biased.
- **Relevance.** This is precisely the structure of about half of the EMT05 P4
  subsets. The failure is therefore **task-relevant, not an artefact**: an EMT05
  run with this instrument could not have established
  H_EMT(P4) = {{30, 33, 35, 37}}, because many stable proper subsets would have
  been UNRESOLVED.

**Stress strata** (reported only):
- ST1 (heavily damped single-channel, like V3 S6): 20/20 resolved.
- ST2 (α ∈ [−0.6, −0.3], outside the domain): 0/15 resolved, with no wrong
  sign. The CI is too wide for strongly damped multichannel modes.

## 6. What was not run (and now never will be, in this line)

| item | status |
|---|---|
| SQ1-3 phasor-TDS holdout T1–T9 | NOT RUN; the traces were never generated |
| G3a/G3b (SG30, SG36, **SG38 blind**) | NOT RUN; SG38 was never simulated |
| G4a/G4b (GFL30, **GFL35, GFL37 blind**) | NOT RUN; the V2 voltage interface was never run in a network |
| EMT04, EMT05, EMT06, EMT07/08, EMT10; the decision gate; EMT09, EMT11–EMT16 | NOT RUN |

**Prepared but never executed** (committed with "never executed" banners, and
guarded by assertions on the preceding gates):
- `sq1/EMTSQ1_refs.py`;
- `sq1/EMTSQ1_tds.py`;
- `sq1/EMTSQ1_gates.py`;
- `sq1/EMTSQ1_science.py`.

## 7. Final state of the EMT evidence

- **Network and operating point (V19/I01):** EMT QUANTITATIVELY REPRODUCED
  (V1 G1). This is the only EMT result.
- **SG and GFL device equivalence in EMT:** not established by any gate.
  - The post hoc diagnostics (V1 D1, V2 DV2-1) indicate exact transcription.
  - They are not gate evidence.
- **All portfolio, policy, design, remediation and nonlinear claims:** no EMT
  evidence either way.
  - They keep their phasor status.
  - They also keep their independent ANDES reproduction: V03, V07 and V21 for
    the custom GFL (Phase 8, GATE 5 PASS).
- **Manuscript.**
  - The existing wording "nonlinear phasor-domain simulation (not EMT)" remains
    correct and required.
  - The network reproduction in ParaEMT may be mentioned.
  - The EMT attempt may be disclosed in the limitations.

  Suggested sentence:

  > *"An electromagnetic-transient corroboration in ParaEMT was preregistered
  > four times; the network and operating point were reproduced, but the
  > campaigns stopped at device-interface and modal-measurement qualification
  > gates before any portfolio result, so no EMT evidence is claimed."*

## 8. Reproducibility

**Chain** (commit 0555f587; frozen estimator f8e0d44f):

```
.venv/tx3-analysis   experiments/paremt_emt/sq1/EMTSQ1_modal_library.py    (SQ1-1 domain)
.venv/xtool-paremt   experiments/paremt_emt/sq1/EMTSQ1_synth.py check      (generator self-check, no estimator)
.venv/xtool-paremt   experiments/paremt_emt/sq1/EMTSQ1_synth.py blind      (SQ1-2)
```

---

## What EMT did and did not validate (final)

The table separates **theorem validity** from **model-fidelity corroboration**.

| Theory / claim | Mathematical status | Phasor evidence | Independent phasor evidence | EMT evidence (V1 + V2 + V3 + SQ1) | Final status |
|---|---|---|---|---|---|
| V01, V02, V08, V22, V26, T04, T06, T08 (theorems) | proved | yes | where applicable | not applicable | PROVED |
| V03 P4 trap, κ = 4 [B02] | — | PCV02; G2 TDS | ANDES Phase 8, 32/32 | none (stopped before EMT05) | VALIDATED (model-specific witness) |
| V06, V07 policy dependence [B01] | structure proved | clean g-only counterfactual | ANDES Phase 8 | none | VALIDATED |
| V10, V11 port derivative, line ranking | proved (V10) | PCV04 | ANDES R3 11/12 | none | PROVED + VALIDATED; VALIDATED (local) |
| V13, V14 remediation, governors | — | PCV04; FC03; G2 | partial ANDES | none | VALIDATED (model-specific); BENCHMARK-SPECIFIC |
| V15–V18; V30, B07, B10 | — / normal form | PCV05; G2 TDS | — | none | as in the ledger |
| V19 network and operating point [I01] | — | reference | ANDES; pandapower after the tap fix | **EMT QUANTITATIVELY REPRODUCED** (V1 G1) | CONDITIONAL (pandapower translation) |
| V20 SG dynamics [I02, I03] | — | internal | ANDES R2/R3, E1 | V1 G3 FAIL; V2 G3a FAIL; V3/SQ1 not reached → EMT UNRESOLVED | VALIDATED (by ANDES) |
| V21 custom GFL [I04] | — | project DAE | **ANDES Phase 8, GATE 5 PASS** | V1 G4 FAIL; later not reached → EMT UNRESOLVED | VALIDATED as a reproduction of computation; EMT reproduction not achieved |
| cumulants, closure, V05, V12, B04, B06, negatives | as in the ledger | as in the ledger | — | not applicable / none | unchanged |
