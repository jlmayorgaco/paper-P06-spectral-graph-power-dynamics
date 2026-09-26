# EMT Scientific Qualification Campaign (EMT-SQ1) — preregistration

Date: 2026-09-12. Branch `research/paremt-emt-validation`, not pushed.

The author approved SQ1 as **one final, separately justified campaign**. It is
committed **before** any SQ1 qualification computation. The only computations
beforehand were the SQ1-1 domain derivation (frozen phasor models, no
estimator) and a generator self-check (dev seed 2222, no estimator).

## 0. Purpose and history

**Purpose.** Qualify the **already-frozen** V3 estimator for the actual
multichannel IEEE-39 scientific measurement domain. Then, **only if it
passes**, run the untouched blind device and portfolio holdouts.

SQ1 tests the existing estimator; it does not develop it.

**History (permanent; nothing rewritten):**
- V1: G3 FAIL, G4 FAIL.
- V2: G3a FAIL (estimator-induced).
- V3: E0 FAIL (criterion A5, the S6 single-channel ringdowns, 14/15).

**About V3 S6.** It remains a stress test that documents a limitation of the
estimator. It is **not** erased or relabelled. SQ1 reports a new stress stratum
alongside, and qualification is decided on the scientific domain only.

**Wording if SQ1 passes (required):**

> "Three earlier preregistered campaigns stopped at infrastructure or
> measurement gates. A final task-specific qualification campaign froze the
> estimator and model interfaces in advance, passed new blind synthetic,
> phasor-TDS and device holdouts, and only then evaluated the EMT scientific
> hypotheses."

## 1. SQ1-0 — frozen instrument and interfaces

**Estimator.** `experiments/paremt_emt/v3/emtv3_estimator.py` at commit
f8e0d44f, sha256 `849472d8804c694eee69701dd251aaa69d506878007a3a4fc88fde1fb3b0817f`
(the V3 spec plus V3-AM1). **No estimator code changes are allowed.** Frozen
elements:
- the fitting model and the V3-AM1 weighting;
- the QR solve;
- the jackknife CI;
- R1–R4;
- the tracker width ±0.15 Hz;
- the sentinel;
- EPS = 0.002;
- the adaptive run length 30 → 60 → 90 s.

Every SQ1 script asserts this hash.

**Frozen interfaces:**
- the v1 SG kernel (`tx4_emt.py`);
- the V2 GFL voltage interface (`tx4_emt_v2.py`, unchanged);
- no numerical damping.

**Frozen holdout tooling:** `EMTV3_tds_traces.py` and `EMTV3_E0_tds.py`
(V3 commit 70d5be8e), used unchanged for T1–T9.

## 2. SQ1-1 — scientific measurement domain (derived from frozen predictions only)

**Target modes.** They come from `results/EMT_PRED/phasor_predictions.csv`:
`band_re` and `band_hz` of every frozen scientific case EMT04–EMT16 (126 cases).
- **Hull:** α ∈ [−0.2290, +0.1990] s⁻¹ and f ∈ [0.5566, 1.0467] Hz.
- **Scientific domain** (hull rounded outward to 0.01 s⁻¹ and 0.05 Hz):
  **α ∈ [−0.23, +0.20] s⁻¹** and **f ∈ [0.55, 1.05] Hz**, which lies inside the
  preregistered 0.2–1.2 Hz band.

**Channels, observability, nuisance structure, duration.** These come from the
**modal library** `results/EMTSQ1/domain/modal_library.npz`, produced by
`experiments/paremt_emt/sq1/EMTSQ1_modal_library.py` (tx3-analysis; no
estimator).
- **Models.** The frozen phasor models of EMT05, EMT08, EMT09, EMT10, EMT11 and
  EMT13 that the frozen G2 builder constructs, **excluding the nine T1–T9
  holdout models** (and EMT04 G_S, which equals T6). That leaves 26 models.
- **Method.** The linear channel residues of every transverse mode (decay
  slower than 5 s⁻¹) to the preregistered pulse:
  `r_cm = (C_c v_m)(w_m·B) ΔP (e^{0.2λ} − 1)/λ`, where C is the v1 §5 channel
  map and B = ∂ẋ/∂P_load20.
- **Validation.** The linear responses reproduce the non-holdout development
  TDS trace P4 BASE to 0.30 % median and 0.83 % max relative error (oscillatory
  part).

| property | frozen-model value |
|---|---|
| channels | 28–32 |
| per-channel target energy share, window [2.2, 30] s | p05 0.0049, p25 0.041, median 0.944, p75 0.992, p95 0.999; **30 % of channels < 0.1** (low observability) |
| nuisance modes ≥ 1 % median energy share | 16 in 26 models: Δf ∈ [−0.37, +0.27] Hz, Δα ∈ [−0.63, −0.003] s⁻¹; in 10 of 26 models a nuisance family (the 0.64-Hz family beside a 0.91-Hz rightmost target) carries **85–92 %** of the median channel energy |
| record duration | 30 s, with preregistered extensions to 60 and 90 s |

**Conclusion.** The minority-target, dominant-nuisance structure is a property
of the actual scientific domain. It is reproduced, not avoided.

## 3. SQ1-2 — fresh blind synthetic qualification

**Generator:** `experiments/paremt_emt/sq1/EMTSQ1_synth.py`, committed with this
document. **Seed 20260913** (never used before).
- **Template.** A random library model. Its channel residues and all its
  non-target modes are kept, with each residue multiplied by U(0.8, 1.25) and
  e^{jU(−0.2, 0.2)}.
- **Target.** The model's target mode is replaced by (α, f) drawn from the
  domain. f is redrawn if any other complex mode of median share ≥ 5 % lies
  within 0.03 Hz.
- **Trend and noise.** Per-channel affine trend `c ~ U(−1, 1)·rms_c` and
  `d ~ U(−0.01, 0.01)·rms_c` per s. White noise with σ = ν·rms_c,
  ν ~ log-U[1e-4, 3e-3].
- **Record.** The response starts at t = 1.2 s; records are 90 s at 1 kHz.
- **Prediction.** f_pred = f_true + U(−0.03, 0.03) Hz.
- **Evaluation.** The frozen targeted tracker, sentinel and adaptive rule.

| stratum | n | content |
|---|---|---|
| SD1 domain | 130 | α ~ U(−0.23, 0.20), f ~ U(0.55, 1.05) |
| SD2 near boundary | 40 | α ~ U(−0.01, +0.01), f ~ U(0.55, 1.05) |
| SD3 low observability | 40 | as SD1, with 30 % of channels having the target residue × U(1e-3, 1e-1) |
| **ST1 stress** (reported only) | 20 | V3-S6-like single-channel, heavily damped: α ~ U(−0.5, −0.3), f ~ U(1.0, 1.6), slow real mode; unit-test settings |
| **ST2 stress** (reported only) | 15 | as SD1, with α ~ U(−0.6, −0.3) (outside the domain) |

**Primary gates on SD1–SD3 (210 cases), fixed here and not altered after
results:**

| gate | rule |
|---|---|
| GA | **zero wrong-sign classifications** (STABLE with α > 0 or UNSTABLE with α < 0) for \|α\| ≥ 0.005 s⁻¹ |
| GB | **≥ 95 % resolved** (R1–R4 all hold) among cases with \|α\| ≥ 0.01 s⁻¹ |
| GC | p95 \|α̂ − α\| ≤ 0.005 s⁻¹ over resolved cases |
| GD | p95 \|f̂ − f\| ≤ 0.005 Hz over resolved cases |
| GE | no exception or crash |

- **Near the boundary** (|α| < 0.01), UNRESOLVED is acceptable.
- **Reported, not gated:** the correct-conclusive fraction; the stress strata.

**If SQ1-2 fails: STOP the ParaEMT line permanently.**

## 4. SQ1-3 — the unused phasor-TDS holdout T1–T9 (only if SQ1-2 passes)

**Traces.** Generated for the first time by the frozen `EMTV3_tds_traces.py
holdout` and evaluated by the frozen `EMTV3_E0_tds.py holdout`. There is no
regeneration or selection based on estimator output. The cases and eigenvalues
are those of prereg V3 §4.2.

**Required:**
- **no wrong-sign verdict** in any case;
- every case with |α_eig| ≥ 0.01 is **resolved**;
- every resolved case: |α̂ − α_eig| ≤ 0.005 s⁻¹ and |f̂ − f_eig| ≤ 0.005 Hz;
- **no UNEXPECTED_MODE flag.** The phasor model has no unexpected mode, so a
  flag is a false alarm.
- T7 (|α| = 0.0013) may be UNRESOLVED.

**If T1–T9 fail: STOP ParaEMT permanently.**

## 5. SQ1-4/5 — synchronous-machine device gates (only if SQ1-3 passes)

**Principle.** Device-equation equivalence is established by **equilibrium,
exact parameters and bases, and state-trajectory agreement**. The ringdown
estimator is **secondary corroboration**.

If the trajectories pass and the estimator is unresolved on a single-channel
ringdown, the result is **DEVICE PASS / MODAL MEASUREMENT UNRESOLVED**. This
separates two questions that V1 and V2 conflated; it does not relax a threshold.

**Cases:**

| case | role |
|---|---|
| SG30 | regression |
| SG36 | seen; diagnostic/regression |
| **SG38** | **blind**: bus 38, S = 7.64783381 + j1.23275808, never simulated |

The unit-test system is the V2 one.

- **G3a (algebraic network, the v1 SG kernel)**
  - finite;
  - equilibrium |ΔP|, |ΔQ| ≤ 1e-4 at 0.9 s;
  - **exact parameters and bases:** every parameter of the EMT machine row
    (ra, xd, xq, x′, T′d0, T′q0, M, D, KA, TA, KS, T4, T5, T6, w) equals the
    canonical `SynchronousMachine` parameters to 1e-12 relative;
  - all 7 states within the unchanged 2 % trajectory rule against the
    quasi-static reference over [1, 10] s (raw 1-ms samples).
  - **Secondary** (reported): the frozen estimator (unit mode) on ω − 1 for EMT
    and reference, labelled
    - MODAL AGREES (both resolved, |Δα|, |Δf| ≤ 0.005);
    - MODAL DISAGREES (both resolved, outside that);
    - MODAL MEASUREMENT UNRESOLVED.
- **G3b (the physical EMT R–L line)**
  - finite, equilibrium, and all 7 states within 2 % against the
    **dynamic-line** reference.
  - **Reported separately:**
    - EMT against the dynamic-line reference: the device-transcription error in
      the EMT embedding;
    - the dynamic-line against the quasi-static reference: the electromagnetic
      line correction, not a transcription error;
    - EMT against the quasi-static reference;
    - 25 µs.
  - **Secondary** modal labels as in G3a.

**On failure:** a G3a or G3b failure stops SQ1.

## 6. SQ1-6 — GFL gates (only if G3a and G3b pass)

**Interface.** The V2 interface, **unchanged**:
- average-value voltage source behind the physical Rf–Lf filter;
- physical filter currents as i_d/i_q;
- 9 integrated controller states;
- linearly extrapolated voltage and current;
- no numerical damping.

**Cases.** R-GFL30 is regression; **B1-GFL35** (S = 2 + j0.5) and
**B2-GFL37** (canonical S) are blind (never run).

- **G4a (algebraic network; the v1 11-state harness)**
  - the v1 EMT03 rules unchanged: finite, equilibrium ≤ 1e-4, 11 states within
    2 %, PLL lock ≤ 1e-3 at T, g = 0 q-tracking ≤ 1e-3;
  - plus an **exact parameter/base check**: the kernel row equals the
    canonical `ConverterParameters` defaults, the leak and w, to 1e-12.
- **G4b (physical EMT embedding)**
  - the V2 G4b checks 1–5 **unchanged**:
    1. finite and no-event stationarity ≤ 1e-4;
    2. no growing Nyquist indicator;
    3. 50 vs 25 µs within 0.01·excursion;
    4. P/Q equilibrium and 15-Hz trajectories within 2 %;
    5. all 11 quantities within 2 % (15-Hz LP), PLL lock, g = 0 q-tracking.
  - **Secondary** (reported): the frozen estimator on the 11 LP'd quantities
    (V3 settings), with modal labels as in G3a.

**On failure.** If G4a or G4b fails: STOP ParaEMT **permanently**. No Norton
interface, no damping workaround, no new interface.

## 7. SQ1-7 — scientific core (only if all previous gates pass)

The core runs exactly as prereg V3 §7:
- the IEEE-39 v1 realization with every GFL on the V2 interface;
- the frozen estimator, with the targeted tracker at the committed band_hz;
- every v1 success rule unchanged.

| block | content | rule |
|---|---|---|
| EMT04 | 100/50/25 µs on base, 30+33+35 and H4 at P4, and H4 at G_S | **G5** (50 vs 25 µs): verdict invariant, \|Δf\| ≤ 0.01 Hz, \|Δα\| ≤ 0.01 where resolved |
| EMT04 holdouts | τ_m 0.5 ms (v1 G2); EMT-native stator (descriptive); no-event stationarity (reported) | — |
| EMT05 | the 16 P4 subsets | **primary:** 15 proper subsets STABLE and H4 UNSTABLE, i.e. **H_EMT(P4) = {{30, 33, 35, 37}}**; secondary MAE ≤ 0.03 and max\|Δf\| ≤ 0.05; amplitude checks 1 % and 4 % |
| EMT06 | the 8 contextual signs | all equal the phasor signs (resolved fits only; otherwise UNRESOLVED) |
| EMT07/08 | the g-only boundary, bracket + ≤ 4 bisections | **g\*_EMT with \|g\*_EMT − 0.20768\| ≤ 0.02**; primary: a resolved UNSTABLE point below and a STABLE point above |
| EMT10 | P4 → 0.08806 → 0.20768 → 0.25 | α̂ decreases monotonically, g = 0.25 is STABLE, and \|g\* − 0.20768\| ≤ 0.02 |

- **Bisection f_pred.** Bisection points have no committed prediction. f_pred
  there is the linear interpolation of the committed band_hz of the neighbouring
  grid points. It is a search neighbourhood only, and it is declared.
- **Stop.** If G5 fails: STOP, with no amendment in SQ1.
- **Interim report** after EMT10.

## 8. SQ1-8 — decision gate after EMT10

- **Corroborated** iff all five hold:
  - EMT04 G5 PASS;
  - EMT05 primary PASS;
  - EMT06 PASS;
  - EMT07/08 primary **and** target PASS;
  - EMT10 PASS.
- **If corroborated,** proceed in this campaign, with the v1 §8 rules and the
  frozen estimator, to:
  - EMT09 (P_inf);
  - EMT11 (k boundary);
  - EMT12 (12 line reinforcements);
  - EMT13 (condenser);
  - EMT14 (governors);
  - EMT15 (nonlinear recovery);
  - EMT16 (robustness).
- **Otherwise (any FAIL or UNRESOLVED): STOP.** The model is not retuned. The
  disagreement or non-resolution is the scientific result.

## 9. Determinism and reporting

**Determinism.**
- SQ1-2 is rerun once, and its summary must be identical.
- EMT05 and the EMT08 grid are rerun once, with byte-identical summary CSVs.

**Reports.**
- `docs/20260912_PAREMT_EMT_SQ1_REPORT.{md,tex,pdf}`, preserving V1 FAIL, V2
  FAIL and V3 E0 FAIL, then SQ1 separately.
- An interim scientific report after EMT10, if it is reached.

No push.
