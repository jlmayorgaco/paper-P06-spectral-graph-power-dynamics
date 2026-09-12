# ParaEMT EMT-validation campaign — preregistration V3 (instrument qualification first)

Date: 2026-09-12. Branch `research/paremt-emt-validation`, not pushed. The author
approved it after the V2 report.

It is committed **before any V3 numerical evaluation**: no V3 estimator run, no
synthetic generation, no TDS trace generation, and no ParaEMT run.
- Prerequisite: V3-0 audit `docs/20260912_EMT_ESTIMATOR_FAILURE_AUDIT.md`
  (commit 4e35164e).

## 0. History (binding; never rewritten)

- **V1:** G3 FAIL, G4 FAIL, EMT04–EMT18 BLOCKED (prereg c2947bd8).
- **V2:** SW PASS, G3a FAIL under its preregistered rule, G3b/G4a/G4b NOT RUN,
  EMT04–EMT18 BLOCKED (prereg 0ab8efd4; results 632374dc).
- **V3.** It was designed after the V1 and V2 failures, the V2 diagnostics and the
  V3-0 audit. Only the holdouts frozen here provide prospective evidence:
  - the E0 blind synthetic set (seed below);
  - the E0 phasor-TDS instrument holdout set;
  - B-SG38;
  - B1-GFL35 and B2-GFL37 (never run in V2).

  V3 is not a repeat of the V1 blind campaign.

**Wording if V3 passes:**

> "V1 and V2 failed their preregistered gates. Subsequent diagnostics
> identified independent interface and estimator defects. V3 addressed those
> defects prospectively and was evaluated on new holdouts."

Never "the previous failures were invalid".

**Frozen files that are not modified:**
- `experiments/paremt_emt/emt_estimator.py` (V1/V2 provenance);
- the v1 kernel `overlay/tx4_emt.py`;
- the V2 kernel `overlay/tx4_emt_v2.py`;
- all V1 and V2 results.

## 1. Order of work and stop rules

1. **V3-1.** This document, committed.
2. **Estimator implementation**, then an **implementation commit**, before any
   blind data exist.
   - During development, code may be exercised only on:
     - the four DV2-2 cases (regression; no longer blind);
     - a **development synthetic set with seed 1111** (generator of §4, strata
       S1–S6, 30 cases). It is never used for qualification.
   - Implementation bugs may be fixed during development.
   - Any change to the **specification** in §2–§3 (model, bounds, window, rules,
     thresholds) would be a numbered amendment, committed before E0 and disclosed.
3. **Gate E0** (instrument qualification). Only after the implementation commit:
   - generate the blind synthetic set (seed 20260912);
   - generate the phasor-TDS holdout traces;
   - run the estimator on both.

   E0 passes only if the synthetic **and** TDS criteria pass. **If E0 fails: STOP
   V3.** No other estimator is tried in this campaign.
4. **Device gates** (ParaEMT, unit tests), only if E0 passes: G3a → G3b → G4a →
   G4b. **The first failure stops V3's ParaEMT scientific validation.** Inside V3
   there is:
   - no numerical damping;
   - no Norton interface;
   - no ideal-current-source interface;
   - no other converter realization.

   There will be no V4 merely to obtain agreement.
5. **Scientific blocks**, only if E0, G3a, G3b, G4a and G4b all pass, in this
   order: EMT04 → EMT05 → EMT06 → EMT07/08 → EMT10. Then **STOP AND REPORT**.

   The k boundary (EMT11), the 12 lines (EMT12), the condenser (EMT13), the
   governors (EMT14), nonlinear recovery (EMT15), the robustness holdout (EMT16),
   P_inf (EMT09) and Kundur (EMT17) wait for the author's decision.
6. **Report:** `docs/20260912_PAREMT_EMT_V3_REPORT.{md,tex,pdf}`. No push.

## 2. The V3 instrument: multi-channel bounded variable-projection damped-mode estimator ("VP")

**Implementation:** a new file, `experiments/paremt_emt/v3/emtv3_estimator.py`.
It does not import `emt_estimator.py`.

### 2.1 Input and preprocessing

**Input.** A uniformly sampled record: t (step h_in; 0.1 s must be an integer
multiple of h_in) and a channel matrix Y (n_ch × N). Also a window [T0, T1], a
search box α ∈ [α_lo, α_hi] and f ∈ [f_lo, f_hi], and a channel-set mode.

**Preprocessing.**
1. **Block-mean decimation to 10 Hz.** Take non-overlapping 0.1-s blocks
   [T0 + kΔ, T0 + (k+1)Δ), Δ = 0.1 s, fully inside the window. Use the sample
   mean of each block, time-stamped at the block centre.
   - Block averaging is linear and time-invariant, so it multiplies every
     exponential e^{st} by a constant complex factor and **leaves s unchanged**.
   - It creates no edge transients.
   - It nulls 60/120-Hz ripple and step-rate content exactly.
2. **Centred time.** `τ = t − (T0 + T1)/2`.
3. **Channel normalization.** For each channel, remove its least-squares affine
   fit (1, τ) over the window, only to measure its oscillatory RMS r_c.
   - Drop the channel if `r_c < 1e-12` or `r_c < 1e-8 · max_c r_c`.
   - Otherwise scale the original channel by 1/r_c. The affine terms remain in
     the model.

### 2.2 Model and variable projection

For channel c:

`y_c(τ) = e^{ατ}[a_c cos 2πfτ + b_c sin 2πfτ] + c_c + d_c τ + ε_c`

**Linear sub-problem.** For fixed (α, f), the design matrix
`B = [e^{ατ}cos, e^{ατ}sin, 1, τ]` (M × 4) is **common to all channels**.
- Each column is scaled to unit 2-norm.
- `(a_c, b_c, c_c, d_c)` for all channels are solved at once by a Householder QR
  of B (numpy `qr`, reduced) with multiple right-hand sides.
- A rank check requires the smallest |R_ii|/|R_11| ≥ 1e-12; otherwise the point
  is infeasible.
- Projected residual: `R(α, f) = Y − Q Qᵀ Y`, and the objective is
  `J(α, f) = ‖R‖_F²`.

**No full-record z^k Vandermonde and no exp(2αT) energy computation** appear
anywhere.

**Bounds.**
- Targeted tracker: α ∈ [−0.8, +0.8] s⁻¹.
- f per the caller: targeted neighbourhood (§3.1), sentinel band (§3.2) or
  unit-test band (§6).

**Optimization.**
1. **Grid:** α step 0.02 s⁻¹ and f step 0.005 Hz over the box, evaluating J.
2. **Refinement:** bounded refinement (`scipy.optimize.least_squares`,
   `method="trf"`, the flattened R as residual vector, `xtol = ftol = gtol = 1e-12`)
   from the 3 lowest grid local minima.
3. **Result:** the lowest refined J.
- Converged means `status > 0` and a finite result.
- Interior means α̂ is more than 0.01 s⁻¹ from its bounds and f̂ is more than
  0.002 Hz from its bounds.

### 2.3 Uncertainty (jackknife; preregistered)

1. **Channel jackknife** (if n_ch ≥ 3): delete one channel and refit locally
   (least_squares started at (α̂, f̂), same bounds).
   `SE_ch = √((n−1)/n · Σ (θ_i − θ̄)²)`.
2. **Time-block jackknife:** split the window into K = 5 contiguous equal
   blocks, delete one block and refit locally. The least-squares fit accepts
   the gap, because τ is explicit.
   `SE_t = √((K−1)/K · Σ (θ_k − θ̄)²)`.
3. **Combined:** `SE_α = max(SE_ch, SE_t, 5e-4 s⁻¹)` and
   `SE_f = max(SE_ch, SE_t, 5e-4 Hz)`.
4. **CI:** the 95 % CI is `θ̂ ± 1.96·SE`.

### 2.4 Resolved flag (every rule must hold)

- **R1.** The refinement converged and every output is finite.
- **R2.** The optimum is interior (§2.2).
- **R3.** The median over retained channels of `‖r_c‖/‖y_c − affine_c‖` is
  ≤ 0.9. This is the per-channel normalized residual: the mode explains at least
  19 % of the oscillatory energy in the median channel. R3 rejects "no-mode"
  fits only.
  - The limitation is known and declared: a single-mode fit cannot separate a
    comparable-amplitude nuisance mode that lies within the effective resolution
    of the record.
  - In that case R4 (the time-block jackknife, which exposes non-stationary
    misfit) is the preregistered guard. Bias should become a wide CI and hence
    UNRESOLVED, not a wrong value.
- **R4.** CI half-widths: `1.96·SE_α ≤ 0.02 s⁻¹` and `1.96·SE_f ≤ 0.02 Hz`.

**No modal value (α̂, f̂) may be used by any rule when the estimator reports
unresolved.**

### 2.5 Verdict (confidence-interval rule; EPS = 0.002 s⁻¹, unchanged from v1)

| verdict | condition |
|---|---|
| **STABLE** | resolved, `α̂ + 1.96·SE_α < −EPS`, and no sentinel flag |
| **UNSTABLE** | resolved, `α̂ − 1.96·SE_α > +EPS`, and no sentinel flag |
| **UNRESOLVED** | otherwise: not resolved (`NOT_RESOLVED`), or the CI overlaps [−EPS, +EPS] (`CI_IN_DECISION_BAND`) |
| **UNEXPECTED_MODE** | the sentinel flags (§3.2); counted as UNRESOLVED by every success rule |

### 2.6 Adaptive duration (scientific IEEE-39 cases, E0 synthetic and E0 TDS; not unit tests)

1. Evaluate at T = 30 s, with window [2.2 s, T].
2. If the verdict is UNRESOLVED **and** R1, R2 and R3 hold and the sentinel did
   not flag, extend **once** to 60 s.
3. If the same condition still holds, a final extension goes to 90 s.

There are no other extensions. A TDS trace that ends earlier (DIVERGED)
uses its available record for the extension step.

## 3. Targeted mode tracker and broad-band sentinel

### 3.1 Targeted tracker (scientific cases)

- **Search neighbourhood.** `f ∈ [f_pred − 0.15, f_pred + 0.15] ∩ [0.2, 1.2]` Hz,
  where f_pred is the committed phasor critical frequency
  (`band_hz`, `results/EMT_PRED/phasor_predictions.csv`).
- **Inputs.** The phasor α (`band_re`) is **never** an input.
- **Channels.** The v1 §5 channel set (`_emt.estimator_channels`): relative
  machine speeds, GFL frequency deviations, bus-angle differences and |V| at the
  ten generator buses.
- **Outputs.** α̂, f̂, their CIs, the normalized residuals (median and global),
  the resolved flag with reasons, the verdict and the duration used.

### 3.2 Broad-band sentinel

- **Input.** The target-removed residual channels R (normalized as in §2.1),
  with no renormalization.
- **Fit.** A VP fit over f ∈ [0.2, 1.2] Hz and α ∈ [−0.8, +0.8], with the §2.3
  CI.
- **Flag UNEXPECTED_MODE** iff all of:
  - the fit is converged and interior;
  - `α_s − 1.96·SE_α,s > +EPS`;
  - it explains ≥ 10 % of the residual energy;
  - `|f_s − f̂| > 0.05` Hz.
- **Limits.**
  - The sentinel can only downgrade a verdict. It never turns UNRESOLVED into
    STABLE or UNSTABLE.
  - A flag is reported for investigation. No further investigation is part of
    V3.

## 4. Gate E0 — instrument qualification

### 4.1 Blind synthetic set (seed 20260912; generated only after the implementation commit)

**Generator:** `experiments/paremt_emt/v3/EMTV3_E0.py`.

**Common properties.**
- Sampling 1 kHz, record 0–90 s.
- The response starts at t_on = 1.2 s:
  `y_c = Σ_m A_cm e^{α_m(t−t_on)} cos(2πf_m(t−t_on) + φ_cm)` for t ≥ t_on,
  and 0 before.
- Amplitudes |A| ~ U(0.2, 1) × random sign; phases ~ U(0, 2π).
- Numerical noise is white, with σ = 1e-4 × the channel's target amplitude,
  unless stated otherwise.
- The neighbourhood centre is `f_pred = f_true + δ`, δ ~ U(−0.03, +0.03) Hz.
  This imitates an imperfect phasor prediction.

**Evaluation.** The targeted tracker plus the sentinel, with the adaptive
duration (the record is truncated at 30, 60 and 90 s).

| stratum | n | content |
|---|---|---|
| S1 single mode | 40 | α ~ U(−0.30, +0.20), f ~ U(0.25, 1.15) Hz, n_ch ~ U{10, …, 40} |
| S2 target + secondary damped mode | 25 | S1 target; secondary mode with f₂ ~ U(0.2, 1.2), \|f₂ − f₁\| ≥ 0.20 Hz, α₂ = α₁ − U(0.15, 0.60), amplitude ratio at t_on ~ U(0.3, 1.5) |
| S3 near boundary | 20 | α ~ U(−0.010, +0.010), f ~ U(0.25, 1.15), n_ch ~ U{10, …, 40} |
| S4 heterogeneous channels | 20 | S1 target; per-channel affine trend c_c ~ U(−1, 1), d_c ~ U(−0.05, 0.05) s⁻¹ (× target amplitude); noise σ = 1e-3; 20 % of channels low-observability (target amplitude × 1e-3) |
| S5 two close modes | 20 | target (α₁, f₁) as S1; second mode f₂ = f₁ ± U(0.04, 0.12) Hz, α₂ = α₁ − U(0.05, 0.30), comparable amplitude (ratio U(0.5, 1.5)) |
| S6 single-channel SG-like ringdown | 15 | 1 channel; f ~ U(1.0, 1.6), α ~ U(−0.5, −0.05); plus a slow real mode B e^{−β(t−t_on)}, β ~ U(0.02, 0.2), B/A ~ U(0.2, 1); window [1.5, 10] s; f search [0.5, 2.0]; no adaptive extension; no sentinel |

Total: 140 cases. For S1–S5 the target is defined as the mode with the larger α.

**Acceptance.** All of A1–A6 (including A3b) must hold.

| id | criterion |
|---|---|
| A1 | no exception or crash in any case |
| A2 | **zero wrong-sign verdicts** (STABLE with α_true > 0, or UNSTABLE with α_true < 0) among all cases with \|α_true\| ≥ 0.005 s⁻¹ |
| A3 | in S1 and S4 (single target), among cases with \|α_true\| ≥ 0.01: ≥ 95 % reach the **correct conclusive** verdict within the adaptive rule |
| A3b | in S2 (target plus nuisance), among cases with \|α_true\| ≥ 0.01: ≥ 80 % reach the correct conclusive verdict. The lower rate is declared because separability is limited by (effective record length) × (frequency separation). Wrong signs remain excluded by A2 and inaccuracy by A4. |
| A4 | over every case whose final target estimate is resolved (all strata): p95 \|α̂ − α_true\| ≤ 0.005 s⁻¹ and p95 \|f̂ − f_true\| ≤ 0.005 Hz |
| A5 | in S6: ≥ 95 % resolved (their accuracy is included in A4) |
| A6 | **DV2-2 regression** (the four historical cases; not blind): no crash; the three cases with \|α\| ≥ 0.01 are conclusive and correct with \|Δα\|, \|Δf\| ≤ 0.005; the +0.0036 case has no wrong sign |

**Reported, not gated:**
- the CI coverage of α_true;
- the S3 and S5 resolution fractions;
- the sentinel flag count.

Near |α| < 0.005, UNRESOLVED is acceptable, and no sign is forced.

### 4.2 Phasor-TDS instrument holdout (fixed here; generated only after the implementation commit)

**Generator:** `experiments/paremt_emt/v3/EMTV3_tds_traces.py`
(`.venv/tx3-analysis`).
- It imports the frozen `experiments/G2_tds.py` without modification and uses
  its `build` and `simulate`.
- Disturbance: the D2 pulse (+2 % bus-20 load for 0.2 s).
- Horizon: 89 s after the pulse, or until the G2 DIVERGED rule fires.
- Times are shifted by +1.0 s so the pulse sits at [1.0, 1.2) s, as in EMT.
- **Channels are the v1 §5 set, computed from the phasor states:**
  - `ω_sg,b − ω_sg39`;
  - GFL `θ′_pll/ω_B − (ω_39 − 1)`, with θ′ from the DAE right-hand side;
  - `∠V_b − ∠V_39`;
  - `|V_b|` at the ten generator buses.
- Sampling: 5 ms.

The independent linear eigenvalues are the committed `band_re` and `band_hz`.

| id | case (k = 1.425, t = 1.5, h = 1) | eigenvalue (α, f) | expected |
|---|---|---|---|
| T1 | P4 H4 (g = 0.03625) | +0.12701, 0.62228 Hz | UNSTABLE |
| T2 | P4 30+33+35 | −0.20469, 0.91602 | STABLE |
| T3 | P4 33 | −0.16963, 0.64977 | STABLE |
| T4 | P4 37 | −0.14800, 0.64049 | STABLE |
| T5 | H4 g = 0.18 | +0.01379, 0.70262 | UNSTABLE |
| T6 | H4 g = 0.25 | −0.01739, 0.71054 | STABLE |
| T7 | H4 g = 0.205 (near boundary) | +0.00125, 0.70610 | not STABLE |
| T8 | P4 H4 + damped condenser 2.0 % (R_0010000) | +0.02636, 0.62090 | UNSTABLE |
| T9 | P4 H4 + damped condenser 3.0 % | −0.02948, 0.62128 (in band; the global α −0.0017 is a real, out-of-band mode) | STABLE |

No frozen governed phasor-TDS traces exist, so no governed case is included.

**Acceptance.** B1–B5 must all hold.

| id | criterion |
|---|---|
| B1 | no crash |
| B2 | every case with \|α_eig\| ≥ 0.01 reaches the correct conclusive verdict within the adaptive rule |
| B3 | T7 is not STABLE (UNRESOLVED or UNSTABLE are acceptable) |
| B4 | every resolved case: \|α̂ − α_eig\| ≤ 0.005 s⁻¹ and \|f̂ − f_eig\| ≤ 0.005 Hz |
| B5 | no UNEXPECTED_MODE flag in any case |

This is qualification of a measurement tool. It is **not** new evidence for the
scientific claims.

**Gate E0 PASS iff A1–A6 (with A3b) and B1–B5 all hold.**

## 5. Device-harness holdouts (ParaEMT unit tests)

**New blind SG holdout B-SG38.** Machine data of bus 38 (Sn 1684.1 MVA,
x′d = x′q = 0.57), with S at its canonical base-PF dispatch
**7.64783381 + j1.23275808** (system pu). P4 settings, the V2 unit-test system,
Pm ×1.02 on [1.0, 1.2) s, 10 s.

**Roles.**

| case | role |
|---|---|
| R-SG30 | regression |
| R-SG36 | regression (seen in V2) |
| B-SG38 | **blind** |
| R-GFL30 | regression |
| B1-GFL35, B2-GFL37 | **blind**, frozen in V2 and never run |

**References.** Canonical classes, Radau; the V2 references for the existing
cases and new references for B-SG38, both quasi-static and dynamic-line.
- Script: `experiments/paremt_emt/v3/EMTV3_refs.py`, which imports the V2
  reference functions.

## 6. Device gates (only after E0 PASS)

**Unit-test estimator settings** (identical for EMT and reference):
- VP on the single channel ω − 1;
- window [1.5, 10] s;
- f ∈ [0.5, 2.0] Hz, α ∈ [−0.8, +0.8];
- time-block jackknife only;
- R1–R4; no adaptive extension; no sentinel.

**Ringdown rule.** EMT and reference must both be **resolved**, with
|Δα| ≤ 0.005 s⁻¹ and |Δf| ≤ 0.005 Hz. If the trajectories pass but the qualified
estimator is unresolved on a clean ringdown (EMT or reference), the gate
**FAILS** as a measurement/device-harness failure, and no α/f are compared.

- **G3a — SG transcription.** Algebraic network, v1 SG kernel. For R-SG30,
  R-SG36 and B-SG38:
  - finite;
  - equilibrium |ΔP|, |ΔQ| ≤ 1e-4 at 0.9 s;
  - every one of the 7 states within 0.02 × excursion against the quasi-static
    reference over [1, 10] s (unchanged 2 % rule, raw 1-ms samples);
  - the ringdown rule.
- **G3b — SG electromagnetic embedding.** The EMT R–L line, against the
  **dynamic-line** reference. Same cases and the same rules (raw 2 %, equilibrium,
  ringdown).
  - **Reported separately:** EMT against the dynamic-line reference; the
    dynamic-line against the quasi-static reference (a measured model-domain
    correction, not an SG transcription error); EMT against the quasi-static
    reference; 25 µs.
- **G4a — GFL transcription.** Exactly V2 G4a: the v1 11-state harness on the
  algebraic network, with the v1 EMT03 rules. Cases R-GFL30, B1-GFL35 and
  B2-GFL37 (27 runs).
- **G4b — GFL physical EMT embedding.** The V2 interface **unchanged**: voltage
  source behind Rf–Lf, measured i_d/i_q, 9 integrated controller states, linear
  extrapolation, no numerical damping. The V2 G4b checks 1–5 are unchanged:
  - no-event stationarity;
  - the Nyquist indicator;
  - 50 vs 25 µs;
  - P/Q;
  - 15-Hz LP trajectories, PLL lock and g = 0 q-tracking.

  **Check 6 is replaced.** The VP estimator runs on the 11 LP'd quantities
  (multi-channel) over [1.05, 4.8] s, with f ∈ [0.2, 4.0] Hz and
  α ∈ [−0.8, +0.8], the same jackknife and R1–R4, applied identically to EMT and
  reference.
  - If the reference is resolved, EMT must be resolved, with
    `|Δα| ≤ 0.005 + 0.02|α_ref|` and `|Δf| ≤ 0.005 + 0.02 f_ref`.
  - If the reference is unresolved, the check is N/A.
  - The use of f above 1.2 Hz lies outside the E0 qualification range. It is used
    only for this identical-processing EMT-vs-reference comparison.

**Stop rules:** §1.4.

## 7. Scientific blocks (only if E0, G3a, G3b, G4a and G4b pass)

The IEEE-39 model is the v1 realization, **except** that every GFL uses the V2
interface (with its own rating). The frozen v1 elements are all unchanged:
- portfolios, policies and line holdouts;
- the frozen predictions `results/EMT_PRED/`;
- the pulse (+2 % bus-20 P, [1.0, 1.2) s);
- Δt 50 µs, 1-ms storage, τ_m = 1 ms;
- the success thresholds.

**Instrument.** Only the broken v1 estimator is replaced: the verdict and α_EMT
come from §2–§3 (targeted tracker, sentinel, adaptive 30 → 60 → 90 s).

**Using α values.**
- A v1 rule that uses "α_EMT even where UNRESOLVED" uses α̂ when the fit is
  **resolved** (R1–R4) but the CI overlaps the decision band.
- If the fit is **not resolved**, no value exists, and the dependent result is
  UNRESOLVED.

The D5 diagnostic drift is unused.

1. **EMT04.**
   - **Cases:** base, 30+33+35 and H4 at P4, and H4 at G_S.
   - **Time steps:** 100, 50 and 25 µs.
   - **G5 acceptance (50 vs 25 µs):** the verdict is invariant, |Δf| ≤ 0.01 Hz
     and |Δα| ≤ 0.01 s⁻¹ where both are resolved. Failure: STOP and amend (v1 §6).
   - **Holdouts** (same acceptance against the primary; H4 at P4 and the base):
     - τ_m 0.5 ms. A failure is the v1 G2 gate: relabel the load realization as a
       model variant.
     - The EMT-native-stator machine variant (descriptive): E′ behind an EMT R–L
       branch R = ra/w, L = x′/(w ω0), realized with the same exact
       ideal-source companion as the V2 GFL filter.
   - **No-event stationarity** of the four cases, 5 s, reported.
2. **EMT05.**
   - **Runs:** the 16 subsets of H4 at P4.
   - **Primary:** all 15 proper subsets STABLE and H4 UNSTABLE.
   - **Secondary:** MAE(α̂ − band_re) ≤ 0.03 and max|f̂ − band_hz| ≤ 0.05 over
     resolved fits.
   - **Amplitude checks:** 1 % and 4 % on 30+33+35 and on H4; |Δα| to 2 % must be
     ≤ 0.01 (resolved fits).
3. **EMT06.**
   - **Effects:** the 8 contextual effects from the EMT05 α̂.
   - **PASS** iff all 8 signs equal the phasor signs. Every α involved must be
     resolved; otherwise the result is UNRESOLVED.
   - The 32 lattice marginal-effect sign agreement is reported.
4. **EMT07/08.**
   - **Points:** g ∈ {0.18, 0.20, 0.205, 0.21, 0.225, 0.25} at H4, k = 1.425,
     plus P4.
   - **Bracket** and at most 4 bisection runs (width ≤ 0.0025), as in v1 §8.
   - **g\*_EMT** is the linear-interpolation zero of α̂ (resolved fits) across
     the final bracket.
   - **Primary:** a resolved UNSTABLE point below g\* and a resolved STABLE point
     above it.
   - **Target:** |g\*_EMT − 0.20768| ≤ 0.02.
5. **EMT10.**
   - **Points:** P4 → g = 0.08806 → g = 0.20768 → g = 0.25.
   - **PASS:** α̂ (resolved) decreases monotonically, g = 0.25 is STABLE, and
     |g\*_EMT − 0.20768| ≤ 0.02.

**Then STOP AND REPORT.**

## 8. Manifests and determinism

- **Manifests.** Every run records the git HEAD, the ParaEMT SHA d79d735a, the
  working-copy diff SHA, the environment hash, the case, parameter and
  operating-point hashes, Δt, duration, disturbance and seed.
- **Determinism.**
  - E0 is rerun once; its summary must be identical.
  - The headline EMT05 and EMT08 grid are rerun once, with byte-identical summary
    CSVs (wall-clock fields excluded).

## 9. What V3 does not claim

- EMT does not prove theorems.
- E0 is instrument qualification, not scientific evidence.
- The phasor-TDS holdout reproduces the phasor model's own dynamics.
- There is no claim of robustness of {30, 33, 35, 37}, no universal four-unit
  claim, and no line-ranking transfer across converter models.
