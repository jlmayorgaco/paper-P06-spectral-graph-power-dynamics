# ParaEMT EMT campaign — preregistration V3 amendments

**Base.** `docs/20260912_PAREMT_EMT_PREREG_V3.md`, commit 6d339006.

Each amendment is committed **before** the runs it affects. No blind E0 data,
no TDS holdout trace, and no V3 ParaEMT run existed when V3-AM1 was committed.

## V3-AM1 — envelope-equalizing weight in the VP fit (estimator specification, §2.2)

**Found during development.** The development runs used only what prereg V3
§1.2 allows: seed 1111 (30 cases), the four non-blind DV2-2 regression cases,
and the non-holdout development TDS case P4 BASE. The first development run of
the implementation as specified in §2 gave:

| item | as specified (first development run) | criterion |
|---|---|---|
| DV2-2 case α = −0.144 (0.95-Hz nuisance mode, α = −0.6) | α̂ = −0.14948, error **0.0055**; CI [−0.1596, −0.1393] covers the truth; verdict STABLE (correct) | A6: \|Δα\| ≤ 0.005 → **fail** |
| development S6 (single-channel SG-like ringdown, 10 s) | 4/5 resolved; the α = −0.428 case failed R4 (CI width) | A5: ≥ 95 % → would fail |
| development S1, S3 and S4 | all resolved and correct; errors ≤ 2e-3 | — |
| development S2 and S5 | resolved; max error 4.4e-3 (S2); one S5 case unresolved | — |

**Cause.** An unweighted least-squares fit of a single damped mode is
dominated by the earliest samples of the window, where a faster-decaying
nuisance mode is largest. This is the declared limitation of prereg V3 §2.4
(R3).

The same early-sample dominance makes the time-block jackknife of a strongly
damped single-channel ringdown depend mostly on its first block. That widens
the CI.

**Amendment.** The only change is to the estimator's fitting procedure. The
model, the parameters optimized (α and f only), the bounds, the grids, the
jackknife definition, R1–R4, the verdict rule, EPS, the sentinel, the adaptive
rule and **every E0 acceptance criterion** are unchanged.

1. **First pass.** An unweighted VP fit, exactly as §2.2, giving α⁽¹⁾.
2. **Weight.** `W(τ) = exp(−α_w τ)`, with `α_w = clip(α⁽¹⁾, −3/L, +3/L)`, where
   L is the window length.
   - The weight's dynamic range over the window is at most e³ ≈ 20.
   - It makes the target envelope approximately flat.
   - It therefore down-weights the early, nuisance-rich samples of decaying
     modes, and the late samples of growing modes, where nonlinearity would
     enter first.
3. **Final pass.** A weighted VP fit: minimize `‖W ⊙ (Y − model)‖_F` with the
   same grid and bounded refinement. W is fixed, not optimized.
4. **Jackknife.** Channel and time-block refits use the same fixed W.
5. **R3.** It uses the weighted residual relative to the weighted affine-removed
   norm, per channel.
6. **Sentinel.** It stays unweighted. Its input is the unweighted residual of the
   final target fit.

**Development result with V3-AM1** (same development data):
- DV2-2: all four cases pass A6. The α = −0.144 case has error 3.6e-4.
- Development set: 30/30 resolved; no wrong sign; p95 |Δα| and p95 |Δf| within
  criteria; S6 5/5.
- Development TDS P4 BASE: α̂ = −0.14469 against −0.14413 (Δ 5.6e-4), and f̂
  against f_eig Δ 1.8e-4 Hz.

**Disclosure.** V3-AM1 was chosen after seeing development results. The blind
qualification set (seed 20260912) and the TDS holdout T1–T9 had not been
generated. Development data are not qualification evidence.
