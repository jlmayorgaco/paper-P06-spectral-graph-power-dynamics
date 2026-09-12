# Frozen v1 EMT modal estimator — failure audit (V3-0)

Date: 2026-09-12. Branch `research/paremt-emt-validation`, not pushed. This is
documentation only.

**Audited file:** `experiments/paremt_emt/emt_estimator.py`.
- Frozen with prereg v1 (c2947bd8) and committed at ff1c4d8d.
- sha256 `38e4f1eb…d5bd7e6`.
- It is **not modified**, because it belongs to the V1 and V2 provenance.

**Evidence.**
- `results/EMTV2/G3a/G3a_diagnostics.json` (V2, commit 632374dc).
- `results/EMTV3/V3-0/estimator_audit.json`. This comes from
  `experiments/paremt_emt/v3/EMTV3_estimator_audit.py`, which re-executes the
  estimator's own arithmetic step by step on the same inputs.

**Permanent history.** Nothing here changes it:
- **V1:** G3 FAIL, G4 FAIL, EMT04–EMT18 BLOCKED.
- **V2:** SW PASS, G3a FAIL under its preregistered rule, G3b/G4a/G4b NOT RUN,
  EMT04–EMT18 BLOCKED.

## 1. What the estimator does (estimator A, `matrix_pencil`)

1. **Prepare.** Take the window [2.2 s, T]. Detrend each channel linearly, scale
   it to unit RMS, and decimate it with a zero-phase FIR to 10 Hz.
2. **Hankel pencil.** Stack the Hankel matrices, with pencil parameter
   L = ⌊N/3⌋, and take their SVD. The **model order** is the number of singular
   values with σ_i/σ_1 ≥ 1e-4, clipped to [2, 30].
3. **Pencil eigenvalues.** Take z_i and set `s_i = ln z_i / Δt`.
4. **Residues.** Solve `lstsq(V, y)` with the **full-record Vandermonde**
   `V[k, i] = z_i^k`, for k = 0 … N−1, and `rcond=None`.
5. **Energy.** Compute `energy_i = Σ_c |res_ic|² · (e^{2 Re(s_i) T} − 1)/(2 Re s_i)`.
   Select in-band modes (0.2–1.2 Hz) with energy ≥ 1 % of the in-band total.
   Report **α = max Re over the selected modes**.
6. **Resolution.** The result is "resolved" iff at least one mode is selected
   and the relative reconstruction residual is ≤ 0.05.

## 2. Findings

### 2.1 `resolved = False` in all four V2 G3a evaluations

The ringdown signal is ω − 1 on the algebraic network.

| evaluation | resolved | residual | order | reported α (s⁻¹) | reported f (Hz) |
|---|---|---|---|---|---|
| R-SG30 EMT | False | 0.99887 | 10 | −0.118368 | 1.177112 |
| R-SG30 reference | False | 0.99887 | 10 | −0.118367 | 1.177112 |
| B-SG36 EMT | False | 0.99875 | 16 | **+6.3134** | 1.1051 |
| B-SG36 reference | False | 0.99874 | 16 | **+6.2976** | 1.1097 |

### 2.2 The reconstruction residual was about 0.999

About 0.999 means the fitted modal model reconstructs essentially **none** of
the signal. The mechanism is §2.6.

### 2.3 The G3a rule nevertheless compared α and f unconditionally

- **The rule.** It was inherited verbatim from v1 EMT02 ("ringdown |Δα| ≤ 0.005,
  |Δf| ≤ 0.005, frozen `matrix_pencil` on ω − 1"). It has no resolution
  condition.
- **v1.** For R-SG30 it compared two unresolved outputs that happened to equal
  the physical mode (§2.6c).
- **V2, bus 36.** It compared two spurious outputs: Δα = 0.0158 > 0.005, so G3a
  FAILED.
- **Contrast.** The V2 prereg did add "where resolved" logic, but only for the new
  G4b estimator R.

### 2.4 The blind B-SG36 EMT run and its reference both return the same spurious ≈ +6.3 s⁻¹ mode

- **Spurious.** A growing mode at +6.3 s⁻¹ is not physical. Both trajectories
  decay, and a fixed-order damped-sinusoid fit gives α = −0.2439564 (EMT) and
  −0.2439563 (reference) at f = 1.250783 Hz.
- **Same in both.** The spurious mode appears identically in the EMT run and in
  the reference. The two differ only in the fourth significant digit of a
  numerically meaningless quantity.
- **The data agree.** The machine trajectories themselves agree to 2.1e-4 of
  their excursion.

### 2.5 DV2-2: three numerical failures and one unresolved result

The inputs are synthetic: 20 channels, 30 s, the target mode plus a damped
0.95-Hz mode (seed 7).

| true α, f | estimator A | reason (audit) |
|---|---|---|
| +0.127, 0.622 Hz | **LinAlgError** (SVD did not converge in lstsq) | spurious \|z\| = 18.9; z^278 ≈ 10^355 overflows to inf |
| −0.144, 0.637 Hz | unresolved (residual 0.9999998) | \|z\| = 4.0; column norms span > 10^300; lstsq keeps rank 1 of 10 |
| −0.0174, 0.71 Hz | **LinAlgError** | \|z\| = 18.7; z^278 ≈ 10^353 overflows |
| +0.0036, 0.705 Hz | **LinAlgError** | \|z\| = 18.7; overflow |

`classify` therefore returns ERROR in three cases and UNRESOLVED in one. The
EMT05–EMT18 portfolio classification depends on it, so it **could not have
produced a single STABLE/UNSTABLE verdict** even if every device gate had
passed.

### 2.6 Why: an ill-conditioned full-record Vandermonde and exponential energy weighting

**(a) Spurious |z| > 1 modes.** The order rule (σ_i/σ_1 ≥ 1e-4) admits noise
and trend-residual directions as "modes". Every audited evaluation has
spurious pencil eigenvalues outside the unit circle:
- up to |z| = 16.7 (R-SG30);
- up to |z| = 2.17 (B-SG36);
- up to |z| = 18.9 (DV2-2).

Detrending and short records (N = 79 samples for the SG ringdown at 10 Hz)
make this worse.

**(b) Full-record Vandermonde.** Each column is `z_i^k` for k up to N − 1. With
|z| = 2.17 and N = 79 the columns span **10^26** in norm, and with |z| = 16.7
they span **10^95**; beyond about 10^308 they overflow (DV2-2).
- `lstsq(..., rcond=None)` truncates singular values relative to the largest,
  which belongs to the exploding column.
- It therefore keeps a numerical rank of **1–3 out of 10–16** and discards the
  physical columns.
- The reconstruction is then ≈ 0, which gives the residual ≈ 1.

**(c) Exponential energy weighting.** `energy_i ∝ |res_i|² · e^{2 Re(s_i) T}`.
- A spurious in-band mode with α = +6.3 s⁻¹ and a residue of only 1e-22 gets a
  weight of **10^43**. Its energy is 0.024, while the physical mode, after
  truncation, has ≈ 1e-92.
- So it is "selected", and `max Re` reports it (B-SG36).
- **R-SG30.** The only in-band mode was the physical one, whose truncated
  residue was exactly 0. The in-band total was then 0, so the threshold
  `energy ≥ 0.01·0` selected it, and the physical value came out **by accident**
  (with a NaN energy fraction).

**(d) No guard.** There is no bound on Re(s) and no rejection of |z| ≫ 1
before the residue solve. The "resolved" flag was computed but, in G3a, was not
consulted.

## 3. Consequences

1. **V1 and V2 histories stand.** The V2 G3a FAIL was produced by the
   preregistered rule acting on this estimator. It is not reinterpreted.
2. **Machine transcription.** On the evidence of trajectories and fixed-order
   fits it was never in doubt: 2e-4 of the excursion, and equal modes to 1e-7.
   This is a post hoc reading, not a gate result.
3. **Estimator unfit.** The v1 estimator A is unfit for the purpose preregistered
   in v1 §5, and so is anything built on `classify`. Every later EMT modal
   verdict needs a **new, independently qualified instrument** (V3).
4. **Estimator B.** The Hilbert envelope did not crash. It resolved the two
   strongly damped or growing DV2-2 cases correctly, but it was unresolved near
   the boundary. It is not the V3 instrument and is not modified.

## 4. Requirements carried into V3 (from this audit)

- **No full-record Vandermonde and no exp(2αT) energy weighting.**
- **Bounded (α, f)**, with the time axis centered on the window so that basis
  functions stay O(e^{|α|·T/2}).
- **A numerically stable linear sub-problem** (QR/SVD on a small, fixed-width
  design matrix).
- **An explicit resolution rule, consulted by every gate.** No modal value may
  be used when it is UNRESOLVED.
- **Qualification on blind synthetic signals and on frozen phasor-TDS traces
  before any EMT use** (gate E0).
