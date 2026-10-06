# ParaEMT EMT-validation campaign — preregistration v1

Date: 2026-09-11.

- **Branch.** `research/paremt-emt-validation`, created from the finalized
  scientific HEAD a70dac94. Not pushed.
- **Timing.** This document, the frozen phasor predictions
  (`results/EMT_PRED/`) and the canonical operating points are committed
  **before any TX4 EMT result**. Nothing below may change after a TX4 EMT
  result has been seen. Any deviation is a numbered amendment, committed
  before the affected rerun and disclosed in the report.

## 0. Purpose and evidence taxonomy

**Question.** Do the principal observable predictions of the frozen phasor
theory survive in an independent electromagnetic-transient (EMT) realization of
the SAME IEEE-39 network and the SAME device/control equations?

EMT provides model-fidelity corroboration of observable consequences. It does
**not** prove theorems.

Every claim is assigned exactly one of:
- THEOREM (no EMT validation required);
- EMT CORROBORATED;
- EMT QUANTITATIVELY REPRODUCED;
- EMT REFUTED;
- EMT UNRESOLVED;
- NOT APPLICABLE TO EMT.

The following are **not** EMT targets: T01–T06 and C01–C06. Connected
cumulants remain secondary and get no EMT figure.

**Wording rule.** Write "EMT corroborates the observable consequence of
theorem X", never "EMT proves X".

## 1. Simulator and environment

- **Official upstream.** `https://github.com/NatLabRockies/ParaEMT_public.git`,
  commit `d79d735a4a587d56c5b88187d1a499195b6b2b84` (2026-01-15), cloned
  2026-09-11 to `external/ParaEMT_upstream`. It is never modified.
- **Working copy.** `external/ParaEMT_tx4`, a clone of the same commit.
  - Its only diff is the documented TX4 network patch
    (`experiments/paremt_emt/tx4_patch_paremt.py`):
    - exact from-side transformer taps in `numba_InitNet`,
      `numba_updateIhis` and `Re_Init`;
    - a switch for ParaEMT's numerical damping elements;
    - stock behaviour stays the default.
  - The custom TX4 models and the time loop are added as new files (§3).
  - The diff SHA is recorded per run.
- **Environment.** `.venv/xtool-paremt`, CPython 3.11.10, with the exact
  upstream requirements (setuptools 72.1.0, xlrd 1.2.0, matplotlib 3.9.2,
  numba 0.60.0, numpy 2.0.1, pandas 2.2.2, scipy 1.14.0) and nothing
  else.
- **Network solver.** Serial; the matrix is factorized once per case. The
  nxmetis/BBD parallel path is not used or repaired.
- **Untouched.** `vendor/paraemt` (the earlier TX3 copy) and
  `.venv/tx3-paraemt`.
- **EMT00 (done before this document; tool-level only).** It is
  deterministic, serializes, snapshots and resumes. The stock IEEE-39
  no-event run is not stationary because upstream never applies
  `pfd.xfmr_k`; with that single tap applied, the start-up jump falls from
  0.62 to 2.2e-4 pu. Tool gate: PASS (`results/EMT00/EMT00_gate.json`).

## 2. Network realization (EMT01)

- **Data source.** The frozen TX4 data, `configs/ias2026/ieee39_network.json`
  (derived from `data/raw/ieee39_full.xlsx`, sha256 9c2048dc…). The frozen
  benchmark is never changed to match ParaEMT. Instead, ParaEMT is given a
  TX4-specific case.
- **Elements.** Built with ParaEMT's own `numba_InitNet` companion models
  (trapezoidal rule):
  - **Lines.** Series R–L and C = b/(2ω0) at each end, balanced three-phase
    with no mutual coupling. Every TX4 line has X > 0.
  - **Transformers.** Rows with `tap != 1` or `trans == 1`, all with
    b = g = 0. Series R–L with an ideal tap t on bus1, the canonical
    convention: y/t² (bus1–bus1), y (bus2–bus2), −y/t (off-diagonals). This
    uses the working-copy tap patch.
  - **Shunts** (buses 4 and 5). Capacitors with `shnt_gb = 100 (g + jb)`,
    g = 0.
- **Numerical damping off** (`net_damping = 0`). ParaEMT's damping resistors
  (Rp = 20/3 · 2L/Δt across R–L branches; Rs = 0.15 Δt/2C in series with
  capacitors) change the 60-Hz admittance by about 1.4e-3. Removing them is a
  documented convention translation; the discretization is then the pure
  trapezoidal rule.
- **Base and units.** 100 MVA, 60 Hz, `bus_basekV = 1` (all per unit);
  instantaneous values in ParaEMT's amplitude-per-unit convention;
  phase-block node ordering.
- **Ybus acceptance.**
  - The 60-Hz positive-sequence nodal admittance assembled from the
    **continuous-time element values actually stored** in the EMT network
    (R, L, C, tap of every branch) must match the canonical TX4 Ybus with
    `max|ΔY| / max|Y| ≤ 1e-8`.
  - The trapezoidal frequency warping at 60 Hz (a relative reactance error of
    about (ω0Δt)²/12, i.e. 3e-5 at 50 µs) is a numerical-method property. It
    is reported separately and is not part of this gate.
- **Equilibrium acceptance** (base portfolio at P4, 2 s no-event run at
  50 µs, positive-sequence phasors averaged over the last 0.5 s):
  - `max|ΔV| ≤ 1e-4` pu and `max|Δθ| ≤ 1e-3` rad against the canonical power
    flow (angles relative to bus 39);
  - P/Q residuals per device reported.
- **Official ParaEMT IEEE-39 vs frozen TX4.** An automated diff of base MVA,
  frequency, bus numbering, base kV, bus types, line and transformer
  endpoints, R/X/B, taps and orientation, shunts, loads, generator P/Q and
  ratings, and |V| and θ. Outputs: `results/EMT01/network_diff.csv`,
  `operating_point_diff.csv`, `docs/20260911_PAREMT_NETWORK_EQUIVALENCE.md`.

## 3. Device realization (frozen equations)

All devices use the same frame convention as the phasor model: a synchronous
frame rotating at ω0 = 2π·60.

**Voltage phasor extraction.** The network-frame voltage phasor of a bus is
extracted from the abc solution by the amplitude-invariant space vector

`V(t) = (2/3)(v_a + a v_b + a² v_c) e^{−jω0 t}`, with `a = e^{j2π/3}`.

This is exact for positive sequence; negative- and zero-sequence content
appears as 2ω0 / ω0 ripple.

**Current injection.** A device current phasor I(t) is injected as

`i_k(t) = Re(I(t) e^{j(ω0 t − 2πk/3)})`, for k = 0, 1, 2.

**Synchronous-frame algebraic Norton element.** A phasor relation
`I = Y (E − V)` with constant complex Y is realized **exactly** as a constant
real 3×3 abc conductance `G_km = (2/3) Re(a^{m−k} Y)`, stamped into the
network matrix, plus the source current `Re(Y E e^{j(ω0 t − 2πk/3)})`.
- Its symmetric part is positive semidefinite for Re Y ≥ 0, so the element is
  passive.
- It reproduces the phasor relation for positive sequence at every frequency,
  with no added electromagnetic state.

### 3.1 Synchronous machine (two-axis, own base, weight w = Sn/100)

- **Equations.** Identical to `SynchronousMachine` (spec
  `docs/20260911_GFL_REPRODUCTION_SPEC.md` §4.1): 7 states
  (δ, ω, e′q, e′d, efd, pss_w, pss_l), first-order AVR, washout-lag power
  PSS, D = 0 and constant Pm (no governor) unless §9 EMT14 applies.
- **Stator.** x′d = x′q for all ten machines (checked), so the stator
  `I = (E′ − V)/(ra + j x′)` with `E′ = (e′q − j e′d) e^{jδ}` is exactly an
  algebraic Norton element with `Y = w/(ra + j x′)` on the system base. The
  EMF source uses the predicted states of the current step.
- **Air-gap power** `Pe = Re(E′ conj(I_dev))`, the same as the phasor
  `vd id + vq iq + ra(id² + iq²)`. `id, iq` come from the phasor I at δ, and
  `|V|` from V(t).
- **Secondary model variant (EMT04 only, descriptive).** An "EMT-native
  stator": E′ behind an R–L branch (R = ra/w, L = x′/(w ω0)) with ParaEMT's
  trapezoidal companion. It adds stator electromagnetic transients and is
  labelled as a variant.

### 3.2 Grid-following converter (11 states)

- **Equations.** Identical to `GridFollowingConverter` (GFL spec §4.2): θ,
  x_pll, p_f, q_f, x_p, x_q, i_d, i_q, x_id, x_iq, x_v; Q/V gain g,
  leak w = 0.05 rad/s, ω_B/xf filter dynamics with nominal cross-coupling.
- **Inputs.** `v_d + j v_q = V(t) e^{−jθ}` and `|V|` from V(t).
- **Output.** The average-value current injection
  `I = w (i_d + j i_q) e^{jθ}`: a state-driven current source, not a Norton
  element.
- **Not included.** Current limits, DC link, PWM delay, saturation,
  ride-through or protection.

### 3.3 Constant-power load (outcome A, declared)

- **Norton part.** A constant-impedance Norton element with
  `Y0 = conj(S0)/|V0|²`, where S0 is the scheduled load and V0 the
  canonical power-flow voltage.
- **Compensation.** A compensating injection

  `I_comp = −(conj(S(t)/V_f) − Y0 V_f)`,

  with V_f a first-order filter of V(t), `τ_m = 1 ms`
  (`dV_f/dt = (V − V_f)/τ_m`, V_f(0) = V0). At equilibrium I_comp = 0.
- **Low-frequency equivalence.** The effective incremental load admittance
  is `Y0 + (Y_cp − Y0)/(1 + sτ_m)`: exact constant power as sτ_m → 0. At
  0.7 Hz the relative deviation of the (Y_cp − Y0) part is ωτ_m = 4.4e-3.
- **Why the filter.** It keeps the explicit coupling numerically stable at
  load buses without capacitance. It is declared, not tuned.
- **Holdout.** τ_m = 0.5 ms (EMT04).
- **Disturbance.** The load pulse acts on S(t) (§5).

### 3.4 Condenser (EMT13) and governor (EMT14)

**Condenser.** The frozen G1 R_0010000 damped condenser at each replaced bus:
- a two-axis machine at rating fraction p of the retired machine;
- inertia scale 0.04, D = 2 (condenser base);
- flux and AVR blends 0, PSS 0, q_share 0 (zero output at equilibrium);
- the GFL carries all P and Q.

It uses the same machine realization as §3.1.

**Governor.** TGOV1N from `configs/ias2026/ieee39_governed_documented_v1.json`
on every machine with Pm ≠ 0 (as `ibr_cycles.models.governed`), 2 states,
valve limits inactive at equilibrium.

## 4. Numerical method and initialization

- **Time loop.** A new jitted loop (`external/ParaEMT_tx4/tx4_emt.py`) using
  ParaEMT's network matrices and `numba_updateIhis`. Per step n:
  1. **predict** the device states:
     `x_p = x_{n−1} + Δt f(x_{n−1}, V_{n−1})`;
  2. assemble the source currents (Norton EMFs, GFL currents, load
     compensation) from x_p;
  3. solve the network `G V_n = I_src + I_his` (dense inverse of the constant
     G);
  4. **correct** the device states:
     `x_n = x_{n−1} + (Δt/2)[f(x_{n−1}, V_{n−1}) + f(x_p, V_n)]`;
  5. update the branch histories (ParaEMT) and the load filter.
- **Time steps.** Δt = 50 µs (primary), 25 µs (holdout), 100 µs (convergence
  table only).
- **Output.** Stored every 1 ms. The raw series go to compressed NPZ under
  `external/paremt_runs/raw/` (hash in the run manifest). Compact summaries
  and the estimator inputs of the headline runs are committed.
- **Initialization.**
  - The canonical power flow (`results/EMT_PRED/tx4_operating_points.json`)
    gives the network phasors.
  - The device states use the phasor `initialize` formulas; the load gets
    V_f = V0.
  - The branch histories are set from the phasor steady state (ParaEMT).
  - One equilibrium serves every portfolio: matched dispatch. The 12 line
    cases use their own canonical power flows.

## 5. Common disturbance and estimator

### Disturbance (all portfolio comparisons)

- **Pulse.** The frozen G2 D2 pulse: +2 % of the bus-20 active load
  (6.8 pu → +0.136 pu) from t = 1.000 s to t = 1.200 s.
  - It is applied to S(t) of the constant-power load: P only, rectangular.
  - It is identical for every case.
  - Bus 20 is a non-candidate load bus, chosen before EMT because the frozen
    phasor TDS uses exactly this disturbance, which keeps EMT directly
    comparable.
- **Amplitude checks.** 1 % and 4 % on EMT05 30+33+35 (stable) and on H4
  (P4).
- **Run length.** 30 s, and 60 s for EMT13. Longer runs are used only by the
  rule in §10.

### Estimator signals (relative, so the neutral common-frequency mode is excluded)

- **Relative machine speeds** `ω_i − ω_39` of every machine present.
- **GFL frequency deviations** `θ'_j/ω_B − (ω_39 − 1)`.
- **Bus-angle differences** `∠V_b − ∠V_39` for the ten generator buses.
- **Voltage magnitudes** |V_b| for the ten generator buses.

### Preprocessing

- The analysis window is [2.2 s, T_end].
- Each channel is detrended (least-squares linear trend) and scaled to unit
  RMS.
- The data are zero-phase FIR anti-aliased and decimated to 50 Hz (Hilbert)
  and to 10 Hz (matrix pencil).

### Estimator A — multi-channel matrix pencil (10 Hz data)

- **Pencil.** A stacked Hankel matrix of all channels, with pencil parameter
  L = ⌊N/3⌋.
- **Model order.** The number of singular values with σ_i/σ_1 ≥ 1e-4,
  between 2 and 30.
- **Modes and energies.** Modes `s_i = ln(z_i)/Δt`. Residues are found by
  least squares per channel. Energy = Σ_channels |residue|² × ∫window |e^{s t}|² dt.
- **Selection.** In-band modes (0.2 ≤ f ≤ 1.2 Hz) whose energy is ≥ 1 % of
  the total in-band energy.
- **Outputs.** `α_A = max Re` over the selected modes; `f_A` is that mode's
  frequency.
- **Resolved iff** at least one mode is selected and the multi-channel
  reconstruction residual is ≤ 0.05 (relative, in window).

### Estimator B — band-limited analytic-signal envelope (50 Hz data)

- **Channel and band.** The channel with the largest in-band Welch power.
  The peak frequency f_pk is taken from the Hann-windowed zero-padded FFT in
  band.
- **Filter.** A 4th-order Butterworth band-pass over
  [f_pk − 0.15, f_pk + 0.15] ∩ [0.2, 1.2] Hz, applied zero-phase.
- **Fit.** The Hilbert analytic signal. `α_B` is the slope of a least-squares
  fit of ln|envelope| on the window trimmed by 3 s at each end; `f_B` is the
  slope of the unwrapped phase divided by 2π.
- **Resolved iff** the fit has R² ≥ 0.90 and spans ≥ 5 cycles.

### Classification (ε = 0.002 s⁻¹, frozen)

- **STABLE** if both estimators are resolved and `α_A, α_B < −ε`.
- **UNSTABLE** if both are resolved and `α_A, α_B > +ε`.
- **UNRESOLVED** otherwise. Near a boundary this is allowed and never forced.
- **Agreement flag** (reported, not required):
  `|α_A − α_B| ≤ 0.02 s⁻¹` and `|f_A − f_B| ≤ 0.03 Hz`.
- **Reported EMT margin and frequency:** `α_EMT := α_A`, `f_EMT := f_A`.
- **Phasor comparator.** The in-band rightmost transverse eigenvalue
  (`band_re`, `band_hz` in `results/EMT_PRED/phasor_predictions.csv`). Where
  the global α⊥ lies outside the band (condenser 3 % and 5 %, some EMT16
  subsets), the EMT verdict concerns the in-band family only. This is stated
  in the result.

## 6. Numerical convergence (EMT04)

- **Cases.**
  - base at P4;
  - 30+33+35 at P4 (stable proper subset);
  - H4 at P4 (unstable);
  - H4 at G_S = (0.25, 1.425, 1.5, 1) (stable near-policy point).

  Each at Δt = 100, 50 and 25 µs.
- **Acceptance (G5).** Between 50 and 25 µs:
  - the verdict is invariant;
  - `|Δf_EMT| ≤ 0.01 Hz`;
  - `|Δα_EMT| ≤ 0.01 s⁻¹` where resolved.
- **Holdouts (descriptive, same acceptance).** τ_m 0.5 ms vs 1 ms, and the
  EMT-native-stator machine variant, both on H4 at P4 and on the base.
- **On failure.** If 50 µs fails, stop and amend before continuing.

## 7. Unit tests (G3, G4)

- **Test system.** The device at bus T, a line of R = 0, X = 0.05 pu (system
  base, EMT R–L), and an infinite bus. The infinite bus is realized as a
  Thevenin source 1.0∠0 behind z = j1e-4 pu, identical on both sides.
- **Operating point.** Device output S = 4.0 + j1.0 pu, the bus-30 device
  (Sn 1040 MVA), P4 settings.
- **Phasor reference.** Integrated from the canonical device classes (Radau,
  rtol 1e-10), sampled at 1 ms.

**EMT02, machine.**
- **Disturbance.** Pm +2 % for 0.2 s at t = 1 s, 10 s run.
- **Compared.** Equilibrium P and Q; δ, ω, e′q, e′d, efd, pss_w, pss_l.
- **Pass.**
  - equilibrium `|ΔP|, |ΔQ| ≤ 1e-4`;
  - every trajectory: `max|x_EMT − x_ph| ≤ 0.02 · max|x_ph − x_ph(0)|` over
    [1 s, 10 s];
  - ringdown `|Δα| ≤ 0.005 s⁻¹`, `|Δf| ≤ 0.005 Hz` (matrix pencil on ω).

**EMT03, GFL.** Cases at g ∈ {0, 0.03625, 0.25}:
- (a) p_ref step of +1 % at 1 s;
- (b) infinite-bus |V| step of −1 % at 1 s (Q/V);
- (c) infinite-bus phase step of +0.02 rad at 1 s (PLL lock).

The run is 5 s.
- **Compared.** All 11 states, P and Q.
- **Pass.**
  - the trajectory rule above for every state;
  - equilibrium `|ΔP|, |ΔQ| ≤ 1e-4`;
  - PLL steady-state `|θ − ∠V_T| ≤ 1e-3` rad;
  - g = 0: q_f returns to q_ref within 1e-3.
- **Gate.** No IEEE-39 result is called equation-equivalent until EMT02 and
  EMT03 pass. The equation map goes to
  `docs/20260911_PAREMT_GFL_EQUATION_MAP.md`.

## 8. Core experiments

In every experiment below: same network, matched P/Q, same ratings, same
disturbance, same estimator, 30 s. The phasor predictions are frozen in
`results/EMT_PRED/phasor_predictions.csv`.

**EMT05 — 16 subsets of H4 at P4**
(`results/EMT05/p4_16_subsets.csv`).
- **Primary PASS.** All 15 proper subsets STABLE and H4 UNSTABLE, so
  `H_EMT = {{30,33,35,37}}` and `κ_EMT = 4`.
- **Secondary targets** (resolved cases; do not override the primary):
  - MAE(α_EMT − band_re) ≤ 0.03 s⁻¹;
  - `max|f_EMT − band_hz|` ≤ 0.05 Hz.
- **Amplitude checks.** For each of the two cases, |α_EMT(1 %) − α_EMT(2 %)|
  and |α_EMT(4 %) − α_EMT(2 %)| must both be ≤ 0.01 s⁻¹.

**EMT06 — contextual effects**, computed from the EMT05 α_EMT.
- `Δ_i α(∅) = α({i}) − α(∅)` and `Δ_i α(H4∖{i}) = α(H4) − α(H4∖{i})`.
- **Hypothesis.** Each replacement is stabilizing alone (< 0) and
  destabilizing last (> 0).
- **PASS** iff the EMT signs equal the phasor signs for all 8 effects.
- All 32 Boolean-lattice marginal effects: the sign-agreement fraction is
  reported.

**EMT07–08 — clean g-only crossing** (H4; k = 1.425, t = 1.5, h = 1)
(`results/EMT08/g_boundary.csv`).
- **Points.** g ∈ {0.18, 0.20, 0.205, 0.21, 0.225, 0.25}, plus P4 from EMT05.
- **Bracket.** The bracket is the consecutive resolved UNSTABLE/STABLE pair.
  Refine by bisection, at most 4 extra runs, stopping when the width is
  ≤ 0.0025.
- **g\*_EMT.** The zero of the linear interpolation of α_EMT across the
  final bracket. α_EMT is used even where the case is UNRESOLVED.
- **Primary PASS.** At least one resolved UNSTABLE point below and one
  resolved STABLE point above g\*_EMT.
- **Target.** `|g*_EMT − 0.20768| ≤ 0.02`.

**EMT09 — P_inf H4.** A multi-coordinate contrast: P4 → P_inf changes both g
and k. PASS: STABLE.

**EMT10 — Newton design.**
- **Points.** P4; the frozen first port-Newton iterate g = 0.08806
  (under-correction); the phasor Newton g\* = 0.20768 (UNRESOLVED allowed);
  and g = 0.25.
- **PASS.** α_EMT decreases monotonically along the sequence, g = 0.25 is
  STABLE, and `|g*_EMT − g*_phasor-Newton| ≤ 0.02`.
- **Not claimed.** Finite-step magnitude.

**EMT11 — k boundary** (g = 0.03625, t = 1.5, h = 1).
- **Points.** k ∈ {1.25, 1.28, 1.30, 1.31, 1.33, 1.35}, bracketed and bisected
  as in EMT08.
- **Target.** `|k*_EMT − 1.3046267| ≤ 0.02`.
- **PASS.** A clear STABLE → UNSTABLE change across k\*_EMT.

**EMT12 — 12 frozen holdout lines × 1.5 at P4 H4.** Each line case is
re-equilibrated on its canonical power flow.
- `Δα_EMT = α_EMT(line) − α_EMT(P4 H4)`.
- **Primary.** Against the exact phasor finite change: sign agreement
  ≥ 10/12 **and** Spearman ≥ 0.80.
- **Secondary.** Against the port derivative and the static proxies
  (ΔgSCR, 1/x, x-weighted betweenness): Spearman and Kendall, plus the
  magnitude error.
- **Not claimed.** Transfer to other converter models.

**EMT13 — damped condenser.** p ∈ {2.0, 2.5, 3.0, 5.0} % at P4 H4, 60 s.
- **Primary.** 2.0 % UNSTABLE and 5.0 % STABLE, in band.
- 2.5 % and 3.0 % are descriptive (2.5 % is expected UNRESOLVED).
- No exact threshold is required.

**EMT14 — documented governors.**
- **Cases.** Governed H4 at P4 (primary: STABLE), and governed H4 at
  g = 0.020, k = 1.425, inside the frozen governed κ = 4 interval (expected
  UNSTABLE, descriptive).
- **Interpretation rule.** The exact witness is model-dependent;
  non-composability and policy dependence are the broader phenomenon.

**EMT15 — near-Hopf finite-disturbance recovery (nonlinear stress, separate).**
- **Points.** H4 at g = 0.22 (phasor −0.0055) and g = 0.30 (−0.0334),
  k = 1.425.
- **Sweep.** The same pulse shape at amplitudes of 2, 5, 10, 20, 40 and 80 %
  of the bus-20 load.
- **Outcome labels** (not equivalent to each other):
  - RETURNS (final 5 s relative-speed spread < 1e-3 pu, and envelope
    decaying);
  - LEAVES_NEIGHBOURHOOD (relative speed > 0.05 pu, or an angle-difference
    change > π/2);
  - FAILS_TO_RECOVER_IN_WINDOW;
  - NUMERICAL_FAILURE.
- **Guard flag, recorded separately** (not applied): GFL terminal |V|
  outside [0.9, 1.1] pu, or PSS output |KS·pss_l| > 0.1.
- **Refinement.** Bisection between the largest RETURNS and the smallest
  non-RETURNS amplitude, 4 steps.
- **Question.** Can a linearly stable near-boundary policy have poor
  finite-disturbance recovery? A negative result is acceptable.

**EMT16 — robustness holdout.**
- **Selection.** Fixed before EMT (`results/EMT_PRED/emt16_selection.json`,
  seed 20260921) from the frozen PCV05 draws of the machine envelopes:
  - A: 6 draws with H = {H4};
  - B: 6 draws with H = {30+33+35};
  - C: 4 draws with H = ∅.
- **Tests.** In A and B: the witness, its immediate proper subsets, and the
  base. In C: H4 and the base.
- **PASS** on each of the following, each in ≥ 80 % of its draws:
  - **A:** the witness UNSTABLE and its subsets STABLE;
  - **B:** {30,33,35} UNSTABLE and its pairs STABLE (the witness is not H4);
  - **C:** H4 STABLE.
- **Rule.** Bus-level frequencies are not promoted.

**EMT17 — Kundur zero-frequency (optional).** Run only if the frozen Kundur
configuration of this repository can be reproduced with the same adapter;
otherwise NOT RUN, keeping the existing evidence.

**EMT18 — high-fidelity model-risk stress.** Not part of this campaign
version. It would need a new versioned model and a new preregistration.

## 9. Software validation

Each run writes a manifest with:
- run_id and git HEAD;
- the ParaEMT upstream SHA and the working-copy diff SHA;
- the environment hash (pip freeze sha256);
- case, model-parameter and operating-point hashes;
- Δt, duration, disturbance and seed (none);
- wall time, and the output path with its sha256.

The deterministic headline experiments (EMT05, and EMT08 at the grid points)
are rerun once. Summary CSVs must be byte-identical, with wall-clock fields
excluded.

## 10. Stopping and failure gates

| gate | condition | action |
|---|---|---|
| G0 | official ParaEMT IEEE-39 cannot run reliably | stop (passed at the tool level, §1) |
| G1 | the EMT01 Ybus or equilibrium tolerance fails | stop comparisons; diagnose; later phases BLOCKED |
| G2 | the constant-power realization materially changes the low-frequency system (the τ_m holdout fails the EMT04 acceptance) | stop the "equation-equivalent" label; relabel as a model variant |
| G3 | EMT02 fails | stop |
| G4 | EMT03 fails | stop |
| G5 | the 50 µs results are not converged against 25 µs | stop, then amend |

- **No retuning** of any physical or estimator parameter to obtain
  agreement.
- **Implementation bugs** found during the campaign may be fixed. Each fix is
  listed as a numbered amendment and every affected run is repeated. Design
  choices in this document are not bugs.
- **Longer-run rule.** If an estimator reports "not resolved" because the
  window contains fewer than 5 cycles or the fit fails, the case may be rerun
  once at 60 s. This is recorded.

## 11. Success gates (decided at the end)

- **A, MODEL.** EMT01, EMT02, EMT03 and EMT04 pass.
- **B, PORTFOLIO.** The EMT05 primary.
- **C, POLICY.** The EMT08 primary.
- **D, DESIGN.** EMT10: the Newton direction restores the EMT case.
- **E, LINE DESIGN.** The EMT12 primary.
- **F, MODEL DEPENDENCE.** The EMT13 and EMT14 primaries.
- **G, NONLINEAR.** EMT15 is explicitly characterized, whether positive or
  negative.

## 12. Deliverables

- **Results.**
  - `results/EMT00`, `results/EMT01` … `results/EMT16`;
  - `results/20260911_EMT_CLAIM_MATRIX.csv` (T01–T10, C01–C11, B01–B10,
    I01–I04, N01–N09);
  - figures F1–F8 (no cumulant figure).
- **Documents.**
  - `docs/20260911_PAREMT_NETWORK_EQUIVALENCE.md`;
  - `docs/20260911_PAREMT_GFL_EQUATION_MAP.md`;
  - `docs/20260911_PAREMT_EMT_FINAL_REPORT.{md,tex,pdf}`.
