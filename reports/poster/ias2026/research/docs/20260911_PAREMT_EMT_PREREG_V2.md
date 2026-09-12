# ParaEMT EMT-validation campaign — preregistration V2

Date: 2026-09-11. Branch `research/paremt-emt-validation`, not pushed. The author
approved it after the v1 report.

It is committed **before any V2 computation**, together with
`docs/20260911_PAREMT_GFL_VOLTAGE_INTERFACE_DERIVATION.md`, and after the
ledger reconciliation V2-0 (commit cf8a0b38).

## 0. Status and history (binding)

- **V2 was designed after observing the V1 failures and diagnostics.**
- **V1 results remain negative:**
  - G3 FAIL;
  - G4 FAIL;
  - EMT04–EMT18 BLOCKED.

  They stay recorded as such (prereg v1, commit c2947bd8; amendment A1,
  9d06436a; results dd1e2ea3).
- **Only newly frozen holdouts can provide prospective evidence.** V2 is not a
  repeat of the original blind campaign.
- **V2 separates three things that V1 conflated:**
  1. exact device-transcription validation;
  2. the electromagnetic-network model discrepancy;
  3. the EMT interface of the grid-following converter.

**Prospective and non-prospective elements.**

| element | status |
|---|---|
| bus-30 SG and GFL unit cases | **regression** (seen in v1) |
| bus-36 SG, bus-35 GFL, bus-37 GFL unit cases | **blind holdouts**, frozen here |
| the voltage-source-behind-Rf–Lf GFL interface (all G4b runs) | new; never run before this commit |
| EMT04–EMT18 on IEEE-39 | never run; the v1 frozen predictions remain the targets |
| v1 diagnostic D5 (IEEE-39 no-event drift) | treated as unseen; used for no threshold, signal or window choice |

## 1. Unchanged from v1

The following carry over from prereg v1 without change:
- the simulator, upstream commit and serial solver;
- the TX4 network patch;
- the network realization (§2; damping off);
- the synchronous-machine realization (§3.1);
- the constant-power load (§3.3);
- the condenser and governor (§3.4);
- the time loop (§4);
- Δt 50 µs primary / 25 µs holdout, 1-ms storage;
- the common disturbance and both estimators (§5, `emt_estimator.py` unchanged);
- ε = 0.002;
- EMT04–EMT18 cases, points, frozen predictions (`results/EMT_PRED/`) and
  success rules (§6, §8, §11);
- manifests (§9);
- the no-retuning rule;
- the unit-test system (§7), with amendment A1.

## 2. Changed: the GFL EMT realization

**Primary V2 realization** (the only interface of this campaign): an
average-value controlled voltage source behind the physical Rf–Lf filter.
- It is specified completely in
  `docs/20260911_PAREMT_GFL_VOLTAGE_INTERFACE_DERIVATION.md`: the command (C),
  the companion (§6), the discretization with linear extrapolation (§7) and the
  initialization (§8).
- The filter currents are the physical i_d and i_q. The controller integrates
  only its 9 other states.
- ParaEMT's numerical damping stays **off**.

**What is retired.** The v1 ideal-current-source interface is retired as an EMT
interface. It is used **only** as the transcription harness of G4a, on an
algebraic network. There, no electromagnetic element sits at the device
terminal, so the harness integrates exactly the 11 TX4 states.

**Everywhere else.** The GFL of every IEEE-39 run (EMT04–EMT18) uses the V2
realization, with its own rating (w = S_n/100). This is the only change to the
v1 IEEE-39 model.

## 3. Unit-test cases (frozen)

**Unit-test system.**
- The v1 system: an EMT R–L line (R = 0, X = 0.05) and an infinite bus 1.0∠0
  behind j1e-4, solidly grounded (A1).
- **Algebraic variant** (G3a, G4a): the T–INF line is replaced by its
  synchronous-frame algebraic Norton element plus the zero-sequence coupling.
  This is exactly the v1 diagnostic D1 construction, `quasi_static_g`.

**Synchronous machine.** P4 settings: k = 1.425, t = 1.5, flux and AVR blends 1,
no governor. Disturbance: Pm ×1.02 on [1.0, 1.2) s; 10-s run.

| id | role | machine data | device output S (system pu) |
|---|---|---|---|
| R-SG30 | regression | bus 30 (Sn 1040) | 4.0 + j1.0 |
| B-SG36 | **blind** | bus 36 (Sn 1025.2) | 5.79999998 + j0.60708066 (canonical base-PF dispatch of bus 36) |

**Grid-following converter.** ConverterParameters defaults; leak 0.05;
g ∈ {0, 0.03625, 0.25}. Disturbances, applied at t = 1 s and held:
- (a) p_ref ×1.01;
- (b) infinite-bus |E| ×0.99;
- (c) E phase +0.02 rad.

Each run is 5 s: 9 runs per operating point.

| id | role | rating | S |
|---|---|---|---|
| R-GFL30 | regression | bus 30, w = 10.4 | 4.0 + j1.0 |
| B1-GFL35 | **blind** | bus 35, w = 10.857 | 2.0 + j0.5 |
| B2-GFL37 | **blind** | bus 37, w = 9.702 | 3.21521338 − j0.27617116 (canonical base-PF dispatch of bus 37) |

**Operating point.** V_T is the fixed point of `V_T = E + (z_s + jX) conj(S/V_T)`,
as in v1.

**References** (`.venv/tx3-analysis`, canonical device classes, Radau,
rtol 1e-10, atol 1e-12, piecewise across events):
- **Quasi-static network** (the frozen phasor model): used for G3a and G4a, and
  for reporting.
- **Dynamic-line network** (for G3b and G4b): the positive-sequence content of
  the EMT line is added as `(X/ω0) dI/dt`.
  - **SG:** I_line is a state, as in v1 D3.
  - **GFL:** the line current equals the device current `w i e^{jθ}`, and V_T
    solves
    `V_T = V_INF + (X/ω0)·w e^{jθ}(f_i + jθ′ i) + jX·w i e^{jθ}` with
    `V_INF = E + z_s I`.
    - It is solved by fixed-point iteration to |ΔV_T| ≤ 1e-14 at every
      right-hand-side evaluation, with at most 200 iterations; otherwise the
      run is an error.
    - The G4b references are output on the 50-µs grid.

## 4. Gates

**Order:** SW → G3a → G3b → G4a → G4b. Each gate passes only if **every** listed
case passes: regression and blind alike.

**Symbols.**
- `exc(x) = max over the window of |x_ref(t) − x_ref(0)|`.
- "2 % rule": `max |x − x_ref| ≤ 0.02 · exc(x)`.
- v1 rules are unchanged wherever v1 wording is cited.

### SW — interface software tests

Tests T1–T4 of the derivation document §9 must pass before any G4b run.

### G3a — SG transcription (algebraic network; v1 SG kernel; 50 µs)

Cases R-SG30 and B-SG36. The v1 EMT02 rules apply unchanged (not relaxed):
- finite;
- equilibrium |ΔP|, |ΔQ| ≤ 1e-4 at t = 0.9 s;
- every one of the 7 states obeys the 2 % rule against the **quasi-static**
  reference over [1, 10] s, on raw 1-ms samples;
- ringdown |Δα| ≤ 0.005 s⁻¹ and |Δf| ≤ 0.005 Hz, using the frozen
  `matrix_pencil` on ω − 1.

### G3b — SG electromagnetic embedding (the v1 EMT R–L line; 50 µs)

Cases R-SG30 and B-SG36, against the **dynamic-line** reference:
- finite;
- equilibrium |ΔP|, |ΔQ| ≤ 1e-4 at t = 0.9 s;
- every state obeys the 2 % rule over [1, 10] s, on raw 1-ms samples;
- ringdown |Δα| ≤ 0.005 s⁻¹ and |Δf| ≤ 0.005 Hz, using the frozen
  `matrix_pencil` on ω − 1.

**Reported, not gated:**
- EMT against the quasi-static reference: the v1-type discrepancy, reported as a
  **model-domain discrepancy** of about f/f0, not a transcription error;
- the dynamic-line against the quasi-static reference;
- the 25 µs runs.

### G4a — GFL transcription (algebraic network; the v1 11-state harness; 50 µs)

Cases R-GFL30, B1-GFL35 and B2-GFL37 (27 runs). The v1 EMT03 rules apply
unchanged, against the **quasi-static** reference:
- finite;
- equilibrium |ΔP|, |ΔQ| ≤ 1e-4 at t = 0.9 s;
- every one of the 11 states obeys the 2 % rule over [1, 5] s, on raw 1-ms
  samples;
- PLL lock |θ − ∠V_T| ≤ 1e-3 rad at T;
- for g = 0, case (b): |q_f − q_ref| ≤ 1e-3 at T.

The v1 D1 bus-30 result is a regression only.

### G4b — GFL physical EMT embedding (V2 interface; the v1 EMT R–L line)

**Runs.**
- The 27 event runs at 50 µs and again at 25 µs.
- 9 no-event runs of 5 s at 50 µs (3 operating points × 3 g).
- Every quantity is recorded at every step: the 9 controller states, the
  measured i_d and i_q, P and Q (system base, from the positive-sequence
  phasors of V_T and the branch current), and V_T.

**Low-frequency comparison operator LP.**
- A zero-phase 4th-order Butterworth low-pass with cutoff 15 Hz (`butter` +
  `sosfiltfilt`), applied at each signal's native rate:
  - EMT: 20 or 40 kHz;
  - reference: 20 kHz.
- It is then sampled on the 1-ms grid. The window is W = [1.0, 4.8] s.
- It removes content far above the phasor model's band: step-rate artifacts, the
  trapezoidal split mode, and the start-of-event doublet of derivation §7.

**Nyquist indicator.** For each 1-ms block k:

`A(k) = max over the steps n in the block of ||V_T,n| − |V_T,n−1||`

**Checks** (all against the **dynamic-line** reference; a floor of 1e-9 is added
to every tolerance of the form 0.02·exc or 0.01·exc):

1. **No-event stationarity** (9 runs). Over all steps of [0, 5] s:
   - `max ||V_T| − |V_T(0)|| ≤ 1e-4`;
   - `max |P − P0| ≤ 1e-4` and `max |Q − Q0| ≤ 1e-4` (system pu);
   - every recorded controller quantity: `max |x − x(0)| ≤ 1e-4`.
2. **No growing step-to-step alternation.**
   - No-event runs:
     - `max A` over [0.1, 5.0] s ≤ 1e-5;
     - `max A` over [4.5, 5.0] ≤ `max(max A over [0.5, 1.0], 1e-9)`.
   - Every event run, at 50 and at 25 µs:
     - `max A` over [4.5, 5.0] ≤ 1e-5;
     - `max A` over [4.5, 5.0] ≤ `max A` over [1.0, 1.5].
   - The log-slope of A over [1.5, 5.0] is reported.
3. **50 µs vs 25 µs** (27 pairs). For all 11 quantities and for P and Q:
   `max_W |LP(x_50) − LP(x_25)| ≤ 0.01 · exc(x) + 1e-9`.
4. **Fundamental-frequency P/Q** (27 runs at 50 µs).
   - Equilibrium: the means of P and Q over [0.8, 0.9] s are within 1e-4 of S.
   - `max_W |LP(P) − LP(P_ref)| ≤ 0.02 · exc(P) + 1e-9`, and the same for Q.
5. **Controller-state trajectories** (27 runs at 50 µs).
   - For all 11 quantities: `max_W |LP(x) − LP(x_ref)| ≤ 0.02 · exc(x) + 1e-9`.
   - PLL lock: `|mean over [4.9, 5.0] of wrap(θ − ∠V_T)| ≤ 1e-3` rad.
   - For g = 0, case (b): `|mean over [4.9, 5.0] of (q_f − q_ref)| ≤ 1e-3`.
6. **Ringdown** (27 runs at 50 µs). Estimator R is applied identically to the
   EMT run and to the reference.
   - **Input:** the 11 quantities, after LP and decimation to 100 Hz, over
     [1.05, 4.8] s. Each channel minus its window mean is scaled to unit RMS;
     channels with RMS < 1e-12 are dropped.
   - **Pencil:** a multi-channel matrix pencil with a stacked Hankel matrix,
     L = ⌊N/3⌋, and order #(σ_i/σ_1 ≥ 1e-4) clipped to [2, 30].
   - **Mode selection:** modes with Im s > 0 and 0.1 ≤ f ≤ 15 Hz whose energy
     is ≥ 1 % of that band's total.
   - **Resolved** iff at least one mode is selected and the reconstruction
     residual is ≤ 0.05.
   - **Output:** α = max Re over the selected modes, and f is that mode's
     frequency.
   - **Pass rule:**
     - If the reference is resolved, the check passes iff the EMT run is
       resolved and `|Δα| ≤ 0.005 + 0.02|α_ref|` and
       `|Δf| ≤ 0.005 + 0.02 f_ref`.
     - If the reference is not resolved, the check is N/A.

**Reported, not gated:**
- the unfiltered comparisons;
- EMT against the quasi-static reference;
- the dynamic-line against the quasi-static reference;
- the check-5 figures at 25 µs.

## 5. Stop rules

1. **Any gate FAIL: STOP.** Diagnose, report, and mark later phases BLOCKED. No
   threshold changes and no retuning.
2. **G4b.** If the voltage-source-behind-Rf–Lf interface cannot pass G4b on the
   blind holdouts: STOP.
   - Norton compensation, numerical damping and any other interface are **not**
     tried in this campaign. They would need a V3.
3. **Numerical-damping sensitivity** (secondary, optional). Only after G4b has
   passed without damping: a small subset may be rerun with the upstream damping
   option.
   - It is reported as a robustness sensitivity, with the ≈ 1.4e-3 change in
     60-Hz admittance.
   - It is never a headline result.
4. **Resuming the portfolio campaign.** Only if SW, G3a, G3b, G4a and G4b all
   pass: EMT04–EMT18 resume **unchanged** from prereg v1:
   - the frozen phasor predictions;
   - the P4 subset set;
   - the g points and k points;
   - the line holdouts;
   - the disturbance;
   - the estimators;
   - ε;
   - the success rules, including the v1 G2 and G5 rules of EMT04.

   The v1 diagnostic D5 drift is not used to modify thresholds, signals or
   windows.

## 6. Reporting

**Report.** `docs/20260911_PAREMT_EMT_V2_REPORT.{md,tex,pdf}`, side by side:
1. the v1 outcome;
2. the post hoc v1 diagnosis;
3. the V2 prospective design;
4. the V2 blind-holdout result;
5. the final EMT result.

**The v1 FAIL is never replaced by a V2 PASS.** If V2 passes, the wording is:

> "The first preregistered EMT realization failed its device-interface gates.
> A second preregistered realization, designed after diagnosing those failures
> and evaluated on new holdouts, passed …"

## 7. Files (to be created after this commit)

| file | role |
|---|---|
| `experiments/paremt_emt/overlay/tx4_emt_v2.py` | V2 kernel: the v1 kernel plus the GFL voltage interface; v1 kernel unchanged |
| `experiments/paremt_emt/v2/EMTV2_refs.py` | references (tx3-analysis) |
| `experiments/paremt_emt/v2/EMTV2_sw_tests.py` | gate SW |
| `experiments/paremt_emt/v2/EMTV2_gates.py` | G3a, G3b, G4a, G4b |
| `experiments/paremt_emt/v2/emtv2_rules.py` | comparison operators and the estimator R, as frozen above |
| `results/EMTV2/` | outputs, gate JSON and manifests |
