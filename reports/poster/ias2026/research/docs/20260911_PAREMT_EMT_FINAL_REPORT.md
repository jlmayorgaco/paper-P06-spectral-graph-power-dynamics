# ParaEMT EMT-validation campaign (IAS2026 / TX4) — final report

**Status: BLOCKED at the device gates.**

- **Network and operating point:** reproduced (gate G1 PASS).
- **Custom devices:** the preregistered unit tests of the custom synchronous
  machine (G3) and the custom grid-following converter (G4) both FAILED.
- **Consequence:** under the preregistered stopping rules, EMT04–EMT18 were not
  run. No portfolio, policy, design, remediation or nonlinear claim receives EMT
  evidence from this campaign.

| item | value |
|---|---|
| branch | `research/paremt-emt-validation` (from a70dac94, not pushed) |
| preregistration | `docs/20260911_PAREMT_EMT_PREREG_V1.md`, commit c2947bd8, before any TX4 EMT result; amendment A1, commit 9d06436a, before the unit-test results |
| simulator | NatLabRockies/ParaEMT_public, commit d79d735a4a587d56c5b88187d1a499195b6b2b84 (2026-01-15), serial LU |
| environment | `.venv/xtool-paremt`, CPython 3.11.10, exact upstream requirements |
| date | 2026-09-11 |

> **Erratum (V2-0 ledger reconciliation;
> `docs/20260911_PAREMT_EMT_LEDGER_RECONCILIATION.md`).** The first version of
> this report took claim statuses from the pre-final ledger
> `results/20260911_CLAIM_MATRIX.csv` instead of the final canonical ledger
> `results/20260911_FINAL_VALIDATION_MATRIX.csv` (V01–V30). The following
> statements are corrected in place:
> 1. **I04 / V21.** "I04 remains NOT YET TESTED" and "single-implementation
>    results" were wrong. The custom-GFL results were independently reproduced
>    in ANDES with the same equations (post-cumulant Phase 8, GATE 5 PASS;
>    V21 VALIDATED as a reproduction of computation). EMT reproduction was not
>    achieved.
> 2. **T08.** The Kundur holdout is 28/29 under the preregistered rule (the miss
>    is a real-pair coalescence), i.e. 28/28 crossings through the origin.
> 3. **§23 claim matrix and closing table.** Regenerated from the final ledger.
>
> **Unchanged:** every EMT measurement, gate and verdict (G1 PASS, G3 FAIL,
> G4 FAIL, EMT04–EMT18 BLOCKED).

---

## 1. Executive summary

**Question.** Do the principal observable predictions of the frozen
phasor-domain theory survive an independent electromagnetic-transient (EMT)
realization of the same IEEE-39 network and the same custom device and control
equations?

**Answer.** This campaign cannot say. Each stage stopped at a preregistered gate:

1. **EMT00 (tool).** The official ParaEMT IEEE-39 case installs, runs
   deterministically, serializes and resumes.
   - Its no-event run is **not** stationary: a 0.62 pu start-up jump.
   - Cause: upstream never applies the one transformer ratio its case carries.
     With that tap applied, the jump falls to 2.2e-4 pu.
   - The tool gate was recorded as "PASS (tool)". The stationarity check as first
     written failed; this is disclosed in §4.
2. **EMT01 (network, G1): PASS.**
   - The frozen TX4 network is built from ParaEMT's own companion elements.
   - Its 60-Hz admittance matches the canonical Ybus to 1.7e-16.
   - The EMT equilibrium matches the canonical power flow to 1.5e-6 pu in |V|
     and 5.9e-6 rad in angle (tolerances 1e-4 and 1e-3).
   - The official ParaEMT IEEE-39 differs from TX4 in 24 of 195 network rows,
     so it is not used.
3. **EMT02 (two-axis machine, G3): FAIL.**
   - Six of seven states track the phasor reference within 0.87 % of their
     excursion, and the ringdown agrees (Δα = −3.1e-4 s⁻¹, Δf = 1.2e-4 Hz).
   - The AVR state efd misses the 2 % trajectory limit: 2.26 %.
   - Diagnostics (not gates) attribute the excess to the electromagnetic
     dynamics of the EMT line, a physical effect of relative size f/f0 ≈ 2 %.
     The machine transcription itself agrees to 2e-4.
4. **EMT03 (11-state GFL, G4): FAIL.**
   - With the preregistered EMT line, an undamped step-to-step (Nyquist) mode
     of the trapezoidal line inductance is amplified by the ideal-current-source
     converter.
   - At g = 0 the three runs go non-finite at 0.76 s. At g > 0 they stay finite
     but saturate at a 0.11 pu alternation, giving state errors 3–1100× the
     excursion.
   - With the line realized algebraically (diagnostic only), all nine cases agree
     within 0.78 %.
5. **EMT04–EMT18: BLOCKED.**

**What changes in the claim ledger.** Nothing is weakened or strengthened,
except that V19/I01 (network and operating point) gains an EMT reproduction.
V21/I04 (the custom-GFL results in a second simulator) keeps its final-ledger
status: VALIDATED by the independent ANDES reproduction (Phase 8). The EMT
reproduction was not achieved. *(Corrected in V2-0; the first version wrongly
said "NOT YET TESTED".)*

**Next step.** A preregistration v2 is proposed for the author's decision
(Appendix A). It is not executed here.

## 2. Scientific questions

**Q-EMT.** Does an EMT realization of the same network and the same device
equations reproduce the following observable consequences of the frozen phasor
theory?
- the P4 portfolio trap (15 proper subsets of {30,33,35,37} stable, the full
  set unstable);
- the controller-only policy boundary (g\* ≈ 0.20768) and the k boundary
  (k\* ≈ 1.3046);
- port-guided design (the Newton direction; line reinforcement);
- the model-dependent remediations (condenser, governors);
- the near-Hopf finite-disturbance behaviour and the robustness holdout.

**Precondition (gate A, MODEL).** The network, the equilibrium and each custom
device must first be shown to be realized equation-equivalently in EMT
(EMT01–EMT04). This precondition failed at EMT02 and EMT03.

## 3. Evidence taxonomy

Every claim receives exactly one label:
- THEOREM (no EMT validation required);
- EMT CORROBORATED;
- EMT QUANTITATIVELY REPRODUCED;
- EMT REFUTED;
- EMT UNRESOLVED;
- NOT APPLICABLE TO EMT.

**Wording.** EMT can corroborate the *observable consequence* of a theorem. It
never proves one.
- T01–T06 and C01–C06 are not EMT targets.
- Connected cumulants remain secondary and get no EMT figure.

In this report, "BLOCKED" means "not run because a preregistered stopping gate
failed". Its taxonomy label is EMT UNRESOLVED.

## 4. ParaEMT provenance and environment

**Upstream.**
- `https://github.com/NatLabRockies/ParaEMT_public.git`, commit d79d735a
  ("Update README.md", 2026-01-15), cloned 2026-09-11.
- `external/ParaEMT_upstream` is pristine (diff sha256 98eb43ea…, empty).
- License: `LICENSE.md` (Alliance for Sustainable Energy, LLC).

**Working copy.** `external/ParaEMT_tx4`, same commit, changed only by
`experiments/paremt_emt/tx4_patch_paremt.py`:
- **`lib_numba.py` and `Lib_BW.py`, +86/−54 lines:**
  - exact from-side transformer taps in `numba_InitNet`, `numba_updateIhis`
    and `Re_Init`;
  - a switch for the numerical-damping resistors;
  - stock behaviour stays the default.
- **Untracked overlay files:** `tx4_emt.py` and `tx4_case.py`, the custom
  device kernel and case builder. They are copied from
  `experiments/paremt_emt/overlay/` by `EMT_setup_tx4.py`.
- **Diff hashes.**
  - Patch only: sha256 190b8611…
  - Patch plus overlay, as recorded in every run manifest: a61698cb…

**Environment.**
- `.venv/xtool-paremt`: CPython 3.11.10, with setuptools 72.1.0, xlrd 1.2.0,
  matplotlib 3.9.2, numba 0.60.0, numpy 2.0.1, pandas 2.2.2 and scipy 1.14.0.
  Only their dependencies were added; nothing was upgraded.
- **Files:** `results/EMT00/environment_manifest.txt`, `pip_freeze.txt`,
  `upstream_provenance.json`.
- **Platform:** Windows 11 Pro 10.0.26200, Intel Core Ultra 9 185H, 22 logical
  CPUs.
- **Solver.** Serial LU; the nxmetis/BBD path was not used or repaired.
- **Reproducibility.** BLAS and numba threads are pinned to 1, which is needed
  for bit-reproducibility.
- **Phasor references.** Computed in the frozen `.venv/tx3-analysis` with the
  canonical device classes, read-only.

**EMT00, stock smoke test** (`results/EMT00/`).
- **Run.** IEEE-39 (systemN = 3):
  - network dimension 117;
  - 39 buses, 34 lines, 12 transformers, 10 generators, 19 loads, 2 shunts;
  - Δt 50 µs, 261 µs per step.
- **Checks passed:**
  - initialization;
  - determinism (bit-identical);
  - serialization and reload (identical);
  - snapshot and resume;
  - a governor-reference step gives a sensible response (+0.0070 pu Pe).
  The exciter step is a no-op in the stock case.
- **Check failed:** stationarity of the no-event run. There is a 0.62 pu
  voltage jump in the first millisecond.
- **Diagnosis** (`tap_diagnostic.json`, rule written before the rerun).
  - The stock case carries one off-nominal ratio (6–31, k = 0.9714), which
    upstream never applies.
  - With only that tap applied, the jump is 2.2e-4 pu and the 10-s drift is at
    the level of the stock power-flow residual.
- **Gate record.** "PASS (tool)", with `smoke_threshold_verdict_as_first_written =
  FAIL` kept in `EMT00_gate.json`.
- **Disclosure.** This reclassification is a judgment: the defect is in the
  upstream case and code, not in the installation. It is recorded, not hidden.
  EMT00 has no TX4 scientific meaning.

## 5. Official IEEE-39 versus frozen TX4 IEEE-39

`results/EMT01/network_diff.csv` and `operating_point_diff.csv` compare the two
cases; details are in `docs/20260911_PAREMT_NETWORK_EQUIVALENCE.md`.

| | official ParaEMT | frozen TX4 |
|---|---|---|
| base, frequency, bus numbering | 100 MVA, 60 Hz, 1–39 | identical |
| 34 lines (R, X, B) and 12 transformer endpoints / R, X | — | identical |
| shunts (buses 4, 5) | — | identical |
| off-nominal taps | one (6–31, 0.9714, to-bus side, not applied upstream) | 11 taps on bus1 (e.g. 31→6 0.9, 19→33 1.07) |
| loads | differ at 8 buses (e.g. bus 39: 11.04 pu; load 7 Q: 8.4 pu) | ANDES `ieee39_full` data (bus 39: 4.0 pu; load 7 Q: 0.84 pu) |
| generator ratings / dispatch | 1000 MVA each (2000 at bus 39) | Sn 836–1684 MVA, frozen dispatch |

**Totals.** 171 of 195 network rows are identical. The official case is
therefore **not** used. A TX4-specific ParaEMT case is built from the frozen
data, which is never modified.

## 6. Equation-equivalent network construction (EMT01)

**Construction.** ParaEMT's own `numba_InitNet` trapezoidal companion models:
- lines as series R–L with C = b/(2ω0) at each end;
- transformers as series R–L with an exact ideal tap on bus1;
- shunts as capacitors (`shnt_gb = 100(g + jb)`).

The numerical damping resistors are switched off (`net_damping = 0`): they
change the 60-Hz admittance by about 1.4e-3.

| quantity | frozen prediction | EMT measurement | preregistered rule | verdict |
|---|---|---|---|---|
| 60-Hz Ybus from stored R/L/C/tap | canonical TX4 Ybus | max\|ΔY\|/max\|Y\| = 1.7e-16 (25/50/100 µs) | ≤ 1e-8 | PASS |
| trapezoidal warp at 60 Hz | (ω0Δt)²/12 | 7.4e-6 / 3.0e-5 / 1.2e-4 at 25/50/100 µs | reported, not gated | — |
| equilibrium \|ΔV\| (P4 base, 2 s no-event, last 0.5 s) | canonical PF | 1.48e-6 pu | ≤ 1e-4 | PASS |
| equilibrium \|Δθ\| rel. bus 39 | canonical PF | 5.85e-6 rad | ≤ 1e-3 | PASS |
| device P / Q residual | 0 | 1.7e-5 / 5.8e-5 pu | reported | — |
| speed deviation | 0 | 8.3e-8 pu | reported | — |

**Gate G1: PASS.** Figure F1 (`figures/20260911_EMT_F1_network_equivalence.pdf`)
shows the per-bus residuals. Two EMT01 CSV reruns were byte-identical.

**Interpretation.** The frozen network and its operating point are reproduced
by an independent EMT network realization. The only changes are documented
convention translations:
- the tap patch fixes an upstream omission;
- damping off;
- the shunt scaling.

This is the only claim-level EMT result of the campaign (I01).

## 7. Constant-power-load treatment

**Outcome A** (declared in prereg §3.3): an EMT-compatible
fundamental-frequency constant-power load.
- A constant-impedance Norton element `Y0 = conj(S0)/|V0|²`.
- A compensating injection `I_comp = −(conj(S(t)/V_f) − Y0 V_f)`, where V_f is a
  1-ms first-order filter of the extracted phasor.
- The incremental admittance is `Y0 + (Y_cp − Y0)/(1 + sτ_m)`: exact constant
  power as sτ_m → 0, with a relative deviation of 4.4e-3 at 0.7 Hz.

There is no silent substitution; ParaEMT's stock RLC load is not used.

| item | frozen prediction | EMT measurement | rule | verdict |
|---|---|---|---|---|
| equilibrium with the CP load | canonical PF | reproduced (§6) | EMT01 | PASS |
| low-frequency equivalence (τ_m 0.5 ms vs 1 ms holdout, EMT04 acceptance) | verdict, α and f invariant | not run | G2 | UNRESOLVED (BLOCKED) |

## 8. SG equation mapping and unit test (EMT02, gate G3)

**Mapping** (prereg §3.1). The two-axis machine is `SynchronousMachine`
(7 states: δ, ω, e′q, e′d, efd, pss_w, pss_l), with a first-order AVR, a
washout-lag PSS, D = 0 and no governor.
- Because x′d = x′q for all ten machines, the stator `I = (E′ − V)/(ra + jx′)`
  is realized **exactly** as a synchronous-frame algebraic Norton element:
  - a constant real 3×3 conductance `G_km = (2/3)Re(a^{m−k}Y)`, with
    Y = w/(ra + jx′);
  - an EMF current source driven by the predicted states.
- The terminal phasor is the space vector
  `V = (2/3)(v_a + a v_b + a² v_c)e^{−jω0t}`.

**Test system** (prereg §7 and amendment A1).
- **Line.** The device at T, then an EMT R–L line (R = 0, X = 0.05 pu), then
  INF.
- **Infinite bus.** Thevenin 1.0∠0 behind j1e-4 pu, solidly grounded through a
  zero-sequence-only conductance (A1: without it the network matrix was
  singular).
- **Operating point.** S = 4 + j1 pu, bus-30 machine data, P4 settings.
- **Disturbance.** Pm +2 % for 0.2 s at 1 s; 10-s run at 50 µs.
- **Reference.** The canonical class integrated by Radau (rtol 1e-10) with a
  quasi-static network.

| quantity | frozen reference | EMT measurement | preregistered rule | verdict |
|---|---|---|---|---|
| equilibrium ΔP, ΔQ | 0 | −8.3e-5, +2.3e-5 pu | ≤ 1e-4 | pass |
| δ, ω, e′q, e′d, pss_w, pss_l: max error / max excursion | 0 | 0.0032, 0.0026, 0.0087, 0.0018, 0.0013, 0.0016 | ≤ 0.02 | pass |
| **efd: max error / max excursion** | 0 | **0.0226** | ≤ 0.02 | **fail** |
| ringdown α, f (matrix pencil on ω) | −0.11837 s⁻¹, 1.17711 Hz | −0.11868 s⁻¹, 1.17723 Hz | \|Δα\| ≤ 0.005, \|Δf\| ≤ 0.005 | pass |

**Gate G3: FAIL.** No threshold was changed.

**Diagnostics (post hoc; they do not lift G3).** Files:
`results/EMT02/EMT02_diagnostics.json` and figure FD1.

| diagnostic | worst state ratio | reading |
|---|---|---|
| D1: the T–INF line as an algebraic synchronous-frame Norton (EMT network algebraic), 50 µs | 1.9e-4 | the machine transcription is exact to integration accuracy |
| D2: the preregistered EMT line at 25 µs | efd 0.0228 | the excess does not shrink with Δt, so it is not a discretization error |
| D3: the preregistered EMT line vs a phasor reference that adds the line's `(X/ω0) dI/dt` | 1.4e-3 | the excess is the electromagnetic dynamics of the line |

**Interpretation.**
- The error oscillates at the 1.18 Hz electromechanical frequency, with
  relative size f/f0 = 1.96 %. This is the known first-order dynamic-phasor
  correction that the quasi-static phasor model omits.
- The preregistered 2 % trajectory rule was therefore marginal for a
  comparison that deliberately crosses modelling domains.
- Because the failure is attributed post hoc, the gate stays FAILED.

## 9. GFL equation mapping and unit test (EMT03, gate G4)

**Mapping.** The equation-by-equation map is in
`docs/20260911_PAREMT_GFL_EQUATION_MAP.md`:
- all 11 states;
- the algebraic quantities, parameters and initialization, identical to
  `GridFollowingConverter` (GFL spec §4.2);
- the average-value interface: dq voltages from the extracted phasor, output
  `I = w(i_d + j i_q)e^{jθ}` injected as an ideal three-phase current source.

No stock converter is used, and no limits, DC link, PWM or protection are
added.

**Test.** The same system as §8, with g ∈ {0, 0.03625, 0.25} and three
disturbances at 1 s:
- (a) p_ref +1 %;
- (b) |E| −1 %;
- (c) E phase +0.02 rad.

Each run is 5 s at 50 µs.

| quantity | rule | EMT measurement | verdict |
|---|---|---|---|
| every state, 9 cases | error ≤ 0.02 × excursion | 3/99 state checks pass; g = 0: non-finite at t = 0.76 s (before the disturbance); g > 0: worst-state ratios 3.2–1103 | fail |
| equilibrium P/Q | ≤ 1e-4 | g = 0: none (non-finite); g = 0.03625: ΔP 2.5e-3, ΔQ −0.31; g = 0.25: ΔP −0.34, ΔQ 0.35 pu | fail |
| PLL steady state | \|θ − ∠V_T\| ≤ 1e-3 | 0.19–0.33 rad (g > 0); non-finite (g = 0) | fail |
| g = 0: q_f → q_ref | ≤ 1e-3 | non-finite | fail |

**Gate G4: FAIL.**

**Diagnostics.** Files: `EMT03_diagnostics.json`,
`EMT03_D1_quasi_static_line.csv`, `EMT03_D5_ieee39_chatter.csv` and figure FD2.
- **D4, mechanism.** With no disturbance, the step-to-step alternation of
  |V_T| grows:
  - at g = 0: exponentially at 26.6 s⁻¹ (factor 1.0013 per step) until the run
    fails;
  - at g = 0.25: at 445 s⁻¹, saturating at about 0.11 pu by 0.1 s.
  - Reading: an ideal current source feeding a lossless trapezoidal inductor
    has an undamped Nyquist mode (the companion voltage alternates sign each
    step). The explicitly coupled converter loop amplifies it.
- **D1, transcription.** With the line algebraic, all nine cases are finite and
  every state is within 0.78 % of its excursion. The transcription is correct;
  the failure is in the EMT interface.
- **D5, IEEE-39.** No-event 30-s runs of the preregistered realization: P4
  base, P4 H4, and H4 at G_S.
  - The start-up alternation at the four GFL buses (5e-6 pu) decays to
    2e-8 (P4 H4) and 7e-9 pu (G_S H4) over 1–10 s. No step-to-step
    instability develops, so the unit-test mechanism does not reproduce here.
  - In P4 H4 over 10–30 s, the consecutive-sample difference rises again to
    2.4e-7 pu, together with a slow |V| drift of 6.8e-5 pu. This smooth growth
    was **not analysed**.
  - No estimator was applied and no phasor prediction was compared.
    **D5 does not lift G4**: the preregistered rule is that no IEEE-39 result
    may be called equation-equivalent until EMT03 passes.

## 10. EMT initialization

**Procedure** (prereg §4).
- The network phasors come from the canonical power flow.
- The device states come from the phasor `initialize` formulas.
- ParaEMT sets the branch histories from the phasor steady state.
- The load filter starts at V_f = V0.

**Verified.**
- In EMT01, the P4 base drifts by at most 1.5e-6 pu.
- In the unit tests, the equilibrium P/Q residual is 8.3e-5 pu (SG).
- The GFL equilibrium was destroyed by the chatter of §9 (g = 0 and g > 0
  alike).

## 11. Numerical convergence (EMT04, gate G5)

| frozen prediction | EMT measurement | rule | verdict |
|---|---|---|---|
| verdicts of base, 30+33+35, H4 (P4) and H4 (G_S) invariant between 50 and 25 µs; \|Δf\| ≤ 0.01 Hz; \|Δα\| ≤ 0.01 s⁻¹ | not run | G5 | UNRESOLVED (BLOCKED by G3/G4) |

**Interpretation.**
- The only Δt evidence is the SG unit test: the efd error ratio is 0.0226 at
  50 µs and 0.0228 at 25 µs (D2). That is consistent with convergence but is
  not an EMT04 result.
- The 60-Hz warp is 3.0e-5 at 50 µs (§6).

## 12. Common disturbance and modal estimator

**Frozen (prereg §5), not exercised on any TX4 portfolio case.**
- **Disturbance.** +2 % of the bus-20 active load (0.136 pu) on [1.0, 1.2) s.
- **Signals.** Relative machine speeds, GFL frequency deviations, and
  bus-angle differences and |V| at the ten generator buses.
- **Estimator A.** Multi-channel matrix pencil at 10 Hz.
- **Estimator B.** Hilbert envelope of the dominant band-passed channel at
  50 Hz.
- **Band and window.** 0.2–1.2 Hz, window [2.2 s, T].
- **Classification.** STABLE / UNSTABLE if both estimators are resolved and
  beyond ±ε = 0.002 s⁻¹; UNRESOLVED otherwise.
- **Implementation.** `experiments/paremt_emt/emt_estimator.py`, used here only
  for the SG ringdown of §8.

## 13–22. Core experiments (all BLOCKED)

Every experiment below was preregistered with frozen phasor predictions
(`results/EMT_PRED/phasor_predictions.csv`, commit c2947bd8). None was run,
because gate A (MODEL) failed at G3 and G4.

**For each experiment:**
- EMT measurement: none;
- verdict: **UNRESOLVED (BLOCKED)**;
- interpretation: no EMT evidence either way.

### 13. EMT05 — P4 portfolio trap

- **Frozen prediction.**
  - The 15 proper subsets of {30,33,35,37} (including ∅) are stable, with
    in-band α from −0.214 to −0.144 s⁻¹ at 0.64–0.92 Hz.
  - H4 is unstable: α = +0.1270 s⁻¹ at 0.622 Hz.
- **Rule.**
  - All 15 proper subsets STABLE and H4 UNSTABLE (κ_EMT = 4).
  - Secondary: MAE(α) ≤ 0.03 s⁻¹ and \|Δf\| ≤ 0.05 Hz.
  - Amplitude checks at 1 % and 4 %.

### 14. EMT06 — contextual intervention sign reversal

- **Prediction.** Each of 30, 33, 35 and 37 is stabilizing alone and
  destabilizing when added last (the signs of 8 effects).
- **Rule.** All 8 EMT signs equal the phasor signs.

### 15. EMT07–08 — g-only controller boundary

- **Prediction.** At H4, k = 1.425, α(g) = +0.0138, +0.0036, +0.0012,
  −0.0011, −0.0076, −0.0174 s⁻¹ at g = 0.18, 0.20, 0.205, 0.21, 0.225, 0.25;
  g\* = 0.20768.
- **Rule.** A resolved UNSTABLE point below and a STABLE point above g\*_EMT.
  Target: \|g\*_EMT − 0.20768\| ≤ 0.02.

**Also blocked.**
- EMT09 (P_inf H4): prediction α = −0.155 s⁻¹; rule STABLE.
- EMT10 (Newton design):
  - prediction: first iterate g = 0.08806 has α = +0.0739; g\* = 0.20768 has
    α = 0; g = 0.25 has α = −0.0174;
  - rule: α_EMT monotone along the sequence, g = 0.25 STABLE, and
    \|g\*_EMT − g\*_Newton\| ≤ 0.02.

### 16. EMT11 — k boundary

- **Prediction.** At g = 0.03625, α(k) = −0.0809, −0.0337, −0.0061,
  +0.0069, +0.0313, +0.0539 at k = 1.25, 1.28, 1.30, 1.31, 1.33, 1.35;
  k\* = 1.3046267.
- **Rule.** \|k\*_EMT − 1.3046267\| ≤ 0.02, with a clear STABLE → UNSTABLE
  change.

### 17. EMT12 — line reinforcement

- **Prediction.** Exact phasor finite changes Δα for 12 frozen holdout lines
  ×1.5 at P4 H4:
  - from −0.0371 (L42), −0.0340 (L35) and −0.0242 (L26);
  - to +0.0013 (L17) and +0.0042 (L31).
- **Rule.** Sign agreement ≥ 10/12 and Spearman ≥ 0.80.
- **Not claimed.** Transfer across converter models.

### 18. EMT13 — condenser

- **Prediction.** At P4 H4, the damped condenser at 2.0 / 2.5 / 3.0 / 5.0 %
  gives in-band α = +0.0264 / −0.0010 / −0.0295 / −0.1569 s⁻¹.
- **Rule.** 2 % UNSTABLE and 5 % STABLE, in band.

### 19. EMT14 — governors

- **Prediction.**
  - Governed H4 at P4: α = −0.0745 s⁻¹ (stable).
  - Governed H4 at g = 0.020: α = +0.0691 s⁻¹ (unstable, descriptive).
- **Rule.** The governed P4 H4 is STABLE.

### 20. EMT15 — nonlinear recovery near the subcritical Hopf

- **Prediction.** H4 at g = 0.22 (α = −0.0055) and g = 0.30 (−0.0334).
- **Protocol.** Pulse amplitudes of 2–80 % with RETURNS / LEAVES /
  FAILS / NUMERICAL labels.
- **Rule.** Characterize the outcome, whether positive or negative.

### 21. EMT16 — robustness holdout

- **Prediction.** From the frozen PCV05 draws (seed 20260921; 74 cases in
  `emt16_selection.json`):
  - A: the witness is H4;
  - B: the witness is {30,33,35};
  - C: H4 is stable.
- **Rule.** Each set passes in ≥ 80 % of its draws.

### 22. EMT17 — optional Kundur zero-frequency check

Not run: it is optional, and the campaign is blocked. T08 rests on its proof
and on the existing phasor holdout: 28/29 real-count changes under the
preregistered rule, i.e. 28/28 crossings through the origin, with 0/445 false
positives.

## 23. Claim-to-EMT matrix

File: `results/20260911_EMT_CLAIM_MATRIX.csv`, generated by
`experiments/paremt_emt/EMT_claim_matrix.py`. *(Regenerated in V2-0 from the
final ledger.)*

**Rows** (43):
- the final canonical ledger `results/20260911_FINAL_VALIDATION_MATRIX.csv`
  (V01–V30);
- the 13 older ids that no V row cross-references (T04, T06, T08, C08, B04, B06,
  B07, B10, N01, N05–N08), from `results/20260911_CLAIM_MATRIX.csv`.

Both ledgers are read only. The column `independent_phasor_reproduction` is kept
separate from `EMT_result`.

| EMT label | claims |
|---|---|
| THEOREM | V01, V02, V08, V22, V26, T04, T06, T08 (8) |
| NOT APPLICABLE TO EMT | V05, V09, V12, V23, V24, V25, V27, V28, C08, N01, N05–N08 (14) |
| EMT QUANTITATIVELY REPRODUCED | V19 (1) |
| EMT UNRESOLVED | V03, V04, V06, V07, V10, V11, V13–V18, V20 (G3 FAIL), V21 (G4 FAIL), V29, V30, B04, B06, B07, B10 (20; BLOCKED except V18 and B06, which were not EMT targets) |
| EMT CORROBORATED / EMT REFUTED | none |

## 24. Negative and unresolved results

**Negative (preregistered gates failed).**
1. **G3.** The SG efd trajectory error is 2.26 % against a 2 % limit.
2. **G4.** The GFL unit test is numerically unstable with the preregistered
   EMT line: 3/9 runs non-finite, 6/9 grossly wrong.
3. **EMT00.** The stock no-event stationarity check failed as first written
   (upstream never applies the transformer tap).

**Unresolved.** EMT04–EMT18 are BLOCKED. Every portfolio, policy, design,
remediation, nonlinear and robustness question remains unanswered by EMT.

**Not done.**
- No retuning.
- No threshold change.
- No stock GFL or machine substitution.
- No constant-power-load substitution.
- No IEEE-39 portfolio comparison.
- The D5 no-event runs were used only to see whether the chatter mechanism
  appears there.

## 25. Limitations

- **Diagnostics are post hoc.** D1–D5 were designed after the gate failures.
  They explain the failures; they do not replace the gates.
- **The unit-test system is part of the preregistration.** A lossless line to
  a stiff source is a demanding topology for an ideal current source. D5 shows
  the IEEE-39 GFL buses behave differently, but the test was chosen in advance
  and failed.
- **The GFL interface is a design choice.** It is an average-value ideal
  current source with explicit predictor–corrector coupling (prereg §4).
  Other interfaces were not preregistered, so none was tried as a gate:
  - a converter voltage behind an EMT filter branch;
  - Norton-compensated injection;
  - ParaEMT's standard numerical damping.
- **A domain-crossing comparison.** The 2 % trajectory rule compares a
  quasi-static phasor reference with an EMT network. The network's own
  electromagnetic dynamics contribute about f/f0 of the excursion at an
  electromechanical frequency f.
- **Scope of EMT01.** It covers the base portfolio at P4 only (SGs and loads).
  The GFL equilibrium in IEEE-39 was not gated.

## 26. Manuscript implications

**No wording changes.**
- The manuscript already states that its simulations are "nonlinear
  phasor-domain simulation (not EMT)" (`docs/FINAL_IAS_SAFE_CLAIMS.md`,
  `docs/FINAL_REJECTED_WORDING.md`). That wording remains correct and required.
- The manuscript must **not** cite this campaign as EMT validation of any
  portfolio, policy or design result.
- The custom-GFL results were reproduced in a second phasor-domain tool (ANDES,
  same equations; V21, Phase 8). They were **not** reproduced in EMT. *(Corrected
  in V2-0; the first version wrongly said "NOT YET TESTED / single-implementation".)*

**Optional sentences** (limitations / reproducibility):
- *"The IEEE-39 network and operating point were independently realized in the
  ParaEMT electromagnetic-transient solver (60-Hz admittance residual 2e-16,
  equilibrium residual 2e-6 pu)."*
- *"An equation-equivalent EMT realization of the converter did not pass its
  preregistered unit test (numerical interface instability); EMT corroboration
  of the portfolio results is therefore not claimed."*

## 27. Reproducibility instructions

Paths are relative to `reports/poster/ias2026/research/`. Python for EMT is
`.venv/xtool-paremt`; for phasor references it is `.venv/tx3-analysis`
(read-only).

1. **Clone and patch.** Clone ParaEMT_public at d79d735a into
   `external/ParaEMT_upstream` and `external/ParaEMT_tx4`, then run
   `experiments/paremt_emt/tx4_patch_paremt.py` and `EMT_setup_tx4.py` (the
   overlay copy).
2. **EMT00.** `EMT00_provenance.py`, `EMT00_stock_smoke.py`,
   `EMT00_tap_diagnostic.py`, `EMT00_gate.py`.
3. **Frozen predictions** (tx3-analysis, already committed):
   `EMT_phasor_predictions.py`, `EMT_export_operating_points.py`.
4. **EMT01.** `EMT01_network_equivalence.py`.
5. **Unit-test references** (tx3-analysis): `EMT02_03_phasor_refs.py`,
   `EMT02_03_dynphasor_ref.py`.
6. **Unit tests and diagnostics** (xtool-paremt): `EMT02_03_unit_tests.py`,
   `EMT02_03_diagnostics.py`, `EMT03_D5_ieee39_chatter.py`. Then
   `EMT02_03_manifest.py` re-executes all three and writes
   `results/EMT02/EMT02_03_manifest.json`.
7. **Claim matrix and figures.** `EMT_claim_matrix.py`, `EMT_figures.py`.
8. **Report.** Build `docs/20260911_PAREMT_EMT_FINAL_REPORT.tex` with pdflatex,
   twice.

**Determinism.**
- Every summary CSV/JSON of EMT01–EMT03 and the diagnostics was regenerated
  twice and is byte-identical. Wall-clock fields sit only in the manifests.
- The phasor references are bit-identical after pinning BLAS to one thread.
  Before pinning they differed by up to 7e-12.
- The claim matrix and the figure-source CSV are byte-identical on rerun.

**Manifests.**
- EMT01: `results/EMT01/EMT01_summary.json` (`manifest`).
- EMT02/EMT03 and diagnostics: `results/EMT02/EMT02_03_manifest.json`, with
  run ids, git HEAD, upstream SHA, working-copy diff SHA, environment hash,
  Δt, duration, disturbance, seed (none), per-script wall time and the output
  sha256.
- Raw EMT01 time series: `external/paremt_runs/raw/` (NPZ, hash in the
  manifest). Unit-test trajectories: `results/EMT0{2,3}/emt_traj.npz`.

**Figures.**
- F1: `figures/20260911_EMT_F1_network_equivalence.{pdf,png}`.
- Diagnostic figures: FD1 `figures/20260911_EMT_FD1_sg_unit_test.{pdf,png}` and
  FD2 `figures/20260911_EMT_FD2_gfl_unit_test.{pdf,png}`.
- F2–F8 need EMT04–EMT16 and are **not produced (BLOCKED)**.
- No cumulant figure.

---

## Answers to the final questions

1. **Is the frozen TX4 IEEE-39 faithfully reproduced in ParaEMT?**
   Yes, for the network and the operating point:
   - Ybus residual 1.7e-16;
   - equilibrium |ΔV| 1.5e-6 pu and |Δθ| 5.9e-6 rad (base portfolio at P4).

   This required a documented tap patch, because upstream never applies
   transformer ratios. The official ParaEMT IEEE-39 is a different case.
2. **Is the custom 10+1-state GFL equation-equivalent in EMT at fundamental
   frequency?** Not established.
   - The equations are transcribed identically, and match numerically within
     0.78 % when the line is algebraic (diagnostic).
   - But the preregistered EMT unit test failed (G4): the ideal-current-source
     interface is numerically unstable with the preregistered EMT line.
3. **Is the two-axis SG equation-equivalent?** Not established by the gate.
   - G3 failed on efd (2.26 % against 2 %); the other states pass and the
     ringdown agrees to 3e-4 s⁻¹.
   - The diagnostics (transcription exact to 2e-4; excess explained by the
     line's dI/dt) suggest the machine is equivalent, but that is a post hoc
     reading.
4. **Proper P4 subsets stable and full H4 unstable?** Not tested (BLOCKED).
5. **κ_EMT at P4?** Unknown (BLOCKED).
6. **α_EMT and f_EMT for the 16 subsets?** Not measured (BLOCKED). The phasor
   predictions are listed in §13.
7. **Contextual sign reversal?** Not tested (BLOCKED).
8. **g-only transition?** Not tested (BLOCKED).
9. **g\*_EMT vs 0.20768?** Not measured (BLOCKED).
10. **k boundary?** Not tested (BLOCKED).
11. **Port-guided retuning?** Not tested (BLOCKED).
12. **Line sensitivities vs finite EMT reinforcement?** Not tested (BLOCKED).
13. **Condenser remediation?** Not tested (BLOCKED).
14. **Governor stabilizes the flagship?** Not tested (BLOCKED).
15. **Finite recovery threshold near the subcritical Hopf?** Not tested
    (BLOCKED). There is no EMT evidence either way.
16. **Robustness holdout?** Not tested (BLOCKED).
17. **Which claims are not EMT-testable?**
    - T01–T07 and C01–C06, which remain proofs.
    - T08, a theorem whose optional EMT17 was not run.
    - The algebraic identities behind T09 and T10: the T10 derivative formula
      is a theorem, and only its observable consequences are EMT-testable.
    - The benchmark cumulant results C07–C09 and the refuted C10–C11 are
      phasor-port constructions with no EMT observable.
18. **Did any manuscript statement become too strong?**
    - No statement became too strong because of EMT: no portfolio-level EMT
      result exists.
    - The existing "phasor-domain, not EMT" wording must stay.
    - The campaign must not be cited as EMT validation.
    - The implicit assumption that a quasi-static network is exact is now
      measured: trajectory deviations of about f/f0 (2 % at 1.2 Hz). This
      is within the manuscript's declared "no claim above a few Hz" scope.
19. **New model-fidelity limitations?** Yes, three:
    - (i) The average-value GFL, as an ideal current source with explicit
      coupling, is numerically fragile in EMT when its terminal sees only
      lossless inductance to a stiff source (Nyquist chatter). A faithful EMT
      transplant needs a deliberate interface design.
    - (ii) The quasi-static network omits the line's dI/dt, which changes
      machine trajectories by about f/f0.
    - (iii) Upstream ParaEMT ignores transformer ratios, and its IEEE-39 case
      is not the TX4 case.
20. **Strongest defensible statement:**
    - **Mathematical proof:** the theorems stand (T01–T08, C01–C06).
    - **Python phasor model:** the frozen phasor results stand.
    - **ANDES:** it independently reproduces, with the same equations, the
      network, the synchronous dynamics and branch sensitivities (V19, V20), and
      the custom-GFL portfolio results: P4 trap 32/32, H and κ equal, P4 vs G_S
      (V21, Phase 8).
    - **Nonlinear phasor TDS:** it confirms the small-signal verdicts (B09).
    - **ParaEMT EMT:** it independently reproduces the network and operating
      point (I01), but it provides **no** corroboration of the portfolio,
      policy, design, remediation or nonlinear claims. Its preregistered device
      gates failed, so those claims remain phasor-domain results. They are
      reproduced in two phasor tools (V21) but not in EMT. *(Corrected in
      V2-0.)*

---

## Appendix A — proposed preregistration v2 (for the author's decision; not executed)

> **Superseded (V2-0).** The author approved a different V2. It separates device
> transcription from EMT embedding (G3a/G3b, G4a/G4b) and preregisters a single
> voltage-source-behind-Rf–Lf GFL interface. Numerical damping and Norton
> compensation are excluded as primary solutions. This appendix is kept only as
> a record and is not executed.

Any v2 is designed **after** seeing the v1 failures and the diagnostics. It must
say so, and it must add blind elements so that it is not a forced match.

**V2-1. Device gates on an algebraic-network unit test.**
- Unit-test system: the line realized as a synchronous-frame algebraic Norton,
  so the test isolates the device transcription. The 2 % trajectory rule and
  the equilibrium, ringdown and PLL rules are unchanged.
- Blind holdout, fixed in v2 before any run: a second operating point and
  device (e.g. bus-35 data, S = 2 + j0.5) and a second disturbance set.
- **Caveat.** v1's D1 already passed on the v1 operating point, so only the
  holdout is blind.

**V2-2. Network electromagnetic effect as a reported quantity, not a gate.**
The EMT-line unit test is compared against the dynamic-phasor reference (D3
type) and reported as a model-fidelity quantity.

**V2-3. One preregistered GFL EMT interface, chosen before any IEEE-39 run.**
Candidates:
- (a) the converter voltage source behind the filter reactor, realized as an EMT
  R–L branch. i_d and i_q come from the branch current; the phasor model's
  nominal cross-coupling differs by the θ′ term, which must be bounded;
- (b) the ideal current source with ParaEMT's standard numerical damping
  switched on. The 1.4e-3 change in 60-Hz admittance is then declared, and
  checked by an EMT04-style holdout with damping off where stable;
- (c) a Norton-compensated current injection.

Acceptance: the preregistered EMT-line GFL unit test must be finite and stable
for all nine cases, and (a) must additionally pass V2-1.

**V2-4. Unchanged.**
- EMT04–EMT16 as in v1, with unchanged thresholds, disturbance, estimator and
  frozen predictions.

**V2-5. Declared history.** v1 failed G3 and G4. That record stays in the
ledger regardless of v2.

---

## What EMT did and did not validate

The table separates **theorem validity** (mathematical status, column 2) from
**model-fidelity corroboration** (columns 3–5).

*(Rebuilt in V2-0 from the final canonical ledger; V-ids, with old ids in
brackets.)*

| Theory / claim | Mathematical status | Phasor evidence | Independent phasor evidence | EMT evidence | Final status |
|---|---|---|---|---|---|
| V01 transverse quotient [T01] | proved | S1 (4.6e-10) | ANDES structural pair 32/32 | not applicable (THEOREM) | PROVED + REPRODUCED |
| V02, V26, T04, T06 (antichain / any-order planning, NP-completeness, Metzler, complex) | proved | used throughout | — | not applicable (THEOREM) | PROVED |
| V08 factorization [T07, T09]; V22 cumulant theorems [C01–C06] | proved | PCV03, FC18; CC02 | — | not applicable (THEOREM) | PROVED (+ benchmark-specific boundary statement for V08) |
| T08 zero-frequency port | proved | Kundur holdout 28/29 (28/28 origin crossings), 0/445 FP | — | EMT17 not run | PROVED |
| V03 P4 trap, κ = 4 [B02] | — | PCV02; G2 TDS | ANDES Phase 8: 32/32, H and κ equal | BLOCKED (EMT05/06) | VALIDATED (model-specific witness) |
| V04 lower-order screening approves H4 | follows from V03 | PCV02/03 | via V03 | BLOCKED (EMT05) | VALIDATED |
| V06, V07 policy-dependent incompatibility [B01] | structure proved | clean g-only counterfactual; F7 | ANDES Phase 8: P4 vs G_S | BLOCKED (EMT07–08/10/11/16) | VALIDATED |
| V10 port derivative, Newton [T10] | proved | PCV04 (g, k Newton) | ANDES R3 11/12, ρ 1.00 | BLOCKED (EMT08/10/11/12) | PROVED + VALIDATED |
| V11 line ranking | — | ρ 0.937; finite sign 9/12 | line sensitivities equation-equivalent in ANDES | BLOCKED (EMT12) | VALIDATED (local ranking) |
| V13 remediation, V14 governors [B03, B05] | — | PCV04; FC03; G2 | ANDES Phase 8 (g side); governor direction | BLOCKED (EMT08/10/11/13/14) | VALIDATED (model-specific); BENCHMARK-SPECIFIC |
| V15–V18 robustness | — | PCV05 | — | BLOCKED (EMT16) / not targeted (V18) | R1 PASS 3/4; R2 NOT ROBUST; R3 descriptive; R4 mixed |
| V30, B07, B10 nonlinear scope [B08, B09] | normal form | G2 TDS 32/32; FC05–07; FC13 | — | BLOCKED (EMT15, EMT05) | NUMERICALLY VALIDATED |
| V19 network and operating point [I01] | — | reference | ANDES 1e-13; pandapower 3e-14 after tap fix | **EMT QUANTITATIVELY REPRODUCED** (Ybus 2e-16, V 1.5e-6) | CONDITIONAL (pandapower translation; memo A) |
| V20 SG dynamics and sensitivities [I02, I03] | — | internal | ANDES R2/R3 12/12; E1 12/12 | G3 FAIL (efd 2.26 %) → EMT UNRESOLVED; EMT12 BLOCKED | VALIDATED |
| V21 custom-GFL results in a second simulator [I04] | — | project DAE | **ANDES Phase 8, GATE 5 PASS** (same equations) | G4 FAIL → EMT UNRESOLVED | VALIDATED as a reproduction of computation; EMT reproduction not achieved |
| V09, V23, V24, V25, C08 (cumulants, closure diagnostics) | V22 | CC campaign; PCV02 | — | not applicable | BENCHMARK-SPECIFIC / DIAGNOSTIC / REFUTED |
| V05 F10 recount; V12 finite-step nonlinearity | — | PCV02; PCV04 | — | not applicable | VALIDATED; NEGATIVE RESULT |
| B04, B06 | — | FC04; FC10 | — | BLOCKED / not targeted | BENCHMARK-SPECIFIC |
| V27, V28, V29, N01, N05–N08 negatives | — | refuted | — | not applicable / BLOCKED (V29) | REFUTED (kept) |
