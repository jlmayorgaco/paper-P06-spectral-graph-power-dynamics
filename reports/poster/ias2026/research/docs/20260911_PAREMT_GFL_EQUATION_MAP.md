# ParaEMT campaign — GFL equation map (prereg §3.2, §7)

**Implementations compared.**
- **Phasor.** `GridFollowingConverter.derivatives / injection / initialize`,
  as specified in `docs/20260911_GFL_REPRODUCTION_SPEC.md` §4.2.
- **EMT.** The jitted kernel `_gfl_eval` in
  `experiments/paremt_emt/overlay/tx4_emt.py` (copied to
  `external/ParaEMT_tx4/tx4_emt.py`), initialized by `gfl_init` in
  `overlay/tx4_case.py`.

No stock ParaEMT converter model is used.

## 1. Interface

| item | phasor model | EMT realization |
|---|---|---|
| terminal voltage | algebraic network phasor V | `V(t) = (2/3)(v_a + a v_b + a² v_c) e^{−jω0 t}` from the abc solution (`_phasors`) |
| PLL-frame voltage | `v_d + j v_q = V e^{−jθ}` | `rot = V(t)·(cos θ − j sin θ)`, identical |
| \|V\| for the Q/V loop | \|V\| | `abs(V(t))`, identical |
| output | `I = w (i_d + j i_q) e^{jθ}` (system base) | the same phasor, injected as `i_k = Re(I e^{j(ω0 t − 2πk/3)})`, k = 0, 1, 2 (`_inject`): an ideal three-phase current source |
| frame | synchronous, ω0 = 2π·60 | the same (the injection and extraction rotations cancel for positive sequence) |

## 2. States and equations (line by line)

In `_gfl_eval` the states are indexed x[0..10]. Every row below is the same
expression on both sides.

| state | spec §4.2 expression | kernel line |
|---|---|---|
| θ (x0) | `θ' = kp_pll v_q + x_pll` | `f[0] = KPPLL*v_q + x[1]` |
| x_pll (x1) | `ki_pll v_q` | `f[1] = KIPLL*v_q` |
| p_f (x2) | `(P − p_f)/τ_p`, with `P = v_d i_d + v_q i_q` | `f[2] = (power − x[2])/TAUP` |
| q_f (x3) | `(Q − q_f)/τ_p`, with `Q = v_q i_d − v_d i_q` | `f[3] = (reactive − x[3])/TAUP` |
| x_p (x4) | `ki_p (p_ref − p_f)` | `f[4] = KIP*(pref_eff − x[2])` |
| x_q (x5) | `ki_q (q_cmd − q_f)` | `f[5] = KIQ*(q_cmd − x[3])` |
| i_d (x6) | `(ω_B/xf)(e_d − v_d − rf i_d + xf i_q)` | `f[6] = (W0/xf)*(e_d − v_d − rf*i_d + xf*i_q)` |
| i_q (x7) | `(ω_B/xf)(e_q − v_q − rf i_q − xf i_d)` | `f[7] = (W0/xf)*(e_q − v_q − rf*i_q − xf*i_d)` |
| x_id (x8) | `ki_i (id_ref − i_d)` | `f[8] = KII*(id_ref − i_d)` |
| x_iq (x9) | `ki_i (iq_ref − i_q)` | `f[9] = KII*(iq_ref − i_q)` |
| x_v (x10) | `ki_v·err − L (x_v − q_ref)` | `f[10] = KIV*err − LEAK*(x[10] − QREF)` |

The algebraic quantities are also identical:
- `err = g (v_ref − |V|)`;
- `q_cmd = kp_v err + x_v`;
- `id_ref = kp_p (p_ref − p_f) + x_p`;
- `iq_ref = −(kp_q (q_cmd − q_f) + x_q)`;
- `e_d = v_d + kp_i (id_ref − i_d) + x_id − xf i_q`;
- `e_q = v_q + kp_i (iq_ref − i_q) + x_iq + xf i_d`.

`ω_B` is the nominal frequency, 2π·60, as in the phasor model (nominal
cross-coupling).

## 3. Parameters and initialization

**Parameters.** `GFL_DEFAULTS` holds the internal `ConverterParameters`
defaults:
- kp_pll 53, ki_pll 1400, τ_p 0.03;
- kp_p 0.20, ki_p 8, kp_q 0.20, ki_q 8;
- kp_i 0.25, ki_i 6;
- xf 0.15, rf 0.01, kp_v 2, ki_v 20.

The leak L is 0.05 rad/s (frozen), and g is the policy coordinate.

**Initialization** (`gfl_init`, same as spec §4.2):
1. `I = conj(S/w/V)`, `θ0 = ∠V`;
2. `i_d + j i_q = I e^{−jθ0}`;
3. `p_ref = |V| i_d`, `q_ref = −|V| i_q`, `v_ref = |V|`;
4. x0 = (θ0, 0, p_ref, q_ref, i_d, −i_q, i_d, i_q, rf i_d, rf i_q, q_ref).

**Disturbances in the kernel:**
- event 3 multiplies p_ref (EMT03 a);
- events 4 and 5 act on the unit-test source (EMT03 b, c).

## 4. What is not the same (declared)

**(i) Coupling to the network.**
- In the phasor model V is algebraic: the network solution at the same
  instant.
- In EMT, the state derivative at step n uses V from the network solve of
  that step. The explicit Heun predictor–corrector of prereg §4 is used, and
  the GFL current is set from the predicted state.
- The GFL is an ideal current source. At a node whose only other connection
  is a lossless inductive branch, the trapezoidal companion of that branch
  has an undamped Nyquist (step-to-step alternating) mode. The
  explicitly coupled converter loop can destabilize it (see §5).

**(ii) Network electromagnetic dynamics.** The EMT line/transformer adds
`(X/ω0) dI/dt` to the quasi-static phasor relation. This is physical, not an
equation difference (EMT02 diagnostic D3).

**(iii) Not modelled on either side:** current limits, DC link, PWM,
saturation, protection and ride-through.

## 5. Verification status

**Preregistered EMT03, gate G4: FAIL** (`results/EMT03/EMT03_summary.json`).
- **g = 0 (cases a, b, c).** The run becomes non-finite at t = 0.76 s,
  before the disturbance.
- **g = 0.03625 and 0.25.** The runs stay finite, but the state errors are
  3 to 1100 times the phasor excursion. The equilibrium P/Q errors reach
  0.35 pu, and the PLL error 0.19–0.33 rad.
- **Cause (D4,** `EMT03_diagnostics.json`**).** A step-to-step alternation of
  the terminal |V| grows from 1.5e-4 pu at 10 ms to about 0.11 pu by 0.1 s,
  where it saturates.

**Transcription check** (diagnostic D1, `EMT03_D1_quasi_static_line.csv`).
Here the T–INF line is realized as a synchronous-frame algebraic Norton, so
the EMT network is algebraic:
- all 9 cases stay finite;
- every state is within 0.78 % of its phasor excursion.

The transcription above is therefore verified. The preregistered EMT03
failure is an **interface / numerical-stability** failure of the current-source
realization in the preregistered test topology. It is not an equation
mismatch.

**IEEE-39 (diagnostic D5,** `EMT03_D5_ieee39_chatter.csv`**; no-event runs
only, with no estimator and no prediction comparison).** The same alternation
does not develop: it decays below 1e-8 pu at the four GFL buses by 1–10 s in
P4 H4 and G_S H4.

This does **not** lift G4. The preregistered unit test failed, and EMT04–EMT18
stay BLOCKED under preregistration v1.
