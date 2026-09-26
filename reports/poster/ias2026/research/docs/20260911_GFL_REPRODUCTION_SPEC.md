# Custom-GFL reproduction in ANDES: specification (Phase 8)

Date: 2026-09-11. Preregistration: `docs/20260911_POST_CUMULANT_VALIDATION_PREREG.md`
(commit 5d0b1986, section 10).

This document is committed **before any comparison**, device-level Gate 0
included. Nothing below may be changed after a comparison has been run.
Deviations are listed in the results document.

## 1. Purpose and what is independent

The aim is to reproduce the internal eigen-verdicts of the 16 subsets of
H4 = {30, 33, 35, 37} at P4 = (0.03625, 1.425, 1.5, 1) and at
G_S = (0.25, 1.425, 1.5, 1). The reproduction uses a separate tool (ANDES 2.0.0),
with the **same device equations** implemented as custom ANDES models.

| layer | internal (tx3-analysis) | reproduction (xtool-andes-gfl) | shared? |
|---|---|---|---|
| network data | `configs/ias2026/ieee39_network.json` (derived from `data/raw/ieee39_full.xlsx`, sha256 9c2048dc…) | the ANDES case `ieee39_full.xlsx`, read by ANDES' own parser (same sha256) | the source file |
| network equations | rectangular current balance `Y v = I(x, v)` | ANDES polar power balance (`Line`, `Shunt`, `PQ`, `Bus`) | no |
| power flow | internal Newton | ANDES `PFlow` | no |
| device equations | `ibr_cycles.models.ieee39_devices` | custom ANDES models `SG2AX`, `GFL11` (sympy strings, section 4) | the equations, by design |
| device parameters | derived in `ieee39_case` | a parameter handoff JSON exported once from the internal side (section 6) | the values |
| initialization | explicit, `initialize()` | ANDES `ConstService` chain from the `PFlow` solution | no |
| Jacobian | central differences (h / 2h) | sympy symbolic → lambdified | no |
| DAE reduction / eigensolver | `reduce_index_one`, numpy | ANDES `EIG` (`As = fx − fy gy⁻¹ gx`) | no |
| symmetry handling | transverse quotient `Z` (R_x, w) | drop the two eigenvalues of smallest modulus (§8.2) | no |
| classifier / hypergraph | FC01 | re-implemented in the reproduction script (§8.3) | no |

The reproduction therefore tests:
- the network model and data parsing;
- the power flow;
- initialization;
- linearization and reduction;
- the eigen pipeline;
- the transcription of the device equations (Gate 0).

It does not test whether those device equations are physically right. That is
the separate question of model adequacy.

## 2. Environment

- **Interpreter.** `.venv/xtool-andes-gfl`, created with `python -m venv` from
  the same base CPython 3.13.14 as `tx3-andes`. The dependencies are pinned to
  the tx3-andes versions: numpy 2.5.2, scipy 1.18.1, sympy 1.14.0, pandas 3.0.5,
  kvxopt 1.3.3.1, matplotlib 3.11.1, openpyxl 3.1.5, xlsxwriter 3.2.9,
  dill 0.4.1, pathos 0.3.5, tqdm 4.70.0, ipywidgets 8.1.9, pyyaml 6.0.3,
  chardet 7.6.0, texttable 1.7.0.
- **ANDES copy.** `git -C vendor/andes archive v2.0.0` (commit eda5163c) is
  extracted to `.venv/xtool-andes-gfl/src/andes`. It is installed editable with
  `SETUPTOOLS_SCM_PRETEND_VERSION=2.0.0`, `--no-deps --no-build-isolation`.
- **Not installed.** The internal package `ibr_cycles`. The reproduction script
  asserts that `importlib.util.find_spec("ibr_cycles") is None`.
- **Generated-code isolation.** Every ANDES process in this env runs with
  `USERPROFILE=.venv/xtool-andes-gfl/home` and an explicit
  `pycode_path=.venv/xtool-andes-gfl/home/pycode`, so `~/.andes` is never read
  or written. ANDES reuses a `pycode` module already in `sys.modules`, so the
  two installs are never loaded in one process.
- **Untouched.** `.venv/tx3-andes`, `vendor/andes` and `~/.andes`. This is
  verified before and after the campaign by `git -C vendor/andes status`, and by
  a sha256 listing of `~/.andes/pycode`.
- **Setup is scripted.** The setup script is
  `experiments/post_cumulant_validation/gfl_repro/setup_xtool_env.py`. The
  model file `pcv_models.py` is copied into the ANDES copy as
  `andes/models/pcv_models.py`. It is registered by adding
  `('pcv_models', ['SG2AX', 'GFL11'])` to `file_classes` right after
  `('static', …)`. This is the only change to the copy.

## 3. Network

- **Case data.** From the ANDES case `ieee39_full.xlsx`, only the sheets `Bus`,
  `PQ`, `PV`, `Slack`, `Shunt`, `Line` and `Area` are kept (when present).
  Every dynamic sheet is removed: `GENROU`, `IEEEX1`, `IEEEST`, `TGOV1N`,
  `Toggler`, and any other.
- **Lines.** ANDES `Line`: the π model with `tap`, `phi`, series `r + jx`, and
  shunt `(g + jb)/2` at each end. This equals the internal Ybus assembly:
  - `y_ff = (y_s + y_c)/|m|²`, `y_tt = y_s + y_c`;
  - `y_ft = −y_s / conj(m)`, `y_tf = −y_s / m`;
  - where `y_s = 1/(r + jx)`, `y_c = (g + jb)/2` and `m = tap·e^{jφ}`.
- **Shunts.** ANDES `Shunt`, equal to the internal `y_ii += g + jb`.
- **Loads.** `PQ` as constant power in the power flow and in TDS/EIG. Set
  `PQ.config.pq2z = 0`, `p2p = q2q = 1`, and `p2z = q2z = p2i = q2i = 0`.
  Internally the load current is `conj(S_L)/conj(v)`.
- **Solution settings.** `PFlow.config.tol = 1e-12`. PV Q limits are not
  enforced (`pv2pq = 0`), as in the internal flow. The system base is 100 MVA
  and the frequency 60 Hz.
- **Replacement.** Every PV / Slack static generator is replaced by exactly one
  dynamic device through `gen` (`replaces=True`):
  - `GFL11` at buses in S;
  - `SG2AX` at every other generator bus.

## 4. Custom device models (exact transcription)

Notation:
- The bus voltage is `V = v e^{ja}` (ANDES polar).
- Each device carries its own base. `w` is its rating weight on the 100 MVA
  system base, `w = Sn/100` (full replacement, IEEE-39 rule).
- The power it withdraws from the bus enters the ANDES bus equations as
  `a: −w·P`, `v: −w·Q` (ANDES sign convention: bus equations sum withdrawals).
- `ω_B = 2π·60`.
- No ANDES per-unit conversion flags (`z`, `power`) are set on any parameter.
  All values are the handoff values, on the device base.

### 4.1 `SG2AX` — two-axis machine + first-order AVR + washout-lag PSS (7 states)

Internal source: `SynchronousMachine.derivatives / injection / initialize`
(`src/ibr_cycles/models/ieee39_devices.py`), with `avr_tb = 0`, no services, and
blends equal to 1.

Terminal voltage in machine axes:
- internal `_rotate_to_dq`: `v_d = Re(V) sin δ − Im(V) cos δ`,
  `v_q = Re(V) cos δ + Im(V) sin δ`;
- polar: `v_d = v sin(δ − a)`, `v_q = v cos(δ − a)`.

| quantity | expression |
|---|---|
| `det` | `ra² + xd1·xq1` |
| `i_d` | `(ra (ed1 − v_d) + xq1 (eq1 − v_q)) / det` |
| `i_q` | `(−xd1 (ed1 − v_d) + ra (eq1 − v_q)) / det` |
| `Pe` (air gap) | `v_d i_d + v_q i_q + ra (i_d² + i_q²)` |
| `P`, `Q` at the terminal | `v_d i_d + v_q i_q`, `v_q i_d − v_d i_q` |
| `δ'` | `ω_B (ω − 1)` |
| `ω'` | `(pm − Pe − D (ω − 1)) / M` |
| `eq1'` | `(efd − eq1 − (xd − xd1) i_d) / Td10` |
| `ed1'` | `(−ed1 + (xq − xq1) i_q) / Tq10` |
| `efd'` | `(KA (vref + KS·xl − v) − efd) / TE` |
| `xw'` | `(Pe − xw) / T6` |
| `xl'` | `(T5 (Pe − xw)/T6 − xl) / T4` |

Here `xw` is `pss_w` and `xl` is `pss_l` (the PSS lag state, not the leakage
reactance). The terminal quantities follow from
`V conj(I) = (v_q − j v_d)(i_q + j i_d)`, where the network phasor of a dq pair is
`(x_q − j x_d) e^{jδ}` (internal `_rotate_to_network`).

The parameter values are:
- `KA = KA_case·k`, `TE = TE_case·t` (the policy);
- `KS`, `T4`, `T5`, `T6` from the case PSS rows;
- `M`, `D`, `ra`, `xd`, `xq`, `xd1`, `xq1`, `Td10`, `Tq10` from the case
  machine rows, on the machine base.

**Initialization** (identical formulas, as ANDES `ConstService`s):
1. `I = conj((p0 + jq0)/w / V)`;
2. `E = V + (ra + j xq) I`, `δ0 = arg E`;
3. `i_d, i_q = Re, Im(I e^{j(π/2 − δ0)})`, and likewise `v_d, v_q`;
4. `eq1 = v_q + ra i_q + xd1 i_d`, `ed1 = v_d + ra i_d − xq1 i_q`;
5. `efd = eq1 + (xd − xd1) i_d`;
6. `pm = Pe`, `vref = v + efd/KA`;
7. `ω = 1`, `xw = Pe`, `xl = 0`.

### 4.2 `GFL11` — grid-following converter with Q/V regulator and leak (11 states)

Internal source: `GridFollowingConverter.derivatives / injection / initialize`,
with `voltage_control = True`, `inertia_emulation = False`, the policy
`voltage_gain = g`, and `voltage_leak = L = 0.05`.

PLL-frame voltage:
- internal: `V e^{−jθ}`;
- polar: `v_d = v cos(a − θ)`, `v_q = v sin(a − θ)`.

| quantity | expression |
|---|---|
| `P`, `Q` | `v_d i_d + v_q i_q`, `v_q i_d − v_d i_q` |
| `err` | `g (v_ref − v)` |
| `q_cmd` | `kp_v·err + x_v` |
| `id_ref` | `kp_p (p_ref − p_f) + x_p` |
| `iq_ref` | `−(kp_q (q_cmd − q_f) + x_q)` |
| `e_d` | `v_d + kp_i (id_ref − i_d) + x_id − xf i_q` |
| `e_q` | `v_q + kp_i (iq_ref − i_q) + x_iq + xf i_d` |
| `θ'` | `kp_pll v_q + x_pll` |
| `x_pll'` | `ki_pll v_q` |
| `p_f'`, `q_f'` | `(P − p_f)/τ_p`, `(Q − q_f)/τ_p` |
| `x_p'`, `x_q'` | `ki_p (p_ref − p_f)`, `ki_q (q_cmd − q_f)` |
| `i_d'` | `(ω_B/xf)(e_d − v_d − rf i_d + xf i_q)` |
| `i_q'` | `(ω_B/xf)(e_q − v_q − rf i_q − xf i_d)` |
| `x_id'`, `x_iq'` | `ki_i (id_ref − i_d)`, `ki_i (iq_ref − i_q)` |
| `x_v'` | `ki_v·err − L (x_v − q_ref)` |

The injection is `w (i_d + j i_q) e^{jθ}` (system base), so the withdrawal terms
are `−w·P` and `−w·Q`.

The parameters are the internal `ConverterParameters` defaults:
- `kp_pll 53`, `ki_pll 1400`, `τ_p 0.03`;
- `kp_p 0.20`, `ki_p 8`, `kp_q 0.20`, `ki_q 8`;
- `kp_i 0.25`, `ki_i 6`;
- `xf 0.15`, `rf 0.01`, `kp_v 2`, `ki_v 20`.

Here `L` is the frozen leak (0.05 rad/s) and `g` is the policy.

**Initialization:**
1. `I = conj((p0 + jq0)/w / V)`, `θ0 = a`;
2. `i_d + j i_q = I e^{−jθ0}`;
3. `p_ref = v i_d`, `q_ref = −v i_q`, `v_ref = v`;
4. the states are
   `(θ0, 0, p_ref, q_ref, i_d, −i_q, i_d, i_q, rf i_d, rf i_q, q_ref)`.

### 4.3 Declared deviation from the preregistration text

§10 anticipated a custom GFL and "the PSS if custom". The ANDES library has no
fourth-order two-axis machine: the registered machines are `GENCLS`, the
sixth-order `GENROU` and `PLBVFU1`. The degraded-GENROU route of F1 is a
singular perturbation, not an exact transcription. The first-order AVR would
also require `SEXS` with a decoupled lead-lag state.

Exact equivalence therefore requires a custom machine. `SG2AX` includes the
AVR and PSS, so the "custom PSS" case of §10 applies to the PSS. No library
converter model is used anywhere.

## 5. Policy mapping

| coordinate | internal | ANDES |
|---|---|---|
| g | `ConverterParameters.voltage_gain` (scales both kp_v and ki_v through `err`) | `GFL11.g` |
| k | `machine_scaling["ka"]` (multiplies the case KA) | `SG2AX.KA = KA_case·k` |
| t | `machine_scaling["ta"]` (multiplies the case TE) | `SG2AX.TE = TE_case·t` |
| h | h = 1 at both points (no per-bus TE scaling) | — |

## 6. Parameter handoff

`gfl_repro/export_internal.py` runs once in tx3-analysis and writes
`results/PCV/PCV06/PCV06_handoff.json`. It contains:
- per generator bus: the `SG2AX` parameters, **before** the policy scaling
  (the case KA and TE), and `w`;
- the GFL defaults, `L`, and the two policy points;
- the Gate 0 vectors (§7);
- the internal reference results for the 32 cases (from
  `results/PCV/PCV02/PCV02_truth_core.csv`, plus the full internal spectra and
  equilibrium voltages, recomputed with the same direct path).

The ANDES-side script reads the network from the xlsx and the device
parameters from the handoff. It never reads internal results before its own
verdicts are written. The comparison is a third step.

## 7. Gate 0 — device-level equivalence

- **Vectors.** Take 20 random states and terminal voltages per instance
  (numpy `default_rng(20260918)`, drawn in the order listed):
  - every SG instance: 10 buses, with the P4 machine parameters
    (k = 1.425, t = 1.5), which are identical at G_S;
  - every GFL instance: 4 buses, at g = 0.03625 and at g = 0.25.
- **Draws.**
  - Each state component is `x_i = x_eq,i + 0.05·max(1, |x_eq,i|)·ε_i`, with
    `ε_i ~ N(0, 1)`.
  - The SG equilibria come from the internal P4 base case (no replacement).
  - The GFL equilibria come from the internal P4 H4 case.
  - `v ~ U(0.9, 1.1)` and `a ~ U(−π, π)`.
  - The setpoints (`pm`, `vref`, `p_ref`, `q_ref`, `v_ref`) are held at
    their equilibrium values.
- **Outputs compared.** All state derivatives, and the withdrawn powers `w·P`,
  `w·Q` (system base).
- **ANDES evaluation.** The ANDES values come from the lambdified `f` and the
  bus-equation contributions of the custom models, evaluated at those inputs
  (`model.f_update` / `g_update` after writing the variable values).
- **Pass.** `max |f_ANDES − f_internal| / max(1, max |f_internal|) ≤ 1e-9`, per
  device and output. Failure stops Phase 8 (S4).

## 8. Validation set and comparison

### 8.1 Cases

The 16 subsets of H4 at P4 and at G_S, 32 cases. P_inf is not run (optional in
§10).

### 8.2 ANDES spectrum

- `PFlow.run()` → `TDS.init()` → `EIG.run()`, giving `EIG.mu` (all
  eigenvalues of the reduced `As`).
- **Initialization residual.** `max |f|, |g|` at the initialized point.
- **Structural pair.** Sort by modulus, and require:
  - the two smallest to have `|λ| < 1e-3`;
  - every other eigenvalue to have `|λ| ≥ 1e-2`.

  Those two are dropped; the rest is the ANDES transverse spectrum.
- **State count.** The retained dimension must equal the internal
  `dim A_perp = n_x − 2`.

### 8.3 Per-case quantities (both sides)

- `α_perp`: the maximum real part of the transverse spectrum.
- `rhp`: the count with `Re > 0`.
- **Critical mode.** The rightmost transverse eigenvalue with
  0.3 ≤ f ≤ 1.5 Hz (Im > 0); its frequency.
- **Verdict.**
  - If `min |Re λ| ≥ 2e-3` over the transverse spectrum: `UNSTABLE` if
    `rhp > 0`, else `STABLE`.
  - Otherwise `BOUNDARY_OR_UNRESOLVED`. This is the sign-count branch of the
    FC01 logic; the internal side keeps its own FC01 classifier labels.
- **H, κ.** The minimal unstable subsets (FC01 `h_of` logic, re-implemented).

### 8.4 Primary tolerances (§10)

- equilibrium: the ANDES initialization residual ≤ 1e-6, and
  `max |V_ANDES − V_internal| ≤ 1e-6` over all buses;
- `|Δα_perp| ≤ 1e-3 s⁻¹`;
- `|Δf_crit| ≤ 1e-3 Hz`;
- RHP counts equal in 32/32;
- verdicts equal in 32/32;
- H and κ equal at P4 and at G_S.

### 8.5 Secondary tolerance

For every internal band eigenvalue (0.3–1.5 Hz, Im > 0), the nearest ANDES
eigenvalue lies within 1e-4.

## 9. Stopping and mismatch diagnosis

- Any verdict or H mismatch, or `|Δα| > 1e-2`, stops the stronger claims.
- Nothing is retuned. The first mismatch is localized in this order:
  1. the power flow (bus voltages);
  2. the initialized device states;
  3. Gate 0 at the equilibrium itself;
  4. the full Jacobian blocks (`fx`, `fy`, `gx`, `gy` versus internal central
     differences, after the variable mapping);
  5. the reduced spectrum.

## 10. Outputs

- `results/PCV/PCV06/`: the handoff, Gate 0 table, ANDES per-case results, spectra
  and comparison.
- `results/20260911_GFL_REPRODUCTION.csv`: one row per case, with internal /
  ANDES α, frequency, rhp, verdict, differences and pass flags.
- `docs/20260911_GFL_REPRODUCTION.md`: environment evidence, Gate 0, tolerances,
  deviations and the conclusion.
