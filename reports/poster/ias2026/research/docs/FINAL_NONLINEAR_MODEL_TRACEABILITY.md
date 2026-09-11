# Final nonlinear model traceability (the equations the repository runs)

No new benchmark is introduced. This is the frozen L0 phasor DAE: the code
that produced F7–F12, G1–G3, Phase I and this campaign. The line-level mapping
of all 33 equations to code and tests is in
`docs/nonlinear_portfolio_v2/EQUATION_CODE_TRACEABILITY.csv`, and the extended
derivation in `docs/nonlinear_portfolio_v2/NONLINEAR_POWER_SYSTEM_MODEL.md`.
Paths below are relative to `src/ibr_cycles/`.

## 0. Common form

    M(x, z, S, q) x' = F(x, z, S, q; theta, xi)
    0                = G(x, z, S, q; theta, xi)

| symbol | meaning in this repository | code |
|---|---|---|
| `x` | device states, stacked slot by slot (machine 7, or 8 with SEXS; GFL 10 + `x_v` [+ `y_vi`]) | `models/ieee39_case.py:144` (`Ieee39Dae`, labels) |
| `z` | bus voltages in interleaved real form `(vx_1, vy_1, …)` | `Ieee39Dae.voltages`, `:160` |
| `S` | portfolio: the replaced buses, each with fraction `rho_b` (`rho = 1` removes the machine) | `ReplacementPlan`, `:53`; `build_dae`, `:352` |
| `theta` | policy and operating point: `g` (leaky Q/V gain), `k` (AVR gain scale), `t` (AVR `TE` scale), `h` (TE heterogeneity); matched dispatch | `experiments/_f7_common.py:90–140` |
| `xi` | uncertain physical parameters. Not varied in this campaign; E35 and E37 varied operating point and machine data | — |
| `q` | controller / limiter regime. **None**: L0 has no limiter, so `q` is empty and `active_set = ABSENT` | `nonlinear/model.py` |
| `M` | identity: every device returns the normalized derivative | `nonlinear/model.py` (`mass`) |
| `F` | `Ieee39Dae.f`, device by device | `models/ieee39_case.py:163` |
| `G` | `Ieee39Dae.g`, KCL | `models/ieee39_case.py:175` |

Theorems are applied to IEEE-39 only where their assumptions are traced to
these equations (see §H).

## A. AC network

| item | equation | code |
|---|---|---|
| Ybus | branch pi model; `y = 1/(r + jx)`; complex tap on the from side; `Yff = (y + jb/2)/abs(t)^2`, `Yft = −y/conj(t)`, `Ytf = −y/t`, `Ytt = y + jb/2`; `Ybus = Cf^T Yf + Ct^T Yt + diag(y_sh)`; no symmetry imposed | L0 `models/ieee39_network.py:95`; L1 `nonlinear/network.py` (tested equal to 1e-12) |
| shunts | `diag(y_sh)` from the source shunt table | `ieee39_network.py:95` |
| current balance | `0 = realify(Ybus V − sum_k I_k(x_k, V) + I_load(V))`, injections positive | `ieee39_case.py:175–197` |
| complex power | `S_k = V conj(I_k)`; branch loss `Re(S_f + S_t) = r abs(y (V_f/t − V_t))^2` | `nonlinear/network.py` (`branch_flows`) |
| power flow | Newton AC power flow, PV/slack from the source; `tol = 1e-12` | `ieee39_network.py:187` |

## B. Synchronous machine: two-axis, 4th order plus controls

Machine base `Sn`; weight `gamma = (1 − rho) Sn/100` on the injection only
(`ieee39_devices.py:154–278`):

    vd = vx sin(delta) − vy cos(delta),  vq = vx cos(delta) + vy sin(delta)      (:69)
    id = (ra (e'd − vd) + x'q (e'q − vq))/det,  iq = (−x'd (e'd − vd) + ra (e'q − vq))/det,  det = ra^2 + x'd x'q   (:176)
    Pe = vd id + vq iq + ra (id^2 + iq^2)                                         (:192)
    delta' = omega_B (omega − 1)
    omega' = (Pm − Pe − D (omega − 1)) / M        (power form; M = 2H on the machine base)   (:202–203)
    e'q' = (efd − e'q − (xd − x'd) id)/T'd0,   e'd' = (−e'd + (xq − x'q) iq)/T'q0
    I = gamma (id + j iq) rotated back to the network frame                       (:236)

- **Order used:** two-axis (transient EMFs `e'q`, `e'd`).
- **Damping:** `D = 0` for every machine, as in the source (GENROU `D = 0`).
- **Governor:** none; `Pm` is constant. The documented TGOV1N appears only in the
  separately versioned `models/governed.py` (FC03).
- **IEEE-68:** the sub-transient model, eqs. 1–6 plus a dummy coil, in torque form
  (`models/ieee68_devices.py:163`).

## C. Exciter / AVR

- **IEEE-39, harmonized first order** (`ieee39_devices.py:217–226`):

      efd' = (KA (vref + vs − abs(V)) − efd)/TE,   KA = 10.1–40, TE = 0.25–0.95 (per machine)

  This merges the IEEEX1 regulator gain `KA` with its exciter time constant
  `TE`. It drops `KE`, the rate feedback `KF1/TF1` and saturation.
  Status: `MODEL_CHANGED` relative to the ANDES source, declared since F1. The
  policy coordinates `k`, `t`, `h` scale `KA` and `TE`.
- **Kundur, SEXS** (`:196–211`): lead-lag `TATB`, `TB` followed by
  `KA/(1 + s TE)`, as in the source, with the limits absent.
- **IEEE-68:** DC4B / ST1A / manual, as documented (`ieee68_devices.py:11–14`).
- **Units:** per unit on the machine base; `vref` is set at initialization.

## D. PSS

- **IEEE-39:** IEEEST mode 3, power input. With `T1 = T2 = T3 = 0` in the
  source, what remains is the washout `T5/T6` and one lag `T4`, with
  `vs = KS x6`, `KS = −2`. Code: `ieee39_devices.py:194–195, 207–208`. The
  documented output limit is ±0.1 pu (LSMAX/LSMIN); it is **absent** in L0 and
  used as a scope guard in FC05.
- **Kundur:** none. **IEEE-68:** speed-input washout plus three lead-lags.

## E. GFL converter (per unit on the converter rating; `ieee39_devices.py:342–445`)

    dq transform     v_dq = V e^{−j theta}                        (:377)
    PLL              theta' = kp_pll vq + x_pll,   x_pll' = ki_pll vq    (SRF, unnormalized)
    P/Q filters      p_f' = (P − p_f)/tau_p,  q_f' = (Q − q_f)/tau_p;  P = vd id + vq iq,  Q = vq id − vd iq
    Q/V (leaky)      q_cmd = kp_v g e + x_v,  e = v_ref − abs(V),  x_v' = ki_v g e − w (x_v − q_ref),  w = 0.05
    outer P, Q       id_ref = kp_p (p_ref − p_f) + x_p;  iq_ref = −(kp_q (q_cmd − q_f) + x_q)
    current control  e_d = vd + kp_i (id_ref − id) + x_id − xf iq;  e_q = vq + kp_i (iq_ref − iq) + x_iq + xf id
    filter           id' = (omega_B/xf)(e_d − vd − rf id + xf iq),  iq' = (omega_B/xf)(e_q − vq − rf iq − xf id)
    current equation injection = gamma (id + j iq) e^{j theta}      (:421)

Assumptions:

- `P_ref` equals the displaced dispatch, and the DC side is ideal and unlimited.
- The feedforward and decoupling are exact. This gives the algebraic identity
  N.34: the current loop sees the grid only through its references and the PLL.
- The rating rule is `Sn` for IEEE-39 and Kundur, and `abs(S_gen)/0.8` for
  IEEE-68.

## F. Loads

| benchmark | load model | code |
|---|---|---|
| IEEE-39, Kundur | constant power, `I_load = conj(S)/conj(V)` | `ieee39_case.py:196` |
| IEEE-68 | constant impedance at the power-flow voltage | `:191` |

ZIP and frequency-dependent loads are absent.

## G. Omitted physics (absent from L0; the consequence for every claim)

| mechanism | status in L0 | consequence |
|---|---|---|
| governor / primary frequency control | absent (source documents TGOV1N; separate governed model FC03) | frozen results are **transverse**; the common frequency is marginal |
| machine damping | `D = 0` (as documented) | exact Jordan partner `A w = omega_B R_x` |
| GFL current saturation and priority | absent | TDS guard `abs(I) <= 1.2 pu` |
| anti-windup | absent | outer PI integrators unbounded |
| DC link, PV power budget | absent (`P_ref` = dispatch, unlimited DC) | D3 restricted to downward dips |
| LVRT / HVRT modes | absent | TDS guard `abs(V_gfl)` in [0.9, 1.1] pu |
| protection (under/over-frequency, out-of-step, UFLS) | absent | TDS guards on speed, PLL and COI frequency |
| exciter limits (IEEEX1 VRMAX/VRMIN) | absent | TDS guard, documented values |
| PSS output limit (LSMAX/LSMIN) | absent | TDS guard, ±0.1 pu |
| switching, PWM, plant-level delay | absent | not an EMT model; no claim above a few Hz |
| line dynamics (RLC network) | absent | quasi-static phasors |
| AGC (ACEc in the source) | absent | no secondary control in either version |

## H. Theorem-to-equation assumptions for IEEE-39

| result used | assumption | traced to |
|---|---|---|
| T1 / transverse quotient | rotation invariance; `D = 0`; constant `Pm`; frequency-independent network | §A (covariant KCL), §B, §E (PLL angle), §G |
| principal-minor factorization | common operating point; exact locality | matched dispatch, BC02 C2 realization |
| TDS labels | guards stand for the absent limiters of §G | `nonlinear/tds.py`; frozen config |
| second-order curvature (FC10) | smooth right-hand side in a fixed active set | L0 has no limiter, so the whole state space is one active set |
