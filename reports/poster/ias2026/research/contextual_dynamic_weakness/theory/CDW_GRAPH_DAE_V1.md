# CDW graph-DAE model V1 (notation and exact model)

This document defines the **exact current model**: the frozen TX4 IEEE-39
phasor DAE, reused unchanged. It keeps that model separate from the **future
generalized model**.

## 1. Sets and variables

- **Buses.** `N = {1, …, 39}` in the frozen order of
  `configs/ias2026/ieee39_network.json`.
- **Branches.** `E = {0, …, 45}`, all in service.
  - Indices 34–45 are transformers with off-nominal taps on the from-bus:
    - 2–30, tap 1.025;
    - 31–6, tap 0.9;
    - 10–32, tap 1.07;
    - 12–11, tap 1.006;
    - 12–13, tap 1.006;
    - 19–20, tap 1.06;
    - 19–33, tap 1.07;
    - 20–34, tap 1.009;
    - 22–35, tap 1.025;
    - 23–36, tap 1.0;
    - 25–37, tap 1.025;
    - 29–38, tap 1.025.
  - Lines carry charging `b`.
- **Shunts.** Bus 4 (b = 1.0 p.u.) and bus 5 (b = 2.0 p.u.).
- **Node classes.**
  - SG: 10 two-axis machines (buses 30–39). Bus 39 is the slack and represents
    the interconnection; it is never replaced.
  - GFL: the frozen grid-following converter. It replaces an SG at a candidate
    bus with matched dispatch and equal rating.
  - Loads: 19 constant-power loads.
  - GFM: **future model only**. None is instantiated in CDW V1.
- **Candidates.** The census is `V9 = {30, …, 38}` and the core is
  `V4 = {30, 33, 35, 37}`.

**Variables.**
- `δ ∈ {0,1}^9`: composition.
- `θ = (g, k, t, h)` plus per-device parameters: control.
  - `g` is the Q/V gain of every GFL.
  - `k` and `t` scale every SG's `K_A` and `T_E`.
  - `h` is the `T_E` heterogeneity exponent.
  - Per-device parameters: `g_i`, the PLL scale `c_i`, `K_{A,i}`.
- `ρ ∈ R^46_{>0}`: branch strengths; `ρ_e = 1` is the frozen network.
- `π`: operating point. It collects the scheduled `P_i` and `V_set,i` at
  generator buses, the loads `S_j`, and the slack `(V_s, θ_s)`.

## 2. Exact branch-terminal network

`Y` is **not** of the form `B diag(y) Bᵀ`, because taps and shunts invalidate
that. The exact construction is

    Y(ρ) = Y_sh + Σ_e ρ_e C_eᵀ Y_e C_e,

- `C_e ∈ {0,1}^{2×39}` selects the from- and to-bus of branch `e`.
- The branch two-port (π-model, tap `m_e = t_e e^{jφ_e}` at the from-bus, all
  `φ_e = 0` here) is

      Y_e = [ (y_e + j b_e/2)/|m_e|²      −y_e/conj(m_e) ]
            [ −y_e/m_e                     y_e + j b_e/2  ],   y_e = 1/(r_e + j x_e).

- `ρ_e` multiplies the **whole** two-port (series and charging) with the tap
  unchanged. This is the TX4 PCV04 convention: `ρ_e = 2` is a parallel identical
  circuit and `ρ_e → 0` is an outage.
- A code identity check requires `Y(1) = load_network().ybus` exactly.

**Susceptance decomposition** (used only in graph-mode studies). With
`b_e^s = Im(−y_e)`,

    L_B = Σ_e ρ_e (b_e^s/t_e) (e_f − e_t)(e_f − e_t)ᵀ,
    Im Y = −L_B − Δ_tap + B_ch + B_sh,
    Δ_tap = Σ_e ρ_e b_e^s [ (1/t_e² − 1/t_e) e_f e_fᵀ + (1 − 1/t_e) e_t e_tᵀ ].

`L_B` is a true weighted Laplacian (with a positive weight when `x_e > 0`). The
residual `Δ_tap + B_ch + B_sh` is diagonal and exact.

## 3. Device equations

These are the frozen TX4 device equations (main.tex eqs. 3–6).
- **SG:** two-axis, `D = 0`, constant `P_m`, first-order AVR
  `T_E e_fd' = K_A (v_ref + v_s − |V|) − e_fd`, and a power-input PSS.
- **GFL:**
  - SRF-PLL `(k_p,pll, k_i,pll)`;
  - power filters with time constant `τ_p`;
  - PI outer loops and PI inner current loops behind `(x_f, r_f)`;
  - leaky Q/V regulator
    `q_cmd = k_p,v g e + x_v`, `x_v' = k_i,v g e − w (x_v − q_ref)`, with
    `e = v_ref − |V|` and `w = 0.05 rad/s`.
- **Loads:** constant power, `I = conj(S/V)`.
- **Network:** `0 = Y z_c − Σ_devices I_dev + Σ_loads I_load`, written in
  rectangular coordinates.

## 4. Equilibrium semantics (CDW_THEORY C4)

- **SPR** (TX4 `solve_case`):
  - the AC power flow is solved with generator buses PV `(P_i, V_set,i)` and
    the slack `(V_s, θ_s)`;
  - every device is initialized at that point (`v_ref`, `p_m` for SGs;
    `p_ref`, `q_ref`, `v_ref` for GFLs).
- **RP:** device references are fixed at their base values. The unknowns are
  `(x, z, p_m,slack)` and the equations are `f = 0`, `g = 0` and
  `∠V_slack = θ_s`.

## 5. Linearization and transverse spectrum

- **Jacobians:** central differences as in TX4 `central_difference_jacobians`,
  then `A = f_x − f_z g_z^{-1} g_x`.
- **Transverse operator:** `A_⊥ = Zᵀ A Z` via TX4
  `rotation_generator`, `frequency_partner` and `transverse_operator`. The
  status uses the TX4 classifier (h/2h error).
- **Eigenvector derivatives:** use the full `A`, whose non-structural spectrum
  equals `σ(A_⊥)`, with left and right eigenvectors from `scipy.linalg.eig`.
- **Bus-voltage mode shape:** `φ = −g_z^{-1} g_x u ∈ C^{78}`, a coordinate
  common to all portfolios.

## 6. Future generalized model (NOT used in CDW V1)

- GFM devices; unmatched or re-dispatched replacement (for this, the PF depends
  on `δ`);
- converter current limits;
- governors in the census;
- frequency-dependent network;
- distribution feeders (for Africano, if it ever becomes available).

Each of these changes the equilibrium and the symmetry assumptions of A2.
