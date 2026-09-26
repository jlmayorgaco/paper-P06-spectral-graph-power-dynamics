# Nonlinear power-system model: L0/L1 equations that actually run

**Scope.** This document writes out the equations that the frozen benchmarks
execute today (level L0) and their matrix rewrite (level L1).

- L1 is tested equal to L0 (`tests/test_nl01_model.py`).
- No fidelity is added. Extensions are L2 (limits, delays, DC link,
  governor) and L3 (dynamic RLC network). Each is a separate, versioned model,
  and none exists yet.
- The v2 note (`TEORIA_NO_LINEAL_MATRICIAL_Y_CERTIFICACION.pdf`, sha256
  `713fe27b…`) is used as a reconciliation specification, not as a source of
  equations. Where it differs from L0, the difference is listed in §9 and L0
  is kept.

Model version: `L1-v2.0` (`src/ibr_cycles/nonlinear/__init__.py`).

## 1. Coordinates, units and bases

- Time in seconds. `omega_B = 2 pi 60` rad/s (`ieee39_devices.OMEGA_B`).
- System base `S_base = 100` MVA. Voltages are in pu of the bus base.
- Bus voltages are algebraic unknowns in **interleaved** real form,
  `z = [vx_1, vy_1, vx_2, vy_2, …]` (`Ieee39Dae.voltages`).
  - The stacked form `[Re V; Im V]` and the permutation between the two forms
    are in `nonlinear/network.py` (`interleave_permutation`, `realify`).
  - The two orders are never mixed.
- Injected current is positive. Loads are subtracted explicitly.
- Every device is written on its **own** base.
  - The weight `gamma = S_device / S_base` multiplies only its network
    injection.
  - Consequences: `P_sys = gamma P_local`, `I_sys = gamma I_local`,
    `X_sys = X_local / gamma`, `M_sys = gamma M_local`.
  - Tested behaviourally in NL00-A: the system-base object gives the same
    current and derivatives to about `1e-15`.
- Mass matrix: **M = I**. Every device returns the normalized derivative.
  Rotor inertia `M = 2H` sits inside the speed equation.

## 2. Network (KCL)

Branch model (`nonlinear/network.py`, L1; `ieee39_network._build_ybus`, L0):

    Cf, Ct                          n_branch x n_bus selectors (from, to)
    y = 1/(r + jx),  b = total charging,  t = |t| e^{j phi} on the FROM side
    Yff = (y + jb/2)/|t|^2, Yft = -y/conj(t), Ytf = -y/t, Ytt = y + jb/2
    Yf = diag(Yff) Cf + diag(Yft) Ct,  Yt = diag(Ytf) Cf + diag(Ytt) Ct
    Ybus = Cf^T Yf + Ct^T Yt + diag(y_sh)

With a phase shifter `Ybus` is not symmetric, and no symmetry is imposed.
The branch loss identity is `Re(S_f + S_t) = r |y (V_f/t - V_t)|^2`. It is
tested on all 46 + 15 + 83 real branches (IEEE-39, Kundur, IEEE-68) and on 40 synthetic phase-shifting
branches.

KCL residual (`Ieee39Dae.g`):

    0 = g(x, z) = realify_interleaved( Ybus V - sum_devices gamma_k I_k(x_k, V_bus(k)) + I_load(V) )

Loads:

- IEEE-39 and Kundur: constant power, `I_load = conj(S_load)/conj(V)`.
- IEEE-68: constant impedance fixed at the power-flow voltage,
  `I_load = conj(S_load) V / |V_0|^2`.
- ZIP and frequency-dependent loads are ABSENT.

## 3. Synchronous machine: two-axis (IEEE-39, Kundur; `SynchronousMachine`)

The rotor frame, with `d` leading `q` by 90° behind the rotor angle
(`_rotate_to_dq`):

    vd = vx sin(delta) - vy cos(delta),     vq = vx cos(delta) + vy sin(delta)
    rd = e'd - vd,  rq = e'q - vq,  det = ra^2 + x'd x'q
    id = (ra rd + x'q rq)/det,              iq = (-x'd rd + ra rq)/det
    P  = vd id + vq iq + ra (id^2 + iq^2)            (air-gap power, machine base)

    delta' = omega_B (omega - 1)
    omega' = (Pm - P - D (omega - 1)) / M            (power form at nominal speed)
    e'q'   = (efd - e'q - (xd - x'd) id) / T'd0
    e'd'   = (-e'd + (xq - x'q) iq) / T'q0

The injection is `gamma (id + j iq)` rotated back:
`I = gamma [id sin(delta) + iq cos(delta), -id cos(delta) + iq sin(delta)]`.

Initialization (`initialize`):

- `delta = arg(V + (ra + j xq) I)`;
- `e'q`, `e'd` and `efd` from the stator and flux equations at steady state;
- `Pm = P`, `vref = |V| + efd/KA`.
- Checked by substitution: `f = g = 0` to about `1e-12`.

**Excitation, IEEE-39.** A harmonized first-order AVR,

    efd' = (KA (vref + vs - |V|) - efd) / TE

with `KA = 10.1` and `TE = 0.25`. This combines the **regulator gain** of
IEEEX1 with its **exciter time constant**. It is a declared harmonization and
not an exact reduction of IEEEX1:

- `KE = -0.05` (self-excited) and the rate feedback are dropped;
- see the frozen IEEE-39 docstring and N18.

The v2 note (§4.3) warns against exactly this kind of block merge. Status:
`MODEL_CHANGED` relative to the ANDES source, and declared.

**Excitation, Kundur.** SEXS as in the source:

    ll' = (u - ll) / TB,  u = vref + vs - |V|
    lead-lag output = TATB u + (1 - TATB) ll
    efd' = (KA (lead-lag output) - efd) / TE        (KA = 20, TE = 0.83, TATB = 0.4, TB = 5)

The source limits EMIN and EMAX are ABSENT.

**PSS, IEEE-39.** IEEEST, mode 3 (power input). The source has
`T1 = T2 = T3 = 0`, so only the washout and one lag remain:

    x5' = (P - x5)/T6,   w = T5 (P - x5)/T6,   x6' = (w - x6)/T4,   vs = KS x6

with `KS = -2`, `T5 = 1`, `T6 = 4.2` and `T4 = 0.75` for machine 1. The limit
`LSMAX` is ABSENT. Kundur has no PSS.

## 4. Synchronous machine: IEEE-68 (`Ieee68Machine`, Singh & Pal eqs. 1–10)

The frame is `V e^{-j delta} = Vq + j Vd`, with q real. There are four rotor
coils plus the dummy coil `E'dc` (`Tc = 0.01` s):

    (Iq + j Id) = [ (E'q kd1 + psi1d kd2 - Vq) + j (E'd kq1 - psi2q kq2 - Vd + E'dc) ] / (Ra + j X''d)
    kd1 = (X''d - Xl)/(X'd - Xl),  kd2 = (X'd - X''d)/(X'd - Xl)   (q analogous)
    Te = E'd Id kq1 + E'q Iq kd1 - Id Iq (X''d - X''q) + psi1d Iq kd2 - psi2q Id kq2

    delta' = omega_B S,   S' = (Tm - Te - D S)/(2H)                   (torque form, slip S)
    T'd0 E'q' = Efd - E'q + (Xd - X'd)[Id + (X'd - X''d)/(X'd - Xl)^2 (psi1d - (X'd - Xl) Id - E'q)]
    T'q0 E'd' = -E'd + (Xq - X'q)[-Iq + (X'q - X''q)/(X'q - Xl)^2 ((X'q - Xl) Iq - E'd - psi2q)]
    T''d0 psi1d' = E'q + (X'd - Xl) Id - psi1d
    T''q0 psi2q' = -E'd + (X'q - Xl) Iq - psi2q
    Tc E'dc' = Iq (X''d - X''q) - E'dc

The exciters and stabilizers:

- DC4B (PID, rate feedback, exponential saturation) and ST1A, as in eqs. 7–8;
- a manual field on G13–G16;
- the speed-input PSS (washout plus three lead-lags) of eq. 9.

All controller limits are ABSENT. The equilibrium is checked to lie inside
them, with the smallest margin 3.0. The model reproduces the report's Table 4
to `5e-4` Hz (G3).

## 5. Grid-following converter (`GridFollowingConverter`, 10 states + x_v [+ y_vi])

The PLL frame is `v_dq = V e^{-j theta}`. It is an SRF-PLL without voltage
normalization, the same form as v2 N.28.

    P = vd id + vq iq,   Q = vq id - vd iq      (dQ/di_q < 0: the sign is explicit)
    theta' = kp_pll vq + x_pll,       x_pll' = ki_pll vq
    p_f' = (P - p_f)/tau_p,           q_f' = (Q - q_f)/tau_p
    q_cmd = kp_v g e + x_v,  e = v_ref - |V|,   x_v' = ki_v g e - w (x_v - q_ref)   (leaky Q/V, w = 0.05)
    id_ref = kp_p (p_ref - p_f) + x_p,        x_p' = ki_p (p_ref - p_f)
    iq_ref = -(kp_q (q_cmd - q_f) + x_q),     x_q' = ki_q (q_cmd - q_f)
    e_d = vd + kp_i (id_ref - id) + x_id - xf iq,   e_q = vq + kp_i (iq_ref - iq) + x_iq + xf id
    id' = (omega_B/xf)(e_d - vd - rf id + xf iq),   iq' = (omega_B/xf)(e_q - vq - rf iq - xf id)
    x_id' = ki_i (id_ref - id),  x_iq' = ki_i (iq_ref - iq)
    injection: gamma (id + j iq) e^{j theta}

Defaults: `kp_pll = 53`, `ki_pll = 1400`, `tau_p = 0.03`, `kp/ki = 0.2/8` (P
and Q), `kp_i/ki_i = 0.25/6`, `xf = 0.15`, `rf = 0.01`, `kp_v/ki_v = 2/20`. All
are per unit on the converter rating.

**Identity N.34 holds in L0.** Substituting `e_d` and `e_q`, the bus voltage
and the rotational terms cancel:

    id' = (omega_B/xf)(kp_i (id_ref - id) + x_id - rf id)

The same holds for `iq'`. The current channel sees the grid only through
`id_ref`, `iq_ref` and the PLL. Status: T (algebraic identity of the code).
The numerical check is part of NL07.

ABSENT in L0:

- current and voltage limits, and anti-windup;
- modulation delay and feedforward filter;
- DC link and the PV power budget.

## 6. Partial replacement and portfolios

A plan assigns `rho_i`. The bus carries the machine at weight
`(1 - rho) Sn/100` and the converter at weight `rho * rating`, with shares of
`P` and `Q` per the policy (matched: proportional).

- IEEE-39 and Kundur: `rating = Sn/100`.
- IEEE-68: `rating = |S_gen| / 0.8`.
- At `rho = 1` the machine slot is not instantiated (`RATING_FLOOR = 1e-6`).
  The endpoint is a **removal**, not a continuous limit.
- A condenser keeps the machine with no active power, and its services
  follow the F8 switches.
- The service "blend = 0" switches create exactly zero rows (dead states).
  These are ledgered in BC00-B.

## 7. Equilibrium, symmetry and the zero ledger

- Initialization is the frozen AC power flow followed by device
  initialization and a DAE Newton solve (`solve_case`).
- The rotation generator `R = (R_x, R_z)` is derived from the frames:
  - `R_x = 1` on every `delta` and `theta_pll`;
  - `R_z = j V`.
- The nonlinear finite-angle test (NL00-C, `phi = 0.37` rad, off-equilibrium)
  gives:
  - `f(x + phi R_x, e^{j phi} V) = f` to `3e-13`;
  - KCL rotates covariantly, to `3e-13`.
- **Jordan partner.** The uniform frequency shift is the direction
  `w = (omega +1, x_pll +omega_B, y_vi +1, 68-bus PSS washout +K)`. It
  satisfies `A w = omega_B R_x` whenever every machine has `D = 0`, which is
  the case in all three benchmarks. After the quotient it is an exact zero of
  `A_q`: a physical **marginal** mode (no primary frequency control),
  reported and never deleted.
- No bus is held fixed in the dynamic model. The power-flow slack carries an
  ordinary dynamic machine (IEEE-39 bus 39 with `M = 100` on its base,
  Kundur machine 1, IEEE-68 machine 16). Declared: that machine is a finite
  inertia, not an infinite bus.

## 8. Software contract (`src/ibr_cycles/nonlinear/model.py`, `PhasorModel`)

| method | returns | L0 status |
|---|---|---|
| `residual_f(x, z, u, p)` | `F = f`, state names and units | wraps `Ieee39Dae.f` |
| `residual_g(x, z, u, p)` | KCL, interleaved, names `KCL_re_busN` | wraps `Ieee39Dae.g` |
| `mass(x, z, p)` | identity | declared (L0 normalized) |
| `outputs(x, z, u)` | `\|V\|` per bus; `P`, `Q` per device (system base) | computed from the same device code |
| `initialize_from_pf()` | frozen equilibrium `(x0, z0)` | the `solve_case` path |
| `symmetry_generators()` | `R_x`, `R_z`, partner `w` | `certification/symmetry.py` |
| `active_set()` | empty | ABSENT (no limiters in L0) |
| `guard_values()` | empty | ABSENT |
| `reset_map()` | `NOT_APPLICABLE` | no discrete events in L0 |

The inputs `u` are the load `dP` and `dQ` at every load bus and `dPm` at every
machine, in system pu. The derivative module (`nonlinear/manifold.py`) and the
port and certification modules call these same residuals. There is one
formulation for TDS, Jacobians and ports.

## 9. Differences from the v2 specification (kept as L0, not changed)

| v2 item | L0 | status |
|---|---|---|
| machine convention N.19 | d leads q by 90° behind delta; IEEE-68 uses q real, d imaginary | declared; transformation in §3–4 |
| mechanics N.20 (power) vs N.22 (torque) | two-axis: power form; IEEE-68: torque form | both valid at nominal speed; not interchangeable off nominal |
| excitation N.25, lead-lag plus lag | IEEE-39: first-order merge of KA and TE (harmonized); Kundur: SEXS exact without limits | `MODEL_CHANGED` for IEEE-39 relative to the source |
| PSS N.26, speed plus two lead-lags | IEEE-39: power input, washout plus lag (the source's nonzero blocks); IEEE-68: speed plus three stages | equal to the source structure, limits absent |
| governor N.27 | absent in all three | L2 extension |
| GFL limits, anti-windup N.36–37 | absent | L2 |
| DC link N.38–40 | absent (`P_ref` = dispatch, unlimited DC) | L2 |
| ZIP loads N.15 | PQ (IEEE-39, Kundur), Z (IEEE-68) | not changed |
| state-dependent mass | `M = I` | declared |
| v2 package code (`code/nonlinear_models.py`, `nonlinear_checks.py`, 625 checks) | **not supplied** with the PDF | the v2 checks could not be reproduced; only the v1 404 checks were |
