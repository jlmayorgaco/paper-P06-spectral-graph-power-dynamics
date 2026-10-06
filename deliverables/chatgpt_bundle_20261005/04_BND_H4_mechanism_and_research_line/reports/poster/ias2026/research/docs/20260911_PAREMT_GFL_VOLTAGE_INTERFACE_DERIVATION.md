# GFL EMT interface v2: average-value voltage source behind the Rf–Lf filter — derivation

Date: 2026-09-11. Written and committed **before any V2 computation**, together
with `docs/20260911_PAREMT_EMT_PREREG_V2.md`. It defines the single, primary GFL
EMT realization of preregistration V2. The analysis in §7 was done on paper to
choose the discretization; no V2 simulation had been run.

Notation follows `docs/20260911_GFL_REPRODUCTION_SPEC.md` §4.2 and
`docs/20260911_PAREMT_GFL_EQUATION_MAP.md`. ω0 = ω_B = 2π·60 rad/s.

## 1. Frames and transforms

- **Instantaneous phase quantities.** x_a, x_b, x_c, in amplitude per unit on the
  100 MVA system base (the ParaEMT convention; EMT01).
- **Complex αβ space vector (amplitude invariant).**

  `x_αβ = (2/3)(x_a + a x_b + a² x_c)`, with `a = e^{j2π/3}`.

  For a balanced three-wire set the inverse is
  `x_k = Re(x_αβ e^{−j2πk/3})`, k = 0, 1, 2. The zero sequence is excluded by
  construction.
- **Synchronous-frame phasor.** `X(t) = x_αβ(t) e^{−jω0 t}`. This is the
  quantity called V(t) or I(t) in prereg v1 §3, and it is constant in
  sinusoidal steady state at 60 Hz.
- **PLL frame.**
  - The frame angle is `φ(t) = ω0 t + θ(t)`, where θ is the TX4 PLL state (its
    angle relative to the synchronous frame).
  - Its rate is `φ̇ = ω0 + θ′`, with `θ′ = kp_pll v_q + x_pll` (TX4 PLL
    equation).
  - `x_dq = x_d + j x_q = x_αβ e^{−jφ} = X e^{−jθ}`.
- **J.** The real matrix [[0, −1], [1, 0]] acting on (x_d, x_q) is
  multiplication by j on x_dq: J i_dq ↔ j i_dq.

## 2. Bases

| quantity | device base (rating S_n) | system base (100 MVA) |
|---|---|---|
| voltage | 1 pu = rated amplitude | same (TX4 case: all bases kV = 1) |
| current | i_dq | `I_sys = w · i_dq e^{jθ}` (synchronous phasor), w = S_n/100 |
| filter resistance | rf = 0.01 | `Rf = rf / w` |
| filter reactance | xf = 0.15 | `Xf = xf / w` |
| filter inductance | xf/ω0 (pu·s) | `Lf = xf / (w ω0)` |

**Values for the V2 cases.** Here w is Sn/100.

| case | w | Rf (system) | Lf (system) |
|---|---|---|---|
| bus-30 regression | 10.4 | 9.615e-4 | 3.826e-5 |
| bus-35 blind | 10.857 | 9.211e-4 | 3.665e-5 |
| bus-37 blind | 9.702 | 1.031e-3 | 4.101e-5 |

In IEEE-39 each GFL uses its own bus rating.

## 3. Physical filter and its dq form

**Physical model.** Per phase, the converter terminal E connects to bus T through
Rf and Lf. The branch current i_k flows from E to T, so it is an **injection
into bus T**:

`Lf di_k/dt = e_k − v_k − Rf i_k`, for k = 0, 1, 2.

**Space vector:** `Lf dI_αβ/dt = E_αβ − V_αβ − Rf I_αβ`.

**Substitution.** With `I_αβ = w i_dq e^{jφ}`, `E_αβ = e_dq e^{jφ}` and
`V_αβ = v_dq e^{jφ}`:

`Lf w (di_dq/dt + j φ̇ i_dq) = e_dq − v_dq − Rf w i_dq`.

Since `Lf w = xf/ω0` and `Rf w = rf`:

**(P)  `di_dq/dt = (ω0/xf)(e_dq − v_dq − rf i_dq) − j φ̇ i_dq`**

In the matrix notation of the campaign request this is
`d(i_dq)/dt = Lf⁻¹(e_dq − v_dq − Rf i_dq) − φ̇ J i_dq`, with device-base
Lf = xf/ω0.

## 4. The frozen TX4 current dynamics

From spec §4.2 (kernel `_gfl_eval`, f[6] and f[7]):

`f_i,TX4 = (ω0/xf)(e_TX4 − v_dq − rf i_dq − j xf i_dq)`

with

`e_TX4 = v_dq + kp_i(i_ref − i_dq) + x_i + j xf i_dq`.

Here `x_i = x_id + j x_iq` and `i_ref = id_ref + j iq_ref`. Written by
component:
- `e_d = v_d + kp_i(id_ref − i_d) + x_id − xf i_q`;
- `e_q = v_q + kp_i(iq_ref − i_q) + x_iq + xf i_d`.

Equivalently,

`f_i,TX4 = (ω0/xf)(kp_i(i_ref − i_dq) + x_i − rf i_dq)`.

**Consequence.** The TX4 model is the physical filter (P) with φ̇ replaced by
the nominal ω0: its "nominal cross-coupling" ignores θ′. It depends on v_dq
only through i_ref, which depends on |V| via q_cmd (Q/V loop with gain g).

## 5. Voltage command (the V2 interface law)

**Command (C):**

`e_cmd = v_dq + rf i_dq + (xf/ω0) [ f_i,TX4(x, v_dq, i_dq) + j φ̇ i_dq ]`

**Proof that (P) then gives exactly the TX4 current dynamics.** Substituting (C)
into (P):

`di/dt = (ω0/xf)(rf i + (xf/ω0)(f_i + jφ̇ i) − rf i) − jφ̇ i = f_i,TX4`

The equality is exact.

**Equivalent form.** Substituting §4:

**`e_cmd = e_TX4 + j (xf/ω0) θ′ i_dq`**

That is, the frozen TX4 controller output plus a small correction for the actual
rotation rate of the PLL frame. At equilibrium θ′ = 0, so `e_cmd = e_TX4`.

**Converter source (system base, synchronous phasor).**
`E_sys = e_cmd e^{jθ}`, applied per phase as
`e_k(t) = Re(E_sys e^{j(ω0 t − 2πk/3)})`.

The source is an ideal three-phase voltage source with grounded neutral, and its
zero-sequence component is zero:
- In balanced operation, which covers every preregistered case, no
  zero-sequence current flows.
- The filter branches do give bus T a zero-sequence path to ground. This is
  declared, and it has no effect in balanced operation.

**No duplicated states.**
- The controller integrates only the 9 states θ, x_pll, p_f, q_f, x_p, x_q,
  x_id, x_iq and x_v.
- i_d and i_q are the **physical filter-branch currents**, measured as
  `i_dq = [(2/3)(i_a + a i_b + a² i_c) e^{−jω0 t}] e^{−jθ} / w`.
- The controller uses the measured i_dq wherever TX4 uses i_d and i_q: P, Q,
  the current loops and the x_id/x_iq integrators.
- f_i,TX4 is evaluated only inside the command (C). It is never integrated.

## 6. Network companion and signs

**Companion.** The filter branch is a trapezoidal R–L companion, with the same
formulas as ParaEMT's R–L branches:
- `Req = Rf + 2Lf/Δt`;
- `i_n = (e_n − v_T,n)/Req + Ihis_n`;
- `Ihis_{n+1} = icf · i_n + Gv1 · (e_n − v_T,n)`, with
  `icf = (2Lf/Δt − Rf)/Req` and `Gv1 = 1/Req`.

**Derivation.** From `L(i_n − i_{n−1})/Δt + R(i_n + i_{n−1})/2 = (u_n + u_{n−1})/2`:

`i_n = u_n/Req + [(2L/Δt − R) i_{n−1} + u_{n−1}]/Req`

**Stamping.** Node E is eliminated exactly, because the source is ideal. At bus
T, for each phase:
- the diagonal gets `+1/Req` in G;
- the right-hand side gets `+ e_k,n/Req + Ihis_k,n`.

G stays constant: it is factorized once per case, as in v1.

**Sign check.** A positive (E − V) drives a positive injection into T, which is
the same sign as the v1 injection `I = w i e^{jθ}`.

## 7. Discretization (fixed before any run) and its analysis

The time loop is the v1 loop (prereg v1 §4): an Euler predictor, one network
solve, then a Heun corrector. For each GFL at step n (t_{n−1} → t_n):

1. **At t_{n−1}.** Measure i_{n−1} using θ_{n−1}, and take V_{n−1} (both
   synchronous phasors). Evaluate the controller derivative
   `F0 = F(x_{n−1}, v_{n−1}, i_{n−1})`.
2. **Predict.** `x_p = x_{n−1} + Δt F0` (9 controller states).
3. **Command at t_n.**
   - Extrapolate linearly: `V̂_n = 2V_{n−1} − V_{n−2}` and
     `Î_n = 2I_{n−1} − I_{n−2}` (synchronous phasors). At n = 1:
     `V̂_1 = V_0`, `Î_1 = I_0`, which is exact at equilibrium.
   - Rotate: `v̂ = V̂_n e^{−jθ_p}`, `î = Î_n e^{−jθ_p}/w`.
   - Frame rate: `θ′_p = kp_pll v̂_q + x_pll,p`.
   - Command: `e_cmd,n = v̂ + rf î + (xf/ω0)[f_i,TX4(x_p, v̂, î) + j(ω0 + θ′_p) î]`.
   - Source: `E_n = e_cmd,n e^{jθ_p}`.
4. **Solve** the network with the filter Norton of §6.
5. **Update** the branch currents and histories (§6).
6. **Correct.** `x_n = x_{n−1} + (Δt/2)(F0 + F(x_p, v_n, i_n))`, where i_n is
   measured using θ_p.

**Why extrapolation, and not the previous sample (analysis only).** Replacing
V_n by V_{n−1} inside (C) leaves an error `(V_{n−1} − V_n)` in the filter
voltage. Its sum over steps **telescopes** to −ΔV_total. Every change of the
terminal voltage therefore leaves a net current impulse of
`(ω0 Δt / xf) · ΔV` (device pu), whatever its speed.
- For the preregistered phase step (ΔV ≈ 0.02 pu), that is 2.5e-3 pu at 50 µs.
- This is larger than the i_q excursion of those unit-test cases (≈ 6e-4 in the
  v1 bus-30 reference).
- It is O(Δt), so it would also show up as non-convergence between 50 and
  25 µs.

With linear extrapolation the error sequence at a jump is a +/− doublet whose
net impulse is zero, and the scheme is second order for smooth V. The same
extrapolation is applied to I for consistency.

**Terminal-node (Nyquist) stability (analysis only).** Take bus T between the
filter Lf and an inductive network L_n, with β = Lf/(Lf + L_n) ∈ (0, 1). Let the
command depend linearly on past terminal voltages, `E_n = E(z)·v`. The
trapezoidal recurrences of the two series inductors then give the
characteristic equation

`(z + 1) [ E(z)(1 − β) − 1 ] = 0`

- **z = −1.** This is the trapezoidal "split" mode of two series inductors. It
  is **neutral** for any command law, because the sources drive only the sum
  mode. Only branch resistance (Rf here) damps it, and it is excited only by
  discontinuities.
- **Previous-sample law, E(z) = 1/z:** `z = 1 − β`, which is stable.
- **Linear extrapolation, E(z) = 2/z − 1/z²:**
  `z² − 2(1 − β)z + (1 − β) = 0`, so `|z|² = 1 − β < 1`, which is stable.

**Contrast with v1.** The v1 ideal current source forced its current into the
inductor. The z = −1 mode of the inductor voltage was then driven by the
controller's instantaneous current response (|V| → q_cmd → i_ref, and the PLL),
and it grew (v1 D4: 1.0013 per step at g = 0).

The controller paths neglected in this analysis enter the command with gains
far below the unit feed-forward gain:
- |V| → e: kp_i·kp_q·kp_v·g ≤ 0.025;
- PLL: (xf/ω0)|i|·kp_pll ≈ 0.01.

## 8. Initialization

1. **Operating point** (unit test or IEEE-39 canonical power flow): V0 at T and
   the device output S.
2. **TX4 initialization** (`gfl_init`): `I = conj(S/(w V0))`, `θ0 = ∠V0`,
   `i0 = I e^{−jθ0}`. The controller states are
   (θ0, 0, p_ref, q_ref, i_d, −i_q, ·, ·, rf i_d, rf i_q, q_ref); the i entries
   are measured, not states.
3. **Equilibrium command:** `e_cmd,0 = e_TX4,0 = v0 + (rf + j xf) i0`, so that
   `E0 = V0 + (Rf + jXf) I_sys,0`.
4. **Filter history.** From the continuous phasor steady state, as ParaEMT does
   for its branches:
   - `i_k,0 = Re(I_sys,0 e^{−j2πk/3})`;
   - `u_k,0 = Re((E0 − V0) e^{−j2πk/3})`;
   - `Ihis_k,1 = icf · i_k,0 + Gv1 · u_k,0`.

   The trapezoidal reactance warp (≈ (ω0Δt)²/12 = 3e-5 at 50 µs) leaves a small
   start-up mismatch. It is measured by the G4b no-event stationarity check.
5. **Extrapolation history:** `V_{−1} := V_0`, `I_{−1} := I_0`.

## 9. Software tests (preregistered; gate SW; before any G4b run)

- **T1. Identity.** Use 2000 random vectors: states around the equilibria of the
  three G4 cases, |v| ∈ [0.9, 1.1], ∠v ∈ [−π, π], θ′ ∈ [−5, 5] rad/s, and
  g ∈ {0, 0.03625, 0.25}.
  - Check that (P) evaluated at `e_cmd` of (C), as produced by the kernel
    function, equals `f_i,TX4` from `_gfl_eval`, with relative error ≤ 1e-12.
  - Check also that `e_cmd = e_TX4 + j(xf/ω0)θ′ i` to 1e-12.
- **T2. abc-level filter.** Prescribe smooth trajectories x(t), v_dq(t) and θ(t)
  over 0.2 s, and integrate:
  - the per-phase physical filter `Lf di_k/dt = e_k − v_k − Rf i_k` with e_k
    from (C), using solve_ivp with rtol 1e-11;
  - the dq ODE `di/dt = f_i,TX4`.

  Pass if the difference of i_dq is ≤ 1e-7 (device pu).
- **T3. Companion.** Compare Req, icf and Gv1 against ParaEMT's `numba_InitNet`
  coefficients for an R–L line of the same R and L (net_damping = 0). Pass at
  1e-12 relative.
- **T4. Bases and signs.** Static, at t = 0, before any network solve. For
  each G4 case, the initialized branch currents must give:
  - measured `i_dq = i0` to 1e-12;
  - `V0 · conj(I_sys,0) = S` (system base) to 1e-12;
  - `e_cmd,0` from the kernel command function equal to
    `v0 + (rf + j xf) i0` to 1e-12.

**Gate SW: PASS iff T1–T4 all pass.**
