# Experiment M derivation and sign conventions

## Per-device power at a mixed bus

The installed Sauer–Pai machine reports `P` and `Q` in **machine pu**. With its initialized nameplate `Sn`, its physical terminal power is

\[
P_{\rm SG}=S_n P_{\rm machine},\qquad Q_{\rm SG}=S_n Q_{\rm machine}.
\]

The mixed builder sets `Sn=(1−ρ)Sn₀`. That scales the current conversion from machine to system base, and the physical stored kinetic energy proxy is `H Sn`. It does **not** fix `P_machine` at its original value. Componentwise initialization can change machine flux, torque, AVR, and governor states. Thus `Sn=(1−ρ)Sn₀` does not imply `P_SG=(1−ρ)P₀`.

For the weighted SimpleGFLDC, the filter's external terminal current is

\[
i_{\rm GFL,ext}=\rho(i_{f,r}+j i_{f,i}).
\]

There is no separate GFL `Sn` parameter in this implementation. Given bus voltage `u=u_r+j u_i`, the directly reconstructed system-base GFL power is

\[
P_{\rm GFL}=S_{\rm base}(u_r i_r+u_i i_i),\qquad
Q_{\rm GFL}=S_{\rm base}(u_i i_r-u_r i_i).
\]

The ZIPLoad observables use negative injection for consumption. Direct initialized observables satisfy

\[
(P,Q)_{\rm net}=(P,Q)_{\rm SG}+(P,Q)_{\rm GFL}+(P,Q)_{\rm ZIP}
\]

to below `1e−10` in every generator bus and audit case. The old analytical design instead sets the port contribution at the **original** SG and full-GFL equilibria to `(1−ρ)i_{SG,0}+ρ i_{GFL,0}`. That construction is a valid mathematical port homotopy, but the current PD initialization does not impose its per-device P/Q constraints.

## ZIP operating point

The installed ZIP model uses `P=Pset(KpZ Vrel²+KpI Vrel+KpC)` and the analogous Q expression, with `Vrel=|u|/Vset`. All frozen loads are impedance-only. The initialized `Vset` is a solved parameter, not necessarily 1 pu. At buses 31 and 39 it is 9.071760316660 and 1.767183595366 pu. This explains why the actual ZIP injections differ from the CSV setpoints, without implying an error in the old initialized-load reconstruction.

## Gauge and reduced dynamics

For each descriptor `M ẋ=A x+B d`, partition differential `x_d` and algebraic `x_a` with `M=diag(I,0)`. If `A_aa` is nonsingular,

\[
A_{\rm red}=A_{dd}-A_{da}A_{aa}^{-1}A_{ad},\quad
B_{\rm red}=B_d-A_{da}A_{aa}^{-1}B_a.
\]

The state alignment in TABLE_M10 is a physical permutation and name conversion only. We remove exactly **one** smallest-magnitude global angular gauge pole before one-to-one finite-pole matching. A second near-zero ExpG value is retained and explicitly audited; magnitude alone does not prove it is gauge.

The all-SG analytical reduced A differs from direct PD by `2.919e−15` relative Frobenius norm. Their raw descriptor A matrices are different algebraic realizations; the reduced dynamics and sampled port transfers are the appropriate physical identity checks in this case.

## PD-exact port closure

The installed `open_loop_linearization` exports bus impedance `Zbus` and static network admittance `Ynw`. The installed NetworkDynamics tests use `feedback(Zbus, Ynw + Yinj; pos=true)`. For these three IEEE-39 audit cases `Yinj` is absent. Writing the bus-port realization as `(M_Z,A_Z,B_Z,C_Z,D_Z)` and the static network admittance as `Y`, define

\[
W=I-YD_Z.
\]

The separate corrected analytical closure is

\[
\begin{aligned}
M_{\rm exact}&=M_Z,\\
A_{\rm exact}&=A_Z+B_ZW^{-1}YC_Z,\\
B_{\rm exact}&=B_ZW^{-1},\\
C_{\rm exact}&=C_Z+D_ZW^{-1}YC_Z,\\
D_{\rm exact}&=D_ZW^{-1}.
\end{aligned}
\]

This is implemented in `src/bnd_model_audit_m/pd_exact_closure.py` from exported component/port matrices. The explicit closure state order groups each bus's physical states; TABLE_M21 compares it against the direct PD descriptor after that permutation. M is exact, A relative error is below `4e−17`, and all finite poles agree within `1e−8 s⁻¹`. The construction has no optimizer and no fitted parameters.

The closure is valid **at the sampled initialized state**. A future analytical design must define a parameterized operating-point map that enforces desired SG/GFL P/Q shares and rejects violated physical bounds before differentiating or optimizing the resulting spectrum.

## Linear and nonlinear trajectory checks

The direct PD load-parameter Jacobian supplies the bus-16 Pset/Qset disturbance channel. For that common input, direct PD and the corrected closure have the same reduced `A`, `B`, and frequency output after physical state permutation. Tiny load-step and SG-speed perturbation trajectories therefore agree numerically; TABLE_M19 records the integrated residuals. The nonlinear pulse uses `Pset,Qset→(1+δ)(Pset,Qset)` on `[1.0,1.1) s`. Comparing with the exact PD linearization tests the validity of the linearization, not the old analytical model. TABLE_M20 reports three amplitudes and its fitted error exponent. A small-amplitude numerical floor must be distinguished from the quadratic asymptotic law.
