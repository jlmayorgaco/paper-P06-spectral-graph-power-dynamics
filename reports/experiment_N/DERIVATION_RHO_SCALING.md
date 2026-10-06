# ExpN device scaling and equilibrium trim

This derivation uses the installed PowerDynamics 5.0.0 equations in `Library/Machines/SauerPaiMachine.jl`, `Library/Controls/{AVRs,Govs}.jl`, and `Library/Renewables/ComposableInverter.jl`, plus the repository's `WeightedSimpleGFLDC` terminal adapter. The numerical identity tests in ExpN determine whether this construction is accepted.

## Contract and bases

At generator bus `i`, let `U_i` be the original initialized complex voltage and `S_i=P_i+jQ_i` the **generator component** power on `Sbase=100 MVA`. The unchanged ZIP injection, including at buses 31 and 39, is separate. The requested external powers are `S_SG=(1−ρ)S_i` and `S_GFL=ρS_i`, with `U=U_i`.

The full-dispatch system-base current is

`I_i = conj((S_i/Sbase)/U_i) = [(P_i/Sbase)U_r+(Q_i/Sbase)U_i]/|U_i|² + j[(P_i/Sbase)U_i−(Q_i/Sbase)U_r]/|U_i|²`.

The SG retains `Sn=(1−ρ)Sn₀`. Its terminal current conversion in the installed machine equations includes `Sbase/Sn`. Keeping the original machine-per-unit equilibrium and references therefore gives system current `(1−ρ)I_i`; `H` in seconds stays fixed and the physical stored kinetic energy proxy `H Sn` scales by `(1−ρ)`. At `ρ=1`, the SG state block is absent.

The GFL filter uses `Rf,Xf` and `Cdc` as per-unit quantities. An aggregate of parallel identical modules with fraction `ρ` has the same per-unit internal ODEs and an external terminal current `ρ i_f`; its physical filter/DC energy and DC input scale with aggregate capacity. The installed `SimpleGFLDC` has no explicit rating parameter. ExpN represents that aggregate using a `PortScale=ρ` terminal adapter, while keeping `Rf,Xf,Cdc`, controller gains and internal per-unit references invariant. At `ρ=0`, the GFL state block is absent. This interpretation is valid only if the port and complete-spectrum identity tests pass against an independently trimmed compiled network.

## Explicit GFL trim at fixed voltage and full-dispatch current

Let `θ=arg(U)`, `V=|U|`, `i_d=P/(Sbase V)` and `i_q=−Q/(Sbase V)` under the installed d-aligned `ri_to_dq` convention. The required inner filter current in global coordinates is `I_i`, independent of `ρ`; the adapter emits `ρ I_i` externally.

The installed PLL equilibrium is `θ=arg(U)`, `Δω_i=Δω=0`. Its gains `Kp,Ki` enter the Jacobian, but do not alter this equilibrium. The DC controller equilibrium is `v_dc_state=V_dc`, `v_dc_i=i_d`. The reactive trim parameter is `iset_q=i_q`.

For `CC1_F=CC1_Fcoupl=0`, the zero-error current controller requires

`V_I,d=V+Rf i_d−Xf i_q`, `V_I,q=Rf i_q+Xf i_d`, and `γ_d=V_I,d/CC1_KI`, `γ_q=V_I,q/CC1_KI`.

The DC input trim is `P_dc=V_I,d i_d+V_I,q i_q=V i_d+Rf(i_d²+i_q²)`. This term is **converter-side** active power, including filter dissipation. It is not the AC terminal `P` parameter. No PLL or current-controller gain is a trim variable.

For a controlled SG, the original AVR `vref` and governor `p_ref` are retained because machine-per-unit `P,Q,V` are unchanged. At bus 39 the original fixed `vf_set,τ_m_set` are retained. SG machine torque and field bounds must be checked after assembly.

The trim is valid only if the compiled network has the specified `P,Q`, all differential/algebraic residuals below the gate, and no declared internal bound violations. A failure is reported as `TRIM_INFEASIBLE`; no reactive redistribution is allowed.
