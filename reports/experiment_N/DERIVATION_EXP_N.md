# Derivation: fixed-support PD-exact model and local co-design certificate

## Component contract and equilibrium

Let `S_i^0=P_i^0+jQ_i^0` be the **initialized generator component** injection on the 100 MVA base and `V_i^0` its frozen complex bus voltage. For `ε_i=1−ρ_i`, the declared devices must inject `S_SG=ε_i S_i^0` and `S_GFL=ρ_i S_i^0`. The initialized ZIP load is a separate component, including at buses 31 and 39. The SG rating is `Sn_SG=ε_i Sn_i^0`; its H in seconds is unchanged, so stored energy `H Sn` scales with ε.

The full-dispatch module current is

`I_i^0 = conj((S_i^0/100)/V_i^0)`.

For the GFL's PLL-aligned dq coordinates, `θ=arg(V_i^0)`, `i_d=(P_i^0/100)/|V_i^0|`, and `i_q=−(Q_i^0/100)/|V_i^0|`. The installed filter equations give `v_{I,d}=|V|+R_f i_d−X_f i_q` and `v_{I,q}=R_f i_q+X_f i_d`. The steady CC integrators are `γ_d=v_{I,d}/K_CC,I` and `γ_q=v_{I,q}/K_CC,I`; the PLL frequency states vanish, `v_dc=V_dc`, `v_dc_i=i_d`, `iset_q=i_q`, and `P_dc=v_{I,d}i_d+v_{I,q}i_q`. The GFL aggregate presents `ρ I_i^0` externally while each module retains this per-unit equilibrium. These are trim references, not design degrees of freedom. See [the scaling derivation](DERIVATION_RHO_SCALING.md) for installed equation details and device-base semantics.

At a fixed architecture and frozen voltage, the normalized internal trim `x_i*` is independent of ρ, Kp, and Ki. Its implicit derivative solves `H_x dx*/dz+H_z=0` and equals zero in these coordinates. The external SG rating and GFL aggregate share still change, so the linearized **port** and network closure change materially. A local numerical correction of three surviving SG states is O(10⁻¹²) and addresses Float64 cancellation at bus 31, not a new physical setpoint or optimization variable.

## Exact physical port closure

For each present device, linearize the installed physical equations:

`M_i δẋ_i=A_i δx_i+B_i δv_i`, `δi_i=C_i δx_i+D_i δv_i`.

At mixed bus i, `Y_inj,i(s)=ε_i Y_SG,i(s)+ρ_i Y_GFL,i(s)` after the common base/sign conversion. Absent devices contribute no hidden states. With the frozen network/ZIP matrix `Y_static`, stack all device states and ports to form

`G_y=Y_static+D`, `A_red=A−B G_y⁻¹ C`.

The rotational gauge vector is removed by an orthonormal quotient `Q`; the physical finite poles are the eigenvalues of `A_q=QᵀA_red Q`. This gives the complete spectrum, including SG internal, GFL internal, collective, and intermodal poles. The all-GFL architecture has a second, **physical** generalized zero chain after quotient, so it fails the −0.05 s⁻¹ margin even though the gauge has been removed.

The port operator at a trial pole is

`T(s)=Y_static+Σ_i [D_i+C_i(sM_i−A_i)⁻¹B_i]`.

The exact network closure and the ordinary reduced matrix are equivalent where the local inverses exist. No fitted black-box transfer function or one-mode constraint is used.

## Total derivatives within an architecture

For an equilibrium equation `H(x*,z)=0`, the general derivative is `dx*/dz=−H_x⁻¹H_z`. Thus `dT/dz=T_z+T_x dx*/dz`, even when the equilibrium term happens to vanish under the present proportional trim. For the reduced matrix,

`dA_red=dA−dB G_y⁻¹C−B G_y⁻¹dC+B G_y⁻¹(dD)G_y⁻¹C`.

Here ρ changes SG/GFL `C,D` with opposite signs, and Kp/Ki change GFL `A,B`. For a simple pole of `A_q`,

`dλ/dz=(lᴴ Qᵀ(dA_red)Q r)/(lᴴr)`.

No derivative crosses `ρ=0` or `ρ=1`: those points remove the corresponding physical state block and require a new architecture. [TABLE_N08](TABLE_N08_derivative_validation.csv) records centered eigenvalue checks; [the matrix-difference addendum](TABLE_N08_matrix_derivative_addendum.csv) verifies the same 90 total derivatives without eigensolver subtraction noise.

## Active set and fixed-support KKT law

The nominal objective is `R=Σ_i P_i^0 ε_i` and the exact constraint is `g(Z)=α(Z)+0.05≤0`, where α is the maximum real part over **every** physical finite pole. At a simple active pole and fixed support, the Lagrangian is `L=R+μg`, `μ≥0`. For the winning support `{38}`, ε38 is the only free retention coordinate; all 20 gain coordinates sit at the frozen box corners. The free stationarity equation is

`P_38^0+μ ∂g/∂ε_38=0`.

At the frozen point `P_38^0=830 MW`, `∂g/∂ε_38≈−106.434487 s⁻¹`, giving `μ≈7.7982243`. The lower Kp bounds require `μ∂g/∂Kp_i≥0`; the upper Ki bounds require `μ∂g/∂Ki_i≤0`. All are strict at this point. LICQ follows from the nonzero free ε spectral derivative together with independent active box normals. Strict box multipliers force all gain components of a critical direction to zero; the active spectral equality then forces its ε component to zero. Hence the critical cone contains only zero and the SOSC is satisfied without estimating a noisy 21×21 Hessian. A 10⁻⁹ s⁻¹ design guard makes the frozen point strictly feasible; its complementarity residual is 8.62×10⁻⁹.

This certificate is local to `{38}`. The deterministic scans of 10 single-SG supports, 45 pairs, and 1,024 architecture budget samples did not certify all continuous connected branches. The global flag is therefore **NO**.

## Active-mode self-energy

Partition `T(s)` at the active mode into a two-coordinate pivot k and the remaining ports r. With `R=T_rr⁻¹`,

`T_eff,k=T_kk+Γ_k`, `Γ_k=−T_kr R T_rk`.

Its total parameter derivative splits as a direct term `dT_kk` plus

`dΓ=−dT_kr R T_rk−T_kr R dT_rk+T_kr R(dT_rr)R T_rk`.

The same split applies to `T_s`; dividing the left/right bilinear forms by the total `s` derivative yields direct and collective contributions to `dλ/dz`. At the active pole the SVD pivot is bus 36. Its direct `ρ38` contribution vanishes by locality; the ≈106.434488 s⁻¹ sensitivity arrives through the collective self-energy. The pairwise pathway matrix `−T_kℓ[T_rr⁻¹]_{ℓm}T_mk` is retained as a diagnostic. `Re(Y_port)` and `Im(Y_port)` are reconstructed from the real 2×2 rectangular blocks before graph analysis. The graph diagnostics do not replace the exact closure and do not imply `K=h(L)`.
