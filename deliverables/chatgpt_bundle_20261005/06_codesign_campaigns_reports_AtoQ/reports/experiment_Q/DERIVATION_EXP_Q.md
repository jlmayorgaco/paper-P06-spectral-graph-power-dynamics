# Derivations and model checks for Experiment Q

## Load-step input and gauge quotient

The bus-16 event is a +100 MW active-load step with Qset held fixed. The installed IEEE-39 ZIP parameters at bus 16 have all active and reactive fractions on the constant-impedance term. For a unit MW increment on the 100 MVA system base and frozen pre-event voltage `V0`, the incremental active-load admittance has magnitude `1/(100 |V0|²)`; its realified current injection is assembled at bus 16 and eliminated through the exact algebraic network Jacobian `Gy`. The state input is therefore `Bload = -B (Gy \ yload)`. The frequency outputs are the declared local SG rotor speed and GFL PLL frequency states. The validated global rotational direction is removed with the exact quotient before computing finite poles or step responses.

The linear state response to a constant step is

`x(t) = A⁻¹ (exp(A t) - I) Bload ΔP`,

with `Δf(t)=Cf x(t)`, `RoCoF(t)=Cf exp(A t) Bload ΔP`, and `Δf∞=-Cf A⁻¹ Bload ΔP`. Extrema use a time mesh only to bracket derivative zeros; each bracket is refined by bisection. The mesh is a cross-check and not a frequency-domain robustness certificate.

## Modal residues and cancellation

For each simple pole of the gauge-free reduced matrix, the scalar load-to-frequency residue is `q_j=(Cf r_j)(l_jᴴ Bload)/(l_jᴴ r_j)`. The scripts reconstruct the frequency and RoCoF trajectories from the residues and compare them with direct matrix-exponential responses. For the frozen ExpP point, the largest local RoCoF channel is bus-38 PLL. The rightmost pole is `-0.0500000011 s⁻¹`, with 99.62% of the right eigenvector norm on GFL states; it contributes about `+40.53 Hz` to the eventual 100 MW linear steady offset. The total is about `+23.10 Hz`, because modal terms cancel. A fast real pole at `-1890.57 s⁻¹` gives the largest single modal contribution to initial RoCoF, about `-6.86 Hz/s` for 100 MW. Summed RoCoF modal terms have cancellation ratio 1.97; the steady-frequency terms have ratio 2.53. Thus the rightmost pole does not explain peak RoCoF.

## Why PLL gains do not supply sustained active power

At synchronized equilibrium the PLL phase detector error is zero. `Kp` and `Ki` multiply this phase error in the installed PLL equations and affect its dynamic poles and transient response, but they do not create a nonzero steady active-power/frequency droop law. Twenty centered gain perturbations at withheld interior mixed points give a largest absolute centered derivative of `F∞` of `1.49e-12 Hz/(MW·gain-unit)`. The largest system matrix condition number in those derivative tests is about `7.29e9`; the small absolute derivative is therefore reported with its conditioning and is not presented as an exact symbolic zero from finite differences alone.

The ExpP one-SG example reinforces the dynamic distinction. Doubling bus-38 Kp moved alpha from `-0.0500000011` to `-0.0499893451 s⁻¹` and doubled unit RoCoF peak from `0.19242` to `0.38485 Hz/s/MW`; halving Ki moved alpha to `-0.03542 s⁻¹`. The observed sampled beta upper witness changed across these gain cases, but those sampled values are not beta-star certificates.

## TGOV1 screening law and failed affine representation

At synchronized steady state, unsaturated TGOV1 has a local governor authority `c_i = Sn_i (1/R_i + DT_i)/60` MW/Hz per retained fraction. This gives the tabulated screening coefficients: 333.33 at bus 30, 233.33 at 31, 266.67 at 32/33/35, 200 at 34, 233.33 at 36/37, 333.33 at 38, and zero at bus 39 (no TGOV1).

The exact full-model test then infers `Kload=-1/F∞,unit - Σ c_i ε_i` for 100 deterministic withheld `(ε,Kp,Ki)` points. Inferred `Kload` ranges from `-3995.97` to `1438.82 MW/Hz`; its relative range normalized by the sample mean magnitude is 43.93%. The complete model is therefore classified `NONAFFINE`. The continuous-knapsack allocation of 57.21 MW at bus 30 is only a screening result. It is not a lower bound and is not used to claim globality.

## Initial local-PLL RoCoF at an architecture boundary

For a local PLL with state `Δω` in rad/s and filtered phase detector error `eθ`, the installed equation is `dΔω/dt=(Δω_i + Kp eθ - Δω)/τ_PLL`. Immediately after the algebraic voltage change caused by the step, the incremental PLL states are still zero, so `RoCoF_PLL(0+) = Kp eθ/(2π τ_PLL)`. This is independent of Ki and scales exactly linearly with Kp at that instant.

The one-sided experiment tests both ends of a single-bus mixed branch, with all other buses all-SG and all four Kp/Ki box corners. At the GFL-insertion boundary `rho=1e-6`, the minimum 100 MW initial PLL RoCoF by bus ranges from 1.013 Hz/s (bus 39) to 11.707 Hz/s (bus 35). At the SG-removal boundary `rho=0.999999` (`epsilon=1e-6`), it ranges from 9.451 Hz/s (bus 37) to 21.919 Hz/s (bus 35). All 80 tested boundary/corner combinations exceed the inherited 0.5 Hz/s limit. The result follows the specified maximum over every present PLL state, including a PLL attached to a very small GFL share. These are branch-priority diagnostics only. They do not prove that interior mixed points fail because the algebraic voltage response changes with the other support choices.

## Robustness evidence

For the all-SG feasibility reference, the exact shifted matrix is tested by the bounded-real CARE/LMI method at `γ=1/beta_req`. Float64 checks give positive `P` minimum eigenvalue `0.0024073`, CARE/LMI maximum residual eigenvalue `-0.08404`, closed-loop abscissa `-0.04807 s⁻¹`, and relative CARE residual `1.86e-6`. This is a numerical bounded-real certificate for `beta_star >= beta_req`; no directed rounding or interval arithmetic was used, so it is not called a formal validated-arithmetic proof.
