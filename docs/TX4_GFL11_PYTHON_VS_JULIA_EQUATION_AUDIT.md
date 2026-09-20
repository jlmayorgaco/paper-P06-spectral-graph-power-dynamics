# TX4 GFL11 Python vs Julia Equation Audit

## Verdict

The device-level parity run passed exactly at the preregistered fixed test point: the maximum absolute difference across all 11 derivatives and both injection components is `0.0` in `TX4_GFL11_DEVICE_PARITY.csv`.

## State order

Both implementations use the same 11-state order:

`[theta_pll, x_pll, p_filt, q_filt, x_p, x_q, i_d, i_q, x_id, x_iq, x_v]`.

This is the frozen Python `CONVERTER_LABELS` tuple plus `x_v`, and the Julia `GFLParams`/`gfl_derivatives` unpacking in `code/tx4/run_tx4_exact_p4_ieee39.jl`.

## Algebraic and power equations

For terminal voltage `v` and PLL angle `theta`, both codes form `v_d + j v_q = v exp(-j theta)`. They then use

- `P = v_d i_d + v_q i_q`,
- `Q = v_q i_d - v_d i_q`,
- `err = g (v_ref - |v|)`,
- `q_cmd = kp_v err + x_v`,
- `i_d_ref = kp_p (p_ref - p_f) + x_p`,
- `i_q_ref = -(kp_q (q_cmd - q_f) + x_q)`.

The inner voltage equations are identical:

- `e_d = v_d + kp_i(i_d_ref-i_d) + x_id - x_f i_q`,
- `e_q = v_q + kp_i(i_q_ref-i_q) + x_iq + x_f i_d`.

The ten original derivatives are unchanged from the frozen GFL10 equations. The eleventh derivative is identical in both languages:

`x_v' = ki_v err - leak (x_v - q_ref)`.

The network injection is identical:

`I = w (i_d + j i_q) exp(j theta)`.

## Initialization

At the frozen AC point, both codes compute device-base current from `conj(S/(w V))`, rotate by the terminal angle, set `p_ref=v_d i_d`, `q_ref=-v_d i_q`, and append `x_v=q_ref`. Consequently the matched reactive policy leaves the equilibrium voltage unchanged across the 16-case census.

## Linearization

Both codes central-difference the same semi-explicit residuals and reduce the index-one DAE as `A = f_x - f_z g_z^{-1} g_x`. Julia additionally exports the reconstructed algebraic voltage modes `dz = -g_z^{-1}g_x dx`; Python exports the same quantity for cross-code MAC matching.

## Source anchors

- Python device: `research/ias2026_last_validation/raw/true_same_model/canonical_source/models/ieee39_devices.py`.
- Python assembly: `research/ias2026_last_validation/raw/true_same_model/canonical_source/models/ieee39_case.py`.
- Julia reproduction: `code/tx4/run_tx4_exact_p4_ieee39.jl`.
- Numerical result: `results/TX4_GFL11_DEVICE_PARITY.csv`.
