# Experiment A — Dynamic Self-Energy Extraction

## 1. Executive result

**EXPERIMENT A: PASS**

The primary bus-33 case uses the installed PowerDynamics.jl 5.0.0 IEEE-39 model, with its generator replaced by the library `SimpleGFLDC`. The descriptor has 189 variables (111 differential and 78 algebraic). The normalized equilibrium residual is 1.6287446460197112e-13; the Schur median/p95 residual is 7.149733680142577e-20/3.757006678717024e-19; the BND bridge is `PASS-EXACT` with primary-band p95 error 4.1861097973902148e-16.


The critical non-gauge pole is `-0.09884700664631951 + -0.43026220070517623im` at 0.06847835606782597 Hz. The exact validated second-order result is the more general `T_q(s)=s²M+sD₀+L₀+Π_q(s)`. Since `T_cc(0)` is singular, no global `Σ(s)` or graph-modal transform is asserted, and Experiment B is not yet ready. Experiment A status is PASS; evidence and scope are detailed below.

## 2. Research question

Whether the detailed mixed SG–GFL IEEE-39 linearization can be reduced by an exact descriptor Schur complement, whether synchronization-visible finite poles are preserved, and whether the actual retained equations admit an exact or quantified angle-level second-order bridge. No controller optimization or placement study is included.

## 3. PowerDynamics model provenance

See `tables/TABLE_A01_model_provenance.csv` and `TABLES_EXP_A.md` for Markdown-rendered tables. The source is the maintained package example at `C:\Users\walla\.julia\packages\PowerDynamics\VzOiZ\docs\examples\ieee39_part1.jl`; the package environment is the repository `Project.toml` and `Manifest.toml`. The reference bus is bus 31. The full state ordering is in `tables/TABLE_A03_state_partition.csv`; baseline and mixed equilibrium vectors with exact names are in `matrices/bus33_baseline_equilibrium.csv` and `matrices/bus33_mixed_equilibrium.csv`; raw matrices and eigenvectors are under `matrices/`.

## 4. Mixed SG–GFL construction

The bus-33 component is `PowerDynamics.Library.ComposableInverter.SimpleGFLDC` using the repository's frozen nominal SimpleGFLDC settings. `replace_bus` reuses that bus's original PowerDynamics power-flow model. No auxiliary controller or gain search was introduced. Its PLL retains the installed `pll₊θ`, `pll₊Δω_rad_s`, and `pll₊Δω_i_rad_s` states; converter/filter, current-control, and DC-link states remain in the full model.

## 5. Equilibrium

Initialization uses `solve_powerflow` followed by `initialize_from_pf!`. Maximum normalized residual is 1.6287446460197112e-13, differential residual 3.0518849357575527e-13, and algebraic residual 1.7175150190951172e-13. Bus-wise PF targets, measured initialized busbar P/Q flowing into the network, mismatch, voltage magnitude and angle are in `tables/TABLE_A02_operating_point.csv`; aggregate P/Q sums are in `TABLE_A02_balance_summary.csv`. PowerDynamics BusBase defines the measured sign convention as P=u_r(-i_r)+u_i(-i_i), Q=u_i(-i_r)-u_r(-i_i).

## 6. Full descriptor linearization

NetworkDynamics `linearize_network` returns `E ẋ = A x`, with `E=diag(1 for differential coordinates, 0 for algebraic coordinates)`. Here `size(E)=(189,189)`, rank(E)=111, and the algebraic dimension is 78. The blocks `F_x,F_y,G_x,G_y` and reduced matrix are preserved as CSV. Algebraic elimination is `A_red=F_x-F_y(G_y\G_x)` after verifying `G_y` condition 2174.668491313648 and σmin 0.4718974651061375. Generalized QZ and reduced finite spectra agree to 3.333913236981773e-9. The rotational/gauge mode count is 1; it is retained in raw files and excluded only from critical engineering classification.

The generalized pencil has 111 finite poles and 78 infinite algebraic eigenvalues; only finite values enter the pole analysis. For `T(λ,z)=λE−A(z)`, `dλ/dz=(lᴴA_zr)/(lᴴEr)` with no extra minus. The synthetic central-perturbation test checks this sensitivity sign; model-specific finite differences check the Jacobian sign. No physical-parameter sensitivity is claimed.

## 7. Retained/condensed state partition

Retained coordinates are 20 differential states: every installed SG rotor angle/speed pair, plus GFL PLL angle `θ` and LPF frequency-deviation state `Δω_rad_s`. The PLL integrator `Δω_i_rad_s` is condensed along with the remaining machine internal, AVR/governor, converter/current-controller/filter/DC-link and network states. The angle-rate map `q̇=Cv` is measured directly from the Jacobian; its relative residual is 0.0. The algebraic variables are eliminated separately through `G_y`.

Pairs are ordered by ascending bus/device: q uses machine δ at buses 30, 31, 32; GFL θ at bus 33; then machine δ at buses 34–39. v uses the matching machine ω or PLL Δω_rad_s states in that order. Angles are rad, synchronous-machine speed deviation is pu on a 60-Hz base, PLL Δω_rad_s is rad/s; bus voltage/current, power-flow targets and terminal power are pu; time is s; eigenvalues are s⁻¹; frequency is |Im(λ)|/(2π) Hz. Exact state names and full 189-coordinate order are listed in TABLE_A03 and the named equilibrium CSVs.

## 8. Exact Schur derivation

After algebraic elimination, `T(s)=sI-A_red`. Reorder the dynamic state indices as retained `r=[q;v]` and condensed `c`; define `T_rr,T_rc,T_cr,T_cc` by those index sets. Then `Π(s)=T_rc(s)T_cc(s)⁻¹T_cr(s)` and `S_r(s)=T_rr(s)-Π(s)`. Production code uses linear solves, not an explicit inverse.

## 9. Numerical validation of exact Schur identity

The imaginary-axis grid spans 0.01–100 Hz with 121 extra points over 0.1–5 Hz; controlled offsets around the least-damped poles are also recorded. For points with `cond(T_cc)<1.0e12`, median/p95/max normalized block reconstruction residuals are 7.149733680142577e-20, 3.757006678717024e-19, and 8.901815061727405e-19. Log-absolute-determinant residual is separately tabulated. Figure `FIG_A03_schur_identity_residual.png` shows the frequency-grid residual against frequency.

## 10. Pole preservation

The full finite spectrum contains 111 poles. For each of the 20 least-damped poles and every pole at or below 5 Hz, `σ_min(S_r(λ))`, `σ_min(T_cc(λ))`, retained participation, classification, and normalized Schur pole residual are in `tables/TABLE_A05_schur_pole_validation.csv`. The accepted check is the singular-value-at-full-pole certificate; no independent nonlinear root search is claimed. Retained-visible pole count is 77; maximum normalized residual is 4.0048975131161274e-18.

## 11. Derivation/test of the second-order bridge

The measured kinematic rows are `q̇=Cv`. Substituting `v=C⁻¹sq` in the retained frequency-state equations after the exact controller-state Schur elimination gives `T_q(s)=s²M+sD₀+L₀+Π_q(s)`, with `M=C⁻¹`, `D₀=-A_vvC⁻¹`, `L₀=-A_vq`, and `Π_q(s)=-A_vc(sI-A_cc)⁻¹(A_cq+sA_cvC⁻¹)`. This form is compared numerically with a separately assembled expression from `S_r(s)`. Classification: `PASS-EXACT`. If `Π_q(0)` is finite, set `L_eff=L₀+Π_q(0)` and `Σ(s)=D₀+[Π_q(s)-Π_q(0)]/s`, whose zero-frequency value is the analytic derivative; otherwise the exact generalized form is kept and no global `Σ(s)` is asserted. The bridge error in 0.1–5 Hz has p95 4.1861097973902148e-16 and max 9.481307414613845e-16.

## 12. Frequency-dependent self-energy

Table `TABLE_A06_frequency_residuals.csv` contains the frequency-dependent norm, leading singular values, Hermitian-part extrema, skew-Hermitian/reactive-part norm, conditioning and bridge residual. The characterized object is `Pi_q`. A global Σ(s) is not defined because the required Π_q(0) evaluation failed: `SingularException(83)`. The Hermitian/skew decomposition of Π_q is only an algebraic diagnostic and is not called damping. `TABLE_A07_sigma_critical_frequency.csv` stores each critical-frequency matrix entry.

## 13. Optional graph-modal diagnostic

The graph-modal transform is not available. Π_q(0) is not finite; no static stiffness was asserted.

## 14. Nonlinear-vs-linear validation

A +0.1% active-load pulse is applied at bus 20 from 1.0 to 1.1 s. The same perturbation enters the full linear descriptor after algebraic reconstruction. `TABLE_A08_TDS_validation.csv` and `TDS_trace.csv` compare bus-30 speed, bus-20 angle, bus-33 PLL angle, and the bus-33 terminal active-power proxy `P=V_r I_r+V_i I_i` (per-unit, current sign as defined by the SimpleGFL filter states). Exact event-time samples at 1.0 and 1.1 s are excluded from summary errors because the nonlinear callback saves both sides of those discontinuities; all samples remain in the raw trace. Maximum NRMSE is 0.00010398809829054601; solver retcode is `Success`.

## 15. Cross-bus validation

Single nominal GFL replacements at buses 30, 35 and 37 were run only after bus 33 passed the exact Schur/pole gates. Summary: `TABLE_A09_cross_bus_validation.csv`; figure `FIG_A11_cross_bus_summary.png`. No parameters were retuned.

| bus | equilibrium_pass | spectral_abscissa | critical_frequency_hz | schur_p95_error | pole_max_error | bnd_classification | bnd_p95_error |
|---|---|---|---|---|---|---|---|
| 33 | true | -0.09884700664631951 | 0.06847835606782597 | 3.757006678717024e-19 | 4.0048975131161274e-18 | PASS-EXACT | 4.1861097973902148e-16 |
| 30 | true | -0.10658545410744248 | 0.07445033432918856 | 3.416795120639614e-19 | 1.0353343592030239e-17 | PASS-EXACT | 2.7551670896934464e-16 |
| 35 | true | -0.09864801168090992 | 0.06840941929179577 | 3.158925430481922e-19 | 3.2555141649438855e-17 | PASS-EXACT | 5.292835677437307e-16 |
| 37 | true | -0.13790490492384927 | 0.0 | 4.275377350980973e-19 | 3.540437278269932e-17 | PASS-EXACT | 3.533998955742157e-16 |


## 16. Gate summary

| gate | observed | status | notes |
|---|---|---|---|
| GATE A0 — MODEL PROVENANCE | PowerDynamics 5.0.0; IEEE-39 package example | PASS | Source path and exact state map saved. |
| GATE A1 — EQUILIBRIUM | 1.6287446460197112e-13 | PASS | F/G residuals both reported. |
| GATE A2 — LINEARIZATION SANITY | 1.5649471803316973e-10; 3.333913236981773e-9 | PASS | Raw descriptor retained; finite-difference probes deterministic. |
| GATE A3 — EXACT SCHUR IDENTITY | 7.149733680142577e-20; 3.757006678717024e-19 | PASS | Conditioned and all-point metrics saved. |
| GATE A4 — RETAINED POLE PRESERVATION | 4.0048975131161274e-18 | PASS | Pole root solving omitted; singular-value test used. |
| GATE A5 — BND SECOND-ORDER BRIDGE | PASS-EXACT; p95=4.1861097973902148e-16 | PASS-EXACT | Π_q is derived from model blocks; global Σ(s) is unavailable because T_cc(0) is singular. |
| GATE A6 — SMALL-SIGNAL/TDS CONSISTENCY | 0.00010398809829054601 | PASS | Small local pulse only. |
| GATE A7 — CROSS-BUS REPRODUCIBILITY | true | PASS | No retuning. |



## 17. What Experiment A establishes

It establishes only the scope supported by the gates: the installed IEEE-39 model's actual DAE structure and state order, exact Schur reconstruction for the declared retained coordinates, pole certificates for visible modes, and the measured angle-level second-order representation/error.

## 18. What Experiment A DOES NOT establish

Experiment A does not establish optimal PLL tuning, optimal SG→GFL placement, global stability, nonlinear transient stability beyond the tested local pulse, EMT validity, field validation, novelty of graph damping, or robustness across operating points.

## 19. Decision for Experiment B

Proceed to `Σ(s) → Σhat(s) → Γ_k(s)`: **NO**. The Schur, pole, exact generalized BND, TDS and cross-bus gates pass, but `T_cc(0)` is singular (`SingularException(83)`), so `Π_q(0)` is unavailable and the prescribed global `Σ(s)` cannot be formed. The next diagnostic is to identify this zero-frequency condensed-block singularity without changing the retained coordinates or fitting parameters; until then characterize `Π_q(s)` only.
