# 01 — Equation ledger (BND / IEEE-39 project)

Date: 2026-10-05. Scope: our own documents in `resources/`, `experiments/`, `reports/`, `temp/`; external papers were used only to challenge novelty.
Confidence codes: **V** verified (proof + independent numeric check), **S** supported, **X** exploratory, **Q** questionable, **R** false / retired.
Kind: ID exact identity, TH theorem, AP approximation, CJ conjecture, EM empirical.
"THEORY" = `experiments/theory_collective_damping_20261003/THEORY.tex` (line numbers). "IEEE9" = `resources/BND_IEEE9_Informe_Validacion.pdf` (5 Oct). "CUAD" = `resources/CUADERNO_MATEMATICO_BND_SPARSE_PLL.md` (5 Oct). "T1–T5" = new checks in `research_gold/checks/`.
Duplicated statements across documents are merged into one entry with cross-references.

## A. Models

| ID | Equation | Source | Kind / conf. | Notes, assumptions, links |
|---|---|---|---|---|
| E01 | Nonlinear DDAE: `E x' = f0(x,z;rho) + sum_i b_i(k_i) e_i(t - tau_i)`, `0 = g(x,z;rho,d)` | CUAD §1; THEORY §2; `src/pd39` | model / V | Lossy AC network, SG with governor, GFL with DC link, P/Q. Delay sits on the PLL detector only. |
| E02 | PLL: `theta' = w_p`, `t_f w_p' = -w_p + xi + k_p e(t-tau)`, `xi' = k_I e(t-tau)` | CUAD §1; IEEE9 §2; poster panel 1 | model / V | Gives loop factor `kappa(s) = (k_p s + k_I) e^{-s tau}` (E20). |
| E03 | Replacement share `P_GFL = sum_i P_i^0 rho_i`; SG inertia `M_i = 2 H_i S_i (1 - rho_i)` | IEEE9 §2; mega prompt | model / V | rho changes equilibrium sharing, inertia and network return together. |
| E04 | Simplified swing network `M d'' + D d' + E Gamma sin(E^T d) = P + u`, `L* = E diag(gamma cos eta*) E^T` | prompt seed A; Compendio p.19; THEORY §freq. limits | reduced model / V in its class | NOT the IEEE-39 model. Used only for E40–E45. |
| E05 | Linearised characteristic `Delta(s) = sE - A0 - sum_i b_i(k_i) c_i^T e^{-s tau_i}` (exact exponentials, no Padé) | CUAD §1; THEORY; `experiments/interaction_decision_20261004/model.py` | ID / V | Exported IEEE-39 model: 204 states (203 physical + gauge), 10 PLL channels, E = I. |

## B. Reductions, low-rank structure, determinants

| ID | Equation | Source | Kind / conf. | Notes |
|---|---|---|---|---|
| E10 | Dynamic Schur: `S_r(s) = sI - A_rr - A_rc (sI - A_cc)^{-1} A_cr`; `det(sI-A) = det(sI-A_cc) det S_r` | Feedback-Dressed Topology p.2; Compendio; THEORY | ID / V | Exact where the eliminated block is regular. Standard algebra. |
| E11 | Self-energy `Sigma_r(s) = A_rc (sI-A_cc)^{-1} A_cr`; return ratio `R_c = G_c A_cr G_r A_rc`; `det(sI-A) = det T_r det T_c det(I - R_c)` | Feedback-Dressed p.2 | ID / V | Document itself calls the algebra standard. |
| E12 | "Master operator" `T(s) = s^2 M + s D + L + Sigma(s)` | Compendio p.19; prompt seed B | interpretation / S | Legitimate only for the retained second-order coordinates of a specific partition; not a general form of the IEEE-39 DDAE. |
| E13 | Physical port impedance `Z(s) = M (Delta_vv - Delta_vh Delta_hh^{-1} Delta_hv)`, `d = Z(s) nu` | THEORY 2027 ff.; CUAD §3; IEEE9 §4 | ID / V | nu = SG speed deviations, h = hidden states. Requires invertible hidden block. |
| E14 | Action-space factorisation `T(s;y,tau) = T_ref + U Theta(y,s,tau) V^H`, rank ≤ 30 (≤ 10 at fixed rho) | mega M2; `ieee39_exact_action_space_codesign_20261002` | ID / V | Reconstruction error 3.2e-16 (21 designs × 3 delay fields). 282-state descriptor. |
| E15 | Determinant lemma `det T = det T_ref det(I_r + Theta V^H T_ref^{-1} U)`; for m PLLs `det Delta_u / det Delta_0 = det(I_m + C_S^T Delta_0^{-1} P_S)` | mega M2; CUAD §4; Who-Moved-Mode p.3 | ID / V | Log-det error 1.2e-12; 120 support tests 6.6e-14 (IEEE9). Classical lemma. |
| E16 | Pair identity `det F det F0 / (det F_i det F_j) = det(I - P)`, `P = K_ij K_ji`, `K_ij = (I+H_ii)^{-1} H_ij` | THEORY 1111 | ID / V | Re-derived and checked on random matrices (this session). Schur complement of the 2×2 block. |
| E17 | Interaction `I = z - z0 - sum_i (z_i - z0) = -(1/(2 pi j m)) oint (s-c) tr[(I-P)^{-1} P_s] ds` | THEORY 1111–1150 | ID / V | Needs equal root counts and root-free contour. First term `(1/(2 pi j m)) oint tr P ds`. |
| E18 | Contour moment `M_1 = (1/2 pi j) oint s Delta'/Delta ds = sum of enclosed roots`; Möbius dividends | Who-Moved-Mode p.3; Synthesis p.3 | ID / V (never simulated in those docs) | Argument principle. Attribution through `log det` is representation-dependent (R, see C18 in file 02). |
| E19 | SCC factorisation `det(I - Q_M) = prod_SCC det(I - Q_c)`; acyclic Q ⇒ `det(I+Q) = 1` | Who-Moved-Mode p.3; Two-Gate p.3–4 | ID / V | Graph-theoretic block triangularity; standard. |

## C. The PLL channel: gains, delay, sensitivity (core of the recommendation)

| ID | Equation | Source | Kind / conf. | Notes |
|---|---|---|---|---|
| E20 | PLL return `F_i(s) = s^2 (1 + t_f s) - (k_p s + k_I) e^{-s tau_i} G_i(s)` | CUAD §7; IEEE9 §3; THEORY 2348 | ID / V | `G_i` = full grid return seen by PLL i with all other loops closed. |
| E21 | Pole assignment: `W = lambda^2 (1+t_f lambda) e^{lambda tau} / G_i(lambda)`, `k_p = Im W / Im lambda`, `k_I = Re W - Re(lambda) k_p` | THEORY 2348, 2668; IEEE9 §3; poster panel | ID / V | Residual 1e-13…1e-16 (IEEE-39), 1.35e-14 (IEEE9). Exact for ONE pole; no dominance guarantee (IEEE9 Fig. 3: assigned pole at −0.30, another pair stays at +0.567). |
| E22 | All-PLL map: prescribe `lambda` and pattern `q`, `y = G_PLL(lambda,rho) q`, `W_i = lambda^2(1+t_f lambda) e^{lambda tau_i} q_i / y_i` | THEORY 2668; `Analytic_PLL_formula_validation_summary.csv` | ID / V | 1000 trials, max 2e-15; invariant to complex scaling of q (1.8e-11). Changes 20 gains. |
| E23 | **Delay transport**: `k_p(h) lambda + k_I(h) = (k_p lambda + k_I) e^{lambda h}`, i.e. `k_p(h) = e^{a h}[k_p cos wh + ((k_I + a k_p)/w) sin wh]`, `k_I(h) = e^{a h}[k_I cos wh - ((|lambda|^2 k_p + a k_I)/w) sin wh]` | THEORY ≈2395–2410 (gain flow); IEEE9 §3 | ID / V | Needs no network model. IEEE9: 45 cases 3.9e-14. T2: protected pole drift 0.0 on IEEE-39 for all ten sites together. |
| E24 | Gain flow `d/dh [k_p; k_I] = [[2a, 1], [-|lambda|^2, 0]] [k_p; k_I]` | THEORY eq. gainflow | ID / V | Matrix has eigenvalues `lambda, conj(lambda)`; E23 is its solution. |
| E25 | **Exact remainder**: `kappa(s;h) - kappa(s;0) = -(s-lambda)(s-conj lambda) e^{-s tau} int_0^h k_p(t) e^{-st} dt` | derived this session from E24 (THEORY states divisibility only) | ID / V (algebra) | The delayed PI interpolates the old loop factor at exactly `lambda, conj(lambda)`; everything else pays a quadratic factor. |
| E26 | **Modal toll (spillover)**: `d mu / dh_i = -k_p,i (mu - lambda)(mu - conj lambda) d mu / d k_I,i` | THEORY 2415 (was "deferred validation") | TH / **V (new: T1)** | IEEE-39: 219 well-scaled mode pairs, max rel. error 3.9e-5, sign agreement 100 %. |
| E27 | Uniqueness: one PLL cannot hold two distinct modes under a delay change (2 real gains, 2 real conditions per complex pole) | THEORY 2452 (corollary) | TH / V | Not an impossibility for real-part-only targets or for other sites. |
| E28 | **Delay budget**: `k_I(h) = 0` at `h* = (1/w) atan2(w k_I, |lambda|^2 k_p + a k_I)`; for `a ≈ 0`, `h* = phi_PI(w)/w` with `phi_PI = atan(k_I/(w k_p))` | THEORY ≈2408; this session | ID / V | Classical phase budget of a PI, here at the protected mode. IEEE-39 @4.91 Hz: 8.99 ms. |
| E29 | Low-frequency (moment) transport: `lambda -> 0` gives `k_p(h) = k_p + k_I h`, `k_I(h) = k_I`, remainder factor `s^2` | this session; consistent with E46 | ID / V | This is the "first-order / Padé" delay compensation. T2/T3 show it destabilises IEEE-39 faster than no retune. |
| E30 | NEP sensitivity `d lambda/dp = - l^H Delta_p v / (l^H Delta_s v)`, `Delta_s = E + sum_i tau_i b_i c_i^T e^{-s tau_i}`; `lambda_kp = lambda * lambda_kI` | CUAD §5; THEORY; mega M6 | ID / V | 42 derivatives worst 1.8e-7; IEEE9 45 cases 2e-5. Standard calculus. |
| E31 | First-order sparse authority: `sum_i h_i^T u_i <= -d`; min-norm single-site fix `u_i* = -d h_i/||h_i||^2`; support bound `sum_{i in S} A_i >= d` | CUAD §6 | TH (linear model) / S | Cauchy–Schwarz / sorted support. Linear model only. |
| E32 | m-PLL impossibility certificate: if `eta^T(g0 - beta) > sum of m largest A_i(eta)` no ≤ m-PLL retune works | CUAD §9 | TH conditional / X | Needs rigorous remainder bounds `beta`; not evaluated on IEEE-39. |

## D. Damping operator (Beyond Nodal Damping)

| ID | Equation | Source | Kind / conf. | Notes |
|---|---|---|---|---|
| E33 | `D_H(w) = (Z(jw) + Z(jw)^H)/2`; harmonic mechanical supply `<d^T nu> = (1/2) nu^H D_H nu` | THEORY 2027; CUAD §3; IEEE9 §4 | ID / V | 60 checks, 1.2e-10. Incremental supply at SG ports; not a Lyapunov function, not a stability test. Independent of M (jwM is skew). |
| E34 | Seed formula `D_eff = D + (Pi(jw) - Pi(jw)^H)/(2jw)` | prompt seed C; mega M11 | AP / Q | Coincides with E33 only for `T = s^2 M + sD + L + Pi` with symmetric M, L. Zero-frequency Schur version is **R** (hidden block condition 3e18). |
| E35 | **Rank law**: one PLL change (gains and/or delay) gives `dZ = a v^H`, `dD_H = (a v^H + v a^H)/2` | THEORY 2027; CUAD Prop. A; IEEE9 §4 | TH / V | Sherman–Morrison inside both Schur complements. IEEE-39 60/60 (err < 5e-8); IEEE9 648/648 (1.8e-9). |
| E36 | Eigenvalues `l± = [Re(v^H a) ± sqrt(|a|^2|v|^2 - Im(v^H a)^2)]/2`, `l+ l- = -(|a|^2|v|^2 - |v^H a|^2)/4 < 0` | CUAD Prop. A; IEEE9 §4 | ID / V | Signature (+,−,0,…): anti-damping created = Cauchy–Schwarz deficit of the two network paths. Degenerate iff a ∥ v. |
| E37 | **Counting law**: m retuned PLLs give `D1 - D0 = (U V^H + V U^H)/2` with m columns ⇒ `n_-(D1) >= n_-(D0) - m`; PSD needs `m >= r` and `rank(V^H Q_-) = r` | IEEE9 §5 | TH / V | Proof: pick q in the negative subspace with `V^H q = 0`. IEEE9 @8 Hz: r = 3, m_PSD = 3, signs interval-certified (55 digits). Gain-size independent. |
| E38 | Nodal/collective split `D_node = diag(diag D_H)`, `D_cross = D_H - D_node` | physical_collective_damping | EM / S | One frozen design: +38.66 vs −737.49. |
| E39 | GSP: `L_M = M^{-1/2} L_P M^{-1/2} = U Lambda U^T`, `D^ = U^T M^{-1/2} D_H M^{-1/2} U`; `||[L,D]||_F^2 = sum (lam_k - lam_l)^2 |D^_kl|^2` | CUAD §11; prompt seed D | ID / V (trivial) | A non-zero commutator is NOT an instability criterion (CUAD §11 says so explicitly). No theorem beyond the identity. |

## E. Frequency response, moments, inertia

| ID | Equation | Source | Kind / conf. | Notes |
|---|---|---|---|---|
| E40 | Moment hierarchy H0, H1, H2 of the full linearised model; third coefficient `K3 = T_tau K_I^{-1} + K_I^{-1} T_f - K_I^{-1} K_P K_I^{-1}` | `network_moment_codesign_20261004`; Moment-Laws PDF | ID / V | 39 outputs × 3 inputs, cubic rel. error 3.6e-6. Delay enters as `k_p -> k_p - tau k_I` (cf. E29). |
| E41 | Area invariance: differences of normalised frequency areas between buses are invariant under PLL-only retuning | same | ID / V (linear) | Signed areas, not nadir/RoCoF. |
| E42 | Apparent inertia `m_app = m + 1^T D1 1 - b^T L^+ b`, `b = D0 1 - (d0/m) M 1`; `int (x_inf - x) dt = P m_app / d0^2`; `m_app < 0` ⇒ overshoot | THEORY 2136, 2173; Moment-Laws pp.24–25 | TH reduced / V | Two-node example m_app = −14. Reduced lossless class. |
| E43 | Fixed-bus area law `A_a = P(m+kappa)/d0^2 + (c_D - a)^T delta - Q/d0`; independent of the distribution of M | `nonlinear_frequency_limits_20261004` NOVELTY review (R1) | TH reduced / V | Finite-amplitude, uses two equilibria. |
| E44 | Inertia floor `m + kappa > d0^2 (e_k - c)^T Lbar_k^+ (e_k - c)` (secant-resistance Laplacian) | THEORY 1398–1463 | TH reduced / V | Necessary, not sufficient. Not validated on IEEE-39 (blocked by trajectory terms). |
| E45 | Two-node sharp infimum `inf m = 2 asin(P/2)/P` | THEORY 1579 | TH reduced / V | Re-checked by nonlinear simulation this session (area −0.005549 at m = 1.025). |
| E46 | Energy balance `E = K + E_f + E_dc`, `E' = P_m + P_dc - load - losses`; finite-time frequency-area identity with governor terms | `physical_frequency_bridge_20261004` | ID / V | IEEE-39 RHS error 3e-10 MW. |

## F. Portfolios, weak elements, older lines

| ID | Equation / statement | Source | Kind / conf. | Notes |
|---|---|---|---|---|
| E50 | Terminal closure `T_S = T_0 + E_S D_S E_S^T`, `K = E^T T_0^{-1} E`, `C_S = I + (I + D_S K_d)^{-1} D_S K_o`; Schur update when adding port j | Theory Companion v2; material_synthesis | ID / V | TX4 model (no governor, D = 0). |
| E51 | Contextual marginal `d_j(S) = alpha(S ∪ j) - alpha(S)`; bottleneck recurrence `beta(S) = min{m(S), max_i beta(S\i)}` | material_synthesis; TX4 | ID / V | Widest-path DP; classical algorithm. |
| E52 | Hidden margin `d_H = min_w sigma_min[I - R_c(jw)]` | Feedback-Dressed p.2 | definition / S | Unsigned; re-expands after crossing (document's own caveat). |
| E53 | Cycle compliance `C^T K_e^{-1} C v = kappa C^T C v`; holonomy commutator | Cycle-Space Geometry | TH lossless / S; holonomy part X | Screening quantity only. |
| E54 | Separation bound `||X||_F <= ||B||_F / sep(A_g, A_c)`; `w_r(alpha) = w_r0/sqrt(alpha)`; `d lambda/dm_i = lambda d lambda/dd_i` | Inertia as Spectral Separator | TH textbook / V | Content is numerical (re-entrant window), not the theorems. |
| E55 | Non-normal bounds `||e^{At}|| <= kappa(V) e^{alpha t}`, resolvent gain; `G_phys = sup_t sigma_max^2(C e^{At} B)` | CCTA draft; Stable-but-Unsafe | textbook / S | CCTA has no numbers. |
| E56 | Reduced design law `K_p* = 2 r L^+`, `K_i* = r^2 L^+` | prompt seed I | CJ / Q | Not found in validated form in the repo; matrix (non-local) gains; consensus-literature form. Not usable for local PLLs. |
| E57 | Scaling `|K^(p)| ∝ alpha^{-p/2}` | Two-Gate p.7 | CJ / **R** | Not reproduced (slope −1.161) in Stable-but-Unsafe p.6. |
| E58 | "lambda_c decreases with inertia loss" (H3), hosting inequality (H1), general dynamic Braess (H2) | Compendio pp.38–41 | — / **R** | Retracted or downgraded by the documents themselves. |
| E59 | Zero-frequency Schur graph damping | mega M11 | — / **R** | Hidden block singular (condition 3.1e18). |

## Structural reading of the ledger

Everything that survives in sections B–D is one fact seen from different sides:

> **A local PLL retune is a rank-one perturbation of the grid characteristic matrix (E05, E14, E35) controlled by two real numbers (E02).**

- Rank one ⇒ determinants reduce to m×m (E15–E17) and the damping change has signature (+,−) (E35–E37).
- Two real numbers ⇒ exactly one complex frequency can be interpolated (E21–E27); the remainder is `(s-lambda)(s-conj lambda) R(s)` (E25).
- The network enters only through residues `d mu/d k_I` (E26, E30) — this is where "beyond nodal" is literal.
