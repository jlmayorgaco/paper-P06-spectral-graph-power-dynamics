# POSTER PLAN (frozen by the reviewer). Agents implement ONE module per task, nothing else.

## Global rules (every agent)
1. Edit ONLY `moduleK.tex` of your task (and new figure scripts/figures under `figures/moduleK/`). Never edit main.tex, refstyle.tex, art.tex, skyline.tex, extra.tex, module1.tex or other modules. Panels 1 and the header/footer are LOCKED.
2. Write LaTeX with the Write/Edit tools only. Never generate LaTeX through a Bash heredoc or an inline Python string (backslashes get corrupted). In Python use chr(92) if you must build a backslash.
3. Coordinates are reference pixels (TikZ x=\PX, y=-\PY, origin page top-left), exactly as in module1.tex. Read module1.tex and refstyle.tex first and reuse their macros: \EquationCallout, \GreenCallout, \WarningNote, \ProvedTag, \ObservedTag, \ProposedTag, \StepNumber, fonts \CondBold \CondSemi \CondReg, colours PGreen PBlue PRed PTeal PText PMuted PGold.
4. Type sizes (physical pt): body 20–22.5, equations 24–32 (hero up to 34), captions/labels >= 17, nothing below 15.5. Do not shrink a panel to fit: shorten text or recompose.
5. Every number printed must come from a repository file; put a LaTeX comment `% source: <path>` next to it. No number from the reference image. No invented data, no screenshots of plots, no AI-generated curves. Plots: matplotlib from CSV, saved as PDF at the exact slot size (1:1), fonts 17–20 pt, palette of refstyle.tex.
6. Compile in your OWN build dir: `lualatex -interaction=nonstopmode -halt-on-error -output-directory=build_<module> main.tex` twice (from the poster folder), render `pdftoppm -r 150 -png -singlefile -x X -y Y -W W -H H build_<module>/main.pdf crop_<module>` for your panel, LOOK at the PNG, fix clipping/overlaps, repeat. Panel pixel boxes at 150 dpi: x_px = 4.651*x_ref, y_px = 4.649*y_ref.
7. Finish with a short report: what is in the panel, every number + source, any doubt. Do not commit, do not push.

## Panel boxes (reference px; body starts 3 px below the bar)
P2 x 405–756, y 119.5–451 (bar to 145.5) | P3 x 765–1149, same y
P4 x 13–396, y 459.5–793.5 (bar to 485.5) | P5 x 405–756 | P6 x 765–1148
P7 x 13–390, y 796.5–1129.5 (bar to 822.5) | P8 x 398–756 | P9 x 764–1148
Lower strip (unnumbered): MAIN CONTRIBUTIONS x 13–548, SCOPE x 555–837, REFERENCES / CONTACT x 844–1148, y 1137.5–1246.5 (bar to 1157.5).

## Module specifications (corrected; do not copy the reference image's math)
### M2 — MODEL: DAE -> TRANSVERSE STABILITY
DAE  E(x,δ)ẋ = f(x,z;δ,K,τ), 0 = g(x,z;δ) + tiny legend (x dynamic states, z network/algebraic, δ portfolio, K gains, τ delays).
Step blocks (StepNumber discs): 1 operating point f(x*,z*)=0, g(x*,z*)=0; 2 linearise [Δẋ;0]=[[f_x,f_z],[g_x,g_z]][Δx;Δz]; 3 eliminate algebra on a regular chart (g_z invertible): A = f_x − f_z g_z^{-1} g_x (Schur complement — do not call it Kron unless network-only); 4 transverse quotient A_⊥ = Zᵀ A Z, Z orthonormal basis of g⊥ where g is the rotational gauge (A g = 0).
Hero box: α_⊥(S,K) = max Re σ(A_⊥);  tags α_⊥<0 stable, α_⊥>0 unstable.
Small note: "With PLL delay: α_⊥ = rightmost root of det Δ_⊥(s), Δ(s) = sI − A_0 − B E(s) C (exact exponentials)."
Visual flow strip: PORTFOLIO → EQUILIBRIUM → DAE JACOBIAN → TRANSVERSE QUOTIENT → RIGHTMOST ROOT.
### M3 — IEEE-39: MINIMAL COLLECTIVE BLOCKER (hero result)
Data: `reports/poster/ias2026/research/bnd_h4_mechanism/results/20260925T204017_549c3c07_h4_crossmode_v3/derived/TX4_PROPER_SUBSET_CLOSURE.csv`, `.../docs/IAS26-010_ROOT_CAUSE.md` (α_⊥(H4)=+0.12700646782830968 s⁻¹, f=0.6222796695036534 Hz), `research/ias2026_last_validation/reports/GATE0_REPORT.md`. Network topology for the one-line: find the IEEE-39 branch list in the repo (e.g. reports/experiment_D/inputs or data/); draw buses as a schematic one-line in TikZ or matplotlib with 30,33,35,37 highlighted.
Content: label "SEPARATE TX4 FIXED-POLICY MODEL CONTRACT (no PLL delay in this witness)"; H4 = {30,33,35,37}; definition α_⊥(R)<0 ∀R ⊊ H4 but α_⊥(H4)=+0.1270 s⁻¹ at f=0.6223 Hz; big callout "15/15 proper subsets stable — FULL COALITION UNSTABLE"; Hasse lattice of the 16 subsets coloured by the real α_⊥ values (CSV) and/or stability-cliff plot; scope note "Canonical nominal witness, not a universal blocker." Do NOT show 40/44 ms here.
### M4 — EXACT NETWORK CLOSURE: LOCAL × COLLECTIVE
T_S(s)=T_0(s)+Σ_{i∈S}E_iΔY_i(s)E_iᵀ; ΔY_i=E_iᵀ[T_{\{i\}}(s)−T_0(s)]E_i; 𝒟=blkdiag(ΔY_i) (device-determined); 𝒦(s)=E^T T_0^{-1}(s) E = 𝒦_d+𝒦_o (network-determined); M=𝒟𝒦; Q=(I+𝒟𝒦_d)^{-1}𝒟𝒦_o.
Hero: det(I+M_S) = Π_{i∈S} det(I+M_ii) · det(I+Q_S). Visual split LOCAL RETURN | COLLECTIVE CLOSURE with arrow diagram ΔY_i × T_0^{-1} → self/cross return → Q.
Audit callout: "canonical stored boundary, f ≈ 0.7064 Hz: min_i σ_min(I+M_ii)=0.4892, σ_min(I+Q_H4)=4.39×10⁻⁸" (source `.../h4_crossmode_v3/claims/tx4_summary.json`; NEVER pair with 0.6223 Hz). Corollary det(I+Q_H)=0 ⇔ −1∈σ(Q_H).
Second evidence line (full delayed model, `research/feedback_cycle_fullmodel_20261005/FINAL_FEEDBACK_CYCLE_REPORT.md`): "Delayed full IEEE-39 at 44 ms: eig(Q) = −1 at every unstable root while all |L_i| ≥ 0.011; closest pair cycle 35–36, |1−p| = 0.008."
Hero sentence: "Nothing failed locally. The network closed the destabilising loop."
### M5 — BEYOND NODAL DAMPING: EFFECTIVE OPERATORS
A exact Schur self-energy: S_r(s)=sI−A_rr−A_rc(sI−A_cc)^{-1}A_cr; T(s)=s²M+sD_0+L+Π(s) (tag EXACT). B even/odd split Π_e, Π_o; L_eff=L+Π_e(0); D_eff=D_0+Π_o(s)/s; M_eff=M+[Π_e(s)−Π_e(0)]/s²; T=s²M_eff+sD_eff+L_eff (tag PROPOSED canonical split); sentence "M_eff is an inertial-like closed-loop operator, not physical rotor inertia." C D_H(Ω)=Herm Z(jΩ); ⟨dᵀν⟩=½ν̂^H D_H(Ω) ν̂; δD_H=½(a v^H + v a^H): one positive and one negative eigenvalue (independent a, v). Plot: real 60-case signed eigenvalues from `reports/poster/ias2026/one_mode_20261005/generated/data/TABLE_02_FINITE_DELAY_RANK_LAW.csv` (check columns), summary `experiments/poster_julia_validation_20261003/RANK_LAW_SUMMARY.json` (60 cases, max rel. error 4.85e-8). Hero: "A local change redistributes collective damping; damping at one frequency is not a stability certificate." Small PROPOSED inset D_eff=h_D(L;s)+R_D, M_eff=h_M(L;s)+R_M.
### M6 — GRID-EMBEDDED PLL GAIN SYNTHESIS
Step 1 open PLL ports: A(s;ρ)z=B(s;ρ)θ, e=C(s;ρ)z+D(s;ρ)θ ⇒ G(s;ρ)=D+CA^{-1}B; target λ=α+jω, pattern q, y=G(λ;ρ)q. Step 2 row condition (λK_{p,i}+K_{I,i})e^{−λτ_i}y_i = λ²(1+t_{f,i}λ)q_i (derivation in one line: sθ=ω, t_fsω=ξ+K_pe_d−ω, sξ=K_Ie_d). Step 3 W_i=λ²(1+t_{f,i}λ)e^{λτ_i}q_i / y_i (denominator y_i!). Step 4 K_{p,i}=Im W_i/ω, K_{I,i}=Re W_i − αK_{p,i}. Hero box K=Φ(ρ,λ,q,τ) — "The controller gain is local; its design is embedded in the full-grid mode." Delay strip: [λK_p(h)+K_I(h)]e^{−λ(τ+h)}=[λK_p+K_I]e^{−λτ}; (dμ/dh)|₀=−K_p(μ−λ)(μ−λ̄)∂μ/∂K_I — "Exact preservation of one mode does not preserve the other modes."
Evidence mini-table (full IEEE-39, exact delay, 44 ms; sources `research_gold/checks/T5_nonlinear_events.csv`, `T5c_nonlinear_events_node1Hz_44ms.csv`, `T9_A_lowfreq_44ms.csv`, poster one_mode_20261005 panel 5): no retune 8 unstable roots, events abort; protect 4.91 Hz → 0 roots, 2/5 events; protect 1.04 Hz → 0 roots, 5/5 events. Verify each value in the CSVs before printing.
### M7–M9 and lower strip
WAIT: they depend on research/nhop_ieee39_20261006 results (reviewer will issue the specs after review).
