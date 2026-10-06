# 02 — Claims / evidence matrix

Statuses: **VERIFIED**, **STRONGLY SUPPORTED**, **SUPPORTED**, **EXPLORATORY**, **UNVERIFIED**, **FALSE / RETIRED**.
Novelty: KT-NA = known theory, new application; KC-NS = known components, possibly new synthesis; NEW? = apparently new, literature check required; NN = not novel.
Equation IDs refer to `01_equation_ledger.md`. New checks T1–T5 are in `research_gold/checks/` (scripts + CSV).

## A. Local PLL retuning under delay (recommended core)

| # | Claim | Basis | Evidence | Counter-evidence / limits | Status | Novelty | Needed before the poster |
|---|---|---|---|---|---|---|---|
| C1 | A PI + delay loop factor can be held fixed at exactly one complex frequency λ when the delay changes; the gains are closed-form and need no network model (E23) | algebra, 2 real unknowns | IEEE9: 45 cases 3.9e-14; IEEE-39 (T2): protected pole drift 0.0 with all 10 PLLs transported together, h up to 10 ms | Positivity/gain bounds limit h (C4) | VERIFIED | KT-NA (lead compensation / dominant-pole placement with delay) | nothing |
| C2 | The remainder is exact: `Δκ(s) = -(s-λ)(s-λ̄) e^{-sτ} ∫_0^h k_p(t)e^{-st}dt`; hence `dμ/dh = -k_p (μ-λ)(μ-λ̄) ∂μ/∂k_I` (E25, E26) | derivative of E23; THEORY 2415 | **T1 (new)**: IEEE-39, 219 well-scaled mode pairs, 3 sites, max rel. error 3.9e-5, sign agreement 100 % | 48 weakly sensitive pairs are below finite-difference resolution (|dμ/dh| < 0.05); excluded, not failed | VERIFIED | KC-NS / NEW? (identity is elementary; not found stated) | cite as "exact identity", not "new theory" |
| C3 | Transport with a common λ composes exactly over any set of sites and any heterogeneous delay changes h_i | each loop factor unchanged at λ | T2: drift 0.0 for 10 sites; contrast with decay-only matching, which does not compose (C12) | Mode shape preserved except PLL integrator states | VERIFIED | KC-NS | nothing |
| C4 | The hideable delay is bounded by the PI's own phase lag at the protected mode: `h* = atan2(ω k_I, |λ|²k_p + α k_I)/ω` (E28) | sign of k_I(h) | IEEE-39 @4.91 Hz: h* = 8.99 ms; T3: 0 roots right of the −0.05 guard at 42, 44, 46, 48, 48.9 ms | With the frozen bound k_I ≥ 61.7 the usable range ends near 46.7 ms | VERIFIED (spectrum) | NN as a formula (phase budget); application is new | state it as a phase budget |
| C5 | **Which mode is protected decides everything.** IEEE-39, uniform delay 40 → 40+h ms at all PLLs: no retune → 4 unstable roots at 42 ms, 8 at 44 ms; protect the 0.05 Hz mode (≡ moment matching `k_p += h k_I`) → 4 at 42 ms, 12 at 44 ms (**worse than nothing**); protect the 4.91 Hz mode → 0 unstable, 0 beyond the guard through 48.9 ms | E25–E29 | **T3 (new)**: banded argument-principle count on the full 204×204 `det Δ`, validated against `eig` at τ = 0 (0 and 1 = gauge, exact) | Floating-point count, \|Im\| < 95 Hz; not interval-certified. One network, one base design | STRONGLY SUPPORTED | KC-NS | **T5** nonlinear events (running) |
| C6 | The toll formula predicts the outcome from one baseline eigen-analysis: for 11 candidate protected modes, predicted vs exact worst root after +4 ms has Spearman 0.99; per-case correlation ≥ 0.998 | E26 summed over sites | **T4 (new)** | First order; 11 catalogued modes | STRONGLY SUPPORTED | KC-NS | nothing |
| C7 | Low-frequency (moment-preserving) compensation and modal protection are the same interpolation at different points (λ → 0 vs λ = λ_c); one PI cannot do both | E27, E29, E40 | T2/T3; independent earlier negative: sitewise moment cancellation → +0.4907 s⁻¹ at 41 ms (`network_moment_codesign`); **T5**: node 4.91 Hz at 44 ms (k_I 246.7 → 137.8) is stable but only 2/5 events meet the guards (0.505 Hz; SG slack ≈ 0) | Tension shown on one system | STRONGLY SUPPORTED | KC-NS | — |
| C7b | The node is a design choice: a node at the heavily damped 1.04 Hz mode keeps the spectrum safe (0 unstable, 0 beyond guard at 42–48 ms, T6) while spending less k_I (213 at 44 ms) and passes **5/5 nonlinear events at 44 ms** (T5c); node 4.91 Hz passes 5/5 at 42 ms (T5b) | E23–E28 | T5b, T5c, T6 | One base design; node chosen among catalogued modes, not optimised; LP-based node selection failed under extrapolation (T7) | STRONGLY SUPPORTED | KC-NS | second base design; 46–48 ms events |
| C7c | Event guards with 0.5 s windows do not detect a slowly growing 4.8 Hz instability: original gains at 42 ms meet the guards for 30 s although 4 roots are unstable (+0.096 s⁻¹) | T3 vs T5b | one event, 30 s | — | VERIFIED (observation) | — | always report root counts with event results |
| C8 | One PLL cannot hold two distinct modes while its delay changes (E27) | 2 real dof | THEORY corollary; IEEE9 Fig. 3 (assigned pole −0.30, another pair +0.567) | Not an impossibility for real-part-only targets or multi-site action | VERIFIED | NN (dof count) | nothing |

## B. Damping operator (supporting panel)

| # | Claim | Basis | Evidence | Limits | Status | Novelty | Needed |
|---|---|---|---|---|---|---|---|
| C9 | One PLL change gives `δD_H` of signature (+,−): never pure damping unless the two network paths are collinear (E35, E36) | Sherman–Morrison in Schur complement | IEEE-39 60/60 (5e-8); IEEE9 648/648 (1.8e-9) | Harmonic supply at one frequency; not a stability statement | VERIFIED | KC-NS (rank-one Hermitian inertia is textbook; port specialisation ours) | nothing |
| C10 | Counting law: r anti-damped directions at ω need ≥ r retuned PLLs, for any gain size (E37) | subspace argument | IEEE9 @8 Hz: r = 3, m_PSD = 3, interval-certified signs | Two PLLs already stabilise that system: PSD at one frequency is a different, stronger contract. Not evaluated on IEEE-39 | VERIFIED (IEEE9) | KC-NS (Weyl/inertia; passivity-enforcement flavour) | optional IEEE-39 evaluation |
| C11 | Nodal vs collective damping terms can have opposite sign (+38.7 vs −737.5) | E38 | one frozen design | single design | SUPPORTED | — | keep as a remark |

## C. Interaction of simultaneous retunes

| # | Claim | Basis | Evidence | Limits | Status | Novelty | Needed |
|---|---|---|---|---|---|---|---|
| C12 | Decay-only matched retunes at buses 30 and 37 do not compose: −0.0657 each, −0.0414 together; I = +0.0243 ± 4e-7 | E16, E17 | interval root boxes; 5/5 events after repair | One pair; 2 of 88 resolved cases violate (317 unresolved); pole stays stable; quadratic predictor also detects it | VERIFIED (one case) | KT-NA | present as witness only |
| C13 | Full-pole matching composes (0/230 violations) | E21/E23 | `COMPENSATED_SUMMARY.json` | exploratory family | SUPPORTED | KT-NA | explained by C3 |
| C14 | Fixed-point repair `t_{k+1} = max(0, I(t_k) - b)` | scalar contraction | 11 updates; root enclosed | sampled slope ≤ 0.385 only | SUPPORTED | KT-NA | — |

## D. Frequency response / inertia

| # | Claim | Basis | Evidence | Limits | Status | Novelty | Needed |
|---|---|---|---|---|---|---|---|
| C15 | Bus frequency areas depend on total inertia only; a network floor `m + κ > d0²(e_k-c)ᵀL̄⁺(e_k-c)` (E43–E45) | area/moment identities | proofs; nonlinear two-node simulation matches to 6 decimals | Reduced lossless class; IEEE-39 transfer blocked; areas are not nadir/RoCoF; prior art: Poolla–Bolognani–Dörfler, Paganini–Mallada, Aalipour–Das | VERIFIED (reduced) | KC-NS | uniform trajectory bounds on IEEE-39 |
| C16 | PLL-only retuning cannot change inter-bus area differences (E41) | moment identity | 39×3 coefficient checks | linear, signed areas | VERIFIED (linear) | KC-NS | — |
| C17 | 94.2 % GFL candidate is unsafe (bus-8 RoCoF 0.513, bus-29 frequency 0.589) | simulation | `physical_frequency_bridge` | — | VERIFIED (negative) | — | keep as retraction |

## E. Older lines

| # | Claim | Evidence | Limits | Status | Novelty |
|---|---|---|---|---|---|
| C18 | Möbius/Shapley/`log det` attributions are physical invariants | counter-example: `T̃ = exp(h(a)/n) T` changes the mixed term without moving poles (`material_synthesis` §2.1) | — | **FALSE / RETIRED** as invariants | NN |
| C19 | H4 = {30,33,35,37}: all 15 proper subsets stable, full set unstable (+0.127 s⁻¹); adding bus 34 restores stability | TX4 census (512 configs), Python/Julia agreement | TX4 model (no governor, D = 0); a witness, mechanism not derived; SG-dominated mode | VERIFIED (witness) | KT-NA |
| C20 | Triple {19–20, 19–33, 16–19} unstable though all 46 singles and 1,035 pairs are stable | Who-Moved-Mode p.4 | reduced surrogate; no artefact in the repo | SUPPORTED | KT-NA |
| C21 | "Stronger grid, less stable": line 19–20 ×1.8 raises λ₂ +0.69 % while damping falls 2.89 % → 0.19 % | Feedback-Dressed p.6 | stylised surrogate tuned to expose the effect; full ANDES counterpart (line 20–34, +1 %) moves ζ by −8e-6 only; 0 of 180 line additions destabilise at nominal PLL gains; Coletta–Jacquod 2016 prior art | SUPPORTED (surrogate) / effect tiny on full model | NN as a phenomenon |
| C22 | Line 26–29 alone +1.45 %, retune alone +5.39 %, both −10.76 % | Two-Gate p.6 | PDF only, no data; margin loss, no crossing | UNVERIFIED in this repo | — |
| C23 | Re-entrant inertia window; sign-opposite inertia steering | Inertia-Separator pp.4–5; ANDES +1 % check (bus 30: −6.4e-4, bus 37: +4.2e-4) | reduced stress study; ANDES effect small | SUPPORTED | KT-NA |
| C24 | Non-normal transient growth explains failures | Stable-but-Unsafe: 4.54× gain gap at matched poles | its own negative controls: no topology-induced non-modal reversal | SUPPORTED (mechanism) | NN (Trefethen) |
| C25 | CDW weak nodes: GOLD-A/B pass, C/D fail informatively; IEEE-68 replication = Case B (reversal 0/16) | CDW campaign | not replicated on second network | SUPPORTED / negative on replication | — |
| C26 | Commutator `[L, D_eff]` predicts instability | — | CUAD §11: explicitly not a criterion | **FALSE / RETIRED** as a criterion | NN |
| C27 | `K_p* = 2rL⁺`, `K_i* = r²L⁺` design law | — | non-local matrix gains; no validated derivation | UNVERIFIED | NN (consensus form) |
| C28 | Maximum SG→GFL replacement / delay frontier / heterogeneous-delay ranking | mega experiment | M9–M17 blocked | UNVERIFIED | — |
| C29 | `α^{-p/2}` scaling; H1–H3 of the Compendio; zero-frequency Schur damping | — | failed replication / retracted | **FALSE / RETIRED** | — |
