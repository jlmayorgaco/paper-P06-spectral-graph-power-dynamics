# 03 — Gold candidates (five, ranked)

Scores 0–10. "Known-risk" is the risk that the result is already in the literature (10 = almost certainly known). "Ready" = probability it is poster-ready now.

| # | Candidate | Novelty | Depth | Physical insight | Surprise | IAS relevance | Defensible | Experiment ready | Fits one poster | Known-risk | Ready |
|---|---|---|---|---|---|---|---|---|---|---|---|
| G1 | **One PLL, one mode**: interpolation law for delay compensation, with exact toll and budget | 6 | 6 | 9 | 8 | 9 | 9 | 9 | 9 | 5 | 9 |
| G2 | Signature and counting laws for the grid damping operator | 6 | 6 | 8 | 7 | 7 | 8 | 7 | 8 | 5 | 7 |
| G3 | Frequency-area conservation and the network inertia floor | 5 | 7 | 7 | 7 | 7 | 6 | 3 | 7 | 7 | 3 |
| G4 | Collective (subset-stable) instability with determinant closure | 4 | 5 | 6 | 7 | 7 | 6 | 6 | 7 | 7 | 6 |
| G5 | Controller-mediated reinforcement reversal ("stronger grid, less stable") | 3 | 5 | 7 | 8 | 7 | 4 | 3 | 7 | 8 | 3 |

---

## G1 — One PLL, one mode (recommended)

**One-sentence claim.** When a PLL's delay changes, its two PI gains can hide that delay from exactly one oscillation mode of the grid; the gains are closed-form and need no network model, every other mode pays an exact toll, and choosing the wrong mode (the low-frequency one that first-order compensation protects) destabilises IEEE-39 sooner than doing nothing.

**Central equations.**

```
k_p(h) λ + k_I(h) = (k_p λ + k_I) e^{λh}                                   (transport, exact)
κ(s;h) − κ(s;0) = −(s−λ)(s−λ̄) e^{−sτ} ∫_0^h k_p(t) e^{−st} dt              (remainder, exact)
dμ/dh = −k_p (μ−λ)(μ−λ̄) ∂μ/∂k_I                                            (toll on any other mode μ)
h* = φ_PI(ω)/ω ,  φ_PI = atan(k_I /(ω k_p))                                 (budget, α ≈ 0)
```

**Why it is non-trivial.** The first line is ordinary lead compensation. The content is in the other three: the error left behind is not "small", it is a known function that vanishes at λ and grows like the squared spectral distance; its effect on the grid is that function times a *collective residue* `∂μ/∂k_I` that only the network supplies; and the limit of what a PI can hide is its own phase lag. Because the rule is local, it composes exactly over all sites and over heterogeneous delays.

**What already exists.** Dominant-pole placement for PI/PID with delay (e.g. [Halder et al.](https://ore.exeter.ac.uk/repository/handle/10871/31228?show=full), [discrete PI-PR designs](https://research.itu.edu.tr/en/publications/design-of-discrete-pi-prsup2sup-controllers-for-time-delayed-syst/)); delay effects on grid-following converters ([Aalborg study](https://vbn.aau.dk/ws/files/492344751/Impact_of_Digital_Control_Delay_on_Stability_of_Grid_Following_Converters.pdf)); PLL retuning guidelines in multi-converter grids ([Huang et al.](https://arxiv.org/pdf/1903.05489)); Smith-predictor compensation of PLL delay ([2025 paper](https://www.citedrive.com/en/discovery/overcoming-pll-delay-constraints-in-virtual-inertia-control-a-smith-predictor-approach/)).

**What is ours.** (i) the exact finite remainder with the quadratic factor and the resulting toll identity; (ii) exact composition over sites and heterogeneous delays without a network model; (iii) the demonstration, with full-spectrum counts on a 204-state IEEE-39 DDE, that moment-preserving compensation is worse than none while cluster-mode protection holds the margin up to the closed-form budget; (iv) prediction of that outcome from one eigen-analysis (Spearman 0.99).

**Current evidence.** T1: toll identity on 219 mode pairs, max rel. error 3.9e-5. T2/T3: no retune → 4 unstable roots at 42 ms, 8 at 44 ms; protect 0.05 Hz mode → 4 and 12; protect 4.91 Hz mode → 0 unstable and 0 beyond the −0.05 guard at 42, 44, 46, 48, 48.9 ms; budget formula gives 48.99 ms. T4: prediction quality. IEEE-9 report: transport on three PLLs gives α = −0.131 and 72/72 sampled robustness scenarios. Earlier independent negative: moment cancellation → +0.49 s⁻¹ at 41 ms.

**Missing proof.** None for the identities. Missing evidence: nonlinear events of the transported design (T5, running); an interval-certified count; a second base design or network.

**Killer experiment.** T5: five frozen events at 44 ms with transported gains (expected pass) versus original gains (expected oscillatory failure).

**Reviewer attack.** "This is phase-lead compensation; the identity is three lines of algebra."

**Response.** Agreed, and the poster says so. The contribution is what the three lines imply and what was measured: the standard low-frequency compensation is the wrong interpolation point for an inverter-dominated grid, the price of any choice is computable in closed form from one eigen-analysis, and the limit is a phase budget. No novelty of the algebraic tools is claimed.

---

## G2 — Signature and counting laws for the damping operator

**Claim.** Retuning one PLL changes the grid's harmonic damping operator by a rank-two term with one positive and one negative eigenvalue; removing r anti-damped directions therefore needs at least r PLLs, whatever the gain size.

**Central equations.** `δD_H = (a vᴴ + v aᴴ)/2`, `ℓ₊ℓ₋ = −(|a|²|v|² − |vᴴa|²)/4 ≤ 0`, `n₋(D¹) ≥ n₋(D⁰) − m`.

**Why non-trivial.** Gain-independent impossibility; "damping" stops being a nodal scalar.

**Exists.** Inertia of rank-one Hermitian parts; Weyl inequalities; passivity enforcement by low-rank perturbation.

**Ours.** Specialisation to physical SG torque–speed ports with exact Schur reduction; interval-certified minimal support on a 9-bus bank (m_PSD = 3).

**Evidence.** IEEE-39 60/60; IEEE-9 648/648; interval signs.

**Missing.** Why an engineer needs `D_H ⪰ 0`: two PLLs already stabilise the 9-bus bank. Evaluation on IEEE-39.

**Killer experiment.** Count negative directions of `D_H` on IEEE-39 across 1–10 Hz and show m < r retunes cannot remove them while robustness tracks PSD.

**Attack.** "Positivity at one frequency is neither necessary nor sufficient for stability." **Response.** True; keep it as the damping-domain face of G1, not as the headline.

---

## G3 — Frequency-area conservation and network inertia floor

**Claim.** The signed frequency area at every bus depends on total inertia only; a localised step forces a crossing unless `m + κ > d0²(e_k−c)ᵀL̄_k⁺(e_k−c)`.

**Why non-trivial.** A conservation law: placement and zero-net-energy control cannot beat it. Sharp two-node infimum `2 asin(P/2)/P`.

**Exists.** Poolla–Bolognani–Dörfler (inertia placement), Paganini–Mallada (aggregate response), Aalipour–Das (integral constraints).

**Ours.** Exact finite-amplitude endpoint form with secant resistance.

**Evidence.** Proofs; nonlinear two-node simulation to six decimals.

**Missing.** IEEE-39 transfer (ZIP loads, losses, governor limits); areas are not nadir or RoCoF.

**Attack.** "First-moment facts are classical and an area is not a grid-code quantity." **Response.** Partly conceded; not ready as a headline.

---

## G4 — Collective instability with determinant closure

**Claim.** Stability of every proper subset of interventions does not imply stability of the set (H4; triple of lines), and the joint factor is `det(I − P)`.

**Exists.** Decentralised-control interaction measures; Coletta–Jacquod.

**Ours.** Exhaustive censuses (512 portfolios; 46 singles + 1,035 pairs), a recovery by adding a fifth GFL, order-dependence (2 of 332 targets).

**Missing.** A mechanism that predicts which sets fail; H4 lives in the TX4 model without governors.

**Attack.** "A witness is not a theorem." **Response.** Conceded; this was the previous poster and should now be background.

---

## G5 — Controller-mediated reinforcement reversal

**Claim.** Strengthening a line can reduce closed-loop damping because the controller self-energy term outweighs the stiffness gain.

**Evidence.** Surrogate: λ₂ +0.69 %, damping 2.89 % → 0.19 %. Full ANDES: line 20–34 +1 % moves ζ by −8e-6; 0 of 180 line additions destabilise at nominal PLL gains.

**Attack.** "Known Braess-type effect; on the full model it is negligible; the surrogate was tuned." **Response.** No good answer today. Do not use.
