# Portfolio decision benchmark: the P4 replacement-portfolio trap (Phases 1, 9 and gates)

Date: 2026-09-11.

- Preregistration: `docs/20260911_POST_CUMULANT_VALIDATION_PREREG.md`
  (5d0b1986).
- Data:
  - `results/20260911_BASELINE_COMPARISON.csv`;
  - `results/20260911_POLICY_COUNTERFACTUAL.csv`;
  - `results/20260911_ROBUSTNESS.csv`;
  - `results/20260911_GFL_REPRODUCTION.csv`;
  - `results/PCV/PCV02–PCV06/`.
- Figure: `figures/20260911_portfolio_decision_case.{pdf,png,svg}`.

Terminology: "lower-order portfolio screening" means screening every single,
pair or triple of candidates (budget-limited subset screening). It is not
described as an industry-standard N-k procedure.

## 1. The planning problem (Phase 1)

A planner proposes to replace the four synchronous units at buses 30, 33, 35
and 37 of the IEEE 39-bus system with matched-dispatch grid-following (GFL)
plants.

- **What is matched.** The scheduled P and Q at every bus are preserved, so
  the power flow, bus voltages and ratings do not change.
- **What the experiment isolates.** Dynamic compatibility, not
  redispatch.
- **Policy.** The converters run at the frozen policy
  P4 = (g, k, t, h) = (0.03625, 1.425, 1.5, 1).
- **Portfolio size.** 2096.6 MW of dispatch and 4270.7 MVA of rating.

| | answer | evidence |
|---|---|---|
| **D1** Is the complete portfolio safe? | **No.** α⊥(H4) = +0.1270 s⁻¹ at 0.622 Hz (2 RHP eigenvalues), recomputed with the FC01 direct path. It agrees with the frozen FC03 value to 4.6e-10 (S1) and is reproduced in ANDES with the transcribed equations (+0.1270077; Phase 8). | PCV02, PCV06 |
| **D2** What is the minimum blocking coalition? | **{30, 33, 35, 37} itself.** All 15 proper subsets are stable (α⊥ from −0.214 to −0.144 s⁻¹), so H = {H4} and κ = 4. There is no smaller blocking set to remove. | PCV02 lattice; Fig. A |
| **D3** Would lower-order portfolio screening approve it incorrectly? | **Yes.** Every single, pair and triple is stable, so screening at any budget of at most 3 units approves H4. The α-space extrapolations also predict H4 stable: B3 (first-order modal sensitivity) −0.203, B4 (additive) −0.213, B5 (pairwise) −0.210, B5b (third-order) −0.184. | PCV02, PCV03; Fig. B |
| **D4** Can aggregate MW/MVA or static strength identify it? | **Not as a verdict, and not across policy.** These metrics have no stability threshold for this model. They are identical between the unstable P4 and the stable G_S / G_S2 / P_inf points: gSCR(H4) = 1.803, min SCR 2.760, max MIIF 0.461, Pg 2096.6 MW, Sn 4270.7 MVA. They do carry ranking information across portfolios (census size-4 ROC-AUC: Sn 0.81, gSCR 0.64, Pg 0.57). | PCV02 Task A/D; PCV03 |
| **D5** Does converter/control policy change the answer? | **Yes, with only g changed.** On the clean line k = 1.425, t = 1.5, h = 1, H4 is unstable for g < g* = 0.2076814045 and every subset is stable above it. At G_S (g = 0.25) α⊥(H4) = −0.0174 and at G_S2 (g = 1.0) α⊥(H4) = −0.1016. | PCV03; Fig. C |
| **D6** Which targeted intervention restores stability? | **Several documented ones, each with a computable boundary:** the converter Q/V gain (g ≥ 0.20768); the SG excitation gain scale (k ≤ 1.30463); a damped condenser (≥ 2.48 % = 106.0 MVA); documented governors (α⊥ −0.0745). No single line reinforcement up to 2× restores stability; the best, L42, gives −0.057 s⁻¹ at 2×. | PCV04; Fig. D |
| **D7** Is the phenomenon robust to reasonable parameter uncertainty? | **The phenomenon mostly survives; the exact witness does not survive machine-parameter uncertainty.** Non-composability: 0.88 (E37 per unit), 0.62 (E37 fleet-wide), 1.00 (converter stress). H = {H4}: 0.69, 0.34, 1.00. These are coverage fractions over deterministic envelopes, not probabilities. | PCV05 |

**The scientific point.** Safe proper subsets do not imply a safe complete
portfolio. Here the complete portfolio is the only unsafe set, so budget-limited
subset screening at any budget below 4 approves it.

## 2. The decision table (Phase 9)

Columns follow the Phase 9 specification. ✔ = answers the question correctly.
✘ = gives a wrong answer. — = not designed for the question (N/A, not a
failure). "enum" = requires enumerating the relevant portfolios.

| question | full eigensolver (B0) | MW/MVA screen (B1) | SCR/gSCR/MIIF (B2) | modal sensitivity (B3) | additive (B4) | pairwise (B5) | exact closure / hypergraph framework |
|---|---|---|---|---|---|---|---|
| Is each individual replacement stable? | ✔ (4 solves) | — | — | ✔ (4/4 stable predicted) | ✔ (exact by construction) | ✔ | ✔ (same answer; from one common realization) |
| Are all pairs stable? | ✔ (6 solves) | — | — | ✔ | ✔ (additive predicts all stable) | ✔ (exact by construction) | ✔ |
| Are all triples stable? | ✔ (4 solves) | — | — | ✔ | ✔ | ✔ (third-order B5b exact by construction) | ✔ |
| Is the complete portfolio stable? | **✔ (1 solve) — correctly UNSTABLE** | — (no threshold) | — (no threshold) | ✘ (−0.203) | ✘ (−0.213) | ✘ (−0.210; B5b −0.184) | ✔ (same verdict; C2 identity residual ≤ 3e-10) |
| What is the minimum blocking coalition? | ✔ with enum (16 solves) | — | — | ✘ (EMPTY) | ✘ (EMPTY) | ✘ (EMPTY) | ✔ H = {H4}, κ = 4, with any-order planning semantics (C1, proved) |
| Does policy change the blocking coalition? | ✔ with enum per policy | ✘ (policy-invariant by construction) | ✘ (policy-invariant by construction) | ✘ exact H at 2/4 prereg points (misses P4; false positives at G_S2) | ✘ exact H at 3/4 (misses P4) | ✘ exact H at 2/4 (misses P4; false positive at G_S, +0.028) | ✔ H(θ) as a policy map; the boundary g* is located by port Newton |
| Does aggregate MW/MVA identify the failure? | — | ✘ as a verdict; ranking AUC Sn 0.81, Pg 0.57 on the census | — | — | — | — | — |
| Does static strength identify the policy change? | — | — | ✘ (identical at P4 and G_S; within-portfolio AUC 0.50) | — | — | — | — |
| Can first-order modal sensitivity identify it? | — | — | — | ✘ at P4; 22/52 exact H on F10 | — | — | — |
| Can pairwise truncation identify it? | — | — | — | — | — | ✘ at P4; 38/52 exact H on F10 (third order: 50/52) | — |
| Can exact full eigenanalysis identify it? | **✔ always, if every portfolio is enumerated** | — | — | — | — | — | (the framework does not claim higher accuracy than B0) |
| What structure beyond brute-force eigenanalysis? | none by itself | — | — | — | — | — | exact reduced representation (C2); minimal incompatibility sets and any-order planning (C1); policy maps; principal-minor closure decomposition; analytical boundary derivatives (C3) |
| Which intervention moves the boundary? | ✔ by re-solving each candidate | — | ✘ (ΔgSCR line ranking Spearman −0.03) | ✔ locally (full-DAE sensitivity = port derivative; Spearman 0.94 on 12 lines) | — | — | ✔ port derivative: sign 4/4 policy coordinates, Spearman 0.94 on lines; boundary by Newton (g: 5 steps, k: 4 steps) |
| Can the boundary shift be predicted analytically? | no (numerical) | — | — | local derivative only | — | — | ✔ local derivative plus Newton to the exact boundary (\|Re s\| ≤ 5e-10). Finite-step magnitude is not linear: one-step g estimate 0.052 vs actual 0.171 |
| Is the claim independently reproduced? | ✔ Phase 8: ANDES with the transcribed equations, 32/32 | n/a | n/a | n/a | n/a | n/a | ✔ same (H, κ at P4 and G_S) — reproduction of computation, not of model adequacy |
| Is it robust to parameter uncertainty? | — | — | — | — | — | — | **phenomenon: mostly** (R1 PASS in 3 of 4 envelopes, PARTIAL fleet-wide); **exact witness: no** under machine uncertainty (R2 0.34–0.69), yes under converter stress (1.00) |

**Credit where it is due.**

- **The full eigensolver.** It determines stability exactly when every
  portfolio is enumerated; at P4 that is 16 solves. The framework's verdicts
  are identical to it by construction, so **the novelty is not "predicting
  stability more accurately than eigenanalysis".**
- **Third-order truncation (B5b).** It is the strongest non-enumerating
  predictor on the frozen F10 set: 50/52 exact H, pooled ROC-AUC 0.995, sign
  accuracy 0.996. But it needs the exact α of 15 of the 16 subsets, and it
  still misses P4.
- **Oracle-assisted closure distance (B8b).** It reaches census size-4 ROC-AUC
  0.94 but requires the true critical frequency. The non-oracle version (B8a)
  reaches 0.86.
- **Static metrics.** They rank portfolio risk moderately across portfolios
  (Sn 0.81) but, by construction, cannot see policy.

## 3. Lower-order hierarchy (Phase 5), P4, H4

| level | α space: predicted α⊥(H4) | determinant space: zero of the truncated characteristic function |
|---|---|---|
| local / first order | B3 −0.2032 (stable) | local factors only: no zero within radius 1 of λ* |
| additive singles | B4 −0.2129 (stable) | — |
| pairwise | B5 −0.2104 (stable) | order ≤ 2: Re s = −0.0148 (stable, misses) |
| third order | B5b −0.1837 (stable) | order ≤ 3: Re s = +0.1267 (flags H4; residual 2.9e-2 in \|s − λ*\|) |
| exact full closure | +0.1270 (B0) | order 4: +0.1270 (residual 2.7e-10) |

**Why the trap occurs.** In α space, every truncation extrapolates the
spectral abscissa, and the subsets' abscissae are all near −0.15 to −0.21. The
failure of the union is not visible in any subset's margin.

In determinant space the order-≤3 truncation already moves its zero to the
right half-plane. This is consistent with the connected-cumulant audit, which
found the flagship **composite**:
- ν = 0.069;
- deleting the connected four-way term moves the zero by only 0.024 s⁻¹;
- the dominant connected cluster is {33, 35}.

**κ = 4 is a minimum destabilizing cardinality; it does not imply an
irreducible four-device interaction.** The full-order Boolean coefficient is
`μ_1234 = χ_1234 + Σ_pairings χ_ab χ_cd`, and here it is dominated by products
of lower-order connected terms that jointly touch all four buses. The Boolean
statement that no truncation below |S| vanishes *at s\** remains exactly true.

## 4. Clean policy counterfactual (Phase 4)

**Identical on P4 → G_S → G_S2.** Only g changes: 0.03625 → 0.25 → 1.0. The
following are unchanged:
- the network Ybus (hash);
- the candidate locations;
- the SG and GFL ratings;
- the scheduled P and Q at 30, 33, 35 and 37;
- k, t and h;
- the equilibrium voltages (within 2.5e-8, the power-flow tolerance);
- every static metric: Pg, Sn, penetration, gSCR, min SCR, max MIIF,
  compactness.

**Changes:**
- α⊥(H4): +0.1270 → −0.0174 → −0.1016;
- the critical frequency: 0.622 → 0.711 → 0.724 Hz;
- the RHP count: 2 → 0 → 0;
- H: {H4} → ∅ → ∅;
- κ: 4 → ∞ → ∞.

The boundary is g* = 0.2076814045 (direct bisection; the frozen value is
0.20768140450).

**Conclusion.** A controller-independent static strength metric cannot, by
construction, encode a controller-policy-dependent change in H or κ. This
does **not** say SCR/gSCR is wrong. They answer a different question (network
strength) and do so consistently.

**P4 vs P_inf** (g = 1, k = 0.5) is reported as a multi-coordinate policy
contrast only; α⊥(H4) there is −0.1554.

## 5. Remediation / design (Phase 6)

No common cost is defined across these unlike interventions, and |χ|
suppression is not a target.

| intervention | initial α⊥ | result | port prediction | Newton / boundary | nonlinear TDS | independent tool |
|---|---|---|---|---|---|---|
| Q/V gain g | +0.1270 | stable for g > 0.2076814 (G_S: −0.0174) | dRe s/dg = −2.451. The one-step estimate 0.052 is 3× short of the actual 0.171, because the derivative weakens to −0.461 at g*. | 5 Newton steps to 0.2076814044 (\|Re s\| 2.4e-10) = bisection | none at G_S; the frozen G2 covers P4 (unstable) and P_inf (stable) | Phase 8 reproduces P4 and G_S verdicts |
| AVR gain scale k | +0.1270 | stable for k < 1.3046267 | dRe s/dk = +0.874. The one-step estimate −0.145 vs actual −0.120. | 4 Newton steps to 1.3046267178 = bisection 1.30462672, inside the frozen F7A bracket | none | no |
| damped condenser | +0.1270 | stable above 2.48 % (106.0 MVA) | N/A (outside the port family) | frozen threshold | frozen G2: 2 % unstable, 3 % stable | no |
| documented governors | +0.1270 | α⊥ −0.0745, H = ∅ | N/A (model change) | — | none | governor damping direction only |
| single-line reinforcement (12 frozen lines, ×1.25/1.5/2) | +0.1270 | no line restores stability up to 2×; best L42 −0.057 s⁻¹ at 2× | port = full-DAE derivative; Spearman 0.937 with the actual 1.5× change; finite-step sign 9/12 | N/A | none | line sensitivities equation-equivalent in ANDES (0a4e18c2) |

Static line proxies rank the same 12 lines with negative Spearman:
- ΔgSCR −0.03;
- 1/x −0.30;
- x-weighted betweenness −0.55.

## 6. Scientific gates

| gate | decision | basis |
|---|---|---|
| **GATE 1** P4 survives a fair comparison | **PASS** | H1 is confirmed by B0. B3, B4, B5 (and B5b) approve H4 at P4. B1, B2 and B7 are unchanged across the clean counterfactual while the verdict changes. |
| **GATE 2** Information beyond enumerating B0 | **PASS** (all five criteria met; (d) and (e) rest on frozen results) | (a) exact reduced representation: C2 identity residual ≤ 3e-10, one common realization. (b) Analytical boundary sensitivity: sign correct for g, k, t, h; Newton reaches the direct boundary within 1e-8 for g and k. (c) Spearman 0.937 ≥ 0.9 on lines. (d) Policy mapping and any-order planning (proved, frozen). (e) Structured attribution (FC18, frozen). Never claimed: higher accuracy than eigenanalysis; polynomial complexity (minimum-cardinality destabilization is NP-complete, T05). |
| **GATE 3** R1 phenomenon robust | **PASS in 3 of 4 envelopes; PARTIAL fleet-wide** | EM-u 0.88, EC 1.00, EMC 0.84, EM-f 0.62 |
| **GATE 4** R2 exact witness robust | **FAIL under E37 fleet-wide machine uncertainty; PARTIAL per unit; PASS under converter stress** | EM-f 0.34, EM-u 0.69, EMC 0.62, EC 1.00. The phenomenon, not the bus set, is robust. |
| **GATE 5** custom GFL independently reproduced | **PASS** | Phase 8, 32/32, as a reproduction of the computation with the same equations |
| **GATE 6** strongest defensible practical claim | see §7 | — |

## 7. Strongest defensible practical claim

In matched-dispatch GFL replacement, a portfolio whose every proper subset is
dynamically safe can itself be unsafe.

On IEEE-39 at the frozen policy P4:
- the four-unit portfolio {30, 33, 35, 37} is the unique minimal blocking
  coalition;
- aggregate MW/MVA and static strength metrics are unchanged when the
  converter Q/V gain alone removes the failure;
- lower-order screening (single/pair/triple, and additive / pairwise /
  third-order extrapolations) approves it.

The exact closure framework returns the same verdict as full eigenanalysis. In
addition it returns:
- the minimal coalition with any-order planning semantics;
- the policy map;
- exact boundary derivatives, which locate the restoring g and k boundaries
  by Newton in 4–5 steps.

The four-unit witness is specific to the machine data. Under a declared E37
machine envelope the phenomenon persists in 62–88 % of draws, but the witness
set changes. In every draw where H4 is unstable, raising g to 1 restores
stability (post hoc).
