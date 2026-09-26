# Canonical theory ledger (2026-09-11)

This is the single canonical claim ledger. It supersedes the status statements in
`docs/TPWRS_REMAINING_GAPS.md`, which predates the cross-tool campaign (see §6).

The per-claim matrix is `results/20260911_CLAIM_MATRIX.csv` (44 claims), generated
by `experiments/connected_cumulants/CC00_claim_matrix.py`. Each claim carries
eight separate evidence fields:

1. mathematical proof;
2. toy / synthetic validation;
3. IEEE-9;
4. IEEE-39;
5. Kundur / IEEE-68;
6. independent tool;
7. nonlinear;
8. final status.

**Status vocabulary.**

| status | meaning |
|---|---|
| PROVED | a proof exists; it may be classical and is credited |
| NUMERICALLY VALIDATED | tested beyond the case it was designed on |
| BENCHMARK-SPECIFIC | observed on one frozen model under its declared conditions |
| CONDITIONAL | holds only under a stated modelling condition, here chiefly the custom GFL or the absence of governors |
| REFUTED | the tested hypothesis does not hold |
| NOT YET TESTED | not yet tested |

**Rule.** No IEEE-39 observation is promoted to a theorem.

## 1. Sources read for this ledger

- `docs/FINAL_TPWRS_CLAIMS.md`;
- `docs/FINAL_TRANSACTION_THEORY_AND_EVIDENCE.md`;
- `results/FINAL_EVIDENCE_TABLE.csv`;
- `docs/FINAL_CROSS_TOOL_VALIDATION.md`;
- `docs/FC18_TARGETED_PORT_CHECKS.md`;
- `theory/TRANSVERSE_STABILITY_QUOTIENT.md`;
- `theory/PRINCIPAL_MINOR_PORTFOLIO_STRUCTURE.md`;
- `theory/FINAL_COMBINATORIAL_THEOREMS.md`;
- `src/ibr_cycles/cycles/cumulants.py` and `tests/test_cumulants.py` (infinitesimal cumulants; unchanged);
- the Monte Carlo summaries `results/paper_mc/MC01–MC04`;
- this campaign: `theory/FINITE_AMPLITUDE_CONNECTED_CUMULANTS.md` and `results/CC/CC02–CC06`.

## 2. Mathematically proved (all ingredients classical unless marked)

| id | statement | where |
|---|---|---|
| T01 | Exact transverse quotient: `C = span{R_x, w}` (Jordan chain; `D = 0`, no governor); `det(sI−A) = s² det(sI−A_perp)`; exact nonlinear transverse system | `TRANSVERSE_STABILITY_QUOTIENT.md` |
| T02 | `H` is an antichain; hyperedge-free ⇒ stable; `H` changes only at transverse axis crossings; band-limited emptiness does not certify | F2C/F2E |
| T03 | Planning hierarchy C ⇒ B ⇒ A; any-order safe iff hyperedge-free | `FINAL_COMBINATORIAL_THEOREMS` §3 |
| T04 | Metzler monotone class ⇒ heredity, A = B = C | ibid. §2 (Perron–Frobenius) |
| T05 | NP-completeness of minimum-cardinality destabilization (sparse-PCA corollary; no priority) | ibid. §1 |
| T06 | Resilience complex: Stanley–Reisner organization (tautological) | ibid. §4 |
| T07 | Principal-minor / Möbius structure and its invariants | `PRINCIPAL_MINOR…` |
| T08 | Symmetry-deflated zero-frequency port closure | `SYMMETRY_DEFLECTED_PORT_CLOSURE.md` |
| C01 | **new:** finite-amplitude connected cumulant, moment–cumulant inversion, uniqueness | `FINITE_AMPLITUDE_CONNECTED_CUMULANTS.md` §1 |
| C02 | **new:** factorization ⇔ vanishing of all mixed cumulants; top cumulant zero does not imply factorization | §2 |
| C03 | **new (elementary):** minimal coalitions are cumulant-connected at their characteristic zero; boundary identity `chi_H = −C_H` | §2 |
| C04 | **new:** scalar ports, `chi_S = (−1)^{|S|−1} Σ spanning cycles` for `|S| ≥ 2` (fails at `|S| = 1`) | §3.1 |
| C05 | **new:** 2×2 ports, block-connected expansion, `chi_ab = −tr(Q_ab Q_ba) + det Q_ab det Q_ba`; the naive holonomy formula is refuted | §3.2 |
| C06 | **new:** `d_alg`, `d_conn` and κ are logically independent (five constructions) | §4 |

The port derivative `ds*/da` (T10) and the Sylvester ratio identity in T09 are
also proved; their IEEE-39 content is ranked below.

## 3. Toy / synthetic validation only

These statements have no benchmark evidence beyond what is listed.

- **MC01 S1–S7 (paper preregistration 7a772808 / f2946257).**
  - S1 transverse identity: 1000/1000.
  - S2 hypergraph axioms: 2000 families, 0 violations.
  - S3 crossings: 1044 located.
  - S4 Metzler: 0 violations in 1000; the signed control violates in 997.
  - S5 clique reduction: 500/500.
  - S6 principal minors: ≤ 3.4e-15.
  - S7 port parity: flip = parity in 817/817 kept.
- **CC02, the new layer.**
  - C1–C3 residuals ≤ 2e-13.
  - T3 same-κ = 4 pair: composite path (ν = 0) against connected ring (ν = 1/2).
  - C4 order table.
  - Connectivity theorem: 0 violations in 184 random minimal coalitions.
- **The only IEEE-9 evidence (E01) is a different question**: a created mode in
  the collective factor under three stress actions. It corroborates the
  collective-factor bookkeeping, not any replacement-portfolio claim.

## 4. Established on IEEE-39, benchmark-specific or conditional

These results hold on the frozen model. The model has:

- the harmonized first-order AVR;
- matched dispatch and the core 30/33/35/37;
- the leaky Q/V converter;
- `D = 0` and no governor.

| id | result | status |
|---|---|---|
| B01 | Policy-dependent minimal incompatibility (30/36/16 distinct H; witness 4 → 3 → 2 → 3 → 4 → ∅) | CONDITIONAL (custom GFL, no governor); Kundur corroborates; IEEE-68 preregistered map empty |
| B02 | Flagship κ = 4 at P4 (α +0.127), empty at P_inf | CONDITIONAL (governors remove it) |
| B03 | Governed replication: policy dependence smaller; P4 coalition stabilized | BENCHMARK-SPECIFIC |
| B04–B06 | Inter-area mechanism; damped condenser 106.0 MVA; 3/327 order-dependent targets | BENCHMARK-SPECIFIC |
| B07–B10 | Nonlinear = spectral inside scope (negative); subcritical Hopf; TDS 32/32; KCL curvature | NUMERICALLY VALIDATED (IEEE-39) |
| T09 | Network-closure anatomy: at 10/10 frozen boundaries the zero comes only from the full-order closure term | BENCHMARK-SPECIFIC (the identity is classical) |
| T10 | Port boundary derivative = DAE eigenvalue sensitivity | NUMERICALLY VALIDATED, and independently reproduced in the equation-equivalent configuration |
| C07 | **Flagship boundary is COMPOSITE** (ν = 0.069 < 0.10; deletion displacement 0.024 s⁻¹) | BENCHMARK-SPECIFIC, threshold-sensitive (CONNECTED at 0.05) |
| C08 | Triple boundaries are CONNECTED at 0.10 (ν 0.14–0.18) and COMPOSITE at 0.20 | BENCHMARK-SPECIFIC, threshold-sensitive |
| C09 | `|chi_P(jω_P)|` separates failing from Pg-matched stable portfolios, AUC 0.952; it equals the closure distance (0.944) and beats every non-definitional diagnostic (best 0.812) | BENCHMARK-SPECIFIC (one policy, n = 35) |

## 5. Independently reproduced

The source is `docs/FINAL_CROSS_TOOL_VALIDATION.md` (commit 0a4e18c2).

| id | what | tool | status |
|---|---|---|---|
| I01 | Network and equilibrium | ANDES 2.0.0 (Ybus 1.1e-13; bus voltages 1.7e-7, fully explained by ANDES's +1e-8 r/x regularization); pandapower 3.4.0 (3e-14) **after the tap-convention translation fix** | CONDITIONAL (the pandapower part depends on accepting the fix) |
| I02 | Equation-equivalent synchronous dynamics (R2/R3, 12 subsets) | ANDES: α 1.5e-6, band 1.6e-6 | NUMERICALLY VALIDATED |
| I03 | Equation-equivalent branch sensitivities of the critical mode, and the port derivative | ANDES E1: 12/12 signs, ρ 1.000 (port 11/12) | NUMERICALLY VALIDATED |
| I04 | Anything involving the custom 10-state GFL: H, κ, P4, P_inf, the GFL line ranking | none | **NOT YET TESTED** independently |

**Conditional on the custom GFL:** B01, B02, B03, B05, B06, C07–C09, the
policy maps, the witness contraction and the GFL branch ranking. The last is
REFUTED as transferable: 7/12 signs across converter models, N09.

## 6. Reconciliation of `docs/TPWRS_REMAINING_GAPS.md`

That document predates the final cross-tool campaign. Its entries now read as
follows.

| gap in the old document | current status | evidence |
|---|---|---|
| 1. "Independent implementation, unresolved (blocking)"; ANDES disagrees on the sign of the machine-removal effect | **No longer blocking.** The old disagreement was between two *different models*: the six-state GENROU and IEEEX1 against the harmonized model. In the equation-equivalent F1 configuration both tools agree: the flagship is +0.329 in both, and 12/12 RHP counts match. Network, equilibrium, synchronous dynamics and branch sensitivities are reproduced. **Still open, not blocking:** an equation-equivalent implementation of the custom GFL (Phase 7, §7). | 0a4e18c2; F1 |
| 2. Machine-parameter sensitivity (E37, 418/1000) | Open; not re-examined in this campaign. It is a reviewer question, not a scientific blocker. The claims are already stated as benchmark-specific. | — |
| 3. Repair characterised at one point | Partly addressed. The port-guided Newton retune reaches the boundary to 2e-10 (FC18); E35 gives 117/240. The adaptive per-point retune is still not run. | FC18 §3 |
| 4. Dispatch dependence map | Open; not a blocker for the policy-dependence claim, which is at fixed dispatch. | — |
| 5. Bus-30 audit underpowered | **Closed, negative:** p = 0.18 after adjustment (UC02). | UC02 |
| 6. Extended precision | Open, low risk. The floating-point bound is 2.8e-10 against the smallest margin. | E40 |
| 7. Closure metric grid-limited | **Superseded:** labels now come from the transverse eigenvalues themselves (FC01 re-audit, 0 label changes), not from the gridded metric. | FC01 |
| 8. One benchmark, one topology | Partly addressed. Kundur corroborates. The IEEE-68 preregistered map is empty. The limitation remains and is stated. | D, G3 |

## 7. Phase 7 (independent GFL reproduction): feasibility, not executed

- **The model is small.** The custom GFL (`ieee39_devices.GridFollowingConverter`)
  has 10 smooth states plus the leaky Q/V state. Its outputs are a current
  injection and PLL, filter, PI and L-filter dynamics. Nothing in it needs
  retuning.
- **The obstacle is the frozen tooling.**
  - ANDES registers models in `andes/models/__init__.py` of the vendored tree.
  - The frozen `tx3-andes` environment runs from that tree (tag v2.0.0).
  - ANDES regenerates its numerical code cache for every model.
  - Adding a model card therefore touches frozen tooling. A correct execution
    needs a **separate** ANDES copy and a separate environment with its own code
    cache.
- **What a correct execution needs:**
  - a device-level equation-equivalence unit test: ANDES residuals against
    `GridFollowingConverter.derivatives/injection` at random states;
  - then the small frozen set: base, singletons, pairs, triples, flagship,
    P4, P_inf and the FC18 switch points.
- **Status: NOT YET TESTED.** Desirable, not a blocker for the present manuscript.

## 8. Negative results preserved verbatim

| id | negative result | evidence |
|---|---|---|
| N01 | Simple-cycle causality **failed** | E18 / UC02. Now also anti-predictive in CC04: the cycle score has AUC 0.084 (stable portfolios have *larger* single-cycle holonomies). The theory explains why: single cycles miss the multi-cycle remainder (39–65 % of `chi` at the boundaries) and the cancellations between cycles. |
| N02 | Static topology or size alone **failed** | UC02 Pg AUC 0.57; CC04 Pg 0.11, MVA 0.72, gSCR 0.68, min SCR 0.47 |
| N03 | Order 4 is **model- and policy-dependent**, not structural | witness contraction; governors; IEEE-68 empty; the flagship boundary is itself composite (C07); C06 proves the orders independent |
| N04 | **Governors change the witness** | FC03: the P4 flagship is stable with TGOV1N |
| N05 | **ANDES does not implement the custom GFL** | §5, §7; cross-model GFL line ranking 7/12 |
| N06–N08 | Rank does not explain κ; the walk expansion is non-executable; certificates failed; no low-degree Boolean representation | FC18, FC08/09, BC02b |
| C10 | Dominant connected support tracks the witness: **refuted** (H5) | CC05 |
| C11 | The most stabilizing intervention is the one that most suppresses `|chi_H|`: **refuted** for lines (H6); the policy agreement is not mechanistic | CC06 |
