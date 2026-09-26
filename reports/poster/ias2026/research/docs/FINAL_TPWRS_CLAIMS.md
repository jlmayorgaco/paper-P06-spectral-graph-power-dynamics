# Final TPWRS claims

Sources:
- the evidence is `docs/FINAL_TRANSACTION_THEORY_AND_EVIDENCE.md`;
- the per-claim fields are `results/FINAL_EVIDENCE_TABLE.csv`;
- the prohibited wording is `docs/FINAL_REJECTED_WORDING.md`.

This document does **not** edit the manuscript.

## 1. The strongest result we can honestly submit

| item | answer | numbers |
|---|---|---|
| **strongest theorem** | Under the exact transverse quotient (A.1), the minimal-incompatibility hypergraph `H_RHP_perp(theta)` is well defined without any small-eigenvalue cutoff. It is locally constant and changes only at transverse imaginary-axis crossings (A.2). A target is safe in **every** implementation order iff it contains no hyperedge (A.4). The most original piece is the symmetry-deflated zero-frequency port (A.9): its sign changes iff the parity of the positive real transverse eigenvalues changes. All ingredients are classical, so these are stated as propositions, not as new mathematics. | re-audit: 0 label changes over 345 229 points; port holdout 28/29, 0/445 FP |
| **strongest IEEE-39 observation** | At fixed network and dispatch, the converter Q/V policy changes *which* minimal SG→IBR replacement coalitions are transversely unstable. For the 2096.6 MW / 4270.7 MVA four-bus flagship, `kappa_perp = 4` at P4 and the hypergraph is empty at `P_inf`. | 30 / 36 / 16 distinct `H` on F7A/B/C; up to 15 `H` per `kappa`; witness contraction 4 → 3 → 2; E14 N6: 25/25 Pg-matched controls leave the inter-area family stable |
| **strongest nonlinear observation** | Within the declared model validity, finite-disturbance composability **coincides** with transverse spectral composability. No transversely stable portfolio fails in scope. The oscillatory incompatibility boundaries are **subcritical** Hopf, a hard loss of stability. | 6/6 (point, family) cases; 84/96 thresholds are guard events; `l1` = +0.0095, +0.025, +0.026 |
| **strongest planning implication** | The hypergraph is exactly the any-order-safe constraint. On the 9-candidate census the MW-optimal plan (3652.3 MW / 6499.9 MVA) is order-independent. Three of 327 stable targets are order-dependent: converting 34 last passes through an unstable state. A damped condenser of 106.0 MVA restores small-signal composability at P4. Keeping a 200 MW disturbance inside the declared converter voltage envelope needs 684–792 MVA; this is an envelope requirement, not a stability requirement. | FC10, FC12 |
| **strongest independent validation** | Out-of-sample: the preregistered port holdout and the synthetic suite. Independent implementation: ANDES reproduces the IEEE-39 network, the base inter-area mode (4.1 %) and the damping effect of the documented governors. ANDES does **not** implement the converter model, so no portfolio result is ANDES-validated. | FC02, FC14, F1/E31 |
| **strongest negative result** | The IEEE-39 four-bus coalition is **not robust to primary frequency control**. With the source's own documented TGOV1N governors (untuned), the P4 flagship is stable. Policy-dependent incompatibility survives only in a smaller region. | `alpha` +0.127 → −0.0745; plane non-empty points 74/399 (frozen 148 + 21 base-unstable); 10 vs 26 distinct `H` |
| **most important limitation** | The frozen benchmarks have no primary frequency restoration and no converter limits. The flagship conclusion depends on the first. Every finite-disturbance question reaches the omitted-limiter envelope before any in-scope instability. | C.5, C.7, G.1–G.2 |

## 2. One primary TPWRS paper story

**Working title:** *Policy-dependent minimal incompatibility in synchronous-to-
inverter replacement portfolios: a transverse spectral hypergraph for
any-order-safe retirement planning.*

1. **Stability notion.** The exact rotation / drift center subspace and the
   transverse quotient. The model-scope sentence is stated up front, and there
   is no eigenvalue cutoff.
2. **Object and propositions.**
   - `H_RHP_perp(theta)` and `kappa`: antichain, local constancy, boundary
     localization.
   - Any-order safety iff the target is hyperedge-free.
   - Classical GN boundary localization, credited as classical.
   - The zero-frequency port closure, validated on the preregistered holdout.
3. **IEEE-39 results.**
   - The policy map and witness contraction.
   - The flagship in correct units.
   - The E14 N6 control.
   - The damped-condenser mechanism.
   - The **governed replication as a conditions section**: the result holds
     without primary frequency control; with documented governors the
     structure persists in a smaller region and the P4 coalition is stabilized.
4. **Corroboration and limits.**
   - Kundur recurs.
   - IEEE-68 reproduces the benchmark, but its preregistered map is empty.
   - The ANDES scope is stated.
5. **Planning.** The any-order-safe constraint, the order-dependent targets,
   and the support rating. The envelope result is presented as a caveat.
6. **Time-domain confirmation.** Small-signal verdicts are confirmed at 32/32.
   Within scope, the finite-disturbance structure equals the spectral one, and
   the boundaries are subcritical. The finite-disturbance margins are
   scope-censored; this is stated as a limitation.

Novelty statement: the object and its measured policy dependence are new, and
so is the exact zero-frequency port closure. The tests, the hypergraph
formalism, the complexity and the topology are credited to the literature.

## 3. Nonlinear follow-up paper

**Not justified by the present data.**
- All four hard criteria (§23 A–D) fail inside model scope.
- The only structure beyond `rho_scope` is the limiter-activation order, which
  the models cannot follow.
- A nonlinear composability paper would need a new, versioned model with
  converter current limiting, PLL / ride-through logic and governors, followed
  by a new preregistered campaign. That is a new exploratory branch, and the
  brief excludes it here.

## 4. Claims for the manuscript (exact allowed wording)

| id | wording | status |
|---|---|---|
| T1 | "Stability is stated modulo the exact rotation (a gauge) and the common-frequency drift (a consequence of absent primary frequency control); no eigenvalue is discarded by magnitude." | PROVED |
| T2 | "Portfolios containing no hyperedge of `H_RHP_perp` are transversely stable, and a target is safe in every implementation order iff it contains no hyperedge." | PROVED |
| T3 | "`H_RHP_perp` changes only where a transverse eigenvalue crosses the imaginary axis; the classical return difference localizes the oscillatory boundaries." | PROVED + classical |
| T4 | "Relocating exactly the two structural center eigenvalues restores the zero-frequency determinant test; on a preregistered Kundur holdout it detected 28 of 29 real-count changes (the miss is a real-pair coalescence) with 0 false positives on 445 negative controls." | NUMERICALLY VALIDATED |
| T5 | "On the frozen IEEE-39 benchmark, at fixed network and dispatch, the converter reactive policy changes which minimal replacement coalitions are transversely unstable (30/36/16 distinct hypergraphs on three policy planes)." | BENCHMARK-SPECIFIC |
| T6 | "Replacing 2096.6 MW of active dispatch (4270.7 MVA of converter rating) at four buses is transversely unstable under one reactive policy although every proper subset is stable; under full voltage-support gain the hypergraph is empty." | BENCHMARK-SPECIFIC |
| T7 | "With the source's documented turbine governors, policy-dependent minimal incompatibility persists in a smaller region, while the four-bus coalition at this operating point is stabilized." | BENCHMARK-SPECIFIC |
| T8 | "A condenser with a damped swing mode restores composability above a point-dependent rating (106.0 MVA at the reference point); an undamped low-inertia condenser is itself unstable." | BENCHMARK-SPECIFIC |
| T9 | "Three of 327 stable targets are final-stable but not any-order safe; the MW-optimal plan is order-independent." | BENCHMARK-SPECIFIC |
| T10 | "Within the declared model validity, finite-disturbance composability coincides with transverse spectral composability for three disturbance families; the examined boundaries are subcritical Hopf bifurcations." | NUMERICALLY VALIDATED (negative) |
| T11 | "Small-signal predictions are confirmed in nonlinear phasor-domain simulation (32/32)." | NUMERICALLY VALIDATED |
| T12 | "An independent tool (ANDES) reproduces the network, the base inter-area mode and the governor damping effect; it does not implement the converter model." | INDEPENDENTLY REPRODUCED (limited) |
| T13 | "Exact minimal-witness search is combinatorial in general, a known consequence of the sparse-PCA reduction; on this benchmark exhaustive search over 16 portfolios is exact." | PROVED (credited) |
| T15 | (optional, methods) "A second-order reduced response is consistent only if the algebraic-manifold (KCL) curvature is retained; dropping it breaks the exact rotation / frequency-drift symmetry at second order and makes the second-order frequency prediction about 10³ times worse than the linear one." | NUMERICALLY VALIDATED (Taylor consistency; not new) |
| T14 | "The frozen IEEE-39 model has no primary frequency restoration and is analysed in relative / transverse coordinates. Absolute common-frequency restoration is outside this benchmark." | mandatory scope sentence |
