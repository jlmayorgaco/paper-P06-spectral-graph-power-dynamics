## 23. Negative results (consolidated)

Every negative result below was preregistered as a possible outcome and is
reported with the same prominence as a positive one.

| id | negative result | where |
|---|---|---|
| HS | α_⊥ is neither submodular nor supermodular at any of 39 policies (median 93 % one-step violations); the tracked-mode margin is no more structured (86 %) | §7.1 |
| — | Shapley/context-averaged attribution is uninformative here — dominated by rare fast-instability contexts, not the modest electromechanical reversals | §7.2 |
| — | Structural remark R2 (frozen partial of a setpoint = 0) narrowly fails its strict 10⁻⁸ bound (1.5·10⁻⁷), most likely a finite-difference artefact, not a violation of the underlying claim | §9 |
| C-H4b | none of four static topology scores predicts the effect of a topology action on α(H4) (median \|ρ\| < 0.6 throughout, several wrong-signed) | §13 |
| Q3 (E12) | witness-changing events are not reliably accompanied by a reproducible modal-support transition (3/9 vs the 7/9-and-≤20 %-controls bar) | §17 |
| GOLD-D | plan-level design never demonstrates an advantage over single-boundary tuning at the H4 scale (single always sufficed); at the V9 census scale neither method converges within budget | §15 |
| GOLD-C | no reduced model achieves both a genuine end-to-end speedup (≥ 3×) and a materially smaller state count (≤ 70 %) — the bottleneck is the nonlinear equilibrium solve, not eigenanalysis | §18–19 |
| H8 | the modal-energy/limiting-mode mismatch is real (16 % of cases) but not "systematic" by the preregistered per-policy majority rule | §20 |
| E17 | the Africano/PV material gate fails outright — no source document, feeder, metric or result table exists in the repository | §21 |
| — | static graph scores are mostly weak predictors of contextual/dynamic node quantities (median \|ρ\| often 0.1–0.6), with one partial exception (effective resistance vs total Q/V sensitivity, 0.80) | §16 |

## 24. Publication candidates

Ranked by evidentiary strength and novelty on this benchmark class.

1. **GOLD-B (dynamic vs static line sensitivity), §10.** The strongest single
   result: a 0.49 Spearman gap over the best static baseline, validated on 32
   held-out conditions with an independently cross-checked (port-form)
   sensitivity engine, and a clean mechanistic account (frozen setpoints have
   zero frozen partial; controller coordinates have frozen = total under
   SPR). Publication-ready on its own.
2. **Contextual sign reversal (GOLD-A, §7).** A clear, preregistered,
   three-tier (global / mode-tracked / electromechanical-only) demonstration
   that a fixed node ranking is provably insufficient on this benchmark, with
   an oracle-fitted ranking baseline that is *maximally* favourable to
   node-only ranking and still falls short. Robust under four uncertainty
   envelopes, including where the failing-set identity itself changes.
3. **Topology alone moves incompatibility (§13).** Both directions
   demonstrated (removal and creation) with a clean negative control (no
   static score predicts it), and a striking census-scale result (κ moves
   across its entire {1,…,5} range under admissible single actions).
4. **Weak corridors beating single-line and pure-spectral partitions (§12).**
   A structurally motivated (transformer-group) corridor construction
   outperforms graph-spectral cutsets, with a clean additivity diagnostic.
5. **The reduced-model negative result (§18–19).** A precise, mechanistic
   account of *why* eigen-stage reduction does not help on this class
   (equilibrium solving dominates), with a certificate that is accurate
   whenever it fires. Useful as a cautionary methodological result for the
   broader reduced-order-modelling literature in this application area.

## 27. Recommended next paper(s)

1. **A focused paper on contextual weakness and dynamic line ranking**
   (GOLD-A + GOLD-B), the two strongest, best-validated results, framed as
   "static and single-context rankings are demonstrably insufficient for
   next-replacement and reinforcement decisions in policy-dependent IBR
   portfolios" — this is the natural, tightly scoped successor to TX4.
2. **A topology/corridor paper** built on §13 and §12: topology alone can
   remove or create incompatibilities, static topology scores do not predict
   which, and structurally motivated (not purely spectral) corridors are the
   ones that generalize.
3. **A short methodological note on the reduced-model negative result**
   (§18–19): a worked demonstration that "reducing the linear stage without
   reducing the equilibrium solve buys nothing" is itself a useful, citable
   caution for reduced-order dynamic security assessment.
4. **Africano/PV** (GOLD-E/F) only if the source material can be obtained; the
   plan in `docs/CDW_AFRICANO_PV_BENCHMARK_PLAN.md` is ready to execute as
   preregistered, unchanged, the moment the gate can be satisfied.
5. **Nonlinear recovery vs α_⊥ improvement** (E22, deferred here) is a natural
   follow-up once a paper is scoped around the plan-level design result,
   since TX4 already established subcritical-Hopf boundaries on this
   benchmark.
