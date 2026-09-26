# CDW novelty boundary (H20)

The boundary is based on the focused 2020–2026 search in
`results/CDW_LITERATURE_GAP_MATRIX.csv`: 57 unique entries, every DOI or arXiv id
retrieved and checked during the search.

- **Scope.** Journal and conference records reachable through Crossref and publisher
  pages, plus arXiv.
- **Limit.** IEEE Xplore full-text search was not available.
- **Consequence.** All absence statements below read "to our knowledge, within this
  search", and "first" is not used.

## 1. Not novel by itself (credit, never claim)

| ingredient | prior work (examples in the matrix) |
|---|---|
| Eigenvalue sensitivity with the operating point re-solved (total, implicit-function derivative) | Smed 1993 (TPWRS), feasible eigenvalue sensitivity; Nam et al. 2000 (TPWRS), IEEE-39, line reactance, frozen vs re-solved load flow in the voltage-stability part; Mendoza-Armenta & Dobson 2016 (TPWRS), redispatch; Li et al. 2019 (EPSR), operational parameters; Yao, Roy, Mathieu 2023 (SEGAN); Joswig-Jones et al. 2026 (arXiv 2607.28764, not peer reviewed), power-flow-chained sensitivity to line admittances with GFL/GFM |
| The same line reinforcement can help or hurt linear stability depending on context | Coletta & Jacquod 2016 (PRE), Braess paradox in oscillator networks; Song, Hill, Liu 2018 (TCNS), sign-changing edge weights; Schäfer et al. 2022 (Nat. Commun.), steady-state loading |
| Topology switching or reconfiguration to change small-signal stability | Saric & Stankovic 2015 (TPWRS); Li, Chiang, Du 2018 (TSG); Khaji & Aghamohammadi 2017 (IJEPES); Oduor et al. 2025 |
| Grid-strength indices (SCR, gSCR, gOSCR, site-dependent SCR, impedance metrics) and their limits | Dong et al. 2019; Zhou et al. 2023; Liu et al. 2024 (gOSCR); Ma et al. 2024; Henderson et al. 2024; Wu et al. 2018; Lamrani et al. 2025 |
| The effect of adding a unit type depends on the existing mix | Li, Green, Gu 2023 (TCAS-I), five hand-picked cases, SGs added into a GFL grid; Markovic et al. 2021 (TPWRS), mix-level; Stanojev et al. 2023 (TPWRS), fixed-order replacement |
| Submodularity / supermodularity theory, and its failure for natural dynamic set functions | Topkis 1978/1998; Lovász 1983; Fujishige 2005; Bach 2013; Summers, Cortesi, Lygeros 2016; Olshevsky 2018 (average control energy not supermodular) |
| The chain (telescoping) identity for marginals and its sub/supermodularity consequences | Classical (see `theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md`) |

## 2. Potentially novel combination (what CDW adds; allowed only as a combination)

Each item is subject to its hardening gate (see the hardened claim matrix):

1. **Systematic, material, nested contextual sign reversal.** The same SG→GFL
   replacement is certified to move the *same tracked electromechanical mode* in
   opposite directions in two *nested* stable portfolios. This holds on held-out
   controller policies, and the direction of the reversal is tied to a certified
   discrete-curvature witness.
2. **Quantified failure of an oracle-optimal fixed node ranking** for next-replacement
   decisions, stratified by modal band. It is reported together with the negative
   finding that this failure is carried mostly by transitions into fast
   converter-control modes.
3. **Portfolio-conditioned total branch sensitivity used as a ranking.** It is
   validated against finite re-equilibrated reinforcements at several magnitudes, on
   held-out policies and envelope draws. It is compared head-to-head with a
   discovery-selected static baseline and with the per-condition best of nine static
   indicators, using paired resampling intervals.
4. **Electromechanical-only topology effects.** Single admissible branch actions both
   remove and create electromechanical incompatibilities, and seven static topology
   scores fail a preregistered prediction rule against the tracked EM effect.
5. **A cross-model check with a frozen, independently validated library GFL chain**
   (TX3 WECC, ANDES). No retuning is allowed. It reports which of the findings transfer
   and which are converter-model-specific.

## 3. Allowed wording

- "We show, on the IEEE 39-bus benchmark with a documented GFL model, that the
  marginal small-signal effect of an SG→GFL replacement is portfolio-dependent: the
  same replacement moves the same electromechanical mode in opposite directions in
  nested stable portfolios."
- "By the classical chain identity, such nested reversals certify discrete curvature of
  both signs; the spectral abscissa is neither submodular nor supermodular on the tested
  class."
- "Coverage across the tested holdout policies is X." "Coverage fraction over the
  declared envelope draws is X."
- "A re-equilibrated (total) first-order sensitivity, a known construction, ranks
  finite branch reinforcements far better than the static strength indicators we
  tested on held-out conditions."
- "To our knowledge, within the scope of our search, no prior work reports …" (only
  for the specific combination, with the search scope stated).
- "Static indicators carry moderate information (median abs(ρ) ≈ 0.5) but do not meet a
  preregistered prediction bar."

## 4. Prohibited wording

- "novel sensitivity formula", "we introduce total eigenvalue sensitivity", "new
  implicit-function sensitivity".
- "first", "first ever", "for the first time" (the search is not exhaustive).
- "topology switching is a new lever" (only the EM-certified combination is ours).
- "weak buses do not exist", "static metrics are useless", "static metrics fail
  everywhere".
- "the probability of reversal is X" (the policies and draws are designed conditions).
- "model-independent", "validated for converters in general", "validated for GFM".
- Any EMT claim. Any Africano/PV claim. Any simple-cycle or connected-cumulant causal
  story.
- "mixing causes reversal", or any causal reading of the two-factor decomposition.
- "robust weak corridor", unless gate H10 passed.
- "the chain identity is new".
