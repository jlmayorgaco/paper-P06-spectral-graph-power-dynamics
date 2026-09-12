# TX4 final author decisions (manuscript closure)

Date: 2026-09-12. Branch `research/paremt-emt-validation`, not pushed.

**Scope.**
- These are the decisions on the three open manuscript choices of
  `docs/20260911_AUTHOR_DECISION_MEMO.md`.
- Manuscript: `reports/papers/tx4_policy_dependent_incompatibility/main.tex`.
- No numerical result, figure, experiment or frozen evidence was changed, and
  no new scientific computation was run.

| item | decision | applied in |
|---|---|---|
| **A. pandapower fix** (equation-preserving tap-convention translation; Ybus identity 1e-13) | **ACCEPTED** | ab9edc5d (through diff B, accept-fix version) |
| **B. Limitations replacement** (stale "within 4.1 %" bullet → memo diff B, accept-fix) | **ACCEPTED** | ab9edc5d |
| **C. C2 clarification** (κ = 4 is the minimum failing-set cardinality, not an irreducible four-device interaction) | **ACCEPTED** | the final manuscript-closure commit |

## A. pandapower fix — ACCEPTED

- **What it is.** The as-supplied `from_ppc` conversion moved the off-nominal
  tap of transformers 35, 37 and 38 to the other winding. The translation
  re-expresses the same two-port exactly: Ybus identical to 1.14e-13, power
  flow 3e-14.
- **Manuscript.** The Limitations states that pandapower needed a
  tap-convention translation of the three transformers whose off-nominal tap is
  on the low-voltage from-bus.
- **Ledger effect.** V19 (I01) was "CONDITIONAL on accepting the pandapower
  translation". That condition is now met. The frozen status text in
  `results/20260911_FINAL_VALIDATION_MATRIX.csv` is left unchanged, and the
  resolution is recorded in the FINAL ledger (§4).

## B. Limitations replacement — ACCEPTED

The "Independent implementation" bullet now states:
- the network and operating point are reproduced in ANDES 2.0.0 and, after
  the translation, in pandapower 3.4.0;
- the equation-equivalent SG eigenvalues and branch sensitivities are
  reproduced (12/12, Spearman 1.00);
- the Phase 8 custom-model reproduction: 32/32 verdicts, H = {{30,33,35,37}}
  and κ = 4 at P4, |Δα⊥| ≤ 1.3e-6 s⁻¹. It is scoped explicitly as reproducing
  the computation, not the adequacy of the converter model; with a library
  converter model the branch ranking changes (7/12).

The EMT-corroboration bullet (ParaEMT closure) was added in the same commit.

## C. C2 clarification — ACCEPTED

**Manuscript changes** (wording as instructed by the author):
1. **C2 contribution (Introduction).** The full-order term closes the
   determinant at the boundary zero, and "this full-order term is a closure
   order, not an irreducible |S|-device interaction."
2. **Section V, C2 subsection.** A new paragraph, "Failing-set cardinality is
   not interaction order", contains the author's text:
   - κ = 4 is the cardinality of the minimum failing replacement portfolio;
   - it is not evidence of an irreducible four-device interaction;
   - the full-support Boolean closure contribution is required at the observed
     boundary but is composite.

   The paragraph then gives:
   - the identity μ₁₂₃₄ = χ₁₂₃₄ + χ₁₂χ₃₄ + χ₁₃χ₂₄ + χ₁₄χ₂₃, checked numerically
     to about 3e-14 (memo C);
   - the deletion displacement, 0.024 s⁻¹;
   - the connected cluster with the largest effect on the zero, {33,35};
   - the P4 lower-order facts: the order-≤3 determinant truncation has its
     zero at +0.1267 s⁻¹ against the exact +0.1270; the spectral-abscissa
     truncations are −0.18 to −0.21.
   - It states that the cumulant decomposition is interpretation only and not
     a contribution.
3. **Fig. 5 caption.** At the four-unit boundary the full-order term contains
   products of pair interactions besides the connected four-device term, so
   the full order is a closure order.
4. **Conclusion.** It adds "(a closure order, not an irreducible interaction
   order)".
5. **"No path attribution" paragraph.** "the order of the interaction" becomes
   "the closure order of the interaction expansion".

**Sources** (all frozen and committed):
- `results/CC/CC03/CC03_boundaries.csv`, row FLAG P4-line:
  - |χ_H| = 0.1521 and |μ_H| = 0.1614;
  - ν_H = 0.069;
  - deletion distance 0.0239 s⁻¹ (converged, to Re −0.0055);
  - dominant δ cluster 33+35.
- `docs/20260911_PORTFOLIO_DECISION_BENCHMARK.md`, lines 94–95: order ≤ 2 gives
  −0.0148 and order ≤ 3 gives +0.1267.
- `results/PCV/PCV03/PCV03_summary.json`: α 0.12701; B3/B4/B5/B5b −0.203,
  −0.213, −0.210, −0.184.
- The 3e-14 identity check is recorded in the author decision memo
  (Facts, "Pair identity").

**One correction to the memo C text.** The memo's proposed Fig. 5 caption said
that the full-order term is "dominated by products of pair interactions". The
frozen CC03 values contradict that. At the flagship boundary |χ₁₂₃₄| = 0.152,
against |μ₁₂₃₄| = 0.161, so the connected four-way term is not small *within
the Boolean coefficient*.
- **What "composite" refers to.** The label comes from the connected share
  ν = 0.069 of all partition terms of the characteristic function at the zero,
  and from the small deletion displacement, 0.024 s⁻¹. It does not come from
  the size of χ₁₂₃₄ relative to μ₁₂₃₄.
- **Label sensitivity.** The label is threshold-sensitive: it would be
  CONNECTED at 0.05.
- **Caption as applied.** It therefore makes only the structural statement:
  the term contains products of pair interactions besides the connected
  four-device term. The paper's quantitative evidence is the deletion
  displacement, not ν.

**Unchanged:**
- the exact network-closure result: Theorem 3, the Boolean statement that no
  truncation below |S| reaches the zero at s\*, Table III, Table IV and Fig. 5;
- all numbers and figures.

**Status of the connected-cumulant layer.** It stays secondary (FINAL ledger
V09 and V22–V25) and is not promoted to a contribution.

## Final checks (on the closure commit)

- **Build.** `main.pdf` was rebuilt from a clean state (aux, bbl, log, out and
  pdf removed; pdflatex, bibtex, pdflatex, pdflatex).
  - It has 12 pages.
  - There are no errors and no undefined references or citations.
  - The only overfull boxes are the two output-routine vboxes that the
    pre-reconciliation draft already had.
- **Statements that must not appear.** None does:

| statement | result |
|---|---|
| irreducible four-device interaction | every occurrence is a negation (C2 contribution, the C2 paragraph, the Fig. 5 caption, the Conclusion) |
| robust identity of {30,33,35,37} | none. The set appears only as the nominal-model P4 witness, conditional on the model and on the absence of governors; "every result is conditional on this model" |
| EMT validation of portfolio claims | none; the Limitations states that none is asserted |
| single-implementation custom GFL | none; the Limitations reports the ANDES custom-model reproduction, scoped to the computation |
| old additive count 30/52 | none; the Discussion gives 22, 31 and 38 of 52 (third order 50) |
