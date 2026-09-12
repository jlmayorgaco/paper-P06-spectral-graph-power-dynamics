# TX4: Policy-Dependent Minimal Incompatibility of Synchronous-to-Inverter Replacement Portfolios: Network-Closure Interactions and Transverse Stability

`main.pdf` is the hardened final draft: 12 pages (the last page holds one reference), IEEEtran journal format. It is
not the IAS poster and not the TX3 paper; neither of those was modified.

## Evidence state

- The final scientific evidence is frozen by tag
  `IAS2026_FINAL_SCIENTIFIC_EVIDENCE_FREEZE` at commit `b9f274e2`.
- The manuscript also cites these post-freeze evidence addenda:
  - `f775db89`: FC18 targeted port checks, i.e. the boundary anatomy in
    Tables III–IV and Figs. 5–7;
  - `7a772808` and `f2946257`: Monte Carlo preregistration and its amendment;
  - `cfe9fe4b`: Monte Carlo results (Table VI, Fig. 9).
- No new experiments were run during manuscript hardening.
- Final reconciliation (2026-09-12, ParaEMT closure):
  - the Limitations now carries the memo B independent-implementation text
    (accept-fix: ANDES, pandapower after the tap-convention translation, and
    the Phase 8 custom-model reproduction scoped to the computation);
  - it also carries the EMT-corroboration statement (no EMT validation of
    the portfolio claims is asserted), with the modal-energy observation
    labelled as outside the claims;
  - the F10 additive count is corrected to 31/52 (FINAL ledger V05);
  - the memo C clarification separates kappa = 4 (the minimum failing-set
    cardinality) from the irreducible connected interaction order;
  - a new "Witness identity" Limitations bullet says {30,33,35,37} is the
    witness of the nominal frozen model, not a robust weak-bus set (retained
    in 34/69/62 % of the fleet-wide/per-unit/combined envelope draws).
- TX4 is frozen after this edit (local tag `TX4_FINAL_MANUSCRIPT_FREEZE`).
- Records, under `reports/poster/ias2026/research/docs/`:
  - `20260912_PAREMT_EMT_CLOSURE.md`;
  - `20260911_FINAL_VALIDATION_LEDGER.md` section 6;
  - `20260912_TX4_FINAL_AUTHOR_DECISIONS.md` (memo A, B and C all
    ACCEPTED).

## Story: three contributions only

- **C1: policy-dependent minimal incompatibility.** H(theta), kappa(theta), the
  transverse definition, and any-order-safe planning. The governor result is
  central: the four-unit witness is conditional on the absence of primary
  frequency control, while policy dependence persists.
- **C2: network-closure anatomy.** Local descriptor-affine actions, exact
  network closure, a principal-minor hierarchy, and a coalition-specific zero
  (Fig. 1). At all ten frozen boundaries, no truncation below the cardinality
  of the changing minimal coalition reaches the zero. That cardinality is a
  closure order, not an irreducible interaction order: kappa = 4 is not
  evidence of an irreducible four-device interaction (memo C,
  connected-cumulant clarification).
- **C3: actionable boundary motion.** Exact boundary normals, iterative Q/V
  retuning, and a symmetry-deflated zero-frequency port (Kundur: 28/28 zero
  crossings, one coalescence correctly not flagged, 0/445 false positives).

Everything else is supporting evidence, robustness analysis, or a limitation:
Monte Carlo validation, the nonlinear scope check, the second-order curvature
audit (Appendix B), the rank negative control, and the non-convergent walk
expansion.

## Explicitly not claimed

- Rank explains kappa (negative control: rank 7–8 throughout 4 → 3 → 2 → 3 → 4).
- A dominant path explains a failure (Neumann walk expansion non-convergent,
  rho ≥ 1.003).
- A simple cycle causes a failure (falsified earlier).
- Nonlinear composability as a contribution.

## Build

`latexmk` needs Perl, which is not installed. Build with:

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

If `main.pdf` is open in a viewer on Windows, close it first; otherwise
pdflatex cannot write the file.

## Figures and sources

Every figure is produced by
`reports/poster/ias2026/research/experiments/paper_tx4/MC03_paper_figures.py`
and has a `*_source.csv` in `figures/`. The only exception is Fig. 1, which is
drawn in TikZ inside `main.tex`.

## Reviewer hardening

`reports/poster/ias2026/research/docs/REVIEWER_1_2_ATTACKS.md` lists 17
objections, each with an evidence-backed response and an explicit concession.

## Deviations disclosed in the paper (Sec. VI)

1. **P1 amendment.** Made before the run; threshold unchanged.
2. **P2.** Failed as implemented. The post-hoc re-bisection is reported
   separately and never counted as a preregistered pass.
3. **S7.** Three coalescence-class draws contained a genuine zero crossing,
   and 183 draws were left undecided by the device-flip rule.
4. **P1 chart.** Two draws needed a redrawn off-equilibrium state.

## Before submission (author action)

- Confirm the affiliation and e-mail line.
- Add funding and acknowledgments, if any.
- Choose the repository or DOI to cite in the Reproducibility statement.
- Check the target journal's page policy (the draft has 12 pages).
