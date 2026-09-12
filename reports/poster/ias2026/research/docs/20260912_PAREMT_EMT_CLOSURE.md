# ParaEMT EMT line — archival closure (TX4)

Date: 2026-09-12. Branch `research/paremt-emt-validation`, not pushed.

**Status: CLOSED PERMANENTLY.**
- No further ParaEMT scientific experiment will be run for TX4.
- This document only archives the record. It adds no new computation.

## 1. One-paragraph summary

Electromagnetic-transient (EMT) corroboration of the TX4 IEEE-39 results was
preregistered four times in ParaEMT (V1, V2, V3, SQ1).
- **What was reproduced.** The IEEE-39 network and operating point (V1, gate
  G1 PASS).
- **Where each campaign stopped.** Each campaign stopped at a pre-specified
  qualification gate:
  - V1: the device-interface gates;
  - V2: the synchronous-machine transcription gate, caused by the frozen v1
    modal estimator;
  - V3 and SQ1: the modal-identification (estimator) qualification gates.
- **What was never run.** No portfolio-level EMT experiment (EMT04 onward) was
  ever run.

Therefore:
- no EMT validation of any portfolio, policy, design, remediation or nonlinear
  claim exists or is asserted;
- the only EMT result is the network/equilibrium reproduction;
- no earlier failure is described as invalid, and none was reinterpreted.

## 2. Campaign record

| campaign | preregistration | stopping gate | outcome (permanent) | result / report commits |
|---|---|---|---|---|
| **V1** | c2947bd8 (+ amendment A1 9d06436a) | G3 (SG), G4 (GFL) | G1 PASS; **G3 FAIL, G4 FAIL**; EMT04–EMT18 BLOCKED | ff1c4d8d, dd1e2ea3, 17832edf, ede6f3eb |
| **V2** | 0ab8efd4 (V2-0 reconciliation cf8a0b38) | G3a | SW PASS; **G3a FAIL**; G3b/G4a/G4b NOT RUN | f8b8fb15, 632374dc, e8ab281e |
| **V3** | 6d339006 (+ V3-AM1 in f8e0d44f; audit 4e35164e) | E0 (instrument) | **E0 FAIL** on criterion A5 only; T1–T9 and all device gates NOT RUN | 229613de, 00c75f52, b9ab930f |
| **SQ1** | 0555f587 | SQ1-2 (blind synthetic qualification) | **SQ1-2 FAIL** (GB, GC); line stopped permanently | fc080445, 651a5ff3, 8f06220e |

**Reports** (`docs/`):
- `20260911_PAREMT_EMT_FINAL_REPORT.*` (V1, with the V2-0 erratum);
- `20260911_PAREMT_EMT_V2_REPORT.*`;
- `20260912_PAREMT_EMT_V3_REPORT.*`;
- `20260912_PAREMT_EMT_SQ1_REPORT.*`.

## 3. Why each campaign stopped

### V1 — device-interface gates

- **G1 PASS.** Network and equilibrium: Ybus 1.7e-16, V 1.5e-6.
- **G3 FAIL.** The SG field voltage deviated 2.26 %, against a 2 % tolerance.
  - Post hoc diagnosis: the EMT line's (X/ω0)dI/dt term.
  - With an algebraic line, the SG transcription was exact to 1.9e-4 (D1,
    diagnostic only).
- **G4 FAIL.** The GFL was an ideal current source into a lossless trapezoidal
  inductor, which produced Nyquist-rate chatter.
- EMT04–EMT18 were blocked.

### V2 — SG transcription gate

- The network was redesigned algebraically, and a GFL voltage source was placed
  behind Rf–Lf. Both were preregistered.
- **SW PASS.**
- **G3a FAIL.**
  - Every SG state trajectory agreed to ≤ 2.1e-4 of its excursion on the blind
    bus-36 case.
  - The preregistered ringdown sub-check still failed. It reused the frozen v1
    matrix-pencil estimator, which returned a spurious mode for both the EMT
    run and its phasor reference.
- The stop rule applied: G3b, G4a and G4b were not run. The GFL voltage
  interface was never run in a network.

### V3 — instrument qualification (E0)

- A new bounded variable-projection damped-mode estimator was frozen before any
  blind data (f8e0d44f).
- **E0 FAIL** on the blind synthetic set (seed 20260912, 140 cases): criterion
  A5 resolved 14 of 15 single-channel ringdowns, against 15/15 required.
- All other criteria passed: zero wrong-sign verdicts; p95 errors 0.0037 s⁻¹
  and 0.00086 Hz.
- The stop rule applied: the TDS holdout T1–T9 and every ParaEMT device gate
  were not run.

### SQ1 — scientific qualification (SQ1-2)

- The V3 estimator was frozen unchanged.
- **Domain.** Derived only from frozen predictions: α ∈ [−0.23, 0.20] s⁻¹,
  f ∈ [0.55, 1.05] Hz, with a modal library of 26 frozen phasor models.
- **SQ1-2 FAIL** on the blind synthetic set (seed 20260913, 210 domain cases):
  - GB: 78.4 % resolved for |α| ≥ 0.01, against ≥ 95 % required;
  - GC: p95 |Δα| = 0.0066 s⁻¹, against ≤ 0.005 required.
- **Passing gates:**
  - GA: zero wrong-sign verdicts;
  - GD: p95 |Δf| of 0.0017 Hz;
  - GE: no crash.
- The determinism rerun was byte-identical.
- The preregistered rule applied: **STOP the ParaEMT line permanently.**

## 4. What was and was not tested

| element | tested in EMT? | result |
|---|---|---|
| IEEE-39 network, power flow, operating point | yes (V1 G1) | **REPRODUCED** |
| SG model in EMT (buses 30, 36; blind 38) | partly (V1 G3; V2 G3a on 30/36) | FAIL at the preregistered gates; blind SG38 never simulated |
| GFL model in EMT (bus 30; blind 35, 37) | V1 G4 only (current-source interface) | FAIL; the V2 voltage interface was never run in a network; blind GFL35/GFL37 never simulated |
| modal estimator, frozen v1 | yes (V2) | unfit as frozen (audit 4e35164e) |
| modal estimator, V3 | synthetic only (V3 E0, SQ1-2) | FAIL (A5; GB/GC); no wrong-sign call |
| phasor-TDS estimator holdout T1–T9 | no | never generated |
| EMT04 base/holdouts; EMT05 P4 lattice; EMT06; EMT07/08 policy boundary; EMT10 Newton | no | never run |
| EMT09, EMT11–EMT18 (P_inf, remediation, lines, governors, robustness, nonlinear) | no | never run |

## 5. Final claim impact

These statuses are recorded in the new `emt_validation` column of
`results/20260911_FINAL_VALIDATION_MATRIX.csv` (FINAL ledger §6). The
mathematical, phasor-TDS and ANDES statuses are unchanged.

| category | claims | EMT status |
|---|---|---|
| EMT network / equilibrium | V19 (I01) | **EMT REPRODUCED** — the only EMT-reproduced claim |
| device models in EMT | V20 (I02, I03, SG), V21 (I04, custom GFL) | **EMT UNRESOLVED**: device gates failed or were not reached |
| portfolio-level claims | V03, V04, V06, V07, V10, V11, V13–V18, V29, V30 | **NEVER EVALUATED** (EMT unresolved) |
| mathematical, diagnostic or refuted | V01, V02, V05, V08, V09, V12, V22–V28 | NOT APPLICABLE |

**Final status of the independence claims** (from the FINAL ledger plus the
EMT column):

| id | claim | independent phasor / static status (unchanged) | EMT status |
|---|---|---|---|
| I01 (V19) | network and operating point implementation-independent | CONDITIONAL on accepting the pandapower tap-convention translation (ANDES Ybus 1e-13, V 1.7e-7); the author accepted the translation on 2026-09-12 (memo B accept-fix), so the condition is met | **REPRODUCED** (ParaEMT V1 G1) |
| I02, I03 (V20) | equation-equivalent SG dynamics and sensitivities independent | VALIDATED (ANDES R2/R3 12/12, α within 1.5e-6 s⁻¹; E1 holdout 12/12 signs, ρ 1.000) | UNRESOLVED (V1 G3 FAIL; V2 G3a FAIL) |
| I04 (V21) | custom GFL results reproduced | VALIDATED as a reproduction of the computation (ANDES Phase 8, GATE 5 PASS, 32/32) | UNRESOLVED (V1 G4 FAIL; later never reached) |

**Manuscript.** The TX4 Limitations section states:

> "EMT corroboration was preregistered in ParaEMT. Although the IEEE-39
> network and operating point were independently reproduced, the campaigns
> terminated at pre-specified device-interface and modal-identification
> qualification gates before any portfolio-level experiment; therefore, no
> EMT validation of the portfolio claims is asserted."

No EMT portfolio claim remains anywhere in the manuscript (see §8).

## 6. Methodological observation (not a claim; not an IEEE-39 EMT result)

> **Dominant observed modal energy need not correspond to the
> stability-limiting eigenmode.**

**Source.** The modal-library and estimator-qualification work of SQ1-1 and
SQ1-2. The evidence comes from **linear modal residues of the frozen phasor
models** and from synthetic signals built from them. It is **not** an EMT
result and **not** an IEEE-39 EMT observation.

**Evidence:**
- In 10 of the 26 frozen phasor models (prereg SQ1 §2):
  - the stability-limiting mode is the rightmost in-band transverse mode, the
    0.91-Hz family;
  - under the preregistered load pulse it carries a minority of the channel
    energy;
  - a neighbouring 0.64-Hz family, which decays only 0.003–0.07 s⁻¹ faster,
    carries 85–92 % of the median channel energy.
- Across the library, 30 % of the channels have a target energy share below
  0.1.
- In SQ1-2, templates of this kind were resolved in 40.7 % of cases, against
  96.5 % for all other templates (post hoc; not a gate).

**Consequence (methodological).**
- A modal-identification step that locks onto the most energetic observed mode
  can report the damping of a mode that does not decide stability.
- Any measurement-based stability verdict must therefore target the mode
  predicted to be stability-limiting and qualify its identifiability in
  advance.
- This is kept outside the TX4 claims. It is mentioned in the manuscript only
  as a labelled remark inside the Limitations item on EMT.

## 7. Archive state

- **Committed:**
  - all preregistrations, code, results and reports of V1–SQ1;
  - the SQ1 scripts prepared but never executed, with banners (fc080445);
  - the V3 prepared scripts (00c75f52).
- **Not committed:** the large raw EMT runs. They live under
  `external/paremt_runs/`, excluded from git through `.git/info/exclude`; their
  manifests and hashes are in the results.
- **Frozen and never modified:**
  - ParaEMT d79d735a (`external/ParaEMT_{upstream,tx4}`);
  - the venv `.venv/xtool-paremt`;
  - the v1 estimator (`experiments/paremt_emt/emt_estimator.py`);
  - the V3 estimator (sha256 `849472d8…3b0817f`).
- **Other untracked files in the research tree** are not part of this line:
  - `docs.zip`, which belongs to others;
  - older `outputs/ias2026/` run directories from the final-closure and
    Monte Carlo campaigns.

## 8. Final manuscript reconciliation and stale-claim search

**Manuscript:**
`reports/papers/tx4_policy_dependent_incompatibility/main.tex`, built as
`main.pdf` (12 pages; pdflatex + bibtex; the only overfull boxes are the two
output-routine boxes that the pre-reconciliation build already had).

**Changes, all in Limitations unless noted:**
1. **Independent implementation.** The stale "within 4.1 %" bullet was replaced
   by memo B, accept-fix version, per the author's decision of 2026-09-12. It
   covers:
   - ANDES, and pandapower after the disclosed tap-convention translation;
   - SG eigenvalues and sensitivities;
   - the Phase 8 custom-model reproduction, scoped to reproducing the
     computation and not the adequacy of the converter model.
2. **EMT corroboration.** A new bullet carries the required sentence verbatim.
   It is followed by the methodological observation, labelled as coming from
   the modal-identification qualification (linear residues of the phasor
   models) and "not an EMT result", and kept outside the claims.
3. **Discussion.** The F10 additive count was changed from 30/52 to 31/52, with
   third order 50/52 (FINAL ledger V05, which supersedes 22/30/38).
4. **Reproducibility.** The addenda now also list 0a4e18c2 and 7773ea8c, and
   the archive of the ParaEMT record.
5. **Bibliography.** Entries added for pandapower (Thurner et al., 2018) and
   ParaEMT (Xiong et al., IEEE Trans. Power Del. 39(2), 2024, as given in the
   ParaEMT upstream README).

**Not applied** (it is still an open author choice): memo C, the C2
closure-order clarification.

**Stale-claim search** (case-insensitive, over `main.tex` and `README.md` after
reconciliation):

| term | hits | assessment |
|---|---|---|
| EMT / ParaEMT | main.tex: the "Model" bullet ("not an electromagnetic transient model"), the EMT-corroboration bullet, and Reproducibility | correct; no EMT portfolio claim |
| single implementation | none | — |
| NOT YET TESTED | none | — |
| custom GFL | none (the text says "transcribed into ANDES as custom models") | consistent with V21 wording |
| ANDES | model source (ANDES IEEE-39 case), governor direction, independent-implementation bullet | consistent with V19–V21; "Andes" also matches the author affiliation (false positive) |
| 28/28 | main.tex Thm 4 holdout list and Table VI; README | consistent: 28/28 origin crossings, with the coalescence disclosed (note b: 29 real-count changes) |
| 28/29 | none | the 28/29 total is expressed as 28/28 plus the disclosed coalescence |
| 4.1% | none (removed) | the stale claim is gone |

The IAS poster (`reports/poster/ias2026/main.tex`) contains no EMT or ParaEMT
mention. It has uncommitted edits by others and was not touched.