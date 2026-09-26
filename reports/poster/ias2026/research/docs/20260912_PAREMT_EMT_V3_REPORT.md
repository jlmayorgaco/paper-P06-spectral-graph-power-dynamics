# ParaEMT EMT-validation campaign (IAS2026 / TX4) — V3 report

**Status: V3 STOPPED at gate E0 (instrument qualification).** Branch
`research/paremt-emt-validation`, not pushed.

**Gate E0, synthetic part (blind, seed 20260912, 140 cases): FAIL on one
preregistered criterion.**
- **A5.** 14 of 15 single-channel SG-like ringdowns (S6) were resolved. The
  criterion is ≥ 95 %, which with n = 15 means 15/15.
- **All other criteria pass:**
  - no crash;
  - **zero wrong-sign verdicts**;
  - 98.3 % correct conclusive in S1/S4 and 100 % in S2;
  - p95 errors of 0.0037 s⁻¹ (α) and 0.00086 Hz (f);
  - all four DV2-2 regression cases.

**Not run (stop rule):**
- the E0 phasor-TDS holdout T1–T9 (still blind);
- every ParaEMT device gate (B-SG38, B1-GFL35 and B2-GFL37 are still blind);
- EMT04–EMT10.

No ParaEMT simulation was run in V3.

## 1. History (binding)

| campaign | preregistration | outcome (permanent) |
|---|---|---|
| V1 | c2947bd8 (+ A1 9d06436a) | G1 PASS; **G3 FAIL, G4 FAIL** (interface / device gates); EMT04–EMT18 BLOCKED |
| V2 | 0ab8efd4 | SW PASS; **G3a FAIL** under its preregistered rule, caused by the frozen v1 estimator; G3b/G4a/G4b NOT RUN; EMT04–EMT18 BLOCKED |
| V3 | 6d339006 (+ AM1 in f8e0d44f) | **E0 FAIL** (synthetic criterion A5); everything after E0 NOT RUN |

**Wording:**

> V1 and V2 failed their preregistered gates. Subsequent diagnostics
> identified independent interface and estimator defects. V3 addressed the
> estimator defect prospectively with a new, independently implemented
> instrument. On its blind qualification set the instrument made no wrong-sign
> call and met the accuracy criteria, but it failed one preregistered
> resolution criterion: single-channel ringdowns, 14/15 resolved against a
> 15/15 requirement. By rule V3 therefore stopped before any ParaEMT device or
> portfolio run.

No V1, V2 or V3 failure is described as invalid.

## 2. V3-0 — audit of the frozen v1 estimator (commit 4e35164e)

Source: `docs/20260912_EMT_ESTIMATOR_FAILURE_AUDIT.md` and
`results/EMTV3/V3-0/estimator_audit.json`. `emt_estimator.py` was not modified.

- **Unresolved everywhere.** `resolved = False` in all four V2 G3a
  evaluations, with a residual of about 0.999. The G3a rule compared α and f
  unconditionally.
- **B-SG36.** The EMT run and its reference share the same spurious
  +6.3 s⁻¹ mode.
- **DV2-2.** Three LinAlgError crashes and one unresolved case.
- **Mechanism:**
  - spurious |z| up to 18.9;
  - the full-record Vandermonde z^k spans 10^26 to 10^95, or overflows (10^355);
  - lstsq keeps rank 1–3 of 10–16;
  - the exp(2αT) energy weight (10^43) selects a spurious mode with a residue
    of 1e-22.

## 3. The V3 instrument (prereg 6d339006; amendment V3-AM1)

**Estimator** (`experiments/paremt_emt/v3/emtv3_estimator.py`): a new
multi-channel bounded variable-projection (VP) damped-mode estimator.
- **Model.** For each channel, one damped mode plus an affine term, with time
  centred on the window.
- **Data.** 0.1-s block means, which leave every exponential's s unchanged.
- **Fit.** One Householder QR of the common, column-scaled 4-column design
  matrix. Only (α, f) are optimized, within bounds (grid, then trf).
- **Uncertainty and decision.**
  - channel and time-block jackknife CI;
  - resolved flag R1–R4;
  - verdict from the CI with EPS = 0.002;
  - targeted tracker (committed f_pred ± 0.15 Hz; α_pred never used);
  - broad-band sentinel;
  - adaptive run length 30 → 60 → 90 s.
- **Removed.** There is no full-record Vandermonde and no exponential energy
  weighting.

**V3-AM1** (`docs/20260912_PAREMT_EMT_PREREG_V3_AMENDMENTS.md`), committed
**before** any blind data existed:
- **Problem.** On the development data (seed 1111, DV2-2, the non-holdout TDS
  P4 BASE), the unweighted fit carried a nuisance-mode bias of 0.0055 on DV2-2.
  It also left one strongly damped S6 development case unresolved.
- **Change.** An envelope-equalizing weight exp(−α_w τ), with
  |α_w| ≤ 3/L, taken from a first unweighted pass and then fixed. The
  optimization is still over (α, f) only, and no criterion was changed.
- **Development result:**
  - DV2-2 bias 0.00036;
  - S6 5/5 on the development set;
  - development TDS P4 BASE Δα = 5.6e-4 against the linear eigenvalue.

## 4. Gate E0 — blind synthetic qualification

**Data:** `results/EMTV3/E0/E0_synthetic.{csv,json}`, 140 cases generated from
seed 20260912 after the implementation commit f8e0d44f.

| stratum | n | resolved | correct conclusive | wrong sign | p95 \|Δα\| (resolved) | p95 \|Δf\| (resolved) |
|---|---|---|---|---|---|---|
| S1 single mode | 40 | 40/40 | 39/40 (the miss is α = +0.0029, UNRESOLVED; not in A3) | 0 | 9.5e-7 | 1.7e-7 |
| S2 target + nuisance | 25 | 25/25 | 25/25 | 0 | 3.4e-3 | 9.5e-4 |
| S3 near boundary (\|α\| ≤ 0.01) | 20 | 20/20 | 12/20 (the rest UNRESOLVED, allowed) | 0 | 3.6e-8 | 7.7e-9 |
| S4 heterogeneous channels | 20 | 19/20 | 19/20 | 0 | 2.2e-3 | 4.7e-4 |
| S5 two close modes | 20 | 17/20 | 17/20 | 0 | 4.9e-3 | 1.3e-3 |
| S6 single-channel SG-like ringdown | 15 | **14/15** | 14/15 | 0 | 5.8e-3 | 8.9e-4 |

| criterion | rule | result | verdict |
|---|---|---|---|
| A1 | no crash | none | pass |
| A2 | 0 wrong-sign verdicts with \|α\| ≥ 0.005 | 0 | pass |
| A3 | S1+S4 with \|α\| ≥ 0.01: ≥ 95 % correct conclusive | 57/58 = 98.3 % | pass |
| A3b | S2 with \|α\| ≥ 0.01: ≥ 80 % | 100 % | pass |
| A4 | p95 \|Δα\| ≤ 0.005 and p95 \|Δf\| ≤ 0.005 (all resolved) | 0.0037 s⁻¹; 0.00086 Hz | pass |
| **A5** | **S6 ≥ 95 % resolved** | **14/15 = 93.3 %** | **FAIL** |
| A6 | DV2-2 regression | 4/4 | pass |

**Reported:** the CI coverage of α_true is 99.3 %; the S3 and S5 resolution
fractions are 100 % and 85 %; the sentinel flagged nothing.

**The failing case.** S6 #137 has α_true = −0.473 s⁻¹ and f = 1.067 Hz, a single
channel with a slow real mode superposed, in the window [1.5, 10] s.
- The estimate was α̂ = −0.455.
- The time-block jackknife CI half-width was 0.027 s⁻¹, above the 0.02 limit
  (R4), with a median residual of 0.26.
- The instrument therefore reported **UNRESOLVED instead of a biased value**.
  This is the preregistered guard working as intended.
- The criterion nevertheless requires every S6 case to resolve, so the gate
  fails.

The unresolved S4 and S5 cases also had R4 CI failures. Their point estimates
were biased (for example α̂ = −0.419 against −0.262), and each was withheld as
UNRESOLVED, never reported as a value.

**Gate E0: FAIL.** V3 stops (prereg V3 §1.3). No other estimator is tried in
this campaign.

## 5. What was not run, and why

| item | status | consequence |
|---|---|---|
| E0 phasor-TDS holdout T1–T9 | NOT RUN (E0 already failed on the synthetic part; stop rule) | the holdout traces were never generated, so the set stays blind |
| G3a, G3b (R-SG30, R-SG36, **B-SG38**) | NOT RUN | B-SG38 has never been simulated |
| G4a, G4b (R-GFL30, **B1-GFL35**, **B2-GFL37**) | NOT RUN | the GFL holdouts are still blind; the V2 voltage interface has still never run in a network |
| EMT04, EMT05, EMT06, EMT07/08, EMT10 | NOT RUN | no IEEE-39 EMT portfolio result exists |

**Prepared but never executed** (committed separately and labelled as such):
- `v3/EMTV3_refs.py`;
- `v3/EMTV3_gates.py`;
- `v3/v3_case.py`;
- `v3/EMTV3_science.py`;
- `overlay/tx4_emt_v3.py` (the descriptive EMT-native-stator variant).

`EMTV3_gates.py` asserts that E0 passed and `EMTV3_science.py` asserts that
every model gate passed, so neither can run under the V3 record.

## 6. Incidents and disclosure

1. **The first blind E0 execution was killed** by the system for low memory,
   mid-way through S5. The cause was system-wide memory pressure from unrelated
   applications. No other process was touched.
   - Its progress log had shown **verdict labels** (not true α values) for
     cases up to S5 #6.
   - The implementation was not changed. E0 was re-executed with the identical
     code and seed.
2. **Determinism.** Prereg §8 requires one E0 rerun. It was executed once more
   from the same commit; the comparison is in
   `results/EMTV3/E0/E0_rerun_identity.json`.
3. **Commit f8e0d44f** (implementation) preceded all blind data. The style-only
   commit 70d5be8e (import order) has a byte-order-mark character in its
   subject line; this is cosmetic.

## 7. Claim ledger

No change. `results/20260911_EMT_CLAIM_MATRIX.csv` (V2-0) stands.
- No claim receives EMT support, and none is EMT-refuted.
- V19 (network and operating point) remains the only EMT QUANTITATIVELY
  REPRODUCED claim.
- V20 (SG) and V21 (custom GFL) keep their final-ledger status, which rests on
  the independent ANDES phasor reproductions.

## 8. For the author (decision memo; nothing here is executed)

**What V3 did establish about the instrument.** The blind evidence (outside
the S6 criterion) is that the new VP estimator does not crash, gives no
wrong-sign verdict, and is accurate to p95 0.004 s⁻¹ in α. Where it cannot
separate modes, it withholds a value (R4) rather than report a biased one.
This is instrument evidence only, not scientific evidence.

**Why V3 stopped.** The preregistered criterion demanded every one of 15
strongly and weakly damped single-channel ringdowns to resolve. One heavily
damped case (α ≈ −0.47 s⁻¹ in an 8.5-s window) did not.

**Options.**
- **(a) Close the EMT line with its present record.** The network and operating
  point are EMT-reproduced (V19). Every other claim keeps its phasor and ANDES
  status. This is consistent with "no V4 merely to obtain agreement".
- **(b) A new, separately justified campaign.** Any continuation would need its
  own preregistration and its own reason, not agreement-seeking. Its blind
  elements would be:
  - the unused TDS holdout T1–T9;
  - B-SG38, B1-GFL35 and B2-GFL37;
  - a fresh synthetic seed.

## 9. Reproducibility

**Chain** (paths relative to `reports/poster/ias2026/research`; commits
6d339006, f8e0d44f and 70d5be8e):

```
.venv/xtool-paremt   experiments/paremt_emt/v3/EMTV3_estimator_audit.py   (V3-0)
.venv/xtool-paremt   experiments/paremt_emt/v3/EMTV3_E0.py dev            (development record)
.venv/tx3-analysis   experiments/paremt_emt/v3/EMTV3_tds_traces.py dev    (development TDS case)
.venv/xtool-paremt   experiments/paremt_emt/v3/EMTV3_E0_tds.py dev
.venv/xtool-paremt   experiments/paremt_emt/v3/EMTV3_E0.py blind          (gate E0, synthetic)
```

---

## What EMT did and did not validate (after V1, V2 and V3)

The table separates **theorem validity** from **model-fidelity corroboration**.
Statuses are those of the final canonical ledger.

| Theory / claim | Mathematical status | Phasor evidence | Independent phasor evidence | EMT evidence (V1 + V2 + V3) | Final status |
|---|---|---|---|---|---|
| V01, V02, V08, V22, V26, T04, T06, T08 (theorems) | proved | yes | where applicable | not applicable | PROVED |
| V03 P4 trap, κ = 4 [B02] | — | PCV02; G2 TDS | ANDES Phase 8, 32/32 | BLOCKED (V1, V2); NOT RUN (V3 stopped at E0) | VALIDATED (model-specific witness) |
| V06, V07 policy dependence [B01] | structure proved | clean g-only counterfactual | ANDES Phase 8, P4 vs G_S | BLOCKED / NOT RUN | VALIDATED |
| V10, V11 port derivative, line ranking | proved (V10) | PCV04 | ANDES R3 11/12, ρ 1.00 | BLOCKED / NOT RUN | PROVED + VALIDATED; VALIDATED (local ranking) |
| V13, V14 remediation, governors | — | PCV04; FC03; G2 | partial ANDES | BLOCKED / NOT RUN | VALIDATED (model-specific); BENCHMARK-SPECIFIC |
| V15–V18 robustness; V30, B07, B10 nonlinear scope | — / normal form | PCV05; G2 TDS | — | BLOCKED / NOT RUN | as in the ledger |
| V19 network and operating point [I01] | — | reference | ANDES; pandapower after the tap fix | **EMT QUANTITATIVELY REPRODUCED** (V1 G1) | CONDITIONAL (pandapower translation) |
| V20 SG dynamics [I02, I03] | — | internal | ANDES R2/R3, E1 | V1 G3 FAIL; V2 G3a FAIL (estimator); V3 NOT RUN → EMT UNRESOLVED | VALIDATED (by ANDES) |
| V21 custom GFL in a second simulator [I04] | — | project DAE | **ANDES Phase 8, GATE 5 PASS** | V1 G4 FAIL; V2 G4 NOT RUN; V3 NOT RUN → EMT UNRESOLVED | VALIDATED as a reproduction of computation; EMT reproduction not achieved |
| cumulant and closure diagnostics, V05, V12, B04, B06, negatives | as in the ledger | as in the ledger | — | not applicable / BLOCKED | unchanged |
