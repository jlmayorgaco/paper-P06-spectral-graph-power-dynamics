# ruff: noqa: E501  -- claim table kept one row per line
"""Claim-to-EMT matrix -> results/20260911_EMT_CLAIM_MATRIX.csv (V2-0 reconciliation).

Source of truth: the FINAL canonical validation ledger results/20260911_FINAL_VALIDATION_MATRIX.csv
(V01-V30, post-cumulant Phase 12, commit a70dac94), read only. The 13 old ids that no V row
cross-references (T04, T06, T08, C08, B04, B06, B07, B10, N01, N05-N08) are appended from the older
ledger results/20260911_CLAIM_MATRIX.csv, which the final ledger keeps in force unchanged.

Two kinds of reproduction are kept apart:
  independent_phasor_reproduction  (from the final ledger; e.g. V21 = ANDES Phase 8, GATE 5 PASS)
  EMT_result                       (this campaign; prereg v1: G1 PASS, G3 FAIL, G4 FAIL, EMT04-18 BLOCKED)
An EMT failure never removes an independent phasor reproduction.
The first version of this file (commit 17832edf) was generated from the pre-final ledger and wrongly
carried I04 as NOT YET TESTED; see docs/20260911_PAREMT_EMT_LEDGER_RECONCILIATION.md.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

RESEARCH = Path(__file__).resolve().parents[2]
RES = RESEARCH / "results"
FINAL = RES / "20260911_FINAL_VALIDATION_MATRIX.csv"
OLD = RES / "20260911_CLAIM_MATRIX.csv"

BLOCK = "BLOCKED (prereg v1 s10: G3 and G4 failed; not run)"
KEEP = "Keep the final-ledger wording; add no EMT support."


def theorem():
    return ("NOT APPLICABLE (mathematical result; prereg v1 s0)", "none", "THEOREM", "n/a", "N/A",
            "Mathematical statement; EMT is not a proof instrument.", "none")


def na(why):
    return (f"NOT APPLICABLE ({why})", "none", "NOT APPLICABLE TO EMT", "n/a", "N/A", "Status unchanged.", "none")


def blocked(appl, exp, interp, impact=KEEP):
    return (appl, exp, "EMT UNRESOLVED", BLOCK, "BLOCKED", interp, impact)


def not_targeted(interp):
    return ("WEAK (not an EMT target in prereg v1)", "none", "EMT UNRESOLVED", "not an EMT target in prereg v1", "N/A", interp, "none")


def old_ids_covered(final_rows):
    cov = set()
    for r in final_rows:
        for tok in re.findall(r"[TCBIN]\d\d(?:-[TCBIN]\d\d)?", r["old_ids"]):
            if "-" in tok:
                a, b = tok.split("-")
                cov |= {a[0] + f"{i:02d}" for i in range(int(a[1:]), int(b[1:]) + 1)}
            else:
                cov.add(tok)
    return cov


def main() -> int:
    s1 = json.loads((RES / "EMT01" / "EMT01_summary.json").read_text())
    s2 = json.loads((RES / "EMT02" / "EMT02_summary.json").read_text())
    s3 = json.loads((RES / "EMT03" / "EMT03_summary.json").read_text())
    d2 = json.loads((RES / "EMT02" / "EMT02_diagnostics.json").read_text())
    d3 = json.loads((RES / "EMT03" / "EMT03_diagnostics.json").read_text())
    with (RES / "EMT02" / "EMT02_trajectory_comparison.csv").open(encoding="utf-8") as f:
        tr2 = {r["state"]: float(r["ratio"]) for r in csv.DictReader(f)}
    worst3 = {}
    with (RES / "EMT03" / "EMT03_trajectory_comparison.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if float(r["max_abs_err"]) < float("inf"):
                worst3[r["test"]] = max(worst3.get(r["test"], 0.0), float(r["ratio"]))
    eq = s1["equilibrium"]
    g3 = (f"G3: efd error {tr2['efd']:.4f} of excursion (> 0.02, FAIL); other 6 states <= {max(v for k, v in tr2.items() if k != 'efd'):.4f}; "
          f"ringdown d_alpha {s2['ringdown']['d_alpha']:.1e} s^-1, d_f {s2['ringdown']['d_f']:.1e} Hz; post hoc: algebraic-line D1 "
          f"{d2['D1_quasi_static_line_50us_vs_quasi_static_phasor']['max_ratio']:.1e}, vs dynamic-phasor D3 {d2['D3_emt_line_50us_vs_dynamic_phasor']['max_ratio']:.1e}")
    g4 = (f"G4: {s3['n_cases_nonfinite']}/9 cases non-finite at {s3['earliest_nonfinite_t_s']} s (g = 0); 6/9 finite with worst-state errors "
          f"{min(worst3.values()):.0f}-{max(worst3.values()):.0f}x the excursion; post hoc algebraic-line D1: all 9 finite, max ratio "
          f"{d3['D1_quasi_static_line_all_9_cases']['max_ratio']:.4f}")

    emt = {
        "V01": theorem(),
        "V02": theorem(),
        "V03": blocked("STRONG", "EMT05, EMT06", "kappa_EMT at P4 unknown. The independent phasor reproduction (ANDES Phase 8, 32/32) stands; it is not an EMT result.",
                       "Cite the ANDES reproduction as a reproduction of computation; do not cite EMT; the witness is not robust (V16)."),
        "V04": blocked("STRONG", "EMT05 (all proper subsets)", "Screening approval of H4 has no EMT evidence."),
        "V05": na("F10 baseline recount on the phasor spectra; not an EMT target"),
        "V06": blocked("STRONG (observable: verdict change at identical static quantities)", "EMT07-08, EMT09",
                       "The by-construction part needs no EMT; the observable g-only transition has no EMT evidence."),
        "V07": blocked("STRONG", "EMT07-08, EMT10, EMT11, EMT16", "Policy-dependent minimal incompatibility has no EMT evidence; ANDES Phase 8 reproduces P4 vs G_S."),
        "V08": theorem(),
        "V09": na("connected cumulants are secondary; no EMT observable, prereg v1 s0"),
        "V10": blocked("STRONG (observable consequences: boundary locations, Newton restoration)", "EMT08, EMT10, EMT11, EMT12",
                       "The derivative identity is a theorem; its EMT observable consequences are untested."),
        "V11": blocked("STRONG", "EMT12", "Line ranking has no EMT evidence; transfer across converter models is not claimed (V28)."),
        "V12": na("negative scope limit of the linear port prediction; EMT10 explicitly does not claim finite-step magnitudes"),
        "V13": blocked("STRONG", "EMT08, EMT10, EMT11, EMT13, EMT14", "No remediation family has EMT evidence."),
        "V14": blocked("STRONG", "EMT14", "Governor dependence has no EMT evidence."),
        "V15": blocked("SECONDARY", "EMT16", "Phenomenon robustness has no EMT evidence."),
        "V16": blocked("SECONDARY", "EMT16", "Witness non-robustness has no EMT evidence; the phasor finding stands."),
        "V17": blocked("SECONDARY", "EMT16", "No EMT evidence."),
        "V18": not_targeted("g* robustness across draws is not part of EMT16 (prereg v1)."),
        "V19": ("STRONG", "EMT01", "EMT QUANTITATIVELY REPRODUCED",
                f"Ybus rel {s1['ybus']['50us']['rel_residual_continuous']:.1e} (<= 1e-8); max|dV| {eq['max_dV_pu']:.2e} pu (<= 1e-4); "
                f"max|dtheta| {eq['max_dAngle_rad']:.2e} rad (<= 1e-3); dP {eq['max_dP_pu']:.1e}, dQ {eq['max_dQ_pu']:.1e} pu",
                "PASS (G1)",
                "A third, EMT, realization of the frozen network and operating point (ParaEMT companion elements, documented tap patch). "
                "The final-ledger status CONDITIONAL (pandapower translation, memo A) is unchanged.",
                "Optional sentence: the network and operating point were also realized in ParaEMT."),
        "V20": ("STRONG", "EMT02 (SG); EMT12 (sensitivities)", "EMT UNRESOLVED", g3 + "; EMT12 BLOCKED", "FAIL (G3); BLOCKED",
                "The preregistered EMT SG gate failed on efd; the post hoc diagnostics attribute the excess to the line's (X/w0)dI/dt, "
                "not to the transcription. The independent phasor reproduction (ANDES R2/R3, E1) is unaffected.",
                "Keep V20 as VALIDATED (by ANDES); no EMT support."),
        "V21": ("STRONG (central EMT target)", "EMT03 -> EMT05-EMT16", "EMT UNRESOLVED", g4, "FAIL (G4) -> BLOCKED",
                "The EMT realization failed its preregistered interface gate (ideal current source + lossless trapezoidal line = Nyquist chatter); "
                "the transcription itself is verified post hoc. This does NOT affect V21: the custom GFL results were independently reproduced "
                "in ANDES with the same equations (Phase 8, GATE 5 PASS).",
                "Keep V21 = VALIDATED as a reproduction of computation in a second phasor-domain tool; EMT reproduction: not achieved in v1."),
        "V22": theorem(),
        "V23": na("|chi| is a phasor-port diagnostic; secondary"),
        "V24": na("refuted claim; EMT is not used to reinterpret failed claims"),
        "V25": na("phasor-port closure diagnostic"),
        "V26": theorem(),
        "V27": na("refuted claim; EMT is not used to reinterpret failed claims"),
        "V28": na("EMT12 uses the same custom GFL, so it cannot test transfer across converter models"),
        "V29": blocked("SECONDARY", "EMT14", "The phasor refutation stands; EMT adds nothing.", "none"),
        "V30": blocked("STRONG/SECONDARY (nonlinear scope)", "EMT15, EMT05", "No EMT evidence for the nonlinear scope."),
    }
    emt_old = {
        "T04": theorem(),
        "T06": theorem(),
        "T08": ("SECONDARY (observable: sign change of the deflated zero-frequency port; Kundur, EMT17 optional)", "EMT17", "THEOREM",
                "not run (optional; campaign BLOCKED)", "N/A",
                "The theorem stands on its proof. Phasor holdout: 28/29 real-count changes under the preregistered rule "
                "(the miss is a real-pair coalescence, not an origin crossing), i.e. 28/28 crossings through the origin; 0/445 false positives.",
                "Quote 28/29 as the preregistered score (28/28 origin crossings as the explained restriction)."),
        "C08": na("connected-cumulant benchmark result; secondary"),
        "B04": blocked("SECONDARY", "EMT05 (f_EMT vs band_hz)", "Mode identity not checked in EMT."),
        "B06": not_targeted("A planning consequence of the verdict table; no EMT evidence."),
        "B07": blocked("STRONG/SECONDARY (nonlinear scope)", "EMT15", "No EMT evidence."),
        "B10": blocked("SECONDARY", "EMT15 (indirect)", "No EMT evidence."),
        "N01": na("refuted claim"),
        "N05": na("refuted claim"),
        "N06": na("refuted claim"),
        "N07": na("refuted claim"),
        "N08": na("refuted claim"),
    }

    with FINAL.open(encoding="utf-8") as f:
        final_rows = list(csv.DictReader(f))
    with OLD.open(encoding="utf-8") as f:
        old_rows = list(csv.DictReader(f))
    assert sorted(r["id"] for r in final_rows) == sorted(emt)
    missing = [r["id"] for r in old_rows if r["id"] not in old_ids_covered(final_rows)]
    assert sorted(missing) == sorted(emt_old), (missing, sorted(emt_old))

    cols = ["claim_id", "old_ids", "ledger_source", "canonical_claim", "proof_status", "existing_phasor_evidence",
            "independent_phasor_reproduction", "EMT_applicability", "EMT_experiment", "EMT_result", "quantitative_error",
            "pass_fail_unresolved", "interpretation", "manuscript_impact"]
    with (RES / "20260911_EMT_CLAIM_MATRIX.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for r in final_rows:
            w.writerow([r["id"], r["old_ids"], "FINAL_VALIDATION_MATRIX", r["claim"], r["status"], r["ieee39_validation"],
                        r["independent_tool_validation"], *emt[r["id"]]])
        for r in old_rows:
            if r["id"] in emt_old:
                w.writerow([r["id"], r["id"], "CLAIM_MATRIX (not cross-referenced by the final ledger; kept in force)", r["claim"],
                            r["final_status"], r["ieee39_validation"], r["independent_tool_validation"], *emt_old[r["id"]]])
    print(f"{len(final_rows)} final-ledger claims + {len(emt_old)} old-ledger claims written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
