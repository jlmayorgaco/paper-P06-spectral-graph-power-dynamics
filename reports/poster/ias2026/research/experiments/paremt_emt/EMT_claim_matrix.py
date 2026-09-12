# ruff: noqa: E501  -- claim table kept one row per line
"""Claim-to-EMT matrix (campaign request section 25) -> results/20260911_EMT_CLAIM_MATRIX.csv.

Claim ids, canonical text and proof status come from the frozen ledger results/20260911_CLAIM_MATRIX.csv
(read, never modified). EMT columns record the preregistered v1 outcome: EMT00/EMT01 PASS, EMT02 (G3)
FAIL, EMT03 (G4) FAIL, hence EMT04-EMT18 BLOCKED. EMT_result uses the six-label taxonomy of prereg s0.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

RESEARCH = Path(__file__).resolve().parents[2]
RES = RESEARCH / "results"

BLOCK = "BLOCKED (prereg s10: G3 and G4 failed; not run)"
THEOREM_NA = ("NOT APPLICABLE (mathematical result; prereg s0)", "none", "THEOREM", "n/a", "N/A",
              "Mathematical statement; EMT is not a proof instrument.", "none")
CUM_NA = ("NOT APPLICABLE (connected cumulants are secondary; no EMT target, prereg s0)", "none", "NOT APPLICABLE TO EMT", "n/a", "N/A",
          "An algebraic decomposition of the phasor port; it has no EMT observable of its own.", "none")
NEG_NA = ("NOT APPLICABLE (negative phasor result; EMT is not used to reinterpret failed claims)", "none", "NOT APPLICABLE TO EMT", "n/a", "N/A",
          "Status unchanged; EMT cannot rescue or re-open a refuted claim.", "none")


def blocked(appl, exp, interp, impact="Keep the current phasor-only wording; add no EMT support."):
    return (appl, exp, "EMT UNRESOLVED", BLOCK, "BLOCKED", interp, impact)


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
    efd = [r for r in d2["D2_emt_line_25us_vs_quasi_static_phasor"]["per_state"] if r["state"] == "efd"][0]["ratio"]
    rows_emt = {}
    for cid in ("T01", "T02", "T03", "T04", "T05", "T06", "T07"):
        rows_emt[cid] = THEOREM_NA
    rows_emt["T08"] = ("SECONDARY (observable: sign change of the deflated zero-frequency port; Kundur, EMT17 optional)", "EMT17",
                       "THEOREM", "not run (optional; campaign BLOCKED at G3/G4)", "N/A",
                       "The theorem stands on its proof; the optional EMT17 corroboration was not attempted.", "none")
    rows_emt["T09"] = ("NOT APPLICABLE (the Sylvester identity is classical; the truncation anatomy is a phasor-port observation)", "none",
                       "NOT APPLICABLE TO EMT", "n/a", "N/A", "Neither the identity nor the port anatomy has an EMT observable.", "none")
    rows_emt["T10"] = blocked("STRONG (observable consequences: boundary locations and finite line changes)", "EMT08, EMT10, EMT11, EMT12",
                              "The derivative formula is a theorem (implicit function theorem); its EMT observable consequences are untested.")
    for cid in ("C01", "C02", "C03", "C04", "C05", "C06"):
        rows_emt[cid] = THEOREM_NA
    for cid in ("C07", "C08", "C09"):
        rows_emt[cid] = CUM_NA
    rows_emt["C10"] = NEG_NA
    rows_emt["C11"] = NEG_NA
    rows_emt["B01"] = blocked("STRONG", "EMT07-08 (g crossing), EMT10, EMT11, EMT16", "Policy-dependent minimal incompatibility has no EMT evidence.")
    rows_emt["B02"] = blocked("STRONG", "EMT05, EMT06, EMT09", "kappa_EMT at P4 is unknown; the flagship remains phasor (and ANDES-variant) evidence only.",
                              "Do not cite EMT for the P4 trap; {30,33,35,37} is not claimed robust.")
    rows_emt["B03"] = blocked("STRONG", "EMT14", "The governed stabilization has no EMT evidence.")
    rows_emt["B04"] = blocked("SECONDARY (f_EMT vs band_hz in EMT05; no dispatch-control EMT test preregistered)", "EMT05", "Mode identity is not checked in EMT.")
    rows_emt["B05"] = blocked("STRONG", "EMT13", "Condenser remediation has no EMT evidence.")
    rows_emt["B06"] = ("WEAK (a planning consequence of the verdict table; not preregistered)", "none", "EMT UNRESOLVED",
                       "not an EMT target in prereg v1", "N/A", "No EMT evidence.", "none")
    for cid in ("B07", "B08", "B09", "B10"):
        rows_emt[cid] = blocked("STRONG/SECONDARY (nonlinear scope)", "EMT15" + (", EMT05" if cid == "B09" else ""),
                                "The nonlinear finite-disturbance scope has no EMT evidence.")
    rows_emt["I01"] = ("STRONG", "EMT01", "EMT QUANTITATIVELY REPRODUCED",
                       f"Ybus rel {s1['ybus']['50us']['rel_residual_continuous']:.1e} (<= 1e-8); max|dV| {eq['max_dV_pu']:.2e} pu (<= 1e-4); "
                       f"max|dtheta| {eq['max_dAngle_rad']:.2e} rad (<= 1e-3); dP {eq['max_dP_pu']:.1e}, dQ {eq['max_dQ_pu']:.1e} pu; "
                       f"trapezoidal 60-Hz warp {s1['ybus']['50us']['rel_warp_discrete_vs_continuous']:.1e} at 50 us (reported separately)",
                       "PASS (G1)",
                       "ParaEMT's own companion elements (with the documented tap patch; upstream never applies xfmr_k) realize the frozen 60-Hz network and equilibrium. "
                       "Official ParaEMT IEEE-39 differs from TX4 in 24/195 network rows (taps, loads, ratings).",
                       "I01 gains a third independent network realization (EMT).")
    rows_emt["I02"] = ("STRONG", "EMT02", "EMT UNRESOLVED",
                       f"efd error {tr2['efd']:.4f} of excursion (> 0.02, FAIL); other 6 states <= {max(v for k, v in tr2.items() if k != 'efd'):.4f}; equilibrium dP {s2['equilibrium']['dP']:.1e}, dQ {s2['equilibrium']['dQ']:.1e}; "
                       f"ringdown d_alpha {s2['ringdown']['d_alpha']:.1e} s^-1, d_f {s2['ringdown']['d_f']:.1e} Hz (both <= 0.005). "
                       f"Diagnostics: algebraic-line D1 {d2['D1_quasi_static_line_50us_vs_quasi_static_phasor']['max_ratio']:.1e}; 25 us D2 efd {efd:.4f}; "
                       f"vs dynamic-phasor reference D3 {d2['D3_emt_line_50us_vs_dynamic_phasor']['max_ratio']:.1e}",
                       "FAIL (G3)",
                       "The preregistered SG gate failed on efd. Post hoc diagnostics (not gates) show that the machine transcription matches to 2e-4, "
                       "and that the excess comes from the EMT line's (X/w0)dI/dt, a physical effect of size f/f0 = 2 % that the quasi-static phasor model omits.",
                       "No EMT support for I02; the ANDES evidence is unchanged.")
    rows_emt["I03"] = blocked("STRONG", "EMT12", "Branch sensitivities have no EMT evidence.")
    rows_emt["I04"] = ("STRONG (central target of the campaign)", "EMT03 -> EMT05-EMT16", "EMT UNRESOLVED",
                       f"G4: {s3['n_cases_nonfinite']}/9 cases non-finite at {s3['earliest_nonfinite_t_s']} s (g = 0); 6/9 finite but with worst-state errors {min(worst3.values()):.0f}-{max(worst3.values()):.0f}x the excursion "
                       f"(step-to-step alternation up to {d3['D4_chatter_preregistered_line']['alternating_amplitude_pu']['t=0.1']:.2f} pu); "
                       f"algebraic-line diagnostic D1: all 9 cases finite, max ratio {d3['D1_quasi_static_line_all_9_cases']['max_ratio']:.4f}",
                       "FAIL (G4) -> BLOCKED",
                       "The GFL transcription is verified (D1), but the preregistered current-source realization is numerically unstable in the preregistered test topology. "
                       "No IEEE-39 portfolio result may be called equation-equivalent.",
                       "I04 stays NOT YET TESTED in any second simulator; the manuscript must keep the custom-GFL results as single-implementation.")
    for cid in ("N01", "N02", "N05", "N06", "N07", "N08"):
        rows_emt[cid] = NEG_NA
    rows_emt["N03"] = blocked("SECONDARY (EMT16 witness variation)", "EMT16", "The phasor refutation stands; EMT adds nothing.", "none")
    rows_emt["N04"] = blocked("SECONDARY (EMT14 governed flagship)", "EMT14", "The phasor refutation stands; EMT adds nothing.", "none")
    rows_emt["N09"] = ("NOT APPLICABLE (EMT12 uses the same custom GFL, so it cannot test transfer across converter models)", "none",
                       "NOT APPLICABLE TO EMT", "n/a", "N/A", "Status unchanged (REFUTED).", "none")

    cols = ["claim_id", "canonical_claim", "proof_status", "existing_phasor_evidence", "EMT_applicability", "EMT_experiment",
            "EMT_result", "quantitative_error", "pass_fail_unresolved", "interpretation", "manuscript_impact"]
    with (RES / "20260911_CLAIM_MATRIX.csv").open(encoding="utf-8") as f:
        ledger = list(csv.DictReader(f))
    ids = [r["id"] for r in ledger]
    assert sorted(ids) == sorted(rows_emt), set(ids) ^ set(rows_emt)
    with (RES / "20260911_EMT_CLAIM_MATRIX.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for r in ledger:
            e = rows_emt[r["id"]]
            w.writerow([r["id"], r["claim"], r["final_status"], r["evidence"], *e])
    print(f"{len(ledger)} claims written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
