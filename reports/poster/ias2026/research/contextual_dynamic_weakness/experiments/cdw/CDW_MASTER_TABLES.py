# ruff: noqa: E501
"""Phase 27: master result table, master claim matrix and run manifest from gate JSONs."""

from __future__ import annotations

import csv
import json

import _infra as I


def g(name):
    p = I.RESULTS / name
    return json.loads(p.read_text()) if p.exists() else {}


def status_of(ok, available=True, blocked=False):
    if blocked:
        return "BLOCKED"
    if not available:
        return "NOT_RUN"
    return "SUPPORTED" if ok else "NOT_SUPPORTED"


def main():
    e1, e2, e34, e6, e7, e8, e9 = (g(f) for f in ("CDW_E1_gates.json", "CDW_E2_gates.json", "CDW_E34_gates.json", "CDW_E6_gates.json",
                                               "CDW_E7_gates.json", "CDW_E8_gates.json", "CDW_E9_gates.json"))
    e11, e12, gc, e16, e17, e23, e5 = (g(f) for f in ("CDW_E11_gates.json", "CDW_E12_gates.json", "CDW_E13_E14_gates.json",
                                                   "CDW_E16_gates.json", "CDW_E17_gates.json", "CDW_E23_gates.json", "CDW_E5_gates.json"))
    code = "contextual_dynamic_weakness/experiments/cdw/"
    rows = [
        ("C-H1", "Contextual sign reversal recurs across holdout policies (stable contexts)", status_of(e1.get("A1_pass"), bool(e1)),
         f"A1={e1.get('A1_frac')}", code + "E01_census.py", "results/CDW_E1_policy_summary.csv", "F1,F2"),
        ("C-H1t", "Reversal survives mode tracking", status_of(e1.get("A2_pass"), bool(e1)), f"A2={e1.get('A2_frac')}", code + "E01_census.py",
         "results/CDW_E1_summary.csv", "F1"),
        ("C-H1b", "Best fixed node-only ranking is materially insufficient", status_of(e1.get("A3_pass"), bool(e1)),
         f"median p*={e1.get('A3_median_pstar')}, median regret frac={e1.get('A3_median_frac_regret')}", code + "E01_census.py",
         "results/CDW_E1_policy_summary.csv", "F1"),
        ("C-H1c", "Contextuality persists under envelopes (A4)", status_of(e2.get("A4_pass"), bool(e2)), f"{e2.get('A4_cov_by_env')}",
         code + "E02_robust.py", "results/CDW_E2_robust_summary.csv", "F5"),
        ("GOLD-A", "Contextual sign reversal invalidates a single intrinsic weak-node ranking",
         status_of(bool(e1.get("A1_pass") and e1.get("A3_pass") and e2.get("A4_pass")), bool(e1 and e2)), "A1 & A3 & A4", code + "E01_census.py; E02_robust.py",
         "results/CDW_E1_gates.json; CDW_E2_gates.json", "F1,F2,F5"),
        ("C-IV", "IFT total derivatives match finite re-equilibrated differences", status_of(e34.get("IV_pass"), bool(e34)),
         f"frac ok={e34.get('IV_frac_ok')}", code + "_sens.py; E34_sens.py", "results/CDW_E3_node_sensitivity.parquet", "F4"),
        ("C-R1", "Structural: frozen = total for controller-only coordinates (SPR); frozen = 0 for setpoints",
         status_of(bool(e34.get("R1_pass") and e34.get("R2_pass")), bool(e34)), f"R1 max rel={e34.get('R1_max_rel_frozen_minus_total')}",
         code + "E34_sens.py", "results/CDW_E34_gates.json", ""),
        ("C-H2n-SPR", "Frozen vs total node sensitivities disagree materially (SPR)", status_of(e34.get("H2_node_SPR_true"), bool(e34)),
         f"{e34.get('H2_node_SPR_frac_material_holdout')}", code + "E34_sens.py", "results/CDW_E34_H2_node_SPR.csv", ""),
        ("C-H2n-RP", "Frozen vs total node sensitivities disagree materially (RP)", status_of(e34.get("H2_node_RP_true"), bool(e34)),
         f"{e34.get('H2_node_RP_frac_material_holdout')}", code + "E34_sens.py", "results/CDW_E34_H2_node_RP.csv", ""),
        ("C-H2l-SPR", "Frozen vs total link sensitivities disagree materially (SPR)", status_of(e34.get("H2_line_SPR_true"), bool(e34)),
         f"{e34.get('H2_line_SPR_frac_material_holdout')}", code + "E34_sens.py", "results/CDW_E34_H2_line_SPR.csv", "F4"),
        ("C-H2l-RP", "Frozen vs total link sensitivities disagree materially (RP)", status_of(e34.get("H2_line_RP_true"), bool(e34)),
         f"{e34.get('H2_line_RP_frac_material_holdout')}", code + "E34_sens.py", "results/CDW_E34_H2_line_RP.csv", ""),
        ("GOLD-B", "Total dynamic link measure beats the best static baseline on held-out finite effects", status_of(e34.get("GB_pass"), bool(e34)),
         f"diff={e34.get('GB_median_diff')}, rho={e34.get('GB_median_rho_Dtotal')}, top5={e34.get('GB_median_top5_Dtotal')}",
         code + "E34_sens.py", "results/CDW_E4_baselines.csv", "F4"),
        ("C-H3", "Robust dynamic weak corridors exist", status_of(e6.get("H3_true"), bool(e6)), f"{e6.get('robust_weak_corridors')}",
         code + "E06_corridors.py", "results/CDW_E6_corridors.csv", "F6"),
        ("C-H4a", "A single topology action removes the P4 incompatibility", status_of(bool(e7.get("removes_H4_at_P4")), bool(e7)),
         f"{e7.get('removes_H4_at_P4')}", code + "E07_topology.py", "results/CDW_E7_topology.parquet", "F7"),
        ("C-H4b", "Static topology scores predict the Δα of topology actions", status_of(any((e7.get("static_predicts") or {}).values()), bool(e7)),
         f"{e7.get('static_median_rho')}", code + "E07_topology.py", "results/CDW_E7_static_prediction.csv", "F7"),
        ("C-EX", "Local stability-equivalent exchange rates are first-order valid (BSTAR, small)", status_of(e8.get("E8_pass_BSTAR_small"), bool(e8)),
         f"{e8.get('BSTAR_small_frac_ratio_le_0.2')}", code + "E08_exchange.py", "results/CDW_E8_exchange.csv", "F8"),
        ("GOLD-D", "Plan-level design safe where single-boundary tuning leaves an unsafe subset", status_of(e9.get("GOLD_D_pass"), bool(e9)),
         f"{len(e9.get('cases', []))} cases", code + "E09_design.py", "results/CDW_E9_design.csv", "F9"),
        ("C-H6s", "Static graph scores predict contextual/dynamic node quantities", status_of(e11.get("static_graph_predicts_any"), bool(e11)),
         "", code + "E11_spectral.py", "results/CDW_E11_spectral.csv", "F10"),
        ("C-H6m1", "Controller policy materially changes graph-modal mixing at fixed network", status_of(e12.get("Q1_material"), bool(e12)),
         f"{e12.get('Q1_rel_range_mu_H4')}", code + "E11_spectral.py", "results/CDW_E12_mixing.csv", "F10"),
        ("C-H6m2", "Modal mixing correlates with contextual reversal", status_of(e12.get("Q2_correlates"), bool(e12)),
         f"rho={e12.get('Q2_spearman')}, p={e12.get('Q2_perm_p')}", code + "E11_spectral.py", "results/CDW_E12_mixing.csv", "F10"),
        ("C-H6m3", "Witness changes correspond to reproducible modal-support transitions", status_of(e12.get("Q3_reproducible"), bool(e12)),
         f"{e12.get('Q3_events_with_support_change')}/{e12.get('Q3_n_events')}", code + "E11_spectral.py", "results/CDW_E12_q3_support.csv", ""),
        ("GOLD-C", "Genuinely reduced model preserves decisions with speedup and certificate/abstention", status_of(gc.get("GOLD_C_pass"), bool(gc)),
         "", code + "E13_reduction.py", "results/CDW_E13_reduction.csv; CDW_E14_summary.csv", "F11,F12"),
        ("C-H8", "Energy-dominant observed mode often differs from the limiting mode (phasor)", status_of(e16.get("systematic"), bool(e16)),
         f"{e16.get('frac_differ_all')}", code + "E16_modal_energy.py", "results/CDW_E16_modal_energy.csv", "F13"),
        ("C-H9", "Findings transfer to the static-injection converter model", status_of(e23.get("H9_transfers"), bool(e23)),
         f"{ {k: e23.get(k) for k in ('reversal_transfers', 'links_transfer', 'corridors_transfer')} }", code + "E23_crossmodel.py",
         "results/CDW_E23_*.csv", ""),
        ("GOLD-E", "Static PV weak-node ranking fails to capture coalition/dynamic hosting (Africano)", "BLOCKED", "E17 gate failed", code + "E17_africano_gate.py",
         "docs/CDW_AFRICANO_MISSING_INPUTS.md", ""),
        ("GOLD-F", "No-storage control/topology co-design increases PV hosting", "BLOCKED", "E17 gate failed", code + "E17_africano_gate.py",
         "docs/CDW_AFRICANO_MISSING_INPUTS.md", ""),
        ("C-NL", "Plan-level tuning that improves alpha may reduce nonlinear recovery margin", "NOT_RUN", "optional E22 deferred", "", "", ""),
    ]
    head = I.git_head()
    with (I.RESULTS / "CDW_MASTER_CLAIM_MATRIX.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["claim_id", "claim", "status", "evidence", "code", "data", "figures", "seed_or_inputs", "commit"])
        for r in rows:
            w.writerow([*r, "prereg_inputs (seeds 20260920/20260921/20260917/20260922/20260923/20260924)", head])
    # master result table: flatten gates
    with (I.RESULTS / "CDW_MASTER_RESULT_TABLE.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["experiment", "metric", "value"])
        for name, d in (("E1", e1), ("E2", e2), ("E34", e34), ("E5", e5), ("E6", e6), ("E7", e7), ("E8", e8), ("E9", e9), ("E11", e11),
                        ("E12", e12), ("E13E14", gc), ("E16", e16), ("E17", e17), ("E23", e23)):
            for k, v in d.items():
                w.writerow([name, k, json.dumps(v, default=str)[:500]])
    st = I.load_status()
    with (I.RESULTS / "CDW_MASTER_RUN_MANIFEST.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["phase", "state", "wall_s", "n_tasks", "n_errors", "raw_files", "log", "started", "updated"])
        for ph, rec in st.get("phases", {}).items():
            raw = I.RAW / ph
            w.writerow([ph, rec.get("state"), rec.get("wall_s"), (rec.get("tasks") or {}).get("n_tasks"), rec.get("n_errors"),
                        len(list(raw.glob("*.json"))) if raw.exists() else 0, f"logs/{ph}.log", rec.get("started"), rec.get("updated")])
    print("master tables written")


if __name__ == "__main__":
    main()
