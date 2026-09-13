# ruff: noqa: E501
"""Emit every number quoted in the paper/report as a LaTeX macro, straight from the result
JSON/CSV files (no transcription). Output: reports/papers/cdw_contextual_dynamic_weakness/cdw_numbers.tex
and results/hardening/CDWH_numbers.json."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json
import re

import pandas as pd

R = HI.RESULTS
PAPER = HI.I.REPO / "reports" / "papers" / "cdw_contextual_dynamic_weakness"


def j(n):
    p = R / n
    return json.loads(p.read_text()) if p.exists() else {}


def f2(x, d=2):
    return "n/a" if x is None or (isinstance(x, float) and x != x) else f"{x:.{d}f}"


def pct(x):
    return "n/a" if x is None or x != x else f"{100 * x:.0f}\\%"


def build():
    N = {}
    h3, h4, h5, h7, h13, h15, h10, h18, h19 = (j("H03_gate.json"), j("H04_gate.json"), j("H05_stats.json"), j("H07_goldb_gate.json"),
                                               j("H13_mixing.json"), j("H15_topology_gate.json"), j("H10_gate.json"), j("H18_gate.json"), j("H19_evidence.json"))
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    for s in ("new", "old", "discovery"):
        q = pol[(pol.split == s) & pol.base_stable]
        tag = {"new": "New", "old": "Old", "discovery": "Disc"}[s]
        N[f"n{tag}Stable"] = str(len(q))
        N[f"n{tag}Unstable"] = str(int(((pol.split == s) & ~pol.base_stable).sum()))
        for lvl in "ABCD":
            N[f"cov{lvl}{tag}"] = f2(h3.get(f"{s}_frac_{lvl}"))
            N[f"cnt{lvl}{tag}"] = f"{int(q[f'{lvl}_has'].sum())}/{len(q)}"
        N[f"covArb{tag}"] = f2(h3.get(f"{s}_frac_arbitrary_global"))
        N[f"cntArb{tag}"] = f"{int((q.n_rev_arbitrary > 0).sum())}/{len(q)}"
        N[f"bothC{tag}"] = f2(h3.get(f"{s}_frac_C_both_dirs"))
        for st in ("FULL", "EM", "SAME"):
            N[f"pstar{st}{tag}"] = f2(h4.get(f"{s}_{st}_median_pstar"))
            N[f"regret{st}{tag}"] = f2(h4.get(f"{s}_{st}_median_frac_regret"))
    N["newCpairs"] = str(h3.get("new_C_pairs"))
    dc = h3.get("new_C_dir_counts", {})
    N["newCsd"], N["newCds"] = str(dc.get("s2d", 0)), str(dc.get("d2s", 0))
    q = h3.get("new_C_mag_quartiles", [None] * 3)
    N["newCmagQone"], N["newCmagMed"], N["newCmagQthree"] = (f2(v, 3) for v in q)
    N["newCemClean"] = pct(h3.get("new_C_frac_em_clean_witness"))
    b = h5.get("new_cov_C", [None] * 3)
    N["covCNewLo"], N["covCNewHi"] = f2(b[1]), f2(b[2])
    b = h5.get("new_FULL_pstar", [None] * 3)
    N["pstarNewLo"], N["pstarNewHi"] = f2(b[1]), f2(b[2])
    b = h5.get("new_FULL_regret", [None] * 3)
    N["regretNewLo"], N["regretNewHi"] = f2(b[1]), f2(b[2])
    sh = h5.get("shift_pstar_new_minus_old", [None] * 3)
    N["shiftPstar"], N["shiftPstarLo"], N["shiftPstarHi"] = f2(sh[0], 3), f2(sh[1], 3), f2(sh[2], 3)
    sh = h5.get("shift_regret_new_minus_old", [None] * 3)
    N["shiftRegret"], N["shiftRegretLo"], N["shiftRegretHi"] = f2(sh[0], 3), f2(sh[1], 3), f2(sh[2], 3)
    t1 = h5.get("T1_sign_test", {})
    N["tOneK"], N["tOneN"] = str(t1.get("k_below_0.90")), str(t1.get("n"))
    l7 = h5.get("L7_cov_by_env_global", {})
    for e, v in l7.items():
        N["lSeven" + re.sub("[^A-Za-z]", "", e)] = f2(v)
    for e, v in h5.get("L7_cov_by_env_EM_samemode", {}).items():
        N["lSevenEM" + re.sub("[^A-Za-z]", "", e)] = f2(v)
    N["transDiagAcc"], N["transOffAcc"] = f2(h4.get("transfer_median_diag_acc")), f2(h4.get("transfer_median_offdiag_acc"))
    N["transDiagReg"], N["transOffReg"] = f2(h4.get("transfer_new_median_diag_regret")), f2(h4.get("transfer_new_median_offdiag_regret_by_test"))
    N["transN"], N["distRho"] = str(h4.get("transfer_n_policies")), f2(h4.get("distance_vs_degradation_spearman"))
    N["nOrders"] = str(h4.get("n_distinct_optimal_orders_full"))
    ra = h4.get("regret_anatomy_new", {})
    N["regretFastShare"] = pct(ra.get("share_out_of_band"))
    N["regretNbad"], N["regretNfast"] = str(ra.get("material")), str(ra.get("choice_out_of_EM_band"))
    # GOLD-B
    for k, v in h7.items():
        if k.startswith("GB_") and isinstance(v, (int, float)) and not isinstance(v, bool):
            kk = (k.replace("GB_", "gb").replace("1.10", "OneTen").replace("1.25", "OneTwentyFive").replace("top5", "topfive")
                  .replace("top3", "topthree").replace("ndcg5", "ndcg").replace("_S1", "_Sone").replace("V9", "Vnine").replace("T2", "Ttwo"))
            name = re.sub("[^A-Za-z]", "", kk)
            assert name not in N, ("macro collision", name, k)
            N[name] = f2(v, 3) if abs(v) < 10 else str(v)
    if "GB_primary_diff_ci95" in h7:
        N["gbPrimaryCiLo"], N["gbPrimaryCiHi"] = f2(h7["GB_primary_diff_ci95"][0], 3), f2(h7["GB_primary_diff_ci95"][1], 3)
    if "GB_prereg_diff_ci95" in h7:
        N["gbPreregCiLo"], N["gbPreregCiHi"] = f2(h7["GB_prereg_diff_ci95"][0], 3), f2(h7["GB_prereg_diff_ci95"][1], 3)
    N["nElig"] = str(h7.get("GB_prereg_n"))
    N["nEligClusters"] = str(h7.get("GB_prereg_n_clusters"))
    m_, e_ = f"{h7.get('R0_max_all', float('nan')):.1e}".split("e")
    N["rZeroMax"] = "$" + m_ + "\\times10^{" + str(int(e_)) + "}$"
    N["ivN"] = str(h7.get("IV_n"))
    rat = h7.get("node_opoint_median_ratio_indirect_over_frozen", {})
    dom = h7.get("node_opoint_frac_indirect_dominates", {})
    N["nodeLoadRatio"], N["nodeLoadDom"] = f2(rat.get("load")), pct(dom.get("load"))
    N["nodeVsetDom"] = pct(dom.get("vset"))
    fh = h19.get("F_conf_holm_p", {})
    fr = h19.get("F_conf_raw_p", {})
    N["tOneP"], N["tOnePholm"] = f"{fr.get('T1_goldA_sign', float('nan')):.1e}", f"{fh.get('T1_goldA_sign', float('nan')):.1e}"
    N["tTwoP"], N["tTwoPholm"] = f2(fr.get("T2_goldB_signflip"), 3), f2(fh.get("T2_goldB_signflip"), 3)
    N["tThreeP"], N["tThreePholm"] = f2(fr.get("T3_mixing"), 3), f2(fh.get("T3_mixing"), 3)
    N["ivLinks"] = pct(h7.get("IV_frac_ok_links"))
    N["heightMaterial"] = pct(h7.get("H8_frac_material_H4"))
    N["heightKendall"] = f2(h7.get("H8_median_kendall_frozen_total_H4"))
    N["heightFlipCond"] = pct(h7.get("H8_conditions_with_sign_flip_H4"))
    N["nodeRone"] = f"{h7.get('node_R1_max_rel_frozen_minus_total_controller', float('nan')):.1e}" if h7 else "n/a"
    # mixing
    N["mixVarC"] = pct(h13.get("var_share_C"))
    N["mixShareC"] = pct(h13.get("share_abs_C_of_abs_C_plus_abs_F_median"))
    t = h13.get("tests", {})
    N["mixOldOrigRho"] = f2(t.get("old|mu11|n_rev", {}).get("rho"))
    N["mixOldFixedRho"] = f2(t.get("old|mu10|n_rev", {}).get("rho"))
    N["mixNewRho"] = f2(h13.get("T3_primary", {}).get("rho"))
    N["mixNewP"] = f2(h13.get("T3_primary", {}).get("p"), 3)
    N["mixNewN"] = str(h13.get("T3_primary", {}).get("n"))
    # topology
    N["topoRemovalPfour"] = ", ".join(h15.get("em_removal_H4_P4", []))
    N["topoNcreate"] = str(h15.get("n_em_creation"))
    N["topoNfastCreate"] = str(h15.get("n_fast_creation"))
    N["topoNpol"] = str(h15.get("n_policies"))
    lc = h15.get("label_counts_H4", {})
    N["topoFastHfour"] = str(lc.get("FAST-MODE DOMINATED", 0))
    for sc, v in h15.get("H16", {}).items():
        key = re.sub("[^A-Za-z]", "", sc.title())
        N[f"stat{key}Rho"] = f2(v.get("median_rho_em"))
        N[f"stat{key}Frac"] = f2(v.get("frac_abs_ge_0.6_em"))
    # corridors
    for k, v in h10.items():
        if isinstance(v, float) and ("TXother" in k or "TXall" in k):
            N["cor" + re.sub("[^A-Za-z]", "", k.replace("TXother", "Oth").replace("TXall", "All").replace("|", ""))] = f2(v)
    N["corRobust"] = ", ".join(h10.get("H10_robust_weak_corridors", [])) or "none"
    # cross model
    for k, v in h18.items():
        if isinstance(v, (float, int, str)) and not isinstance(v, bool):
            N["alt" + re.sub("[^A-Za-z]", "", k.title())] = f2(v) if isinstance(v, float) else str(v)
    N["layersPos"] = str(h19.get("n_positive"))
    PAPER.mkdir(parents=True, exist_ok=True)
    lines = ["% generated by experiments/cdw_hardening/CDWH_NUMBERS.py - do not edit"]
    for k, v in sorted(N.items()):
        if re.fullmatch("[A-Za-z]+", k):
            lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")
    (PAPER / "cdw_numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    HI.write_json("CDWH_numbers.json", N)
    return N


if __name__ == "__main__":
    n = build()
    print(len(n), "macros")
