# ruff: noqa: E501
"""R23 (after the R19 decision, CASE B): LaTeX macros and the IEEE-68 table for the CDW manuscript, read from results/.
Writes reports/papers/cdw_contextual_dynamic_weakness/{cdw68_numbers.tex, tab_cdw68.tex} and copies the used figures."""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
import shutil

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

PAPER = R.RESEARCH.parents[2] / "papers" / "cdw_contextual_dynamic_weakness"
RES = R.RESULTS


def j(n):
    return json.loads((RES / n).read_text())


def f2(x, d=2):
    return f"{x:.{d}f}"


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def ab_transfer():
    a = pd.read_csv(RES / "CDW68_R10_A_branch_table.csv")
    b = pd.read_csv(RES / "CDW68_R10_B_branch_table.csv")
    a["pid"] = a.cond.str.split("|").str[-1]
    b["pid"] = b.cond.str.split("|").str[-1].str.replace("B68", "P68")
    rho, kt, top = [], [], []
    for pid, ga in a.groupby("pid"):
        m = ga.merge(b[b.pid == pid], on="e", suffixes=("_a", "_b"))
        ta, tb = m["T_1.5_a"], m["T_1.5_b"]
        rho.append(spearmanr(ta, tb).statistic)
        kt.append(kendalltau(ta, tb).statistic)
        ra, rb = set(m.e[ta.rank(ascending=False) <= 5]), set(m.e[tb.rank(ascending=False) <= 5])
        top.append(len(ra & rb) / 5)
    sd = [(np.sign(g.P_tot) != np.sign(g.P_fro)).mean() for _, g in a.groupby("pid")]
    return float(np.median(rho)), float(np.median(kt)), float(np.median(top)), float(np.median(sd)), float(np.nanmax(np.abs(a[["T_1.1", "T_1.25", "T_1.5"]].to_numpy()))), \
        float(np.nanmax(np.abs(b[["T_1.1", "T_1.25", "T_1.5"]].to_numpy()))), float(np.nanmax(np.abs(a["Tem_1.5"])))


def marg(name):
    m = pd.read_parquet(RES / name)
    t = m[m.lvT].t_delta.dropna()
    return float(t.min()), float(t.max()), int((t <= -R.TAU_MAT).sum()), int((t >= R.TAU_MAT).sum()), float(t.median())


def main():
    r7a, r7b, r12, r15, r14 = j("CDW68_R07_A_gate.json"), j("CDW68_R07_B_gate.json"), j("CDW68_R12_gate.json"), j("CDW68_R15_uncertainty.json"), j("CDW68_R14_tds_posthoc.json")
    A, B = r12["A"]["R12"], r12["B"]["R12"]
    rho_ab, kt_ab, top_ab, sd_tf, tmax_a, tmax_b, temmax = ab_transfer()
    ma, mb = marg("CDW68_R07_A_REAL_marginals.parquet"), marg("CDW68_R07_B_REAL_marginals.parquet")
    mng, msp = marg("CDW68_R13_NOGOV_marginals.parquet"), marg("CDW68_R13_SP33_marginals.parquet")
    nA = sum(v["n_draws"] for k, v in r15.items() if k.startswith("A_E"))
    nB = sum(v["n_draws"] for k, v in r15.items() if k.startswith("B_E"))
    kA = sum(round(v["frac_draws_level_D"] * v["n_draws"]) + round(v["frac_draws_level_T"] * v["n_draws"]) for k, v in r15.items() if k.startswith("A_E"))
    kB = sum(round(v["frac_draws_level_D"] * v["n_draws"]) + round(v["frac_draws_level_T"] * v["n_draws"]) for k, v in r15.items() if k.startswith("B_E"))
    rk = r15.get("A_ranking_draws", {})
    mac = {
        "sxRevA": f"{r7a['A_REAL']['k_D']}/{r7a['A_REAL']['n_eligible_base_stable']}",
        "sxRevB": f"{r7b['B_REAL']['k_D']}/{r7b['B_REAL']['n_eligible_base_stable']}",
        "sxRevTA": f"{r7a['A_REAL']['k_T']}/{r7a['A_REAL']['n_eligible_base_stable']}",
        "sxRevTB": f"{r7b['B_REAL']['k_T']}/{r7b['B_REAL']['n_eligible_base_stable']}",
        "sxRevNOGOV": f"{r7a['ablation_NOGOV']['k_T']}/{r7a['ablation_NOGOV']['n_eligible_base_stable']}",
        "sxRevSP": f"{r7a['ablation_SP33']['k_T']}/{r7a['ablation_SP33']['n_eligible_base_stable']}",
        "sxEmMinA": f2(ma[0], 4), "sxEmMaxA": f2(ma[1], 3), "sxEmDestA": str(ma[3]), "sxEmMedA": f2(ma[4], 4),
        "sxEmMinB": f2(mb[0], 4), "sxEmMaxB": f2(mb[1], 3), "sxEmDestB": str(mb[3]), "sxEmMedB": f2(mb[4], 4),
        "sxStabNOGOV": str(mng[2]), "sxStabSP": str(msp[2]), "sxEmMinSP": f2(msp[0], 4),
        "sxRhoA": f2(A["median_rho_tot"], 3), "sxRhoB": f2(B["median_rho_tot"], 3),
        "sxRhoFroA": f2(A["median_rho_P_fro"], 2), "sxRhoConvA": f2(A["median_rho_P_conv"], 2),
        "sxAdvA": f2(A["median_advantage"], 2), "sxAdvB": f2(B["median_advantage"], 2),
        "sxTopA": f2(A["median_top5_tot"], 1), "sxTopB": f2(B["median_top5_tot"], 1),
        "sxStatA": f2(A["median_rho_static_sel"], 2), "sxStatB": f2(B["median_rho_static_sel"], 2),
        "sxTtwoA": f2(r12["F68"]["raw"]["T2A"], 4), "sxTtwoAholm": f2(r12["F68"]["holm"]["T2A"], 3), "sxTtwoB": f2(r12["F68"]["raw"]["T2B"], 3),
        "sxRoneA": sci(r12["A"]["R11"]["median_rel_err"]), "sxRoneRhoA": f2(r12["A"]["R11"]["spearman"], 5),
        "sxTruthMedA": sci(A["truth_abs_max_median"]), "sxTruthMedB": sci(B["truth_abs_max_median"]),
        "sxTruthMaxA": sci(tmax_a), "sxTruthMaxB": sci(tmax_b),
        "sxEmRho": f2(A["posthoc_EM_tracked"]["median_rho"], 3), "sxEmTop": f2(A["posthoc_EM_tracked"]["median_top5"], 1),
        "sxEmTruthMed": f2(A["posthoc_EM_tracked"]["truth_abs_max_median"], 3), "sxEmTruthMax": f2(temmax, 3), "sxEmNmat": f"{A['posthoc_EM_tracked']['n_material_median']:.0f}",
        "sxAbRho": f2(rho_ab, 2), "sxAbKendall": f2(kt_ab, 2), "sxAbTop": f2(top_ab, 1), "sxSignDisTF": f"{100 * sd_tf:.0f}",
        "sxDrawsA": f"{kA}/{nA}", "sxDrawsB": f"{kB}/{nB}", "sxNdrawsA": str(nA), "sxNdrawsB": str(nB),
        "sxDrawRho": f2(rk["rho_tot_quantiles"][1], 3) if rk else "n/a", "sxDrawRhoEm": f2(rk["rho_em_posthoc_quantiles"][1], 3) if rk else "n/a",
        "sxTdsPairs": str(len(r14["pairs"])), "sxTdsAgree": str(sum(p["sign_agree_1"] + p["sign_agree_2"] for p in r14["pairs"])),
        "sxTdsMaxDiff": sci(max(max(abs(p["tds_d1"] - p["lin_d1"]), abs(p["tds_d2"] - p["lin_d2"])) for p in r14["pairs"])),
    }
    out = ["% generated by ieee68_replication/code/R23_paper_numbers.py from results/; do not edit"]
    out += [rf"\newcommand{{\{k}}}{{{v}}}" for k, v in mac.items()]
    (PAPER / "cdw68_numbers.tex").write_text("\n".join(out) + "\n", encoding="utf-8")
    tab = r"""\begin{table}[!t]
\centering
\caption{Replication on the IEEE 68-bus system with PST governors and damping
(16 policies per converter model; branch ranking on 12 holdout conditions at
$\gamma_e=1.5$)}\label{tab:cdw68}
\footnotesize
\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{@{}lcc@{}}
\toprule
 & A: custom GFL & B: WECC GFL\\
\midrule
level-D nested reversal & \sxRevA & \sxRevB\\
level-T (EM-tracked) reversal & \sxRevTA & \sxRevTB\\
EM-tracked marginals (s$^{-1}$) & $[\sxEmMinA,\,\sxEmMaxA]$ & $[\sxEmMinB,\,\sxEmMaxB]$\\
materially stabilizing marginals & 0 & 0\\
fixed unit ranking, regret & 0 & 0\\
\midrule
median $\rho(\Dtot,\text{finite})$ & \sxRhoA & \sxRhoB\\
selected static index, $\rho$ & \sxStatA{} (el.\ dist.) & \sxStatB{} ($\Delta$gSCR)\\
median paired advantage & \sxAdvA & \sxAdvB\\
median top-5 precision & \sxTopA & \sxTopB\\
largest effect on $\aperp$ (s$^{-1}$) & $\sxTruthMaxA$ & $\sxTruthMaxB$\\
\bottomrule
\end{tabular}
\end{table}
"""
    (PAPER / "tab_cdw68.tex").write_text(tab, encoding="utf-8")
    for f in ("CDW68_F6_ranking_methods.pdf", "CDW68_F12_crossbench.pdf", "CDW68_F9_ablation.pdf"):
        shutil.copy2(R.FIGS / f, PAPER / "figures" / f)
    print(json.dumps(mac, indent=1))


if __name__ == "__main__":
    main()
