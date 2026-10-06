# ruff: noqa: E501
"""E5 static vs dynamic weakness table (descriptive; prereg E5). No new solves."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C


def aggregate():
    from ibr_cycles.diagnosis.baselines import nodal_metrics

    import E34_sens as E

    net = C.load_network()
    nm = nodal_metrics(net, C.pf(net))
    graph = pd.read_csv(I.RESULTS / "CDW_E11_node_graph_scores.csv").set_index("bus")
    e1 = pd.read_csv(I.RESULTS / "CDW_E1_summary.csv")
    e1h = e1[e1.pid.str.startswith("H")]
    pol = pd.read_csv(I.RESULTS / "CDW_E1_policy_summary.csv")
    stable_h = set(pol[(pol.split == "holdout") & pol.base_stable].pid)
    e1h = e1h[e1h.pid.isin(stable_h)]
    e10 = pd.read_csv(I.RESULTS / "CDW_E10_shapley.csv")
    e10h = e10[e10.pid.isin(stable_h)]
    e3 = pd.read_parquet(I.RESULTS / "CDW_E3_node_sensitivity.parquet")
    e3h = e3[(e3.split == "holdout") & (e3.target == "V9") & (e3.semantics == "SPR")]
    rows = []
    for b in C.V9:
        ng = e3h[(e3h.kind == "g") & (e3h.idx == b)].d_total
        npl = e3h[(e3h.kind == "pll") & (e3h.idx == b)].d_total
        g1 = e1h[e1h.i == b]
        rows.append({"bus": b, "SCR": nm[b].scr, "thevenin_absZ": nm[b].thevenin_magnitude,
                     "fiedler_abs": graph.loc[b, "abs_fiedler"], "reff_to_39": graph.loc[b, "reff_to_39"],
                     "NG_total_median": float(ng.median()) if len(ng) else np.nan, "NP_total_median": float(npl.median()) if len(npl) else np.nan,
                     "NG_norm_abs_median": float((ng.abs() * E.RANGE["g"]).median()) if len(ng) else np.nan,
                     "frac_destab_median": float(g1.frac_destab.median()), "frac_neutral_median": float(g1.frac_neutral.median()),
                     "reversal_policy_frac": float(g1.rev_stable.mean()), "entropy_median": float(g1.entropy.median()),
                     "shapley_median": float(e10h[e10h.i == b].shapley.median()) if len(e10h) else np.nan})
    nd = pd.DataFrame(rows).set_index("bus")
    static_weak = set(nd.SCR.nsmallest(3).index)
    dyn_crit = set(nd.NG_norm_abs_median.nlargest(3).index)
    cat = []
    for b in nd.index:
        c = []
        if b in static_weak and b in dyn_crit:
            c.append("STATIC+DYNAMIC")
        elif b in static_weak:
            c.append("STATIC ONLY")
        elif b in dyn_crit:
            c.append("DYNAMIC ONLY")
        if nd.loc[b, "reversal_policy_frac"] >= 0.5:
            c.append("CONTEXT-DEPENDENT")
        if nd.loc[b, "frac_neutral_median"] >= 0.8:
            c.append("ROBUSTLY NEUTRAL")
        cat.append(";".join(c) or "-")
    nd["category"] = cat
    nd.to_csv(I.RESULTS / "CDW_E5_static_vs_dynamic_nodes.csv")
    # branches
    bl = E.static_link_baselines(C.V4).set_index("e")
    lk = pd.read_parquet(I.RESULTS / "CDW_E4_link_sensitivity.parquet")
    h = lk[lk.split.str.startswith("holdout") & (lk.target == "H4") & (lk.semantics == "SPR")]
    agg = h.groupby("idx").agg(Dtotal_median=("d_total", "median"), Dfrozen_median=("d_frozen", "median"),
                               finite_x15_median=("large_x1.5_alpha", "median"), alpha0_median=("alpha0", "median"))
    agg["finite_effect_median"] = agg.finite_x15_median - agg.alpha0_median
    rank_std = h.assign(eff=-(h["large_x1.5_alpha"] - h.alpha0)).groupby("cond").eff.rank(ascending=False)
    h = h.assign(rank=rank_std)
    agg["rank_iqr"] = h.groupby("idx")["rank"].quantile(0.75) - h.groupby("idx")["rank"].quantile(0.25)
    br = bl.join(agg)
    sw = set(br.S4_elecdist.nlargest(3).index)
    dc = set((-br.Dtotal_median).nlargest(3).index)
    br["category"] = ["STATIC+DYNAMIC" if e in sw and e in dc else "STATIC ONLY" if e in sw else "DYNAMIC ONLY" if e in dc else "-" for e in br.index]
    br.to_csv(I.RESULTS / "CDW_E5_static_vs_dynamic_branches.csv")
    out = {"node_kendall_SCR_vs_NG": AN.kendall(-nd.SCR, nd.NG_norm_abs_median),
           "node_kendall_SCR_vs_fracdestab": AN.kendall(-nd.SCR, nd.frac_destab_median),
           "node_categories": nd.category.to_dict(),
           "branch_kendall_S4_vs_Dtotal": AN.kendall(br.S4_elecdist, -br.Dtotal_median),
           "branch_kendall_S1_vs_Dtotal": AN.kendall(br.S1_absP, -br.Dtotal_median),
           "branch_categories": {int(k): v for k, v in br.category.items() if v != "-"}}
    I.atomic_write_json(I.RESULTS / "CDW_E5_gates.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
