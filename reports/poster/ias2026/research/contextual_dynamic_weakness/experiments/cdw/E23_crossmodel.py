# ruff: noqa: E501
"""E23 cross-model holdout: TX4 StaticInjection (constant-power converter limit), prereg E23."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _sens as S
import E01_census as E1
import E06_corridors as E6

PHASE_C = "E23c"
PHASE_L = "E23l"
POLICIES = ("D01", "D03", "D11", "D14", "H01", "H02", "H03", "H04", "H05", "H06")


def tasks_c():
    return [{**t, "phase": PHASE_C} for t in E1.tasks(policies=list(POLICIES), device="static")]


def run_task_c(task):
    return E1.run_task(task)


def tasks_l():
    cor = E6.corridors()
    return [{"phase": PHASE_L, "pid": p, "corridors": cor} for p in POLICIES]


def run_task_l(task):
    theta = C.policy(task["pid"])
    eng = S.Engine(C.V4, theta, device="static")
    specs = [{"kind": "line", "idx": e, "value": 1.0} for e in range(46)]
    d = eng.derivatives(specs, "SPR")
    base = C.eval_portfolio(C.V4, theta, modes=False, device="static")["alpha"]
    items = []
    for it in d["items"]:
        r = C.eval_portfolio(C.V4, theta, modes=False, device="static", network=C.network_with(scale={it["idx"]: 1.5}))
        items.append({"idx": it["idx"], "d_total": it["d_total"][0], "d_frozen": it["d_frozen"][0], "alpha_x1.5": r["alpha"]})
    cors = []
    for name, cut in task["corridors"].items():
        r = C.eval_portfolio(C.V4, theta, modes=False, device="static", network=C.network_with(scale={e: 1.5 for e in cut}))
        cors.append({"corridor": name, "alpha_x1.5": r["alpha"]})
    return {"alpha0": base, "items": items, "corridors": cors, "lam": d["lam"]}


def aggregate():
    df, modes = E1.load(PHASE_C, device="static")
    mg = E1.marginals(df, modes)
    pi = E1.per_intervention(mg)
    pol_s = pi.groupby("pid").rev_stable.any()
    base_s = df[df.S == "BASE"].set_index("pid").alpha
    pol_g = pd.read_csv(I.RESULTS / "CDW_E1_policy_summary.csv").set_index("pid")
    rows = []
    for pid in POLICIES:
        rows.append({"pid": pid, "static_rev_stable": bool(pol_s.get(pid, False)), "static_base_alpha": float(base_s.get(pid, np.nan)),
                     "gfl_rev_stable": bool(pol_g.rev_stable.get(pid, False)) if pid in pol_g.index else None})
    rev = pd.DataFrame(rows)
    link = pd.read_parquet(I.RESULTS / "CDW_E4_link_sensitivity.parquet")
    link = link[(link.semantics == "SPR") & (link.target == "H4")]
    e6 = pd.read_csv(I.RESULTS / "CDW_E6_corridors.csv")
    lk, co = [], []
    for r in I.Store(PHASE_L).all():
        if not r.get("ok"):
            continue
        pid = r["task"]["pid"]
        st = pd.DataFrame(r["items"]).set_index("idx")
        st_eff = -(st["alpha_x1.5"] - r["alpha0"])
        g = link[(link.pid == pid) & link.env.isna()].set_index("idx")
        if len(g) == 46:
            g_eff = -(g["large_x1.5_alpha"] - g["alpha0"])
            lk.append({"pid": pid, "kendall_finite": AN.kendall(g_eff.sort_index(), st_eff.sort_index()),
                       "kendall_dtotal": AN.kendall(g.d_total.sort_index(), st.d_total.sort_index()),
                       "static_rho_dtotal_vs_finite": AN.spearman(-st.d_total, st_eff)})
        sc = pd.DataFrame(r["corridors"])
        sc["eff"] = -(sc["alpha_x1.5"] - r["alpha0"])
        gc = e6[(e6.pid == pid) & e6.env.isna()]
        if len(gc):
            t3s = set(sc.sort_values("eff", ascending=False).corridor.head(3))
            t3g = set(gc.sort_values("delta_finite").corridor.head(3))
            co.append({"pid": pid, "top3_overlap": len(t3s & t3g) / 3})
    lk, co = pd.DataFrame(lk), pd.DataFrame(co)
    rev.to_csv(I.RESULTS / "CDW_E23_crossmodel_reversal.csv", index=False)
    lk.to_csv(I.RESULTS / "CDW_E23_crossmodel_links.csv", index=False)
    co.to_csv(I.RESULTS / "CDW_E23_crossmodel_corridors.csv", index=False)
    bs = rev[rev.static_base_alpha < 0]
    gate = {"static_frac_policies_with_reversal": float(bs.static_rev_stable.mean()) if len(bs) else np.nan,
            "n_static_base_stable": int(len(bs)),
            "median_link_kendall_finite": float(lk.kendall_finite.median()) if len(lk) else np.nan,
            "median_link_kendall_dtotal": float(lk.kendall_dtotal.median()) if len(lk) else np.nan,
            "median_corridor_top3_overlap": float(co.top3_overlap.median()) if len(co) else np.nan}
    gate["reversal_transfers"] = bool(gate["static_frac_policies_with_reversal"] >= 0.5)
    gate["links_transfer"] = bool(gate["median_link_kendall_finite"] >= 0.5)
    gate["corridors_transfer"] = bool(gate["median_corridor_top3_overlap"] >= 2 / 3)
    gate["H9_transfers"] = bool(gate["reversal_transfers"] and gate["links_transfer"] and gate["corridors_transfer"])
    I.atomic_write_json(I.RESULTS / "CDW_E23_gates.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
