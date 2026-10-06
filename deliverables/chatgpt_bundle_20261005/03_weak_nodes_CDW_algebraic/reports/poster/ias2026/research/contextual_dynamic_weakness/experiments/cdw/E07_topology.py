# ruff: noqa: E501
"""E7 topology reconfiguration: single-branch outages (operational) and doublings (planning)."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C

PHASE = "E07"
POLICIES = ("D01", "D03", "D11", "D14", "H01", "H02", "H03", "H04", "H05", "H06")
VGUARD = (0.90, 1.10)


def actions():
    out = [{"action": f"out{e}", "outage": [e], "scale": {}} for e in range(46)]
    out += [{"action": f"dbl{e}", "outage": [], "scale": {str(e): 2.0}} for e in range(46)]
    return out


def network_of(act):
    return C.network_with(scale={int(k): v for k, v in act["scale"].items()} or None, outage=set(act["outage"]) or None)


def admissibility():
    rows = []
    base = C.pf(C.load_network())
    for act in actions():
        rec = {"action": act["action"]}
        if act["outage"] and not C.connected(set(act["outage"])):
            rec.update(admissible=False, reason="islanding")
            rows.append(rec)
            continue
        net = network_of(act)
        fl = C.pf(net)
        vm = np.abs(fl.voltages)
        rec.update(pf_converged=bool(fl.converged), vmin=float(vm.min()), vmax=float(vm.max()),
                   in_095_105=bool((vm >= 0.95).all() and (vm <= 1.05).all()),
                   max_dv=float(np.abs(vm - np.abs(base.voltages)).max()))
        rec["admissible"] = bool(fl.converged and vm.min() >= VGUARD[0] and vm.max() <= VGUARD[1])
        rec["reason"] = "" if rec["admissible"] else ("pf" if not fl.converged else "voltage")
        # generator Q against limits (reported, not enforced)
        viol = 0
        for b, spec in net.pv.items():
            q = (fl.injection(b, net.ybus) + net.loads.get(b, 0j)).imag
            viol += (q > spec["qmax"] + 1e-6) or (q < spec["qmin"] - 1e-6)
        rec["n_q_limit_violations"] = int(viol)
        rows.append(rec)
    return pd.DataFrame(rows)


def tasks():
    adm = admissibility()
    adm.to_csv(I.RESULTS / "CDW_E7_admissibility.csv", index=False)
    ok = set(adm[adm.admissible].action)
    out = []
    acts = [a for a in actions() if a["action"] in ok]
    for pid in POLICIES:
        out.append({"phase": PHASE, "kind": "core", "pid": pid, "action": "NOMINAL", "act": {"outage": [], "scale": {}}})
        for a in acts:
            out.append({"phase": PHASE, "kind": "core", "pid": pid, "action": a["action"], "act": a})
    subs = C.subsets(C.V9)
    for a in [{"action": "NOMINAL", "outage": [], "scale": {}}] + acts:
        for c in range(0, 512, 64):
            out.append({"phase": PHASE, "kind": "census", "pid": "D01", "action": a["action"], "act": a,
                        "subsets": [list(s) for s in subs[c: c + 64]]})
    return out


def run_task(task):
    theta = C.policy(task["pid"])
    net = network_of(task["act"]) if (task["act"]["outage"] or task["act"]["scale"]) else None
    subs = C.subsets(C.V4) if task["kind"] == "core" else [tuple(s) for s in task["subsets"]]
    recs = []
    for s in subs:
        r = C.eval_portfolio(s, theta, modes=False, network=net)
        recs.append({"S": C.label(s), **{k: r.get(k) for k in ("status", "alpha", "lam_hz")}})
    return {"records": recs}


def hyper(status, cands):
    st = {s: status.get(C.label(s)) for s in C.subsets(cands)}
    if st[()] != "STABLE":
        return "BASE_" + str(st[()]), np.nan, []
    unsafe = [s for s, v in st.items() if v == "UNSTABLE"]
    minimal = sorted((s for s in unsafe if not any(set(r) < set(s) for r in unsafe)), key=lambda s: (len(s), s))
    return "|".join(C.label(e) for e in minimal) or "EMPTY", (min(len(e) for e in minimal) if minimal else np.inf), minimal


def static_scores(act):
    from ibr_cycles.diagnosis.baselines import generalized_scr, nodal_metrics

    out = set(act["outage"])
    scale = {int(k): v for k, v in act["scale"].items()}
    L, _ = C.laplacian_B(scale or None, out or None)
    w = np.linalg.eigvalsh(L)
    net = network_of(act) if (out or scale) else C.load_network()
    fl = C.pf(net)
    nm = nodal_metrics(net, fl)
    return {"fiedler": float(w[1]), "kirchhoff": float(len(w) * np.sum(1 / w[1:])),
            "gscr_H4": float(generalized_scr(net, C.V4, fl)), "min_scr_core": float(min(nm[b].scr for b in C.V4))}


def aggregate():
    recs = [r for r in I.Store(PHASE).all() if r.get("ok")]
    adm = pd.read_csv(I.RESULTS / "CDW_E7_admissibility.csv")
    core_rows = []
    for r in recs:
        t = r["task"]
        if t["kind"] != "core":
            continue
        status = {x["S"]: x["status"] for x in r["records"]}
        alpha = {x["S"]: x["alpha"] for x in r["records"]}
        H, kap, _ = hyper(status, C.V4)
        core_rows.append({"pid": t["pid"], "action": t["action"], "H": H, "kappa": kap, "alpha_H4": alpha.get("30+33+35+37"),
                          "status_H4": status.get("30+33+35+37")})
    core = pd.DataFrame(core_rows)
    nom = core[core.action == "NOMINAL"].set_index("pid")
    core["d_alpha_H4"] = core.alpha_H4 - core.pid.map(nom.alpha_H4)
    core["H_changed"] = core.H != core.pid.map(nom.H)
    core["removes_H4_P4"] = (core.pid == "D01") & (core.status_H4 == "STABLE")
    nomk = core.pid.map(nom.kappa)
    core["creates_small_edge"] = (core.kappa <= 3) & ~(nomk <= 3)
    stat = {a: static_scores(next(x for x in actions() if x["action"] == a) if a != "NOMINAL" else {"outage": [], "scale": {}})
            for a in core.action.unique()}
    for k in ("fiedler", "kirchhoff", "gscr_H4", "min_scr_core"):
        core[f"d_{k}"] = core.action.map(lambda a, k=k: stat[a][k] - stat["NOMINAL"][k])
    core.to_parquet(I.RESULTS / "CDW_E7_topology.parquet", index=False)
    # census at P4
    by = {}
    for r in recs:
        if r["task"]["kind"] == "census":
            by.setdefault(r["task"]["action"], []).extend(r["records"])
    cen = []
    for a, lst in by.items():
        status = {x["S"]: x["status"] for x in lst}
        if len(status) < 512:
            continue
        H, kap, minimal = hyper(status, C.V9)
        cen.append({"action": a, "H_V9": H, "kappa_V9": kap, "n_edges": len(minimal),
                    "n_unstable": int(sum(v == "UNSTABLE" for v in status.values()))})
    cen = pd.DataFrame(cen)
    cen.to_csv(I.RESULTS / "CDW_E7_census_P4.csv", index=False)
    pred = []
    for pid, g in core[core.action != "NOMINAL"].groupby("pid"):
        for k in ("fiedler", "kirchhoff", "gscr_H4", "min_scr_core"):
            pred.append({"pid": pid, "score": k, "spearman": AN.spearman(g[f"d_{k}"], g.d_alpha_H4)})
    pred = pd.DataFrame(pred)
    pred.to_csv(I.RESULTS / "CDW_E7_static_prediction.csv", index=False)
    gate = {"n_admissible_outages": int(adm[adm.action.str.startswith("out")].admissible.sum()),
            "n_admissible_doublings": int(adm[adm.action.str.startswith("dbl")].admissible.sum()),
            "removes_H4_at_P4": core[core.removes_H4_P4].action.tolist(),
            "creates_small_edge": core[core.creates_small_edge][["pid", "action", "H"]].to_dict("records"),
            "n_H_changes": int(core.H_changed.sum()),
            "static_predicts": {k: bool((pred[pred.score == k].spearman.abs() >= 0.6).mean() >= 0.75) for k in pred.score.unique()},
            "static_median_rho": pred.groupby("score").spearman.median().round(4).to_dict()}
    if len(cen):
        n0 = cen[cen.action == "NOMINAL"]
        gate["census_nominal_H_V9"] = n0.H_V9.iloc[0] if len(n0) else None
        gate["census_n_actions_changing_H_V9"] = int((cen.H_V9 != (n0.H_V9.iloc[0] if len(n0) else None)).sum())
    I.atomic_write_json(I.RESULTS / "CDW_E7_gates.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
