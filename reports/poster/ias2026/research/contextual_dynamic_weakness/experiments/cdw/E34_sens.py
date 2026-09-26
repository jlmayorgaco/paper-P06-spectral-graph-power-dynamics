# ruff: noqa: E501
"""E3 (node) and E4 (link) total re-equilibrated sensitivity, with baselines (prereg E3/E4)."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _sens as S

PHASE = "E34"
DISC = ("D01", "D03", "D08", "D11", "D14")
HOLD = tuple(f"H{i:02d}" for i in range(1, 13))
ENVS = ("EM-f", "EM-u", "EC", "EMC")
TARGETS = {"H4": C.V4, "V9": C.V9}
RANGE = {"g": 1.0, "pll": 0.4, "ka": 0.3, "vset": 0.04, "load": 0.2, "line": 1.0}
CHUNK = 12


def load_buses():
    return sorted(C.load_network().loads)


def node_specs(T, theta):
    g = theta[0]
    sp = [{"kind": "g", "idx": i, "value": g} for i in T]
    sp += [{"kind": "pll", "idx": i, "value": 1.0} for i in T]
    sp += [{"kind": "ka", "idx": i, "value": 1.0} for i in C.SG_BUSES if i not in T]
    sp += [{"kind": "vset", "idx": i, "value": 0.0} for i in C.SG_BUSES]
    sp += [{"kind": "load", "idx": j, "value": 1.0} for j in load_buses()]
    return sp


def link_specs():
    return [{"kind": "line", "idx": e, "value": 1.0} for e in range(46)]


def conditions():
    out = []
    for pid in DISC:
        for tg in TARGETS:
            out.append({"pid": pid, "target": tg, "split": "discovery", "env": None, "draw": None, "sems": ["SPR", "RP"]})
    for pid in HOLD:
        for tg in TARGETS:
            out.append({"pid": pid, "target": tg, "split": "holdout", "env": None, "draw": None, "sems": ["SPR", "RP"]})
    for env in ENVS:
        for k in range(5):
            out.append({"pid": "D01", "target": "H4", "split": "holdout_env", "env": env, "draw": k, "sems": ["SPR"]})
    return out


def tasks():
    out = []
    for c in conditions():
        theta = C.policy(c["pid"])
        for fam, specs in (("node", node_specs(TARGETS[c["target"]], theta)), ("link", link_specs())):
            for sem in c["sems"]:
                for k in range(0, len(specs), CHUNK):
                    out.append({"phase": PHASE, **{kk: c[kk] for kk in ("pid", "target", "split", "env", "draw")},
                                "family": fam, "semantics": sem, "specs": specs[k: k + CHUNK]})
    # conventional (base-portfolio) eigenvalue sensitivity, SPR derivative only
    seen = set()
    for c in conditions():
        key = (c["pid"], c["env"], c["draw"])
        if key in seen:
            continue
        seen.add(key)
        out.append({"phase": PHASE, "pid": c["pid"], "target": "BASE", "split": c["split"], "env": c["env"], "draw": c["draw"],
                    "family": "conv", "semantics": "SPR", "specs": link_specs()})
    out.append({"phase": PHASE, "pid": "D01", "target": "H4", "split": "discovery", "env": None, "draw": None,
                "family": "port", "semantics": "SPR", "specs": link_specs()})
    return out


def _draw(t):
    if t.get("env") is None:
        return None
    for d in C.cdw_draws():
        if d["envelope"] == t["env"] and d["draw"] == t["draw"]:
            return d
    raise KeyError(t)


def steps(sp, semantics):
    k, v = sp["kind"], sp["value"]
    if k == "g":
        small = (v + 1e-3, v - 1e-3)
        large = {"plus": min(v + 0.1, 1.0), "minus": max(v - 0.1, 0.0)}
    elif k == "vset":
        small = (v + 1e-4, v - 1e-4)
        large = {"plus": v + 0.01, "minus": v - 0.01}
    elif k == "line":
        small = (v * (1 + 1e-3), v * (1 - 1e-3))
        large = {"x1.5": 1.5, "x2": 2.0} if semantics == "SPR" else {"x1.5": 1.5}
    else:
        small = (v * (1 + 1e-3), v * (1 - 1e-3))
        large = {"plus": v * 1.1, "minus": v * 0.9}
    return small, large


def run_task(task):
    theta = C.policy(task["pid"])
    draw = _draw(task)
    members = () if task["target"] == "BASE" else TARGETS[task["target"]]
    eng = S.Engine(members, theta, draw=draw)
    sem = task["semantics"]
    if task["family"] == "port":
        return {"items": port_items(eng, task["specs"])}
    d = eng.derivatives(task["specs"], sem)
    lam = complex(*d["lam"])
    items = []
    for it in d["items"]:
        rec = dict(it)
        if task["family"] != "conv":
            small, large = steps(it, sem)
            fp, fm = eng.finite(it, small[0], sem, lam), eng.finite(it, small[1], sem, lam)
            rec["fd_small"] = (fp["tracked_re"] - fm["tracked_re"]) / (small[0] - small[1]) if (fp and fm) else None
            A00, _ = eng.A_of(eng.w0(sem), None, sem)
            rec["alpha0"] = float(eng.transverse(A00).real.max())
            for name, val in large.items():
                f = eng.finite(it, val, sem, lam)
                rec[f"large_{name}_value"] = val
                rec[f"large_{name}_alpha"] = f["alpha"] if f else None
                rec[f"large_{name}_tracked"] = f["tracked_re"] if f else None
        items.append(rec)
    return {"lam": d["lam"], "gap2": d["gap2"], "R0": d["R0"], "match_err": d.get("match_err"), "items": items,
            "alpha0_transverse": float(C.transverse_eval(eng.case, modes=False)["alpha"])}


def port_items(eng, specs):
    """T-form (port) total derivative via re-solved Jacobians at fixed s* (independent formula)."""

    a0, _, jac0 = C.physical_matrices(eng.case)
    lam, *_ = eng.critical_T(a0, eng.case.dae, eng.case.equilibrium.z)
    out = []

    def tmat(jac, s):
        n = jac.fx.shape[0]
        return jac.gz + jac.gx @ np.linalg.solve(s * np.eye(n) - jac.fx, jac.fz)

    T0 = tmat(jac0, lam)
    U, sv, Vh = np.linalg.svd(T0)
    p, q = U[:, -1], Vh[-1].conj()
    n = jac0.fx.shape[0]
    Rinv = np.linalg.inv(lam * np.eye(n) - jac0.fx)
    dTs = -jac0.gx @ Rinv @ Rinv @ jac0.fz
    den = np.vdot(p, dTs @ q)
    for sp in specs:
        h = 1e-3
        cp = S.spr_case(eng, sp, 1.0 + h)
        cm = S.spr_case(eng, sp, 1.0 - h)
        if cp is None or cm is None:
            out.append({**sp, "d_port": None})
            continue
        jp = C.physical_matrices(cp)[2]
        jm = C.physical_matrices(cm)[2]
        dTa = (tmat(jp, lam) - tmat(jm, lam)) / (2 * h)
        ds = -np.vdot(p, dTa @ q) / den
        out.append({**sp, "d_port": [float(ds.real), float(ds.imag)], "sigma_min": float(sv[-1])})
    return out


# ------------------------------------------------------------------ baselines --
def static_link_baselines(T=None):
    import networkx as nx

    from ibr_cycles.diagnosis.baselines import generalized_scr, short_circuit_ybus

    net = C.load_network()
    flow = C.pf(net)
    V = flow.voltages
    order = list(net.bus_idx)
    idx = {b: i for i, b in enumerate(order)}
    br = C.branches()
    Z = np.linalg.inv(short_circuit_ybus(net))
    L, _ = C.laplacian_B()
    Lp = np.linalg.pinv(L)
    w, U = np.linalg.eigh(L)
    u2 = U[:, 1]
    G = nx.Graph()
    for b in br:
        G.add_edge(b["f"], b["t"], weight=abs(b["x"]), e=b["e"])
    eb = nx.edge_betweenness_centrality(G, weight="weight")
    dvdq = qv_sensitivity(net, flow)
    rows = []
    gscr0 = generalized_scr(net, tuple(T), flow) if T else np.nan
    for b in br:
        f, t = idx[b["f"]], idx[b["t"]]
        y = 1 / complex(b["r"], b["x"])
        m = b["tap"]
        yff = (y + 1j * b["b"] / 2) / m**2
        yft = -y / m
        i_f = yff * V[f] + yft * V[t]
        s_f = V[f] * np.conj(i_f)
        e = np.zeros(len(order))
        e[f], e[t] = 1, -1
        rec = {"e": b["e"], "f": b["f"], "t": b["t"], "trans": b["trans"],
               "S1_absP": abs(s_f.real), "S2_absS": abs(s_f), "S3_absz": abs(complex(b["r"], b["x"])),
               "S4_elecdist": abs(Z[f, f] + Z[t, t] - Z[f, t] - Z[t, f]),
               "S5_reff": float(e @ Lp @ e), "S6_fiedler": float((u2[f] - u2[t]) ** 2),
               "S7_betweenness": eb.get((b["f"], b["t"]), eb.get((b["t"], b["f"]), np.nan)),
               "S8_dvdq": max(dvdq.get(b["f"], 0.0), dvdq.get(b["t"], 0.0))}
        if T:
            net2 = C.network_with(scale={b["e"]: 1.5})
            fl2 = C.pf(net2)
            rec["S9_dgscr"] = generalized_scr(net2, tuple(T), fl2) - gscr0
        rows.append(rec)
    return pd.DataFrame(rows)


def qv_sensitivity(net, flow):
    """diag of the reduced QV sensitivity (dV/dQ) at PQ buses (0 at voltage-controlled buses)."""

    order = list(net.bus_idx)
    n = len(order)
    V = flow.voltages
    slack = order.index(net.slack_bus)
    pvset = {order.index(b) for b in net.pv}
    ang = [i for i in range(n) if i != slack]
    pq = [i for i in range(n) if i != slack and i not in pvset]
    Y = net.ybus

    def mism(th, vm):
        v = vm * np.exp(1j * th)
        s = v * np.conj(Y @ v)
        return s.real, s.imag

    th0, vm0 = np.angle(V), np.abs(V)
    h = 1e-7
    J = np.zeros((len(ang) + len(pq), len(ang) + len(pq)))
    cols = [("t", i) for i in ang] + [("v", i) for i in pq]
    for c, (kind, i) in enumerate(cols):
        thp, vmp, thm, vmm = th0.copy(), vm0.copy(), th0.copy(), vm0.copy()
        if kind == "t":
            thp[i] += h
            thm[i] -= h
        else:
            vmp[i] += h
            vmm[i] -= h
        pp, qp = mism(thp, vmp)
        pm_, qm = mism(thm, vmm)
        J[:, c] = np.concatenate([(pp - pm_)[ang], (qp - qm)[pq]]) / (2 * h)
    na = len(ang)
    Jpt, Jpv, Jqt, Jqv = J[:na, :na], J[:na, na:], J[na:, :na], J[na:, na:]
    JR = Jqv - Jqt @ np.linalg.solve(Jpt, Jpv)
    sens = np.linalg.inv(JR)
    return {order[i]: float(sens[k, k]) for k, i in enumerate(pq)}


# ------------------------------------------------------------------ aggregate --
def aggregate():
    recs = I.Store(PHASE).all()
    bad = [r for r in recs if not r.get("ok")]
    rows = []
    for r in recs:
        if not r.get("ok"):
            continue
        t = r["task"]
        for it in r["items"]:
            row = {k: t[k] for k in ("pid", "target", "split", "env", "draw", "family", "semantics")}
            row.update({k: v for k, v in it.items() if k not in ("d_frozen", "d_total", "d_port")})
            for k in ("d_frozen", "d_total", "d_port"):
                if it.get(k) is not None:
                    row[k] = it[k][0]
                    row[k + "_im"] = it[k][1]
            row["gap2"] = r.get("gap2")
            row["alpha0_T"] = r.get("alpha0_transverse")
            rows.append(row)
    df = pd.DataFrame(rows)
    df["pkey"] = df.pid + "|" + df.env.fillna("-").astype(str) + "|" + df.draw.fillna(-1).astype(int).astype(str)
    df["cond"] = df.pkey + "|" + df.target
    node = df[df.family == "node"].copy()
    link = df[df.family == "link"].copy()
    node.to_parquet(I.RESULTS / "CDW_E3_node_sensitivity.parquet", index=False)
    link.to_parquet(I.RESULTS / "CDW_E4_link_sensitivity.parquet", index=False)
    df[df.family.isin(["conv", "port"])].to_parquet(I.RESULTS / "CDW_E4_conv_port.parquet", index=False)
    gates = {"n_failed_tasks": len(bad)}
    gates.update(iv_gate(pd.concat([node, link])))
    gates.update(structural_checks(node))
    gates.update(h2(node, "node"))
    gates.update(h2(link, "line"))
    gates.update(gold_b(link, df))
    gates.update(port_check(df))
    I.atomic_write_json(I.RESULTS / "CDW_E34_gates.json", gates)
    return gates


def iv_gate(d):
    d = d[d.fd_small.notna() & d.d_total.notna()].copy()
    rng = d.kind.map(RANGE)
    err = (d.fd_small - d.d_total).abs()
    ok = (err <= 0.05 * d.d_total.abs()) | (err * rng <= 1e-4)
    ex = d[(d.gap2 < 1e-3)]
    return {"IV_n": int(len(d)), "IV_frac_ok": float(ok.mean()), "IV_pass": bool(ok.mean() >= 0.95),
            "IV_n_near_tie": int(len(ex)), "IV_frac_ok_by_sem": d.assign(ok=ok).groupby("semantics").ok.mean().round(4).to_dict(),
            "IV_frac_ok_by_kind": d.assign(ok=ok).groupby("kind").ok.mean().round(4).to_dict()}


def structural_checks(node):
    spr = node[node.semantics == "SPR"]
    dyn = spr[spr.kind.isin(["g", "pll", "ka"])]
    rel = ((dyn.d_frozen - dyn.d_total).abs() / dyn.d_total.abs().clip(lower=1e-12))
    vs = node[node.kind == "vset"]
    return {"R1_n": int(len(dyn)), "R1_max_rel_frozen_minus_total": float(rel.max()) if len(rel) else np.nan,
            "R1_pass": bool((rel <= 1e-6).all()) if len(rel) else False,
            "R2_max_abs_frozen_vset": float(vs.d_frozen.abs().max()) if len(vs) else np.nan,
            "R2_pass": bool((vs.d_frozen.abs() <= 1e-8).all()) if len(vs) else False}


def h2(d, name):
    out = {}
    for sem, g0 in d.groupby("semantics"):
        per = []
        for cond, g in g0.groupby("cond"):
            g = g[(g.gap2 >= 1e-3)]
            rng = g.kind.map(RANGE)
            sel = g[(g.d_total.abs() * rng) >= C.TAU_MAT]
            if len(sel) < 3:
                continue
            dis = float((np.sign(sel.d_frozen) != np.sign(sel.d_total)).mean())
            kt = AN.kendall(g.d_frozen, g.d_total)
            split = g.split.iloc[0]
            per.append({"cond": cond, "split": split, "sign_disagree": dis, "kendall": kt, "material": (dis >= 0.10) or (kt <= 0.6)})
        p = pd.DataFrame(per)
        if len(p):
            p.to_csv(I.RESULTS / f"CDW_E34_H2_{name}_{sem}.csv", index=False)
            h = p[p.split.str.startswith("holdout")]
            frac = float(h.material.mean()) if len(h) else np.nan
            out[f"H2_{name}_{sem}_frac_material_holdout"] = frac
            out[f"H2_{name}_{sem}_true"] = bool(frac >= 0.25)
            out[f"H2_{name}_{sem}_median_sign_disagree"] = float(h.sign_disagree.median()) if len(h) else np.nan
            out[f"H2_{name}_{sem}_median_kendall"] = float(h.kendall.median()) if len(h) else np.nan
    return out


def link_metrics(link, conv):
    """Per condition: rank metrics of every predictor against finite x1.5 (SPR)."""

    base_cache = {}
    rows = []
    for cond, g in link[link.semantics == "SPR"].groupby("cond"):
        g = g.sort_values("idx")
        T = TARGETS.get(g.target.iloc[0])
        key = g.target.iloc[0]
        if key not in base_cache:
            base_cache[key] = static_link_baselines(T)
        st = base_cache[key].set_index("e").loc[g.idx.to_numpy()]
        truth = -(g["large_x1.5_alpha"].to_numpy(float) - g.alpha0.to_numpy(float))
        truth2 = -(g["large_x2_alpha"].to_numpy(float) - g.alpha0.to_numpy(float))
        preds = {c: st[c].to_numpy(float) for c in st.columns if c.startswith("S")}
        preds["Dfrozen"] = -g.d_frozen.to_numpy(float)
        preds["Dtotal"] = -g.d_total.to_numpy(float)
        rp = link[(link.semantics == "RP") & (link.cond == cond)].sort_values("idx")
        if len(rp) == len(g):
            preds["DtotalRP"] = -rp.d_total.to_numpy(float)
        cv = conv[conv.pkey == g.pkey.iloc[0]].sort_values("idx")
        if len(cv) == len(g):
            preds["Dconv"] = -cv.d_total.to_numpy(float)
        for name, p in preds.items():
            mat = np.abs(truth) >= C.TAU_MAT
            sa = float((np.sign(p[mat]) == np.sign(truth[mat])).mean()) if mat.any() else np.nan
            rows.append({"cond": cond, "split": g.split.iloc[0], "target": key, "pid": g.pid.iloc[0], "predictor": name,
                         "spearman": AN.spearman(p, truth), "kendall": AN.kendall(p, truth), "spearman_x2": AN.spearman(p, truth2),
                         "top5": AN.topk_precision(p, truth, 5), "sign_acc": sa, "n_material": int(mat.sum()),
                         "gap2": float(g.gap2.iloc[0]), "alpha0": float(g.alpha0.iloc[0])})
    return pd.DataFrame(rows)


def gold_b(link, df):
    conv = df[df.family == "conv"]
    m = link_metrics(link, conv)
    m.to_csv(I.RESULTS / "CDW_E4_baselines.csv", index=False)
    statics = [p for p in m.predictor.unique() if p.startswith("S")]
    disc = m[(m.split == "discovery") & (m.target == "H4")]
    sel = disc[disc.predictor.isin(statics)].groupby("predictor").spearman.median().idxmax()
    h = m[m.split.str.startswith("holdout") & (m.target == "H4")]
    piv = h.pivot_table(index="cond", columns="predictor", values="spearman")
    pivt = h.pivot_table(index="cond", columns="predictor", values="top5")
    diff = (piv["Dtotal"] - piv[sel]).median()
    oracle_static = piv[statics].max(axis=1)
    out = {"GB_static_selected_on_discovery": sel, "GB_n_holdout_conditions": int(len(piv)),
           "GB_median_rho_Dtotal": float(piv["Dtotal"].median()), "GB_median_rho_static_selected": float(piv[sel].median()),
           "GB_median_diff": float(diff), "GB_median_top5_Dtotal": float(pivt["Dtotal"].median()),
           "GB_median_rho_oracle_static": float(oracle_static.median()),
           "GB_median_diff_vs_oracle_static": float((piv["Dtotal"] - oracle_static).median()),
           "GB_median_rho_Dfrozen": float(piv["Dfrozen"].median()),
           "GB_median_rho_Dconv": float(piv["Dconv"].median()) if "Dconv" in piv else np.nan}
    out["GB_pass"] = bool(out["GB_median_diff"] >= 0.20 and out["GB_median_rho_Dtotal"] >= 0.70 and out["GB_median_top5_Dtotal"] >= 0.6)
    hv = m[m.split.str.startswith("holdout") & (m.target == "V9")].pivot_table(index="cond", columns="predictor", values="spearman")
    if len(hv):
        out["GB_V9_median_rho_Dtotal"] = float(hv["Dtotal"].median())
        out["GB_V9_median_rho_static_selected"] = float(hv[sel].median()) if sel in hv else np.nan
    return out


def port_check(df):
    p = df[df.family == "port"].sort_values("idx")
    a = df[(df.family == "link") & (df.semantics == "SPR") & (df.pid == "D01") & (df.target == "H4") & (df.split == "discovery")].sort_values("idx")
    if len(p) != 46 or len(a) != 46:
        return {"port_check": "incomplete"}
    rel = (p.d_port.to_numpy(float) - a.d_total.to_numpy(float))
    return {"port_vs_total_max_abs": float(np.abs(rel).max()), "port_vs_total_spearman": AN.spearman(p.d_port, a.d_total)}


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
