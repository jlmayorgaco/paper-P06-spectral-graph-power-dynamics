# ruff: noqa: E501
"""H14 EM tracking / classification, H15 H_EM vs H under topology, H16 static scores vs EM truth.

Prereg H14-H16. Uses raw/H_H14 (core lattice with modes at 16 policies) and, for the P4 V9
census re-analysis, raw/E07 (stored rightmost frequency)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _infra as I
import E01_census as E1
import E07_topology as E7

TAU = C.TAU_MAT


def em_modes(ms):
    return [m for m in ms if HI.in_band(m["hz"])]


def em_rightmost(ms):
    e = em_modes(ms)
    return max(e, key=lambda m: m["re"]) if e else None


def hyper(status: dict, hz: dict, em=False):
    st = {s: status[C.label(s)] for s in C.subsets(C.V4)}
    if st[()] != "STABLE":
        return "BASE_" + str(st[()]), np.nan
    if em:
        bad = [s for s, v in st.items() if v == "UNSTABLE" and HI.in_band(hz[C.label(s)])]
    else:
        bad = [s for s, v in st.items() if v == "UNSTABLE"]
    minimal = sorted((s for s in bad if not any(set(r) < set(s) for r in bad)), key=lambda s: (len(s), s))
    return "|".join(C.label(e) for e in minimal) or "EMPTY", (min(len(e) for e in minimal) if minimal else np.inf)


def load():
    data = {}
    for r in I.Store("H_H14").all():
        if not r.get("ok"):
            continue
        t = r["task"]
        data[(t["pid"], t["action"])] = {x["S"]: {**{k: x.get(k) for k in ("status", "alpha", "lam_hz", "lam_re")},
                                                  "modes": [E1.compact(m) for m in x.get("modes", [])]} for x in r["records"]}
    return data


def classify(data):
    rows = []
    for (pid, act), recs in data.items():
        if act == "NOMINAL":
            continue
        nom = data.get((pid, "NOMINAL"))
        if nom is None:
            continue
        for S, post in recs.items():
            pre = nom[S]
            rec = {"pid": pid, "action": act, "S": S, "status_pre": pre["status"], "status_post": post["status"],
                   "alpha_pre": pre["alpha"], "alpha_post": post["alpha"], "hz_pre": pre["lam_hz"], "hz_post": post["lam_hz"]}
            rec["d_full"] = (post["alpha"] - pre["alpha"]) if (pre["alpha"] is not None and post["alpha"] is not None) else np.nan
            m0 = em_rightmost(pre["modes"])
            unresolved = (pre["status"] in ("INFEASIBLE", "BOUNDARY_OR_UNRESOLVED") or post["status"] in ("INFEASIBLE", "BOUNDARY_OR_UNRESOLVED")
                          or m0 is None)
            m1, mac = (None, np.nan)
            if not unresolved:
                m1, mac = E1.match(m0, em_modes(post["modes"]))
            rec["mac"] = mac
            rec["d_em"] = (m1["re"] - m0["re"]) if m1 is not None else np.nan
            rec["em_hz_pre"] = m0["hz"] if m0 else np.nan
            if unresolved:
                cls = "UNRESOLVED"
            elif m1 is None:
                cls = "MODE-SWITCH"
            elif rec["d_em"] <= -TAU:
                cls = "EM-STABILIZING"
            elif rec["d_em"] >= TAU:
                cls = "EM-DESTABILIZING"
            else:
                cls = "EM-NEUTRAL"
            fast = bool(m1 is not None and ((not HI.in_band(pre["lam_hz"])) or (not HI.in_band(post["lam_hz"])))
                        and abs(rec["d_full"] - rec["d_em"]) >= TAU)
            rec["em_class"] = cls
            rec["fast_flag"] = fast
            rec["label"] = cls if cls in ("UNRESOLVED", "MODE-SWITCH") else ("FAST-MODE DOMINATED" if fast else cls)
            rows.append(rec)
    return pd.DataFrame(rows)


def hypergraphs(data):
    rows = []
    for (pid, act), recs in data.items():
        status = {S: v["status"] for S, v in recs.items()}
        hz = {S: v["lam_hz"] for S, v in recs.items()}
        H, k = hyper(status, hz)
        He, ke = hyper(status, hz, em=True)
        rows.append({"pid": pid, "action": act, "H": H, "kappa": k, "H_EM": He, "kappa_EM": ke,
                     "status_H4": status["30+33+35+37"], "hz_H4": hz["30+33+35+37"], "alpha_H4": recs["30+33+35+37"]["alpha"]})
    hg = pd.DataFrame(rows)
    nom = hg[hg.action == "NOMINAL"].set_index("pid")
    hg["kappa_EM_nom"] = hg.pid.map(nom.kappa_EM)
    hg["kappa_nom"] = hg.pid.map(nom.kappa)
    hg["H_EM_nom"] = hg.pid.map(nom.H_EM)
    hg["status_H4_nom"] = hg.pid.map(nom.status_H4)
    hg["hz_H4_nom"] = hg.pid.map(nom.hz_H4)
    em_fail_nom = (hg.status_H4_nom == "UNSTABLE") & hg.hz_H4_nom.between(*HI.EM_BAND)
    hg["em_removal_H4"] = em_fail_nom & (hg.status_H4 == "STABLE") & (hg.action != "NOMINAL")
    hg["em_creation"] = (hg.action != "NOMINAL") & (hg.kappa_EM <= 3) & ~(hg.kappa_EM_nom <= 3) & ~hg.H_EM.str.startswith("BASE")
    hg["fast_creation"] = (hg.action != "NOMINAL") & (hg.kappa <= 3) & ~(hg.kappa_nom <= 3) & ~hg.em_creation & ~hg.H.str.startswith("BASE")
    return hg


def census_hem():
    """Re-analysis of the old E7 P4 V9 census: H vs H_EM (stored rightmost frequency)."""

    by = {}
    for r in I.Store("E07").all():
        if r.get("ok") and r["task"]["kind"] == "census":
            by.setdefault(r["task"]["action"], []).extend(r["records"])
    rows = []
    for a, lst in by.items():
        st = {x["S"]: x["status"] for x in lst}
        hz = {x["S"]: x.get("lam_hz") for x in lst}
        if len(st) < 512:
            continue
        subs = C.subsets(C.V9)
        if st["BASE"] != "STABLE":
            continue
        out = {}
        for em in (False, True):
            bad = [s for s in subs if st[C.label(s)] == "UNSTABLE" and (not em or HI.in_band(hz[C.label(s)]))]
            mini = [s for s in bad if not any(set(r) < set(s) for r in bad)]
            out["em" if em else "full"] = (min((len(e) for e in mini), default=np.inf), len(mini))
        rows.append({"action": a, "kappa_V9": out["full"][0], "n_edges_V9": out["full"][1], "kappa_EM_V9": out["em"][0], "n_edges_EM_V9": out["em"][1]})
    return pd.DataFrame(rows)


def static_scores_all(actions):
    from ibr_cycles.diagnosis.baselines import generalized_scr, nodal_metrics

    acts = {a["action"]: a for a in E7.actions()}
    acts["NOMINAL"] = {"action": "NOMINAL", "outage": [], "scale": {}}
    L0, order = C.laplacian_B()
    idx = {b: i for i, b in enumerate(order)}
    rows = {}
    for a in ["NOMINAL"] + list(actions):
        act = acts[a]
        out, scale = set(act["outage"]), {int(k): v for k, v in act["scale"].items()}
        L, _ = C.laplacian_B(scale or None, out or None)
        w = np.linalg.eigvalsh(L)
        net = E7.network_of(act) if (out or scale) else C.load_network()
        fl = C.pf(net)
        nm = nodal_metrics(net, fl)
        V = fl.voltages
        smax, loss = 0.0, 0.0
        for br in C.branches():
            if br["e"] in out:
                continue
            rho = scale.get(br["e"], 1.0)
            f, t = idx[br["f"]], idx[br["t"]]
            y = rho / complex(br["r"], br["x"])
            m = br["tap"]
            ysh = rho * 1j * br["b"] / 2
            i_f = (y + ysh) / m**2 * V[f] - y / m * V[t]
            i_t = -y / m * V[f] + (y + ysh) * V[t]
            sf, st = V[f] * np.conj(i_f), V[t] * np.conj(i_t)
            smax = max(smax, abs(sf), abs(st))
            loss += (sf + st).real
        Lp = np.linalg.pinv(L)
        core = [idx[b] for b in (30, 33, 35, 37, 39)]
        reff = [Lp[a_, a_] + Lp[b_, b_] - 2 * Lp[a_, b_] for k, a_ in enumerate(core) for b_ in core[k + 1:]]
        rows[a] = {"fiedler": float(w[1]), "kirchhoff": float(len(w) * np.sum(1 / w[1:])), "gscr_H4": float(generalized_scr(net, C.V4, fl)),
                   "min_scr_core": float(min(nm[b].scr for b in C.V4)), "max_flow": float(smax), "losses": float(loss), "mean_reff_core": float(np.mean(reff))}
    df = pd.DataFrame(rows).T
    return df - df.loc["NOMINAL"]


def run():
    data = load()
    cl = classify(data)
    cl.to_parquet(HI.RESULTS / "H14_topology_em.parquet", index=False)
    hg = hypergraphs(data)
    hg.to_csv(HI.RESULTS / "H15_hypergraphs.csv", index=False)
    cen = census_hem()
    cen.to_csv(HI.RESULTS / "H15_census_P4_H_vs_HEM.csv", index=False)
    h4 = cl[cl.S == "30+33+35+37"]
    out = {"n_actions": int(cl.action.nunique()), "n_policies": int(cl.pid.nunique()),
           "label_counts_all": cl.label.value_counts().to_dict(), "label_counts_H4": h4.label.value_counts().to_dict()}
    out["em_removal_H4_P4"] = hg[hg.em_removal_H4 & (hg.pid == "D01")].action.tolist()
    out["em_removal_H4_any_policy"] = hg[hg.em_removal_H4][["pid", "action"]].to_dict("records")
    out["em_creation"] = hg[hg.em_creation][["pid", "action", "H_EM"]].to_dict("records")
    out["n_em_creation"] = int(hg.em_creation.sum())
    out["n_fast_creation"] = int(hg.fast_creation.sum())
    out["Q13_yes"] = bool(len(out["em_removal_H4_P4"]) > 0 and out["n_em_creation"] > 0)
    out["n_H_changes"] = int(((hg.H != hg.pid.map(hg[hg.action == "NOMINAL"].set_index("pid").H)) & (hg.action != "NOMINAL")).sum())
    out["n_HEM_changes"] = int(((hg.H_EM != hg.H_EM_nom) & (hg.action != "NOMINAL")).sum())
    if len(cen):
        n0 = cen[cen.action == "NOMINAL"].iloc[0]
        out["census_P4_nominal"] = n0.to_dict()
        out["census_P4_kappa_EM_counts"] = cen[cen.action != "NOMINAL"].kappa_EM_V9.value_counts().to_dict()
        out["census_P4_kappa_counts"] = cen[cen.action != "NOMINAL"].kappa_V9.value_counts().to_dict()
    # ------------------------------------------------------------- H16 statics --
    st = static_scores_all(sorted(cl.action.unique()))
    st.to_csv(HI.RESULTS / "H16_static_scores.csv")
    pred = []
    for pid, g in h4.groupby("pid"):
        ok = g[g.em_class.isin(["EM-STABILIZING", "EM-DESTABILIZING", "EM-NEUTRAL"]) & ~g.fast_flag]
        for sc in st.columns:
            x_em = st.loc[ok.action, sc].to_numpy(float)
            rho_em = AN.spearman(x_em, ok.d_em) if len(ok) >= 10 else np.nan
            gg = g[np.isfinite(g.d_full)]
            rho_full = AN.spearman(st.loc[gg.action, sc].to_numpy(float), gg.d_full)
            pred.append({"pid": pid, "score": sc, "n_em": int(len(ok)), "rho_em": rho_em, "rho_full": rho_full})
    pr = pd.DataFrame(pred)
    pr.to_csv(HI.RESULTS / "H16_static_prediction.csv", index=False)
    summ = {}
    for sc, g in pr.groupby("score"):
        e = g.rho_em.dropna()
        summ[sc] = {"n_pol_em": int(len(e)), "median_rho_em": float(e.median()) if len(e) else np.nan,
                    "frac_abs_ge_0.6_em": float((e.abs() >= 0.6).mean()) if len(e) else np.nan,
                    "predicts_em": bool(len(e) and (e.abs() >= 0.6).mean() >= 0.75),
                    "median_rho_full": float(g.rho_full.median()), "predicts_full": bool((g.rho_full.abs() >= 0.6).mean() >= 0.75)}
    out["H16"] = summ
    out["H16_any_static_predicts_em"] = bool(any(v["predicts_em"] for v in summ.values()))
    HI.write_json("H15_topology_gate.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
