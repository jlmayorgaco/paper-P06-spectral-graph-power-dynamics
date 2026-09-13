# ruff: noqa: E501
"""H07 GOLD-B on the new holdout (prereg H6, H7) and H08 frozen-vs-total decomposition."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001  (first: sets sys.path for the CDW modules)

import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _hstats as ST
import _infra as I
import E34_sens as E34

GAMMAS = ("1.10", "1.25", "1.50")
SEL = "S1_absP"  # selected on discovery in the old E4, never reselected
RANGE = E34.RANGE


def load_h06():
    rows = []
    for r in I.Store("H_H06").all():
        if not r.get("ok"):
            continue
        t = r["task"]
        for it in r["items"]:
            row = {k: t.get(k) for k in ("pid", "target", "source", "env", "draw", "family")}
            row.update({k: v for k, v in it.items() if k not in ("d_frozen", "d_total", "d_port", "semantics")})
            for k in ("d_frozen", "d_total", "d_port"):
                if it.get(k) is not None:
                    row[k] = it[k][0]
            row["lam_re"], row["lam_im"] = (r.get("lam") or [np.nan, np.nan])
            row["gap2"] = r.get("gap2")
            row["R0"] = r.get("R0")
            rows.append(row)
    df = pd.DataFrame(rows)
    df["pkey"] = df.pid + "|" + df.env.fillna("-").astype(str) + "|" + df.draw.fillna(-1).astype(int).astype(str)
    df["cond"] = df.pkey + "|" + df.target
    df["cluster"] = np.where(df.env.notna(), "ENV_" + df.env.fillna(""), df.pid)
    return df


def pair_concordance(p, t, tau=C.TAU_RES):
    p, t = np.asarray(p, float), np.asarray(t, float)
    ok = np.isfinite(p) & np.isfinite(t)
    p, t = p[ok], t[ok]
    good = tot = 0.0
    for a in range(p.size):
        for b in range(a + 1, p.size):
            if abs(t[a] - t[b]) < tau:
                continue
            tot += 1
            s = np.sign(p[a] - p[b]) * np.sign(t[a] - t[b])
            good += 1.0 if s > 0 else (0.5 if s == 0 else 0.0)
    return good / tot if tot else np.nan


def ndcg(p, t, k=5):
    p, t = np.asarray(p, float), np.asarray(t, float)
    ok = np.isfinite(p) & np.isfinite(t)
    p, t = p[ok], t[ok]
    if p.size < k:
        return np.nan
    rel = t - t.min()
    disc = 1 / np.log2(np.arange(2, k + 2))
    dcg = (rel[np.argsort(-p)[:k]] * disc).sum()
    idcg = (np.sort(rel)[::-1][:k] * disc).sum()
    return float(dcg / idcg) if idcg > 0 else np.nan


def metrics(p, t, dynamic):
    mat = np.isfinite(t) & (np.abs(t) >= C.TAU_MAT) & np.isfinite(p)
    return {"spearman": AN.spearman(p, t), "kendall": AN.kendall(p, t), "concord": pair_concordance(p, t),
            "top3": AN.topk_precision(p, t, 3), "top5": AN.topk_precision(p, t, 5), "ndcg5": ndcg(p, t),
            "sign_acc": float((np.sign(p[mat]) == np.sign(t[mat])).mean()) if (dynamic and mat.any()) else np.nan,
            "n_finite": int(np.isfinite(t).sum()), "n_material": int(mat.sum())}


def link_table(df):
    link = df[df.family == "link"].copy()
    conv = df[df.family == "conv"]
    port = df[df.family == "port"]
    statics = {T: E34.static_link_baselines(TT).set_index("e") for T, TT in (("H4", C.V4), ("V9", C.V9))}
    rows = []
    for cond, g in link.groupby("cond"):
        g = g.sort_values("idx")
        if len(g) != 46:
            continue
        T = g.target.iloc[0]
        st = statics[T].loc[g.idx.to_numpy()]
        lam_hz = abs(g.lam_im.iloc[0]) / (2 * np.pi)
        preds = {c: (st[c].to_numpy(float), False) for c in st.columns if c.startswith("S")}
        preds["Dfrozen"] = (-g.d_frozen.to_numpy(float), True)
        preds["Dtotal"] = (-g.d_total.to_numpy(float), True)
        cv = conv[conv.pkey == g.pkey.iloc[0]].sort_values("idx")
        if len(cv) == 46:
            preds["Dconv"] = (-cv.d_total.to_numpy(float), True)
        pt = port[port.cond == cond].sort_values("idx")
        if len(pt) == 46 and pt.d_port.notna().all():
            preds["Dport"] = (-pt.d_port.to_numpy(float), True)
        for gam in GAMMAS:
            ok = g[f"g{gam}_ok"].fillna(False).astype(bool).to_numpy()
            for truth_kind in ("FULL", "SAME"):
                if truth_kind == "FULL":
                    t = -(g[f"g{gam}_alpha"].to_numpy(float) - g.alpha0.to_numpy(float))
                else:
                    t = -(g[f"g{gam}_tracked_re"].to_numpy(float) - g.lam_re.to_numpy(float))
                t = np.where(ok, t, np.nan)
                if np.isfinite(t).sum() < 20:
                    continue
                for name, (p, dyn) in preds.items():
                    rows.append({"cond": cond, "pid": g.pid.iloc[0], "target": T, "cluster": g.cluster.iloc[0],
                                 "set": "fresh_draws" if g.env.notna().iloc[0] else "new", "gamma": gam, "truth": truth_kind,
                                 "predictor": name, "em_cond": bool(HI.in_band(lam_hz)), "lam_hz": lam_hz,
                                 "near_tie": bool(g.gap2.iloc[0] < 1e-3), "alpha0": float(g.alpha0.iloc[0]),
                                 "n_inadmissible": int((~g[f"g{gam}_admissible"].fillna(False).astype(bool)).sum()),
                                 **metrics(p, t, dyn)})
    return pd.DataFrame(rows)


def gate_on(m, label):
    """GOLD-B gate on a metrics subset (H4 targets, gamma 1.50, FULL unless filtered before)."""

    piv = m.pivot_table(index=["cond", "cluster"], columns="predictor", values="spearman").reset_index()
    pv5 = m.pivot_table(index="cond", columns="predictor", values="top5")
    if not len(piv) or "Dtotal" not in piv:
        return {f"{label}_n": 0}
    statics = [c for c in piv.columns if isinstance(c, str) and c.startswith("S")]
    diff = (piv["Dtotal"] - piv[SEL]).to_numpy(float)
    pt, lo, hi, lo5 = ST.cluster_boot(diff, piv.cluster.to_numpy())
    oracle = piv[statics].max(axis=1).to_numpy(float)
    out = {f"{label}_n": int(len(piv)), f"{label}_n_clusters": int(piv.cluster.nunique()),
           f"{label}_median_rho_Dtotal": float(np.nanmedian(piv["Dtotal"])),
           f"{label}_median_rho_Dfrozen": float(np.nanmedian(piv["Dfrozen"])),
           f"{label}_median_rho_Dconv": float(np.nanmedian(piv["Dconv"])) if "Dconv" in piv else np.nan,
           f"{label}_median_rho_S1": float(np.nanmedian(piv[SEL])),
           f"{label}_median_rho_oracle_static": float(np.nanmedian(oracle)),
           f"{label}_median_diff": pt, f"{label}_diff_ci95": [lo, hi], f"{label}_diff_lo5": lo5,
           f"{label}_mean_diff": float(np.nanmean(diff)),
           f"{label}_median_diff_vs_oracle_static": float(np.nanmedian(piv["Dtotal"].to_numpy(float) - oracle)),
           f"{label}_median_top5_Dtotal": float(np.nanmedian(pv5["Dtotal"])),
           f"{label}_median_top5_S1": float(np.nanmedian(pv5[SEL])),
           f"{label}_frac_cond_Dtotal_beats_every_static": float(np.mean(piv["Dtotal"].to_numpy(float) > oracle))}
    out[f"{label}_pass"] = bool(lo > 0.20 and out[f"{label}_median_rho_Dtotal"] >= 0.80 and out[f"{label}_median_top5_Dtotal"] >= 0.70)
    obs, p = ST.cluster_signflip_median(diff, piv.cluster.to_numpy(), shift=0.20)
    out[f"{label}_T2_signflip_p"] = p
    for k in ("kendall", "concord", "top3", "ndcg5", "sign_acc"):
        pk = m.pivot_table(index="cond", columns="predictor", values=k)
        out[f"{label}_median_{k}_Dtotal"] = float(np.nanmedian(pk["Dtotal"])) if "Dtotal" in pk else np.nan
        out[f"{label}_median_{k}_S1"] = float(np.nanmedian(pk[SEL])) if SEL in pk else np.nan
    return out


def decomposition(df):
    link = df[df.family == "link"].copy()
    link["indirect"] = link.d_total - link.d_frozen
    link["ratio"] = link.indirect.abs() / link.d_frozen.abs().clip(lower=1e-12)
    link["cancel"] = np.sign(link.indirect) != np.sign(link.d_frozen)
    link["flip"] = np.sign(link.d_total) != np.sign(link.d_frozen)
    per, ex = [], []
    for cond, g in link.groupby("cond"):
        g = g.copy()
        g["rk_t"] = (-g.d_total).rank(ascending=False)
        g["rk_f"] = (-g.d_frozen).rank(ascending=False)
        g["disp"] = (g.rk_t - g.rk_f).abs()
        sel = g[g.d_total.abs() * RANGE["line"] >= C.TAU_MAT]
        kt = AN.kendall(g.d_frozen, g.d_total)
        dis = float(sel.flip.mean()) if len(sel) >= 3 else np.nan
        per.append({"cond": cond, "target": g.target.iloc[0], "set": "fresh_draws" if g.env.notna().iloc[0] else "new",
                    "kendall_frozen_total": kt, "sign_disagree_material": dis,
                    "material": bool((dis >= 0.10) or (kt <= 0.6)) if np.isfinite(dis) else bool(kt <= 0.6),
                    "median_ratio": float(g.ratio.median()), "frac_ratio_gt1": float((g.ratio > 1).mean()),
                    "frac_cancel": float(g.cancel.mean()), "n_flip": int(g.flip.sum()), "max_disp": float(g.disp.max())})
        top5 = g.reindex(g.d_total.abs().sort_values(ascending=False).index[:5])
        a = top5.sort_values("disp").iloc[0]
        b = g.sort_values("disp").iloc[-1]
        fl = g[g.flip]
        c = fl.reindex(fl.d_total.abs().sort_values(ascending=False).index[:1]) if len(fl) else None
        ex.append({"cond": cond, "agree": int(a.idx), "differ": int(b.idx), "flip": int(c.idx.iloc[0]) if c is not None and len(c) else -1})
    per, ex = pd.DataFrame(per), pd.DataFrame(ex)
    link[["cond", "pid", "target", "env", "draw", "idx", "d_frozen", "d_total", "indirect", "ratio", "cancel", "flip"]].to_parquet(
        HI.RESULTS / "H08_link_decomposition.parquet", index=False)
    per.to_csv(HI.RESULTS / "H08_link_decomposition_per_condition.csv", index=False)
    ex.to_csv(HI.RESULTS / "H08_examples_per_condition.csv", index=False)
    return per, ex


def node_decomposition():
    rows = []
    for r in I.Store("H_H08n").all():
        if not r.get("ok"):
            continue
        t = r["task"]
        for it in r["items"]:
            rows.append({"pid": t["pid"], "kind": it["kind"], "idx": it["idx"], "d_frozen": it["d_frozen"][0], "d_total": it["d_total"][0],
                         "fd_small": it.get("fd_small")})
    nd = pd.DataFrame(rows)
    nd["indirect"] = nd.d_total - nd.d_frozen
    nd["rel_frozen_total"] = (nd.d_frozen - nd.d_total).abs() / nd.d_total.abs().clip(lower=1e-12)
    nd["ratio"] = nd.indirect.abs() / nd.d_frozen.abs().clip(lower=1e-12)
    nd["iv_ok"] = ((nd.fd_small - nd.d_total).abs() <= 0.05 * nd.d_total.abs()) | ((nd.fd_small - nd.d_total).abs() * nd.kind.map(RANGE) <= 1e-4)
    nd.to_parquet(HI.RESULTS / "H08_node_decomposition.parquet", index=False)
    ctrl = nd[nd.kind.isin(["g", "pll", "ka"])]
    net = nd[nd.kind.isin(["vset", "load"])]
    return {"node_R1_max_rel_frozen_minus_total_controller": float(ctrl.rel_frozen_total.max()),
            "node_R1_pass": bool((ctrl.rel_frozen_total <= 1e-6).all()),
            "node_vset_max_abs_frozen": float(nd[nd.kind == "vset"].d_frozen.abs().max()),
            "node_opoint_median_ratio_indirect_over_frozen": {k: float(v.ratio.median()) for k, v in net.groupby("kind")},
            "node_opoint_frac_indirect_dominates": {k: float((v.indirect.abs() > v.d_frozen.abs()).mean()) for k, v in net.groupby("kind")},
            "node_controller_median_ratio": {k: float(v.ratio.median()) for k, v in ctrl.groupby("kind")},
            "node_IV_frac_ok": float(nd.iv_ok.mean())}


def run():
    df = load_h06()
    df.to_parquet(HI.RESULTS / "H06_links_raw.parquet", index=False)
    lk = df[df.family == "link"]
    iv = lk[lk.fd_small.notna()]
    ivok = ((iv.fd_small - iv.d_total).abs() <= 0.05 * iv.d_total.abs()) | ((iv.fd_small - iv.d_total).abs() <= 1e-4)
    m = link_table(df)
    m.to_csv(HI.RESULTS / "H06_metrics.csv", index=False)
    gate = {"IV_frac_ok_links": float(ivok.mean()), "IV_n": int(len(iv)),
            "n_conditions": int(lk.cond.nunique()), "n_conditions_failed_R0": int((lk.groupby("cond").R0.max() > 1e-8).sum()),
            "frac_steps_converged_g1.50": float(lk["g1.50_ok"].fillna(False).astype(bool).mean()),
            "frac_steps_admissible_g1.50": float(lk["g1.50_admissible"].fillna(False).astype(bool).mean())}
    r0 = lk.groupby("cond").R0.max()
    elig = set(r0[r0 <= 1e-8].index)  # prereg H6 eligibility (R0 <= 1e-8), applied literally
    gate["n_conditions_R0_eligible"] = int(len(elig))
    gate["R0_max_all"] = float(r0.max())
    allh4 = m[(m.target == "H4") & (m.gamma == "1.50") & (m.truth == "FULL")]
    gate.update(gate_on(allh4[allh4.cond.isin(elig)], "GB_prereg"))
    gate["GB_hardened_pass"] = gate.get("GB_prereg_pass", False)
    prim = allh4  # sensitivity analysis: every solved condition (R0 <= 2.2e-7; IV check 100 %)
    gate.update(gate_on(prim, "GB_primary"))
    gate.update(gate_on(prim[prim.em_cond], "GB_EM"))
    gate.update(gate_on(m[(m.target == "H4") & (m.gamma == "1.50") & (m.truth == "SAME")], "GB_SAME"))
    gate.update(gate_on(prim[~prim.near_tie], "GB_no_near_tie"))
    gate.update(gate_on(prim[prim.set == "new"], "GB_new_policies"))
    gate.update(gate_on(prim[prim.set == "fresh_draws"], "GB_fresh_draws"))
    for gam in ("1.10", "1.25"):
        gate.update(gate_on(m[(m.target == "H4") & (m.gamma == gam) & (m.truth == "FULL")], f"GB_g{gam}"))
    gate.update(gate_on(m[(m.target == "V9") & (m.gamma == "1.50") & (m.truth == "FULL")], "GB_V9"))
    per, ex = decomposition(df)
    h4 = per[per.target == "H4"]
    gate["H8_frac_material_H4"] = float(h4.material.mean())
    gate["H8_median_kendall_frozen_total_H4"] = float(h4.kendall_frozen_total.median())
    gate["H8_median_frac_ratio_gt1_H4"] = float(h4.frac_ratio_gt1.median())
    gate["H8_conditions_with_sign_flip_H4"] = float((h4.n_flip > 0).mean())
    gate["H8_example_modes"] = {k: int(ex[k].mode().iloc[0]) for k in ("agree", "differ", "flip") if len(ex)}
    gate.update(node_decomposition())
    HI.write_json("H07_goldb_gate.json", gate)
    return gate


if __name__ == "__main__":
    import json

    print(json.dumps(run(), indent=1, default=str))
