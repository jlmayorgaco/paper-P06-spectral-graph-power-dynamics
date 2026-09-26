# ruff: noqa: E501
"""R11 derivative validation and R12 line-ranking gate on IEEE-68 (docs/CDW68_PREREG_V1.md section 11), per model.

Truth: finite re-equilibrated -[alpha(gamma) - alpha(1)] of the target V68. Predictors: -D_tot, -D_fro, -D_conv (A);
the numerical total derivative -(alpha(1.001) - alpha(0.999))/0.002 (B); static indices (raw scores). The static
comparator is the index with the highest median Spearman on the 4 discovery conditions (gamma = 1.5), per model.
Post hoc (exploratory, Model A only): the same metrics for the rightmost EM-band mode.
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import itertools
import json

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr, wilcoxon

STATIC = ["S1_absP", "S2_absS", "S3_absz", "S4_elecdist", "S5_reff", "S6_fiedler", "S7_betweenness", "S8_dvdq", "S9_dgscr"]
GAMMAS = (1.10, 1.25, 1.50)
SEED_BOOT = 20260921


def metrics(score, truth, tau=R.TAU_MAT):
    s, t = np.asarray(score, float), np.asarray(truth, float)
    ok = np.isfinite(s) & np.isfinite(t)
    s, t = s[ok], t[ok]
    if s.size < 5 or np.all(s == s[0]) or np.all(t == t[0]):
        return {"spearman": np.nan, "kendall": np.nan, "top5": np.nan, "ndcg5": np.nan, "sign_acc": np.nan, "concord": np.nan, "n": int(s.size)}
    rho = spearmanr(s, t).statistic
    kt = kendalltau(s, t).statistic
    top_s, top_t = set(np.argsort(-s)[:5]), set(np.argsort(-t)[:5])
    rel = np.argsort(np.argsort(t)) / (len(t) - 1)  # truth rank in [0, 1]
    order = np.argsort(-s)[:5]
    ideal = np.sort(rel)[::-1][:5]
    disc = 1 / np.log2(np.arange(2, 7))
    ndcg = float((rel[order] * disc).sum() / (ideal * disc).sum())
    mat = np.abs(t) >= tau
    sign_acc = float(np.mean(np.sign(s[mat]) == np.sign(t[mat]))) if mat.any() else np.nan
    i, j = np.triu_indices(len(t), 1)
    conc = float(np.mean(np.sign(s[i] - s[j]) == np.sign(t[i] - t[j])))
    return {"spearman": float(rho), "kendall": float(kt), "top5": len(top_s & top_t) / 5, "ndcg5": ndcg, "sign_acc": sign_acc, "concord": conc, "n": int(s.size)}


def signflip_p(x, thr=0.20):
    d = np.asarray(x, float) - thr
    d = d[np.isfinite(d)]
    n = d.size
    if n == 0:
        return np.nan
    obs = np.median(d)
    cnt = 0
    tot = 0
    for signs in itertools.product((1, -1), repeat=n):
        tot += 1
        cnt += np.median(d * np.array(signs)) >= obs
    return cnt / tot


def boot_ci(x, n=10000, seed=SEED_BOOT):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return [np.nan, np.nan]
    rng = np.random.default_rng(seed)
    b = np.median(rng.choice(x, size=(n, x.size), replace=True), axis=1)
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def load_A():
    design = json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())
    split = {p["id"]: p["split"] for p in design["policies_A"]}
    der, fin, base = {}, {}, {}
    for r in R.Store("R10A").all():
        assert r["ok"], r["record"].get("error")
        t, rec = r["task"], r["record"]
        d = rec["der"]
        for it in d["items"]:
            der[(t["pid"], it["e"])] = {**it, "lam": d["lam"][0], "gap2": d["gap2"], "R0": d["R0"], "eig_cond": d["eig_cond"],
                                        "lam_em": (d["lam_em"] or [np.nan])[0]}
        for e, fg in rec["finite"].items():
            fin[(t["pid"], int(e))] = fg
    for r in R.Store("R10C").all():
        assert r["ok"], r["record"].get("error")
        for it in r["record"]["der"]["items"]:
            base[(r["task"]["pid"], it["e"])] = it["d_total"]
    rows = []
    for (pid, e), it in der.items():
        fg = fin.get((pid, e), {})
        row = {"model": "A", "cond": pid, "split": split[pid], "e": e, "lam": it["lam"], "gap2": it["gap2"], "R0": it["R0"], "eig_cond": it["eig_cond"],
               "P_tot": -it["d_total"], "P_fro": -it["d_frozen"], "P_conv": -base.get((pid, e), np.nan),
               "P_em_tot": -it.get("em_d_total", np.nan), "P_em_fro": -it.get("em_d_frozen", np.nan), "lam_em": it["lam_em"]}
        for gm in GAMMAS:
            f = fg.get(str(gm), {})
            row[f"T_{gm}"] = -(f.get("alpha", np.nan) - it["lam"]) if f.get("ok") else np.nan
            row[f"Ttr_{gm}"] = -(f.get("tracked_re", np.nan) - it["lam"]) if f.get("ok") else np.nan
            row[f"Tem_{gm}"] = -(f.get("em_tracked_re", np.nan) - it["lam_em"]) if f.get("ok") else np.nan
            row[f"nvout_{gm}"] = f.get("n_v_out", np.nan)
        fp, fm = fg.get("1.0001"), fg.get("0.9999")
        if fp and fm and fp.get("ok") and fm.get("ok"):
            row["FD_small"] = (fp["tracked_re"] - fm["tracked_re"]) / 2e-4
        rows.append(row)
    return pd.DataFrame(rows)


def load_B():
    design = json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())
    split = {p["id"]: p["split"] for p in design["policies_B"]}
    recs = [json.loads(p.read_text()) for p in sorted((R.RAW / "B68L").glob("*.json"))]
    base, val = {}, {}
    for r in recs:
        t, rec = r["task"], r["record"]
        a = rec.get("alpha") if rec.get("status") in ("STABLE", "UNSTABLE", "BOUNDARY_OR_UNRESOLVED") else np.nan
        if t["e"] == -1:
            base[t["pid"]] = a
        else:
            val[(t["pid"], t["e"], float(t["gamma"]))] = a
    rows = []
    for pid in split:
        a0 = base.get(pid, np.nan)
        for e in design["eligible_branches"]:
            row = {"model": "B", "cond": pid, "split": split[pid], "e": e, "lam": a0}
            ap, am = val.get((pid, e, 1.001), np.nan), val.get((pid, e, 0.999), np.nan)
            row["P_tot"] = -(ap - am) / 0.002
            sp, sm = val.get((pid, e, 1.0001)), val.get((pid, e, 0.9999))
            if sp is not None and sm is not None:
                row["FD_small"] = (sp - sm) / 2e-4
                row["D_h1e3"] = (ap - am) / 0.002
            for gm in GAMMAS:
                row[f"T_{gm}"] = -(val.get((pid, e, gm), np.nan) - a0)
            rows.append(row)
    return pd.DataFrame(rows)


def r11(df, model):
    q = df[df.FD_small.notna()]
    if model == "A":
        d, fd = -q.P_tot.to_numpy(), q.FD_small.to_numpy()
    else:
        d, fd = q.D_h1e3.to_numpy(), q.FD_small.to_numpy()
    if not len(q):
        return {"n": 0}
    sign = float(np.mean(np.sign(d) == np.sign(fd)))
    rel = np.abs(d - fd) / np.maximum(np.abs(fd), 1e-4)
    rho = float(spearmanr(d, fd).statistic)
    return {"n": int(len(q)), "sign_agreement": sign, "median_rel_err": float(np.median(rel)), "spearman": rho,
            "abs_err_median": float(np.median(np.abs(d - fd))), "fd_abs_median": float(np.median(np.abs(fd))),
            "pass": bool(sign >= 0.95 and np.median(rel) <= 0.05 and rho >= 0.95)}


def gate(df, model, static):
    df = df.merge(static, on="e", how="left")
    preds = ["P_tot"] + (["P_fro", "P_conv", "P_em_tot"] if model == "A" else []) + STATIC
    rows = []
    for cond, g in df.groupby("cond"):
        for gm in GAMMAS:
            for p in preds:
                truth = g[f"Tem_{gm}"] if p == "P_em_tot" else g[f"T_{gm}"]
                m = metrics(g[p], truth)
                rows.append({"model": model, "cond": cond, "split": g.split.iloc[0], "gamma": gm, "predictor": p, **m,
                             "truth_abs_max": float(np.nanmax(np.abs(truth))) if truth.notna().any() else np.nan,
                             "n_material": int((np.abs(truth) >= R.TAU_MAT).sum())})
    mt = pd.DataFrame(rows)
    disc = mt[(mt.split == "discovery") & (mt.gamma == 1.5) & mt.predictor.isin(STATIC)]
    sel = disc.groupby("predictor").spearman.median().idxmax()
    hold = mt[(mt.split == "holdout") & (mt.gamma == 1.5)]
    tot = hold[hold.predictor == "P_tot"].set_index("cond")
    st = hold[hold.predictor == sel].set_index("cond")
    adv = (tot.spearman - st.spearman).dropna()
    out = {"selected_static": sel, "discovery_static_medians": disc.groupby("predictor").spearman.median().round(3).to_dict(),
           "holdout_n": int(len(tot)), "median_rho_tot": float(tot.spearman.median()), "median_rho_static_sel": float(st.spearman.median()),
           "median_advantage": float(adv.median()), "advantage_ci95": boot_ci(adv.to_numpy()), "median_top5_tot": float(tot.top5.median()),
           "median_kendall_tot": float(tot.kendall.median()), "median_ndcg5_tot": float(tot.ndcg5.median()),
           "best_static_per_cond_median": float(hold[hold.predictor.isin(STATIC)].groupby("cond").spearman.max().median()),
           "truth_abs_max_median": float(tot.truth_abs_max.median()), "n_material_branches_median": float(tot.n_material.median()),
           "T2_signflip_p": signflip_p(adv.to_numpy()), "wilcoxon_p": float(wilcoxon(adv.to_numpy() - 0.20, alternative="greater").pvalue) if len(adv) > 5 else np.nan}
    for gm in GAMMAS:
        h = mt[(mt.split == "holdout") & (mt.gamma == gm) & (mt.predictor == "P_tot")]
        out[f"median_rho_tot_g{gm}"] = float(h.spearman.median())
    if model == "A":
        for p in ("P_fro", "P_conv", "P_em_tot"):
            out[f"median_rho_{p}"] = float(hold[hold.predictor == p].spearman.median())
        em = hold[hold.predictor == "P_em_tot"]
        out["posthoc_EM_tracked"] = {"median_rho": float(em.spearman.median()), "median_top5": float(em.top5.median()),
                                     "truth_abs_max_median": float(em.truth_abs_max.median()), "n_material_median": float(em.n_material.median())}
    out["R12_pass"] = bool(out["median_rho_tot"] >= 0.70 and out["median_advantage"] >= 0.20 and out["median_top5_tot"] >= 0.60)
    return out, mt


def main():
    static = pd.read_csv(R.RESULTS / "CDW68_R10_static.csv")[["e"] + STATIC]
    res, mts = {}, []
    for model, loader in (("A", load_A), ("B", load_B)):
        try:
            df = loader()
        except (FileNotFoundError, AssertionError) as ex:
            res[model] = {"error": repr(ex)[:200]}
            continue
        if not len(df):
            continue
        df.to_csv(R.RESULTS / f"CDW68_R10_{model}_branch_table.csv", index=False)
        res[model] = {"R11": r11(df, model)}
        g, mt = gate(df, model, static)
        res[model]["R12"] = g
        mts.append(mt)
    if mts:
        pd.concat(mts).to_csv(R.RESULTS / "CDW68_R12_metrics.csv", index=False)
    # Holm over F68 = {T1, T2A, T2B}; T1 comes from R09
    r9 = json.loads((R.RESULTS / "CDW68_R09_gate.json").read_text()) if (R.RESULTS / "CDW68_R09_gate.json").exists() else {}
    ps = {"T1": (r9.get("A", {}).get("T1_sign_test", {}) or {}).get("p_one_sided"),
          "T2A": res.get("A", {}).get("R12", {}).get("T2_signflip_p"), "T2B": res.get("B", {}).get("R12", {}).get("T2_signflip_p")}
    valid = {k: v for k, v in ps.items() if v is not None and np.isfinite(v)}
    order = sorted(valid, key=lambda k: valid[k])
    holm, run = {}, 0.0
    for r_, k in enumerate(order):
        run = max(run, min(1.0, (len(order) - r_) * valid[k]))
        holm[k] = run
    res["F68"] = {"raw": ps, "holm": holm, "note": "T1 not computable when no material preference exists (FULL stratum empty)"}
    R.write_json(R.RESULTS / "CDW68_R12_gate.json", res)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
