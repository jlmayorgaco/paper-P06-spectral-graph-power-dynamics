# ruff: noqa: E501
"""H10 size-matched corridor null (gate) and H11 corridor cross-validation table (prereg H9-H11).

Effect = -(alpha(H4 after) - alpha(H4 base)) (positive = stabilizing). Frozen corridor names that
share an identical branch set (K3_01=K4_03, K3_02=K4_02, K3_12=K4_23) are ranked as one set."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _cdw as C
import _infra as I
import H09_corridors as H9

SETS_NEW = ("new", "fresh_draws")
SETS_OLD = ("old", "old_draws")


def load_h09():
    rows = []
    for r in I.Store("H_H09").all():
        if not r.get("ok"):
            continue
        t = r["task"]
        key = H9.cond_key(t)
        for it in r["items"]:
            rows.append({"cond": key, "set": t["set"], "pid": t["pid"], "corridor": it["corridor"], "size": it["size"],
                         "budget": it["budget"], "alpha": it["alpha"], "alpha0": r["alpha0"], "lam_hz0": r["lam_hz0"],
                         "status0": r["status0"], "eff": -(it["alpha"] - r["alpha0"])})
    return pd.DataFrame(rows)


def load_h10():
    groups = {}
    eff = {}
    meta = {}
    for r in I.Store("H_H10").all():
        if not r.get("ok"):
            continue
        t = r["task"]
        key = H9.cond_key(t)
        meta[key] = {"set": t["set"], "lam_hz0": r["lam_hz0"], "alpha0": r["alpha0"]}
        for g, res in zip(t["groups"], r["res"], strict=True):
            eff.setdefault(key, {})[tuple(g)] = -(res[0] - r["alpha0"]) if res[0] is not None else np.nan
            groups[tuple(g)] = True
    return eff, meta


def percentile(x, null):
    null = np.asarray([v for v in null if np.isfinite(v)])
    if not null.size or not np.isfinite(x):
        return np.nan, 0
    return float(((null < x).sum() + 0.5 * (null == x).sum()) / null.size), int(null.size)


def run():
    h9 = load_h09()
    h9.to_csv(HI.RESULTS / "H09_corridors.csv", index=False)
    cor = H9.corridor_defs()
    uniq = {}
    for name, cut in cor.items():
        uniq.setdefault(tuple(sorted(cut)), []).append(name)
    canon = {name: "=".join(names) for names in uniq.values() for name in names}
    # ---------------------------------------------------------------- H11 table --
    prim = h9[h9.budget == "L1_0.50"].copy()
    prim["uset"] = prim.corridor.map(canon)
    rows = []
    for bud, g0 in h9.assign(uset=h9.corridor.map(canon)).groupby("budget"):
        u = g0.drop_duplicates(["cond", "uset"])
        u = u.assign(rank=u.groupby("cond").eff.rank(ascending=False, method="min"))
        for (st, us), g in u.groupby(["set", "uset"]):
            rows.append({"budget": bud, "set": st, "corridor_set": us, "size": int(g["size"].iloc[0]), "n_cond": int(len(g)),
                         "stab_frac": float((g.eff >= C.TAU_MAT).mean()), "top3_frac": float((g["rank"] <= 3).mean()),
                         "median_eff": float(g.eff.median()), "median_rank": float(g["rank"].median())})
    tab = pd.DataFrame(rows)
    # ---------------------------------------------------------------- H10 null ---
    eff, meta = load_h10()
    ng = HI.null_groups()
    recs = []
    for key, e in eff.items():
        pc = prim[prim.cond == key].drop_duplicates("uset")
        for r in pc.itertuples():
            k = int(r.size)
            for fam in ("A", "B"):
                null = [e.get(tuple(sorted(g)), np.nan) for g in ng[fam][str(k)]]
                p, n = percentile(r.eff, null)
                recs.append({"cond": key, "set": meta[key]["set"], "corridor_set": r.uset, "size": k, "family": fam,
                             "percentile": p, "n_null": n, "eff": r.eff, "em_cond": HI.in_band(meta[key]["lam_hz0"])})
    nl = pd.DataFrame(recs)
    nl.to_csv(HI.RESULTS / "H10_null_percentiles.csv", index=False)
    # consistency: k=1 null (H10) vs single-branch x1.5 in H06 (new policies, H4)
    try:
        lk = pd.read_parquet(HI.RESULTS / "H06_links_raw.parquet")
        lk = lk[(lk.family == "link") & (lk.target == "H4")]
        rows_c = []
        for r in lk.to_dict("records"):
            fresh = isinstance(r["env"], str)
            key = f"{r['pid']}|{'FRESH' if fresh else '-'}|{r['env'] if fresh else '-'}|{int(r['draw']) if pd.notna(r['draw']) else -1}"
            if key in eff and (int(r["idx"]),) in eff[key] and r.get("g1.50_ok"):
                rows_c.append({"cond": key, "e": int(r["idx"]), "h10": eff[key][(int(r["idx"]),)], "h06": -(r["g1.50_alpha"] - r["alpha0"])})
        cc = pd.DataFrame(rows_c)
        cc.to_csv(HI.RESULTS / "H10_k1_consistency.csv", index=False)
        chk = float((cc.h10 - cc.h06).abs().max()) if len(cc) else None
    except FileNotFoundError:
        chk = "H06 not available"
    gate = {"k1_null_vs_H06_max_abs_diff": chk}
    for us in sorted(nl.corridor_set.unique()):
        for fam in ("A", "B"):
            for label, sets in (("new", SETS_NEW), ("old", SETS_OLD)):
                q = nl[(nl.corridor_set == us) & (nl.family == fam) & nl.set.isin(sets)]
                gate[f"{us}|{fam}|{label}_frac_ge95"] = float((q.percentile >= 0.95).mean()) if len(q) else np.nan
                gate[f"{us}|{fam}|{label}_median_pct"] = float(q.percentile.median()) if len(q) else np.nan
                qe = q[q.em_cond]
                gate[f"{us}|{fam}|{label}_EM_frac_ge95"] = float((qe.percentile >= 0.95).mean()) if len(qe) else np.nan
    tested = [canon[n] for n in ("TXother", "TXall", "K2_01", "K3_01", "K3_02", "K3_12", "K4_13", "TXcore")]
    robust = sorted({us for us in tested if gate.get(f"{us}|A|new_frac_ge95", 0) >= 0.75})
    gate["H10_robust_weak_corridors"] = robust
    gate["H10_pass"] = bool(robust)
    for us in tested:
        m = nl[(nl.corridor_set == us) & (nl.family == "A") & nl.set.isin(SETS_NEW)].percentile
        tab.loc[tab.corridor_set == us, "new_median_null_pct_A"] = float(m.median()) if len(m) else np.nan
    tab.to_csv(HI.RESULTS / "H11_corridor_table.csv", index=False)
    HI.write_json("H10_gate.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
