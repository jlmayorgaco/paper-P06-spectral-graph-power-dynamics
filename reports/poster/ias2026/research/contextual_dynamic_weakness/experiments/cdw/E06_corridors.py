# ruff: noqa: E501
"""E6 dynamic weak corridors from structure-only constructions (prereg E6)."""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C

PHASE = "E06"
SEED = 20260922
TXCORE = (34, 40, 42, 44)


def kmeans(X, k, seed=SEED, n_init=50, iters=300):
    rng = np.random.default_rng(seed)
    best, best_inertia = None, np.inf
    for _ in range(n_init):
        # k-means++ seeding
        cent = [X[rng.integers(len(X))]]
        for _ in range(1, k):
            d2 = np.min([((X - c) ** 2).sum(1) for c in cent], axis=0)
            cent.append(X[rng.choice(len(X), p=d2 / d2.sum())])
        cent = np.array(cent)
        for _ in range(iters):
            lab = np.argmin(((X[:, None, :] - cent[None]) ** 2).sum(2), axis=1)
            new = np.array([X[lab == j].mean(0) if (lab == j).any() else cent[j] for j in range(k)])
            if np.allclose(new, cent):
                break
            cent = new
        inertia = ((X - cent[lab]) ** 2).sum()
        if inertia < best_inertia - 1e-12:
            best, best_inertia = lab.copy(), inertia
    return best


def corridors() -> dict:
    """Structure-only corridor construction (frozen rule)."""

    L, order = C.laplacian_B()
    d = np.diag(L)
    Dm = np.diag(1 / np.sqrt(d))
    Ls = Dm @ L @ Dm
    w, U = np.linalg.eigh(Ls)
    br = C.branches()
    out = {}
    for k in (2, 3, 4):
        X = U[:, 1:k]
        X = X / np.linalg.norm(X, axis=1, keepdims=True).clip(1e-12)
        lab = kmeans(X, k)
        comm = {order[i]: int(lab[i]) for i in range(len(order))}
        for a in range(k):
            for b in range(a + 1, k):
                cut = [e["e"] for e in br if {comm[e["f"]], comm[e["t"]]} == {a, b}]
                if cut:
                    out[f"K{k}_{a}{b}"] = cut
    trans = [e["e"] for e in br if e["trans"] == 1]
    out["TXcore"] = list(TXCORE)
    out["TXother"] = [e for e in trans if e not in TXCORE]
    out["TXall"] = trans
    return out


def conditions():
    import E34_sens as E

    return [c for c in E.conditions() if c["target"] == "H4"]


def tasks():
    cor = corridors()
    I.atomic_write_json(I.RESULTS / "CDW_E6_corridor_definitions.json", cor)
    out = []
    for c in conditions():
        out.append({"phase": PHASE, **{k: c[k] for k in ("pid", "target", "split", "env", "draw")}, "corridors": cor})
    return out


def run_task(task):
    import E34_sens as E

    theta = C.policy(task["pid"])
    draw = E._draw(task)
    base = C.eval_portfolio(C.V4, theta, modes=False, draw=draw)
    res = {"alpha0": base["alpha"], "items": []}
    for name, cut in task["corridors"].items():
        r = C.eval_portfolio(C.V4, theta, modes=False, draw=draw, network=C.network_with(scale={e: 1.5 for e in cut}))
        res["items"].append({"corridor": name, "n_branches": len(cut), "alpha_x1.5": r["alpha"], "status": r["status"]})
    return res


def aggregate():
    link = pd.read_parquet(I.RESULTS / "CDW_E4_link_sensitivity.parquet")
    link = link[(link.semantics == "SPR") & (link.target == "H4")]
    base = pd.read_csv(I.RESULTS / "CDW_E4_baselines.csv") if (I.RESULTS / "CDW_E4_baselines.csv").exists() else None
    import E34_sens as E

    st = E.static_link_baselines(C.V4).set_index("e")
    rows = []
    for rec in I.Store(PHASE).all():
        if not rec.get("ok"):
            continue
        t = rec["task"]
        pkey = f"{t['pid']}|{t['env'] or '-'}|{t['draw'] if t['draw'] is not None else -1}"
        g = link[link.pkey == pkey].set_index("idx")
        for it in rec["items"]:
            cut = t["corridors"][it["corridor"]]
            dC = it["alpha_x1.5"] - rec["alpha0"]
            singles = (g.loc[cut, "large_x1.5_alpha"] - g.loc[cut, "alpha0"]).sum() if set(cut) <= set(g.index) else np.nan
            rows.append({"pid": t["pid"], "env": t["env"], "draw": t["draw"], "split": t["split"], "corridor": it["corridor"],
                         "n_branches": it["n_branches"], "delta_finite": dC, "sum_single_finite": singles,
                         "NA": abs(dC - singles) / max(abs(singles), C.TAU_RES) if np.isfinite(singles) else np.nan,
                         "sum_dtotal": g.loc[cut, "d_total"].sum() * 0.5 if set(cut) <= set(g.index) else np.nan,
                         "sum_static_S4": st.loc[cut, "S4_elecdist"].sum(), "sum_static_S1": st.loc[cut, "S1_absP"].sum()})
    df = pd.DataFrame(rows)
    df.to_csv(I.RESULTS / "CDW_E6_corridors.csv", index=False)
    del base
    # robust weak corridor rule on holdout conditions
    h = df[df.split.str.startswith("holdout")].copy()
    h["cond"] = h.pid + "|" + h.env.fillna("-").astype(str) + "|" + h.draw.fillna(-1).astype(int).astype(str)
    top3, improves = {}, {}
    rhos, rhos_single, rhos_static = [], [], []
    for cond, g in h.groupby("cond"):
        g = g.sort_values("delta_finite")
        for name in g.corridor.head(3):
            top3[name] = top3.get(name, 0) + 1
        for r in g.itertuples():
            improves[r.corridor] = improves.get(r.corridor, 0) + (r.delta_finite <= -C.TAU_MAT)
        rhos.append(AN.spearman(-g.sum_dtotal, -g.delta_finite))
        rhos_single.append(AN.spearman(-g.sum_single_finite, -g.delta_finite))
        rhos_static.append(AN.spearman(g.sum_static_S4, -g.delta_finite))
    n = h.cond.nunique()
    robust = [c for c in top3 if top3[c] / n >= 0.75 and improves.get(c, 0) / n >= 0.75]
    gate = {"n_holdout_conditions": int(n), "top3_frequency": {k: v / n for k, v in sorted(top3.items(), key=lambda x: -x[1])},
            "robust_weak_corridors": robust, "H3_true": bool(robust),
            "median_rho_sum_dtotal": float(np.nanmedian(rhos)), "median_rho_sum_single_finite": float(np.nanmedian(rhos_single)),
            "median_rho_sum_static_S4": float(np.nanmedian(rhos_static)),
            "median_NA": float(h.NA.median()), "max_NA": float(h.NA.max())}
    gate["corridor_ranking_predictive"] = bool(gate["median_rho_sum_dtotal"] >= 0.70)
    I.atomic_write_json(I.RESULTS / "CDW_E6_gates.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=1, default=str))
