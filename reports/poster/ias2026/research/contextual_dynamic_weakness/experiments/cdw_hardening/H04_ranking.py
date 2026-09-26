# ruff: noqa: E501
"""H04 fixed-ranking insufficiency (oracle DP ranking), stratified FULL / EM / same-mode, and the
train-policy x test-policy transfer matrix (prereg H4, L5, L6). Uses H03_marginals.parquet."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001  (first: sets sys.path for the CDW modules)

import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C

V9 = C.V9
IDX = {u: k for k, u in enumerate(V9)}
STRATA = {"FULL": "lvA", "EM": None, "SAME": "lvC"}


def stratum(mg, s):
    if s == "FULL":
        return mg[mg.lvA]
    if s == "EM":
        return mg[mg.lvA & mg.em]
    return mg[mg.lvC]


def contexts(g):
    return {S: dict(zip(gg.i, gg.delta, strict=True)) for S, gg in g.groupby("mask")}


def wmatrix(ctx, tau=C.TAU_MAT):
    W = np.zeros((9, 9))
    for d in ctx.values():
        for a in d:
            for b in d:
                if a != b and d[a] < d[b] - tau:
                    W[IDX[a], IDX[b]] += 1
    return W


def accuracy(order, W):
    pos = {u: r for r, u in enumerate(order)}
    tot = W.sum()
    if tot == 0:
        return np.nan
    good = sum(W[i, j] for i in range(9) for j in range(9) if pos[i] < pos[j])
    return good / tot


def regret(order, ctx, tau=C.TAU_MAT):
    pos = {V9[u]: r for r, u in enumerate(order)}
    regs = []
    for d in ctx.values():
        if len(d) < 2:
            continue
        choice = min(d, key=lambda u: pos[u])
        regs.append(d[choice] - min(d.values()))
    regs = np.array(regs)
    if not regs.size:
        return np.nan, np.nan, np.nan, 0
    return float((regs >= tau).mean()), float(np.median(regs)), float(regs.max()), int(regs.size)


def run():
    mg = pd.read_parquet(HI.RESULTS / "H03_marginals.parquet")
    pol = pd.read_csv(HI.RESULTS / "H03_policy_summary.csv").set_index("pid")
    rows, orders, ctxs, Ws = [], {}, {}, {}
    for pid, g in mg.groupby("pid"):
        for s in STRATA:
            ctx = contexts(stratum(g, s))
            W = wmatrix(ctx)
            order, val = AN.optimal_linear_order(W)
            pstar = val / W.sum() if W.sum() else np.nan
            fr, med, mx, n = regret(order, ctx)
            rows.append({"pid": pid, "split": HI.split_of(pid), "base_stable": bool(pol.loc[pid, "base_stable"]), "stratum": s,
                         "p_star": pstar, "n_triples": int(W.sum()), "n_ctx": n, "frac_regret": fr, "median_regret": med,
                         "max_regret": mx, "order": "-".join(str(V9[j]) for j in order)})
            if s == "FULL":
                orders[pid], ctxs[pid], Ws[pid] = order, ctx, W
    rk = pd.DataFrame(rows)
    rk.to_csv(HI.RESULTS / "H04_ranking.csv", index=False)
    # ----------------------------------------------------------------- transfer --
    P = [p for p in sorted(orders) if pol.loc[p, "base_stable"]]
    acc = pd.DataFrame(index=P, columns=P, dtype=float)
    reg = pd.DataFrame(index=P, columns=P, dtype=float)
    for a in P:
        for b in P:
            acc.loc[a, b] = accuracy(orders[a], Ws[b])
            reg.loc[a, b] = regret(orders[a], ctxs[b])[0]
    acc.to_csv(HI.RESULTS / "H04_transfer_accuracy.csv")
    reg.to_csv(HI.RESULTS / "H04_transfer_regret.csv")
    th = {p: HI.theta_of(p) for p in P}
    X = np.array([[np.sqrt(th[p][0]), (th[p][1] - 0.5) / 1.8, (th[p][2] - 0.5) / 2.5, th[p][3] / 2] for p in P])
    Dm = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    A = acc.to_numpy(float)
    deg = np.diag(A)[None, :] - A  # p*(b) - acc(a->b), rows a (train), cols b (test)
    off = ~np.eye(len(P), dtype=bool)
    rho = AN.spearman(Dm[off], deg[off])
    rng = np.random.default_rng(HI.SEED_PERM)
    perm = []
    for _ in range(10000):
        q = rng.permutation(len(P))
        perm.append(AN.spearman(Dm[np.ix_(q, q)][off], deg[off]))
    mantel_p = float((np.abs(perm) >= abs(rho)).mean())
    R = reg.to_numpy(float)
    off_reg_by_test = np.array([np.nanmedian(R[off[:, j], j]) for j in range(len(P))])
    diag_reg = np.diag(R)
    splits = np.array([HI.split_of(p) for p in P])
    newmask = splits == "new"
    full = rk[rk.stratum == "FULL"]
    gate = {}
    for split in ("new", "old", "discovery"):
        for s in STRATA:
            q = rk[(rk.split == split) & rk.base_stable & (rk.stratum == s)]
            gate[f"{split}_{s}_n"] = int(len(q))
            gate[f"{split}_{s}_median_pstar"] = float(q.p_star.median()) if len(q) else np.nan
            gate[f"{split}_{s}_median_frac_regret"] = float(q.frac_regret.median()) if len(q) else np.nan
            gate[f"{split}_{s}_range_pstar"] = [float(q.p_star.min()), float(q.p_star.max())] if len(q) else []
    gate["H4_new_pass"] = bool(gate["new_FULL_median_pstar"] <= 0.90 and gate["new_FULL_median_frac_regret"] >= 0.10)
    gate["L5_pass"] = gate["H4_new_pass"]
    gate["transfer_n_policies"] = len(P)
    gate["transfer_median_diag_acc"] = float(np.nanmedian(np.diag(A)))
    gate["transfer_median_offdiag_acc"] = float(np.nanmedian(A[off]))
    gate["transfer_median_diag_regret"] = float(np.nanmedian(diag_reg))
    gate["transfer_median_offdiag_regret_by_test"] = float(np.nanmedian(off_reg_by_test))
    gate["transfer_new_median_offdiag_regret_by_test"] = float(np.nanmedian(off_reg_by_test[newmask])) if newmask.any() else np.nan
    gate["transfer_new_median_diag_regret"] = float(np.nanmedian(diag_reg[newmask])) if newmask.any() else np.nan
    gate["L6_pass"] = bool(gate["transfer_new_median_offdiag_regret_by_test"] >= 0.10
                           and gate["transfer_new_median_offdiag_regret_by_test"] > gate["transfer_new_median_diag_regret"])
    gate["distance_vs_degradation_spearman"] = rho
    gate["distance_vs_degradation_mantel_p_descriptive"] = mantel_p
    gate["n_distinct_optimal_orders_full"] = int(full[full.base_stable].order.nunique())
    # anatomy of material regret (new holdout, FULL): is the ranked choice a jump out of the EM band?
    tot = bad = fast = fastbig = 0
    for pid in [p for p in P if HI.split_of(p) == "new"]:
        pos = {V9[u]: r for r, u in enumerate(orders[pid])}
        g = mg[(mg.pid == pid) & mg.lvA]
        for _, gg in g.groupby("mask"):
            if len(gg) < 2:
                continue
            tot += 1
            d, hz = dict(zip(gg.i, gg.delta, strict=True)), dict(zip(gg.i, gg.hz_Si, strict=True))
            ch = min(d, key=lambda u: pos[u])
            if d[ch] - min(d.values()) >= C.TAU_MAT:
                bad += 1
                out = not (HI.EM_BAND[0] <= hz[ch] <= HI.EM_BAND[1])
                fast += out
                fastbig += out and d[ch] >= 1.0
    gate["regret_anatomy_new"] = {"decisions": tot, "material": bad, "choice_out_of_EM_band": fast, "choice_out_of_band_and_delta_ge_1": fastbig,
                                  "share_out_of_band": fast / bad if bad else float("nan")}
    HI.write_json("H04_gate.json", gate)
    return gate


if __name__ == "__main__":
    import json

    print(json.dumps(run(), indent=1, default=str))
