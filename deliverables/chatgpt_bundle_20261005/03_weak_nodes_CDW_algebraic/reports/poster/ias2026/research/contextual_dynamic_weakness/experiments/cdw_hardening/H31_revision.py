# ruff: noqa: E501
"""H31 revision analyses requested by the internal reviewers. POST-HOC and EXPLORATORY: they use
existing data (plus cheap static PF/gSCR evaluations) and change no preregistered verdict.
Output: results/hardening/H31_revision.json (+ CSVs)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _hdata as HD
import H03_nested as H3

R = HI.RESULTS
V9 = C.V9
IDX = {u: k for k, u in enumerate(V9)}


def tau_curve(mg, pairs, pol):
    out = []
    stab_new = pol[(pol.split == "new") & pol.base_stable].pid
    stab_old = pol[(pol.split == "old") & pol.base_stable].pid
    for tau in (0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05):
        for split, pids in (("new", stab_new), ("old", stab_old)):
            nC = nCb = nD = nDb = 0
            for pid in pids:
                pc = pairs[(pairs.pid == pid) & (pairs.level == "C") & (pairs.mag >= tau)]
                pdd = pc[pc.lvD.fillna(False).astype(bool)]
                nC += len(pc) > 0
                nCb += (pc.dir == "s2d").any() and (pc.dir == "d2s").any()
                nD += len(pdd) > 0
                nDb += (pdd.dir == "s2d").any() and (pdd.dir == "d2s").any()
            out.append({"tau": tau, "split": split, "n": len(pids), "C": nC, "C_both": nCb, "D": nD, "D_both": nDb})
    return pd.DataFrame(out)


def min_regret_ranking(ctx, tau=C.TAU_MAT):
    """Exact regret-minimizing fixed ranking: DP over placed-prefix sets. Placing u after prefix P makes
    u the choice in every context S with P subset S and u not in S."""

    contexts = [(sum(1 << IDX[u] for u in V9 if u not in d), d) for d in ctx.values() if len(d) >= 2]  # mask of units IN S
    n = len(contexts)
    best = np.full(512, np.inf)
    best[0] = 0.0
    arg = {}
    for P in range(512):
        if not np.isfinite(best[P]):
            continue
        for u in V9:
            b = 1 << IDX[u]
            if P & b:
                continue
            cost = 0
            for smask, d in contexts:
                if (P & smask) == P and not (smask & b) and u in d:
                    cost += (d[u] - min(d.values())) >= tau
            Q = P | b
            if best[P] + cost < best[Q]:
                best[Q] = best[P] + cost
                arg[Q] = (P, u)
    order, Q = [], 511
    while Q:
        P, u = arg[Q]
        order.append(u)
        Q = P
    order.reverse()
    return best[511] / n if n else np.nan, order


def gscr_table():
    from ibr_cycles.diagnosis.baselines import generalized_scr, nodal_metrics

    net = C.load_network()
    fl = C.pf(net)
    nm = nodal_metrics(net, fl)
    g = {}
    for S in C.subsets(V9):
        g[sum(1 << IDX[u] for u in S)] = (generalized_scr(net, tuple(S), fl) if S else np.inf,
                                                       min(nm[b].scr for b in S) if S else np.inf)
    return g


def static_sequencing(mg, pol, gs):
    """Portfolio-conditioned static screens for the next replacement: choose the unit whose resulting
    portfolio has the largest gSCR (or largest minimum SCR). Regret vs the context-aware argmin."""

    rows = []
    for pid in pol[pol.base_stable].pid:
        g = mg[(mg.pid == pid) & mg.lvA]
        for strat, sel in (("FULL", g), ("EM", g[g.em])):
            for S, gg in sel.groupby("mask"):
                if len(gg) < 2:
                    continue
                d = dict(zip(gg.i, gg.delta, strict=True))
                best = min(d.values())
                for name, k in (("gscr", 0), ("minscr", 1)):
                    ch = max(d, key=lambda u: gs[S | (1 << IDX[u])][k])
                    rows.append({"pid": pid, "split": HI.split_of(pid), "stratum": strat, "screen": name, "mask": S,
                                 "material": (d[ch] - best) >= C.TAU_MAT})
    df = pd.DataFrame(rows)
    return df.groupby(["pid", "split", "stratum", "screen"]).material.mean().rename("frac_regret").reset_index()


def em_clean_squares(cen, pol):
    rows = []
    for pid in pol[pol.base_stable].pid:
        d = cen[pid]
        a, st, hz = d["alpha"], d["status"], d["hz"]
        pos = neg = tot = 0
        for S in range(512):
            for bi in range(9):
                for bj in range(bi + 1, 9):
                    if (S >> bi & 1) or (S >> bj & 1):
                        continue
                    four = (S, S | 1 << bi, S | 1 << bj, S | 1 << bi | 1 << bj)
                    if not all(st[x] == "STABLE" and HD.in_em(hz[x]) for x in four):
                        continue
                    dd = a[four[3]] - a[four[1]] - a[four[2]] + a[four[0]]
                    tot += 1
                    pos += dd > C.TAU_RES
                    neg += dd < -C.TAU_RES
        rows.append({"pid": pid, "split": HI.split_of(pid), "n_em_clean_squares": tot, "n_pos": pos, "n_neg": neg,
                     "both_signs": bool(pos and neg)})
    return pd.DataFrame(rows)


def fast_mode_facts(cen, pol):
    rows = []
    for pid in pol[(pol.split == "new") & pol.base_stable].pid:
        d = cen[pid]
        for m in range(512):
            rows.append({"pid": pid, "alpha": d["alpha"][m], "hz": d["hz"][m], "status": d["status"][m], "size": bin(m).count("1")})
    df = pd.DataFrame(rows)
    fast = df[df.alpha >= 1]
    return {"n_portfolios": int(len(df)), "n_fast": int(len(fast)), "frac_fast_real": float((fast.hz < 1e-3).mean()) if len(fast) else np.nan,
            "fast_alpha_quantiles": fast.alpha.quantile([0.05, 0.5, 0.95]).round(1).tolist() if len(fast) else [],
            "fast_size_counts": fast["size"].value_counts().sort_index().to_dict()}


def gz_cond(pol):
    import _infra as I

    rows = []
    for rec in I.Store("H_H01").all():
        pid = rec["task"]["pid"]
        for r in rec["records"]:
            rows.append({"pid": pid, "alpha": r.get("alpha"), "hz": r.get("lam_hz"), "gz_cond": r.get("gz_cond")})
    df = pd.DataFrame(rows)
    df["class"] = np.where(df.alpha >= 1, "fast (alpha>=1)", np.where(df.alpha > 0, "unstable EM/slow", "stable"))
    return df.groupby("class").gz_cond.describe(percentiles=[0.5, 0.95]).round(1).to_dict("index")


def dtot_vs_topology():
    """D_tot (H06, SPR, H4) vs finite doubling (gamma 2) and outage (gamma 0) effects from H14 at HARDENING_H01-H06."""

    lk = pd.read_parquet(R / "H06_links_raw.parquet")
    lk = lk[(lk.family == "link") & (lk.target == "H4") & lk.env.isna()]
    cl = pd.read_parquet(R / "H14_topology_em.parquet")
    h4 = cl[cl.S == "30+33+35+37"]
    rows = []
    for pid in HI.HPOL[:6]:
        g = lk[lk.pid == pid].set_index("idx")
        t = h4[h4.pid == pid]
        for kind, dg in (("dbl", 1.0), ("out", -1.0)):
            q = t[t.action.str.startswith(kind)].copy()
            q["e"] = q.action.str[3:].astype(int)
            q = q[np.isfinite(q.d_full)]
            pred = g.loc[q.e, "d_total"].to_numpy(float) * dg
            rows.append({"pid": pid, "action": kind, "n": int(len(q)), "rho_full": AN.spearman(pred, q.d_full.to_numpy(float)),
                         "rho_em": AN.spearman(pred[np.isfinite(q.d_em)], q.d_em.to_numpy(float)[np.isfinite(q.d_em)]),
                         "top5_full": AN.topk_precision(-pred, -q.d_full.to_numpy(float), 5)})
    return pd.DataFrame(rows)


def run():
    mg = pd.read_parquet(R / "H03_marginals.parquet")
    pairs = pd.read_parquet(R / "H03_nested_pairs.parquet")
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    out = {"label": "POST-HOC reviewer-requested analyses; exploratory; no preregistered verdict changes"}
    stable_new = set(pol[(pol.split == "new") & pol.base_stable].pid)
    pc = pairs[(pairs.level == "C") & pairs.pid.isin(stable_new)]
    pdd = pc[pc.lvD.fillna(False).astype(bool)]
    out["levelC_base_stable"] = {"pairs": int(len(pc)), "s2d": int((pc.dir == "s2d").sum()), "d2s": int((pc.dir == "d2s").sum()),
                                 "em_clean": float(pc.em_clean.mean()), "distinct_policy_unit": int(pc[["pid", "i"]].drop_duplicates().shape[0]),
                                 "distinct_S2_marginals": int(pc[["pid", "i", "S2"]].drop_duplicates().shape[0]), "m1_pairs": int((pc.m == 1).sum()),
                                 "max_mag": float(pc.mag.max())}
    out["levelD_base_stable"] = {"pairs": int(len(pdd)), "s2d": int((pdd.dir == "s2d").sum()), "d2s": int((pdd.dir == "d2s").sum()),
                                 "policies": int(pdd.pid.nunique()),
                                 "policies_both_dirs": int(sum(((pdd.pid == p) & (pdd.dir == "s2d")).any() and ((pdd.pid == p) & (pdd.dir == "d2s")).any() for p in stable_new)),
                                 "max_mag": float(pdd.mag.max()), "mag_quartiles": pdd.mag.quantile([0.25, 0.5, 0.75]).round(4).tolist()}
    zeta = pc.mag / (2 * np.pi * pc[["hz1", "hz2"]].min(axis=1))
    out["levelC_damping_ratio_quartiles_pct"] = (100 * zeta.quantile([0.25, 0.5, 0.75])).round(3).tolist()
    ex = pdd.sort_values("mag", ascending=False)
    exs, seen = [], set()
    for r in ex.itertuples():
        if r.i in seen:
            continue
        seen.add(r.i)
        exs.append({"pid": r.pid, "i": r.i, "S1": r.S1, "S2": r.S2, "d1": r.d1, "d2": r.d2, "hz1": r.hz1, "hz2": r.hz2, "mac12": r.mac12, "dir": r.dir})
        if len(exs) == 3:
            break
    out["levelD_examples"] = exs
    tc = tau_curve(mg, pairs, pol)
    tc.to_csv(R / "H31_tau_curve.csv", index=False)
    out["tau_curve"] = tc.to_dict("records")
    # regret-optimal fixed ranking
    import H04_ranking as H4

    rr = []
    for pid in sorted(pol[pol.base_stable].pid):
        g = mg[mg.pid == pid]
        for s in ("FULL", "EM", "SAME"):
            ctx = H4.contexts(H4.stratum(g, s))
            fr, order = min_regret_ranking(ctx)
            rr.append({"pid": pid, "split": HI.split_of(pid), "stratum": s, "min_regret_rate": fr, "order": "-".join(map(str, order))})
    rr = pd.DataFrame(rr)
    rr.to_csv(R / "H31_min_regret_ranking.csv", index=False)
    out["min_regret_ranking_median"] = rr.groupby(["split", "stratum"]).min_regret_rate.median().round(4).to_dict()
    out["min_regret_ranking_median"] = {f"{k[0]}|{k[1]}": v for k, v in out["min_regret_ranking_median"].items()}
    # portfolio-conditioned static sequencing screens
    gs = gscr_table()
    ss = static_sequencing(mg, pol, gs)
    ss.to_csv(R / "H31_static_sequencing.csv", index=False)
    out["static_sequencing_median_regret"] = {f"{a}|{b}|{c}": v for (a, b, c), v in ss.groupby(["split", "stratum", "screen"]).frac_regret.median().round(4).to_dict().items()}
    # EM-clean one-step squares
    cen = HD.load_census("H_H01", modes=False)
    cen.update(HD.load_census("E01", pids=set(HI.OLD_HOLD), modes=False))
    sq = em_clean_squares(cen, pol[pol.split.isin(["new", "old"])])
    sq.to_csv(R / "H31_em_clean_squares.csv", index=False)
    out["em_clean_squares"] = {s: {"policies": int((sq.split == s).sum()), "both_signs": int(sq[sq.split == s].both_signs.sum()),
                                   "median_squares": float(sq[sq.split == s].n_em_clean_squares.median())} for s in ("new", "old")}
    out["fast_modes_new"] = fast_mode_facts(cen, pol)
    out["gz_cond_by_class"] = gz_cond(pol)
    dt = dtot_vs_topology()
    dt.to_csv(R / "H31_dtot_vs_topology.csv", index=False)
    out["dtot_vs_topology_median"] = {k: {"rho_full": float(v.rho_full.median()), "rho_em": float(v.rho_em.median()), "top5": float(v.top5_full.median())}
                                      for k, v in dt.groupby("action")}
    h15 = json.loads((R / "H15_topology_gate.json").read_text())
    cr = pd.DataFrame(h15["em_creation"])
    out["topology_creation_by_split"] = cr.pid.map(HI.split_of).value_counts().to_dict()
    out["topology_creation_all_outages"] = bool(cr.action.str.startswith("out").all())
    HI.write_json("H31_revision.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
