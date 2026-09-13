# ruff: noqa: E501
"""R9 fixed node ranking on IEEE-68 (docs/CDW68_PREREG_V1.md section 10), Models A and B (REAL).

Port of H04: resolved preference = material difference (Delta_a < Delta_b - tau); p* = agreement of the best fixed
order (all 720 orders of the six units enumerated); material regret of the p*-order; regret-optimal order;
stability-screened order; strata FULL / EM / LEVELD (level-C marginals) / TRACKED (level-T objective); transfer.
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import itertools
import json

import numpy as np
import pandas as pd
from scipy.stats import binom

V = R.V68
IDX = {u: k for k, u in enumerate(V)}
ORDERS = list(itertools.permutations(V))
TAU = R.TAU_MAT


def stratum(mg, s):
    if s == "FULL":
        g = mg[mg.lvA]
        return g.assign(obj=g.delta)
    if s == "EM":
        g = mg[mg.lvA & mg.em]
        return g.assign(obj=g.delta)
    if s == "LEVELD":
        g = mg[mg.lvC]
        return g.assign(obj=g.delta)
    g = mg[mg.lvT]
    return g.assign(obj=g.t_delta)


def contexts(g):
    return {S: dict(zip(gg.i, gg.obj, strict=True)) for S, gg in g.groupby("mask")}


def wmatrix(ctx):
    W = np.zeros((len(V), len(V)))
    for d in ctx.values():
        for a in d:
            for b in d:
                if a != b and d[a] < d[b] - TAU:
                    W[IDX[a], IDX[b]] += 1
    return W


def accuracy(order, W):
    tot = W.sum()
    if tot == 0:
        return np.nan
    pos = {u: r for r, u in enumerate(order)}
    return sum(W[IDX[a], IDX[b]] for a in V for b in V if a != b and pos[a] < pos[b]) / tot


def regret(order, ctx, stable=None):
    pos = {u: r for r, u in enumerate(order)}
    regs = []
    for S, d in ctx.items():
        if len(d) < 2:
            continue
        cand = sorted(d, key=lambda u: pos[u])
        choice = cand[0]
        if stable is not None:
            ok = [u for u in cand if stable.get((S, u), False)]
            if not ok:
                continue
            choice = ok[0]
        regs.append(d[choice] - min(d.values()))
    regs = np.array(regs)
    if not regs.size:
        return np.nan, 0
    return float((regs >= TAU).mean()), int(regs.size)


def analyse(tag, mg, pol):
    rows, orders, ctxs = [], {}, {}
    el = set(pol[pol.base_stable].cond)
    for cond, g in mg.groupby("cond"):
        stab = {(r.mask, r.i): r.status_Si == "STABLE" for r in g.itertuples()}
        for s in ("FULL", "EM", "LEVELD", "TRACKED"):
            ctx = contexts(stratum(g, s))
            W = wmatrix(ctx)
            accs = [accuracy(o, W) for o in ORDERS]
            if np.all(np.isnan(accs)):
                best = ORDERS[0]
                pstar = np.nan
            else:
                k = int(np.nanargmax(accs))
                best, pstar = ORDERS[k], float(accs[k])
            fr, n = regret(best, ctx)
            regs = [regret(o, ctx)[0] for o in ORDERS]
            ropt = float(np.nanmin(regs)) if not np.all(np.isnan(regs)) else np.nan
            scr, n_scr = regret(best, ctx, stable=stab)
            rows.append({"model": tag, "cond": cond, "eligible": cond in el, "stratum": s, "n_decisions": n, "n_pref": int(W.sum()),
                         "p_star": pstar, "frac_regret": fr, "regret_optimal": ropt, "screened_regret": scr, "order": "-".join(map(str, best))})
            if s in ("FULL", "TRACKED"):
                orders[(cond, s)], ctxs[(cond, s)] = best, ctx
    rk = pd.DataFrame(rows)
    tr = []
    for s in ("FULL", "TRACKED"):
        conds = sorted(c for c in el)
        for a in conds:
            for b in conds:
                tr.append({"model": tag, "stratum": s, "fit": a, "use": b, "regret": regret(orders[(a, s)], ctxs[(b, s)])[0]})
    return rk, pd.DataFrame(tr)


def main():
    out, rks, trs = {}, [], []
    for tag, pref in (("A", "CDW68_R07_A_REAL"), ("B", "CDW68_R07_B_REAL")):
        p = R.RESULTS / f"{pref}_marginals.parquet"
        if not p.exists():
            continue
        mg = pd.read_parquet(p)
        pol = pd.read_csv(R.RESULTS / f"{pref}_policy_summary.csv")
        rk, tr = analyse(tag, mg, pol)
        rks.append(rk)
        trs.append(tr)
        g = {}
        for s in ("FULL", "EM", "LEVELD", "TRACKED"):
            q = rk[(rk.stratum == s) & rk.eligible]
            med_n = float(q.n_decisions.median()) if len(q) else 0.0
            g[s] = {"median_n_decisions": med_n, "median_n_material_pref": float(q.n_pref.median()) if len(q) else 0.0,
                    "median_p_star": float(q.p_star.median()) if q.p_star.notna().any() else None,
                    "median_frac_regret": float(q.frac_regret.median()) if q.frac_regret.notna().any() else None,
                    "median_regret_optimal": float(q.regret_optimal.median()) if q.regret_optimal.notna().any() else None,
                    "median_screened_regret": float(q.screened_regret.median()) if q.screened_regret.notna().any() else None,
                    "evaluable": bool(med_n >= 20), "n_policies": int(len(q)),
                    "insufficient_bar_met": bool(q.p_star.notna().any() and q.p_star.median() <= 0.90 and q.frac_regret.median() >= 0.10)}
            t = tr[tr.stratum == s] if s in ("FULL", "TRACKED") else None
            if t is not None and len(t):
                diag = t[t.fit == t.use].regret
                off = t[t.fit != t.use].regret
                g[s]["transfer_diag_median_regret"] = float(diag.median())
                g[s]["transfer_offdiag_median_regret"] = float(off.median())
        q = rk[(rk.stratum == "FULL") & rk.eligible & rk.p_star.notna()]
        kk, nn = int((q.p_star < 0.90).sum()), int(len(q))
        g["T1_sign_test"] = {"k_below_090": kk, "n": nn, "p_one_sided": float(binom.sf(kk - 1, nn, 0.5)) if nn else None}
        g["headline_allowed"] = bool(g["FULL"]["insufficient_bar_met"] and ((g["LEVELD"]["evaluable"] and g["LEVELD"]["insufficient_bar_met"]) or (g["TRACKED"]["evaluable"] and g["TRACKED"]["insufficient_bar_met"])))
        out[tag] = g
    pd.concat(rks).to_csv(R.RESULTS / "CDW68_R09_ranking.csv", index=False)
    pd.concat(trs).to_csv(R.RESULTS / "CDW68_R09_transfer.csv", index=False)
    R.write_json(R.RESULTS / "CDW68_R09_gate.json", out)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
