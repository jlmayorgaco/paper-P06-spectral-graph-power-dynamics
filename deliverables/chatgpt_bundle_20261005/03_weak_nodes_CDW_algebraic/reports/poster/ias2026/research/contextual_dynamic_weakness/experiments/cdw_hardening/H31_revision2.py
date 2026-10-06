# ruff: noqa: E501
"""H31 round-2 revision analyses (POST-HOC, EXPLORATORY; existing data only). Extends
results/hardening/H31_revision.json. Requested by the internal reviewers (docs/CDW_REVIEWER{1,2}_FINAL.md)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _infra as I
import E01_census as E1

R = HI.RESULTS
V9 = C.V9


def screened_ranking(mg, pol):
    """Agreement-optimal (DP) order screened for feasibility: first ranked unit whose S+i is STABLE."""

    rk = pd.read_csv(R / "H04_ranking.csv")
    rk = rk[rk.stratum == "FULL"].set_index("pid")
    rows = []
    for pid in pol[pol.base_stable].pid:
        pos = {int(u): k for k, u in enumerate(rk.loc[pid, "order"].split("-"))}
        g = mg[(mg.pid == pid) & mg.lvA]
        n_all = n_nostable = 0
        reg_plain, reg_screen = [], []
        for _, gg in g.groupby("mask"):
            if len(gg) < 2:
                continue
            n_all += 1
            d = dict(zip(gg.i, gg.delta, strict=True))
            st = dict(zip(gg.i, gg.status_Si == "STABLE", strict=True))
            best = min(d.values())
            if not any(st.values()):
                n_nostable += 1
                continue
            plain = min(d, key=lambda u: pos[u])
            screened = min((u for u in d if st[u]), key=lambda u: pos[u])
            reg_plain.append(d[plain] - best >= C.TAU_MAT)
            reg_screen.append(d[screened] - best >= C.TAU_MAT)
        rows.append({"pid": pid, "split": HI.split_of(pid), "frac_no_stable_option": n_nostable / n_all if n_all else np.nan,
                     "regret_plain_on_feasible": float(np.mean(reg_plain)) if reg_plain else np.nan,
                     "regret_screened_on_feasible": float(np.mean(reg_screen)) if reg_screen else np.nan})
    return pd.DataFrame(rows)


def loco_branch_list():
    lk = pd.read_parquet(R / "H06_links_raw.parquet")
    lk = lk[(lk.family == "link") & (lk.target == "H4")].copy()
    lk["truth"] = -(lk["g1.50_alpha"] - lk.alpha0)
    lk["rank"] = lk.groupby("cond").truth.rank(ascending=False)
    rows = []
    for cond, g in lk.groupby("cond"):
        cl = g.cluster.iloc[0]
        other = lk[lk.cluster != cl].groupby("idx")["rank"].mean()
        g = g.set_index("idx").sort_index()
        fixed = -other.reindex(g.index).to_numpy(float)
        t = g.truth.to_numpy(float)
        rows.append({"cond": cond, "set": "fresh_draws" if isinstance(g.env.iloc[0], str) else "new", "rho_fixed_list": AN.spearman(fixed, t),
                     "rho_Dtotal": AN.spearman(-g.d_total.to_numpy(float), t), "rho_Dfrozen": AN.spearman(-g.d_frozen.to_numpy(float), t),
                     "top5_fixed_list": AN.topk_precision(fixed, t, 5), "top5_Dtotal": AN.topk_precision(-g.d_total.to_numpy(float), t, 5)})
    df = pd.DataFrame(rows)
    out = {}
    for s, q in (("all", df), ("new", df[df.set == "new"]), ("fresh_draws", df[df.set == "fresh_draws"])):
        out[s] = {k: float(q[k].median()) for k in ("rho_fixed_list", "rho_Dtotal", "rho_Dfrozen", "top5_fixed_list", "top5_Dtotal")}
        out[s]["frac_Dtotal_beats_fixed_list"] = float((q.rho_Dtotal > q.rho_fixed_list).mean())
    return df, out


def alt_em_tracked(case_dir=None):
    rows = [json.loads(p.read_text()) for p in (case_dir or (R / "alt" / "cases")).glob("*.json")]
    core = {(r["pid"], r["S"]): r for r in rows if r["net"] == "NOMINAL"}
    V4 = C.V4
    out = {}
    for pid in sorted({k[0] for k in core}):
        recs = {S: core[(p, S)] for (p, S) in core if p == pid}
        if recs.get("BASE", {}).get("status") != "STABLE":
            continue
        em = {}
        for lab, r in recs.items():
            ms = [E1.compact(m) for m in r.get("modes", []) if 0.1 <= m["hz"] <= 2.0]
            em[lab] = (max(ms, key=lambda m: m["re"]) if ms else None, ms, r["status"])
        marg = {}
        for S in C.subsets(V4):
            m0, _, st0 = em[C.label(S)]
            if st0 != "STABLE" or m0 is None:
                continue
            for i in V4:
                if i in S:
                    continue
                m1r, ms1, _ = em[C.label(S + (i,))]
                if m1r is None:
                    continue
                m1, _ = E1.match(m0, ms1)
                if m1 is None or abs(m1["re"] - m1r["re"]) > 1e-12:
                    continue
                marg[(sum(1 << V4.index(u) for u in S), i)] = m1["re"] - m0["re"]
        nested = 0
        for i in V4:
            ks = [(k[0], v) for k, v in marg.items() if k[1] == i]
            for a, va in ks:
                for b, vb in ks:
                    if a != b and (a & b) == a and ((va <= -C.TAU_MAT <= C.TAU_MAT <= vb) or (va >= C.TAU_MAT and vb <= -C.TAU_MAT)):
                        nested += 1
        stab = [r for r in recs.values() if r["status"] == "STABLE"]
        out[pid] = {"n_em_tracked_marginals": len(marg), "n_stabilizing": int(sum(v <= -C.TAU_MAT for v in marg.values())),
                    "n_destabilizing": int(sum(v >= C.TAU_MAT for v in marg.values())), "n_nested_reversals": int(nested),
                    "n_stable_portfolios": len(stab),
                    "frac_stable_alpha_pinned_at_minus_0.1": float(np.mean([abs(r["alpha"] + 0.1) < 1e-6 for r in stab])) if stab else np.nan}
    return out


def h18_uniform():
    d = pd.read_csv(R / "H18_crossmodel.csv")
    t = d[d.alt_base_status == "STABLE"]
    agree = float((t.E_agree * t.E_n).sum() / t.E_n.sum())

    def v3(x, hi, lo, strict_lo=False):
        return "TRANSFERS" if x >= hi else ("PARTIAL" if (x > lo if strict_lo else x >= lo) else "MODEL-SPECIFIC")
    return {"n_tested": int(len(t)), "C_median_kendall": float(t.C_kendall.median()), "C2_median_diff": float(t.C2_diff.median()),
            "C2_n_below_0.20": int((t.C2_diff < 0.20).sum()), "D_frac_top1_agree": float(t.D_top_agree.mean()), "E_pooled_sign_agree": agree,
            "E_n_pairs": int(t.E_n.sum()),
            "verdicts": {"C": v3(t.C_kendall.median(), 0.5, 0.2), "C2": v3(t.C2_diff.median(), 0.2, 0.0), "D": v3(t.D_top_agree.mean(), 0.75, 0.5), "E": v3(agree, 0.75, 0.5)}}


def gap2_sensitivity(pairs, pol):
    gap = {}
    for rec in I.Store("H_H01").all():
        for r in rec["records"]:
            gap[(rec["task"]["pid"], r["S"])] = r.get("gap2", np.nan)
    stable_new = list(pol[(pol.split == "new") & pol.base_stable].pid)
    pc = pairs[(pairs.level == "C") & pairs.pid.isin(stable_new)].copy()

    def add(lab, i):
        s = [] if lab == "BASE" else [int(u) for u in lab.split("+")]
        return C.label(tuple(sorted(s + [int(i)])))
    pc["gap4"] = [min(gap.get((r.pid, r.S1), np.nan), gap.get((r.pid, add(r.S1, r.i)), np.nan), gap.get((r.pid, r.S2), np.nan),
                      gap.get((r.pid, add(r.S2, r.i)), np.nan)) for r in pc.itertuples()]
    rob = pc[pc.gap4 >= C.TAU_MAT]
    return {"frac_pairs_gap_below_mag": float((pc.gap4 < pc.mag).mean()), "frac_pairs_gap_below_tau": float((pc.gap4 < C.TAU_MAT).mean()),
            "coverage_C_gap_ge_tau": f"{rob.pid.nunique()}/{len(stable_new)}",
            "coverage_D_gap_ge_tau": f"{rob[rob.lvD.fillna(False).astype(bool)].pid.nunique()}/{len(stable_new)}"}


def direction_counts(pairs, pol):
    stable_new = set(pol[(pol.split == "new") & pol.base_stable].pid)
    out = {}
    for lvl in ("A", "B", "C", "D"):
        q = pairs[(pairs.level == ("C" if lvl == "D" else lvl)) & pairs.pid.isin(stable_new)]
        if lvl == "D":
            q = q[q.lvD.fillna(False).astype(bool)]
        s2d, d2s = set(q[q.dir == "s2d"].pid), set(q[q.dir == "d2s"].pid)
        out[lvl] = {"s2d_policies": len(s2d), "d2s_policies": len(d2s), "both_policies": len(s2d & d2s), "n": len(stable_new)}
    return out


def mixing_terms():
    d = pd.read_csv(R / "H13_mixing_decomposition.csv")
    vt = (d.C + d.F).var()
    return {"var_C_over_var_total": float(d.C.var() / vt), "var_F_over_var_total": float(d.F.var() / vt),
            "2cov_over_var_total": float(2 * np.cov(d.C, d.F)[0, 1] / vt)}


def topology_within_type():
    cl = pd.read_parquet(R / "H14_topology_em.parquet")
    st = pd.read_csv(R / "H16_static_scores.csv", index_col=0)
    h4 = cl[(cl.S == "30+33+35+37") & cl.em_class.isin(["EM-STABILIZING", "EM-DESTABILIZING", "EM-NEUTRAL"]) & ~cl.fast_flag]
    out = {}
    for kind in ("dbl", "out"):
        for sc in ("fiedler", "kirchhoff", "gscr_H4", "min_scr_core", "losses"):
            rh = np.array([AN.spearman(st.loc[g.action, sc].to_numpy(float), g.d_em) for _, g in h4[h4.action.str.startswith(kind)].groupby("pid") if len(g) >= 10], float)
            out[f"{kind}|{sc}"] = {"median_rho": float(np.nanmedian(rh)), "frac_abs_ge_0.6": float(np.nanmean(np.abs(rh) >= 0.6))}
    return out


def order_distance():
    rk = pd.read_csv(R / "H04_ranking.csv")
    rk = rk[(rk.stratum == "FULL") & rk.base_stable]
    orders = [[int(u) for u in o.split("-")] for o in rk.order]
    taus = []
    for a in range(len(orders)):
        for b in range(a + 1, len(orders)):
            ra, rb = ({u: k for k, u in enumerate(o)} for o in (orders[a], orders[b]))
            taus.append(AN.kendall([ra[u] for u in V9], [rb[u] for u in V9]))
    return {"median_pairwise_kendall": float(np.median(taus)), "q25": float(np.percentile(taus, 25)), "q75": float(np.percentile(taus, 75))}


def run():
    mg = pd.read_parquet(R / "H03_marginals.parquet")
    pairs = pd.read_parquet(R / "H03_nested_pairs.parquet")
    pol = pd.read_csv(R / "H03_policy_summary.csv")
    out = json.loads((R / "H31_revision.json").read_text())
    sc = screened_ranking(mg, pol)
    sc.to_csv(R / "H31_screened_ranking.csv", index=False)
    out["screened_ranking_median"] = {s: {k: float(sc[sc.split == s][k].median()) for k in ("frac_no_stable_option", "regret_plain_on_feasible", "regret_screened_on_feasible")}
                                      for s in ("new", "old", "discovery")}
    lc, lo = loco_branch_list()
    lc.to_csv(R / "H31_loco_branch_list.csv", index=False)
    out["loco_fixed_branch_list"] = lo
    out["alt_em_tracked"] = alt_em_tracked()
    out["h18_uniform_tested_policies"] = h18_uniform()
    out["gap2_sensitivity"] = gap2_sensitivity(pairs, pol)
    out["direction_counts_new"] = direction_counts(pairs, pol)
    out["mixing_variance_terms"] = mixing_terms()
    out["topology_within_type"] = topology_within_type()
    out["optimal_order_distance"] = order_distance()
    HI.write_json("H31_revision.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
