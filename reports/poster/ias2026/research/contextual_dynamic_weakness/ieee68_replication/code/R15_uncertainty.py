# ruff: noqa: E501
"""R15 uncertainty holdout (deterministic stress envelopes E05 = +-5 %, E10 = +-10 %; not probabilities).

Per draw: all 64 portfolios at the reference policy (P68_10 / B68_10); level-D and level-T reversal presence,
strongest witness, share of materially stabilizing EM-tracked marginals; Model A branch-ranking Spearman on the
preregistered subset (first 5 draws of each envelope). Writes results/CDW68_R15_uncertainty.json and per-draw CSV.
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import _rev68 as V


def per_draw(records, model):
    by = {}
    for draw, S, rec in records:
        by.setdefault(draw, []).append((draw, S, rec))
    rows = []
    for draw, recs in sorted(by.items()):
        cen = V.census_from_records(recs)
        mg, pr, pol = V.analyse(cen)
        env = draw.split("_")[1]
        t = mg[mg.lvT & np.isfinite(mg.t_delta)]
        wit = None
        for lvl, sel in (("D", "lvD"), ("T", "lvT2")):
            if len(pr) and pr[sel].any():
                q = pr[pr[sel]].sort_values("mag", ascending=False).iloc[0]
                wit = wit or f"{lvl}:{q.i}|{R.label(R.members_of(int(q.S1)))}|{R.label(R.members_of(int(q.S2)))}"
        rows.append({"model": model, "draw": draw, "envelope": env, "n_portfolios": int(len(recs)), "n_stable": int(pol.n_stable.iloc[0]),
                     "base_stable": bool(pol.base_stable.iloc[0]), "alpha_max": float(np.nanmax(cen[draw]["alpha"])), "alpha_min": float(np.nanmin(cen[draw]["alpha"])),
                     "crit_em_frac": float(pol.frac_crit_em.iloc[0]) if np.isfinite(pol.frac_crit_em.iloc[0]) else 0.0,
                     "D_has": bool(pol.D_has.iloc[0]), "T_has": bool(pol.T_has.iloc[0]), "A_has": bool(pol.A_has.iloc[0]),
                     "T_n_stab_material": int((t.t_delta <= -R.TAU_MAT).sum()), "T_n_destab_material": int((t.t_delta >= R.TAU_MAT).sum()),
                     "T_min": float(t.t_delta.min()) if len(t) else np.nan, "T_max": float(t.t_delta.max()) if len(t) else np.nan,
                     "global_max_abs_delta": float(mg[mg.lvA].delta.abs().max()) if mg.lvA.any() else np.nan, "strongest_witness": wit})
    return pd.DataFrame(rows)


def ranking_subset():
    rows = []
    for r in R.Store("R15L").all():
        if not r["ok"]:
            continue
        t, rec = r["task"], r["record"]
        d = rec["der"]
        for it in d["items"]:
            f = rec["finite"].get(str(it["e"]), {}).get("1.5", {})
            if f.get("ok"):
                rows.append({"draw": t["draw"], "e": it["e"], "P_tot": -it["d_total"], "T": -(f["alpha"] - d["lam"][0]),
                             "P_em": -it.get("em_d_total", np.nan), "Tem": -(f.get("em_tracked_re", np.nan) - (d["lam_em"] or [np.nan])[0])})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    out = []
    for draw, g in df.groupby("draw"):
        out.append({"draw": draw, "rho_tot": float(spearmanr(g.P_tot, g["T"]).statistic), "rho_em_posthoc": float(spearmanr(g.P_em, g.Tem, nan_policy="omit").statistic),
                    "truth_abs_max": float(g["T"].abs().max())})
    return pd.DataFrame(out)


def main():
    res = {"label": "deterministic stress envelopes E05 (+-5 %) and E10 (+-10 %); coverages are fractions of designed draws, not probabilities"}
    parts = []
    a = [(r["task"]["draw"], r["task"]["S"], r["record"]) for r in R.Store("R15A").all() if r["ok"]]
    if a:
        parts.append(per_draw(a, "A"))
    bdir = R.RAW / "B68D"
    if bdir.exists():
        b = [(r["task"]["draw"], r["task"]["S"], r["record"]) for r in (json.loads(p.read_text()) for p in sorted(bdir.glob("*.json")))]
        if b:
            parts.append(per_draw(b, "B"))
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    df.to_csv(R.RESULTS / "CDW68_R15_draws.csv", index=False)
    for (model, env), g in df.groupby(["model", "envelope"]):
        res[f"{model}_{env}"] = {"n_draws": int(len(g)), "frac_draws_level_D": float(g.D_has.mean()), "frac_draws_level_T": float(g.T_has.mean()),
                                 "frac_draws_level_A": float(g.A_has.mean()), "frac_draws_any_stabilizing_EM_material": float((g.T_n_stab_material > 0).mean()),
                                 "median_T_destab_material": float(g.T_n_destab_material.median()), "T_min_over_draws": float(g.T_min.min()),
                                 "all_portfolios_stable_frac": float((g.n_stable == g.n_portfolios).mean()), "alpha_max_over_draws": float(g.alpha_max.max()),
                                 "witnesses": g.strongest_witness.dropna().value_counts().to_dict()}
    rk = ranking_subset() if (R.RAW / "R15L").exists() else pd.DataFrame()
    if len(rk):
        rk.to_csv(R.RESULTS / "CDW68_R15_ranking_draws.csv", index=False)
        res["A_ranking_draws"] = {"n": int(len(rk)), "rho_tot_quantiles": rk.rho_tot.quantile([0.1, 0.5, 0.9]).tolist(),
                                  "rho_em_posthoc_quantiles": rk.rho_em_posthoc.quantile([0.1, 0.5, 0.9]).tolist(),
                                  "truth_abs_max_median": float(rk.truth_abs_max.median())}
    res["witness_identity_note"] = "the nominal reference policy has no level-D or level-T witness, so witness identity is undefined unless a draw creates one"
    R.write_json(R.RESULTS / "CDW68_R15_uncertainty.json", res)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
