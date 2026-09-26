# ruff: noqa: E501
"""Descriptive analytics for the report (no new model evaluation): reads result tables only.

Writes results/CDW_REPORT_ANALYTICS.json. Every number quoted in the report that is not a gate
value comes from here.
"""

from __future__ import annotations

import json

import _infra as I
import numpy as np
import pandas as pd

import _cdw as C

R = I.RESULTS


def e1_analytics():
    mg = pd.read_parquet(R / "CDW_E1_contextual_marginals.parquet")
    pol = pd.read_csv(R / "CDW_E1_policy_summary.csv")
    p = pd.read_parquet(R / "CDW_E1_portfolios.parquet")
    hs = set(pol[(pol.split == "holdout") & pol.base_stable].pid)
    out = {}
    # alpha heavy tail: fast instabilities
    out["portfolio_status_counts"] = p.status.value_counts().to_dict()
    out["frac_unstable_alpha_gt_1"] = float(((p.alpha > 1) & (p.status == "UNSTABLE")).sum() / max((p.status == "UNSTABLE").sum(), 1))
    out["alpha_max"] = float(p.alpha.max())
    st = mg[mg.stable_S & mg.pid.isin(hs)]
    out["n_stable_context_marginals_holdout"] = int(len(st))
    out["frac_same_mode_holdout"] = float(st.same_mode.mean())
    # electromechanical-only reversal: both S and S+i stay slow (|alpha(S+i)| < 1) and crit in 0.1-2 Hz
    em = st[(st.alpha_Si.abs() < 1.0) & st.crit_hz_Si.between(0.1, 2.0)]
    rows = []
    for (pid, i), g in em.groupby(["pid", "i"]):
        neg, pos = g[g.delta <= -C.TAU_MAT], g[g.delta >= C.TAU_MAT]
        if len(neg) and len(pos):
            rows.append({"pid": pid, "i": i, "strength": float(min(-neg.delta.min(), pos.delta.max()))})
    er = pd.DataFrame(rows)
    out["em_only_reversal_policy_frac"] = float(er.pid.nunique() / len(hs)) if len(er) else 0.0
    out["em_only_reversal_pairs"] = int(len(er))
    out["em_only_strength_median"] = float(er.strength.median()) if len(er) else np.nan
    # tracked reversal strengths
    t = st[st.same_mode]
    rows = []
    for (pid, i), g in t.groupby(["pid", "i"]):
        neg, pos = g[g.delta <= -C.TAU_MAT], g[g.delta >= C.TAU_MAT]
        if len(neg) and len(pos):
            rows.append({"pid": pid, "i": int(i), "strength": float(min(-neg.delta.min(), pos.delta.max())),
                         "S_stab": neg.loc[neg.delta.idxmin()].S, "d_stab": float(neg.delta.min()),
                         "S_destab": pos.loc[pos.delta.idxmax()].S, "d_destab": float(pos.delta.max())})
    tr = pd.DataFrame(rows)
    tr.to_csv(R / "CDW_E1_tracked_reversals_holdout.csv", index=False)
    out["tracked_reversal_pairs_holdout"] = int(len(tr))
    out["tracked_strength_quantiles"] = tr.strength.quantile([0.25, 0.5, 0.75, 1.0]).round(4).to_dict() if len(tr) else {}
    out["tracked_reversal_per_unit"] = tr.groupby("i").pid.nunique().to_dict() if len(tr) else {}
    # per-unit sign classes in holdout stable contexts at the material threshold
    s = pd.read_csv(R / "CDW_E1_summary.csv")
    sh = s[s.pid.isin(hs)]
    out["per_unit_median_fracs"] = sh.groupby("i")[["frac_stab", "frac_neutral", "frac_destab"]].median().round(3).to_dict()
    out["node_only_pstar_holdout"] = pol[pol.pid.isin(hs)].p_star.describe().round(3).to_dict()
    out["node_only_regret_holdout"] = pol[pol.pid.isin(hs)].frac_regret.describe().round(3).to_dict()
    # tau_res sensitivity
    sr = pd.read_csv(R / "CDW_E1_summary_tau_res.csv")
    out["tau_res_rev_stable_policy_frac_holdout"] = float(sr[sr.pid.isin(hs)].groupby("pid").rev_stable.any().mean())
    return out


def main():
    out = {"E1": e1_analytics()}
    I.atomic_write_json(R / "CDW_REPORT_ANALYTICS.json", out)
    print(json.dumps(out, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()
