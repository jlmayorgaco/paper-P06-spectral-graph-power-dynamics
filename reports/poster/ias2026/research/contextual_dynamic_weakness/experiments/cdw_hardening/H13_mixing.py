# ruff: noqa: E501
"""H13 two-factor mixing decomposition and fixed-frequency correlations (prereg H12/H13).

C = 0.5[(mu10-mu00)+(mu11-mu01)], F = 0.5[(mu01-mu00)+(mu11-mu10)] (decomposition only)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _hstats as ST
import _infra as I


def run():
    rows = {r["task"]["pid"]: r for r in I.Store("H_H12").all() if r.get("ok")}
    b = [r for r in I.Store("H_H12b").all() if r.get("ok")][0]
    mu00 = b["mu00"]
    df = pd.DataFrame([{"pid": p, "split": HI.split_of(p), "status_H4": r["status"], "alpha_H4": r["alpha"], "wc_hz": r["wc_hz"],
                        "mu10": r["mu_ref"], "mu_060": r["mu_060"], "mu_070": r["mu_070"], "mu_080": r["mu_080"],
                        "mu11": r["mu_wc"], "mu01": b["mu_ref_at_wc"][p]} for p, r in rows.items()])
    df["mu00"] = mu00
    df["C"] = 0.5 * ((df.mu10 - df.mu00) + (df.mu11 - df.mu01))
    df["F"] = 0.5 * ((df.mu01 - df.mu00) + (df.mu11 - df.mu10))
    df["check"] = (df.C + df.F) - (df.mu11 - df.mu00)
    pol = pd.read_csv(HI.RESULTS / "H03_policy_summary.csv").set_index("pid")
    rk = pd.read_csv(HI.RESULTS / "H04_ranking.csv")
    rk = rk[rk.stratum == "FULL"].set_index("pid")
    df["base_stable"] = df.pid.map(pol.base_stable)
    df["n_rev"] = df.pid.map(pol.n_rev_arbitrary)
    df["n_nested"] = df.pid.map(pol.A_n_units)
    df["regret"] = df.pid.map(rk.frac_regret)
    df["wc_in_EM"] = df.wc_hz.between(*HI.EM_BAND)
    df.to_csv(HI.RESULTS / "H13_mixing_decomposition.csv", index=False)
    out = {"max_abs_decomposition_check": float(df.check.abs().max()), "mu00": mu00}
    tot = (df.mu11 - df.mu00)
    out["share_abs_C_of_abs_C_plus_abs_F_median"] = float(np.median(df.C.abs() / (df.C.abs() + df.F.abs()).clip(lower=1e-15)))
    out["var_share_C"] = float(df.C.var() / (df.C.var() + df.F.var()))
    out["cov_C_F"] = float(np.cov(df.C, df.F)[0, 1])
    out["var_total_change"] = float(tot.var())
    out["rel_range_mu10"] = float((df.mu10.max() - df.mu10.min()) / df.mu10.median())
    out["rel_range_mu11"] = float((df.mu11.max() - df.mu11.min()) / df.mu11.median())
    out["rel_range_mu01"] = float((df.mu01.max() - df.mu01.min()) / df.mu01.median())
    for s in ("new", "old", "discovery"):
        q = df[df.split == s]
        out[f"{s}_median_C"], out[f"{s}_median_F"] = float(q.C.median()), float(q.F.median())
        out[f"{s}_IQR_C"], out[f"{s}_IQR_F"] = q.C.quantile([0.25, 0.75]).tolist(), q.F.quantile([0.25, 0.75]).tolist()
    # ---------------------------------------------------------- correlations ----
    tests = {}
    for split in ("new", "old"):
        q = df[(df.split == split) & df.base_stable]
        for pred in ("mu10", "mu_060", "mu_070", "mu_080", "C", "mu11"):
            for tgt in ("n_rev", "n_nested", "regret"):
                rho, p, n = ST.perm_spearman(q[pred], q[tgt])
                tests[f"{split}|{pred}|{tgt}"] = {"rho": rho, "p": p, "n": n}
    prim = tests["new|mu10|n_rev"]
    out["T3_primary"] = prim
    sec_keys = [k for k in tests if k != "new|mu10|n_rev" and not k.split("|")[1] == "mu11"]
    hs = ST.holm({k: tests[k]["p"] for k in sec_keys})
    for k in sec_keys:
        tests[k]["p_holm_Fmix"] = hs[k]
    out["tests"] = tests
    HI.write_json("H13_mixing.json", out)
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
