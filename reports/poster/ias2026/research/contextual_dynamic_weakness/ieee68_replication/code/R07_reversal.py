# ruff: noqa: E501
"""R6 portfolio table + R7 primary replication (level D) with the preregistered secondary (level T) + R8 chain
identity, for Model A (REAL 16 policies; ablations NOGOV/SP33 on 6 policies) and, when present, Model B.

    python R07_reversal.py A     -> results/CDW68_PORTFOLIOS.parquet (A rows), CDW68_R07_A_*.{parquet,csv,json}
    python R07_reversal.py B     -> Model B rows appended, CDW68_R07_B_*
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import beta

import _rev68 as V

SEED_BOOT = 20260920


def boot_cov(flags, n_boot=10000, seed=SEED_BOOT):
    flags = np.asarray(flags, float)
    if not flags.size:
        return [np.nan, np.nan, np.nan]
    rng = np.random.default_rng(seed)
    bs = rng.choice(flags, size=(n_boot, flags.size), replace=True).mean(1)
    return [float(flags.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def clopper(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def flat_rows(model, store, cond_key):
    rows = []
    for r in store:
        t, rec = r["task"], r["record"]
        crit = next((m for m in rec.get("modes", []) if m.get("crit")), {})
        emt = next((m for m in rec.get("modes", []) if m.get("em_top")), {})
        rows.append({"model": model, "variant": t.get("variant", "REAL"), "cond": cond_key(t), "pid": t["pid"], "draw": t.get("draw"), "S": t["S"],
                     "size": 0 if t["S"] == "BASE" else len(t["S"].split("+")), "status": rec.get("status"), "alpha": rec.get("alpha"),
                     "lam_hz": rec.get("lam_hz"), "rhp": rec.get("rhp"), "gap2": rec.get("gap2"), "em_top_re": rec.get("em_top_re"),
                     "em_top_hz": rec.get("em_top_hz"), "crit_part_em": crit.get("part_em"), "emtop_part_em": emt.get("part_em"),
                     "eq_residual": rec.get("eq_residual", rec.get("init_residual")), "gz_cond": rec.get("gz_cond"), "n_x": rec.get("n_x"),
                     "min_abs_re": rec.get("min_abs_re")})
    return rows


def summarize(tag, cen, conds, out_prefix):
    mg, pr, pol = V.analyse({c: cen[c] for c in conds})
    mg.to_parquet(R.RESULTS / f"{out_prefix}_marginals.parquet", index=False)
    if len(pr):
        pr.to_parquet(R.RESULTS / f"{out_prefix}_pairs.parquet", index=False)
    pol.to_csv(R.RESULTS / f"{out_prefix}_policy_summary.csv", index=False)
    el = pol[pol.base_stable]
    n = int(len(el))
    g = {"tag": tag, "n_conditions": int(len(pol)), "n_eligible_base_stable": n, "n_base_unstable": int((~pol.base_stable).sum())}
    for lvl in ("A", "B", "C", "D", "T"):
        k = int(el[f"{lvl}_has"].sum()) if n else 0
        g[f"cov_{lvl}"] = k / n if n else np.nan
        g[f"k_{lvl}"] = k
        g[f"cov_{lvl}_both_dirs"] = float(el[f"{lvl}_both_dirs"].mean()) if n else np.nan
        g[f"cov_{lvl}_boot"] = boot_cov(el[f"{lvl}_has"].to_numpy()) if n else None
        g[f"cov_{lvl}_clopper"] = clopper(k, n) if n else None
    if len(pr):
        for lvl, sel in (("D", pr.lvD), ("T", pr.lvT2), ("C", pr.level == "C"), ("A", pr.level == "A")):
            q = pr[sel & pr.cond.isin(el.cond)]
            g[f"{lvl}_pairs"] = int(len(q))
            g[f"{lvl}_dir_counts"] = q.dir.value_counts().to_dict() if len(q) else {}
            g[f"{lvl}_mag_quartiles"] = q.mag.quantile([0.25, 0.5, 0.75]).tolist() if len(q) else []
            g[f"{lvl}_mag_max"] = float(q.mag.max()) if len(q) else np.nan
            g[f"{lvl}_distinct_cond_unit"] = int(q[["cond", "i"]].drop_duplicates().shape[0]) if len(q) else 0
            g[f"{lvl}_gap_robust_cov"] = (int(q[q.gap_min >= R.TAU_MAT].cond.nunique()) / n) if (len(q) and n) else 0.0
            if lvl in ("D", "T") and len(q):
                hz = q.hz1.to_numpy(float)
                g[f"{lvl}_zeta_quartiles_pct"] = (100 * q.mag / (2 * np.pi * np.where(hz > 0, hz, np.nan))).quantile([0.25, 0.5, 0.75]).tolist()
        allrev = pr[pr.cond.isin(el.cond)]
        g["chain_identity_max_abs_resid"] = float(np.nanmax(np.abs(allrev.chain_resid))) if len(allrev) else np.nan
        g["chain_identity_n_checked"] = int(np.isfinite(allrev.chain_resid).sum()) if len(allrev) else 0
    g["primary_gate_pass"] = bool(n and g["cov_D"] >= 0.5 and (g.get("D_mag_quartiles") or [0, 0])[1] >= R.TAU_MAT)
    g["frac_stable_crit_em_median"] = float(pol.frac_crit_em.median())
    g["base_crit_hz_median"] = float(pol.base_hz.median())
    g["kappa_values"] = pol.kappa_RHP.value_counts(dropna=False).to_dict()
    return g, pol, pr


def census_A():
    st = list(R.Store("R06A").all())
    assert all(r["ok"] for r in st)
    recs = [(f"{r['task']['variant']}|{r['task']['pid']}", r["task"]["S"], r["record"]) for r in st]
    return V.census_from_records(recs), st


def main(mode):
    if mode == "A":
        cen, st = census_A()
        rows = flat_rows("A", st, lambda t: f"{t['variant']}|{t['pid']}")
        pd.DataFrame(rows).to_parquet(R.RESULTS / "CDW68_PORTFOLIOS_A.parquet", index=False)
        design = json.loads((R.INPUTS / "cdw68_design_v1.json").read_text())
        real = [f"REAL|{p['id']}" for p in design["policies_A"]]
        out = {}
        gA, polA, prA = summarize("A-REAL", cen, real, "CDW68_R07_A_REAL")
        out["A_REAL"] = gA
        curve = V.coverage_curve(cen, real)
        curve.to_csv(R.RESULTS / "CDW68_R07_A_REAL_tau_curve.csv", index=False)
        out["A_REAL"]["tau_curve"] = curve.to_dict("records")
        for v in ("REAL", "NOGOV", "SP33"):
            conds = [f"{v}|{p}" for p in design["ablation_policies"]]
            g, _, _ = summarize(f"A-{v}-ablation", cen, conds, f"CDW68_R13_{v}")
            out[f"ablation_{v}"] = g
        R.write_json(R.RESULTS / "CDW68_R07_A_gate.json", out)
        print(json.dumps({k: {kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict)) or kk in ("D_mag_quartiles", "T_mag_quartiles", "D_dir_counts", "T_dir_counts", "cov_D_clopper", "cov_T_clopper")} for k, v in out.items()}, indent=1, default=str))
    elif mode == "B":
        st = [r for r in (json.loads(p.read_text()) for p in sorted((R.RAW / "B68C").glob("*.json")))]
        recs = [(f"REAL|{r['task']['pid']}", r["task"]["S"], r["record"]) for r in st]
        cen = V.census_from_records(recs)
        rows = flat_rows("B", st, lambda t: f"REAL|{t['pid']}")
        pd.DataFrame(rows).to_parquet(R.RESULTS / "CDW68_PORTFOLIOS_B.parquet", index=False)
        conds = sorted(cen)
        g, pol, pr = summarize("B-REAL", cen, conds, "CDW68_R07_B_REAL")
        curve = V.coverage_curve(cen, conds)
        curve.to_csv(R.RESULTS / "CDW68_R07_B_REAL_tau_curve.csv", index=False)
        g["tau_curve"] = curve.to_dict("records")
        g["status_counts"] = pd.Series([r["record"].get("status") for r in st]).value_counts().to_dict()
        R.write_json(R.RESULTS / "CDW68_R07_B_gate.json", {"B_REAL": g})
        print(json.dumps({kk: vv for kk, vv in g.items() if not isinstance(vv, (list, dict)) or kk in ("D_mag_quartiles", "T_mag_quartiles", "D_dir_counts", "T_dir_counts", "status_counts")}, indent=1, default=str))
    # combined portfolio table
    parts = [pd.read_parquet(p) for p in (R.RESULTS / "CDW68_PORTFOLIOS_A.parquet", R.RESULTS / "CDW68_PORTFOLIOS_B.parquet") if p.exists()]
    pd.concat(parts, ignore_index=True).to_parquet(R.RESULTS / "CDW68_PORTFOLIOS.parquet", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
