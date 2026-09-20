"""Derive preregistered TX4 robustness summaries from the master table."""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.contingency_tables import mcnemar

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
H4 = "30+33+35+37"
PARAMS = ("g", "k", "t", "h", "epsilon", "damping", "inertia")
BOUNDS = {
    "g": (0.020, 0.250), "k": (0.750, 1.750), "t": (0.750, 1.750),
    "h": (0.500, 1.500), "epsilon": (0.900, 1.100),
    "damping": (0.000, 0.200), "inertia": (0.800, 1.200),
}
THRESHOLD = 1e-8


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return float("nan"), float("nan")
    p = success / total
    den = 1.0 + z * z / total
    centre = (p + z * z / (2 * total)) / den
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
    return centre - half, centre + half


def read_master() -> pd.DataFrame:
    path = RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv"
    if not path.exists():
        path = RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz"
    df = pd.read_csv(path)
    for col in (*PARAMS, "alpha_all", "alpha_EM", "critical_frequency_hz", "ell_H4", "q_H4", "delta_H4", "eta_H4"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ("H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE", "stable_all", "stable_EM"):
        df[col] = df[col].astype(str).str.lower().eq("true")
    return df


def condition_table(h4: pd.DataFrame) -> pd.DataFrame:
    cols = ["campaign", "condition_id", "seed", *PARAMS, "alpha_all", "alpha_EM", "critical_frequency_hz", "ell_H4", "q_H4", "delta_H4", "eta_H4", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE", "evaluator_tier"]
    out = h4[cols].drop_duplicates("condition_id").copy()
    out["stable_full"] = out["alpha_all"] <= THRESHOLD
    out["stable_mode"] = out["alpha_EM"].notna() & (out["alpha_EM"] <= THRESHOLD)
    return out


def write_table(df: pd.DataFrame, name: str) -> None:
    df.to_csv(RESULTS / name, index=False)


def main() -> None:
    df = read_master()
    h4 = df[df["portfolio"] == H4].copy()
    cond = condition_table(h4)
    # Required master formats.
    df.to_parquet(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet", index=False)
    df.to_csv(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz", index=False, compression="gzip")
    df[df["campaign"] == "1D"].to_csv(RESULTS / "TX4_1D_SWEEPS.csv", index=False)
    h4[df["campaign"] == "2D" if False else h4["campaign"] == "2D"].to_parquet(RESULTS / "TX4_2D_PHASE_MAPS.parquet", index=False)

    summary_rows = []
    for campaign in ("QMC", "MC"):
        sub = cond[cond["campaign"] == campaign]
        for endpoint in ("H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"):
            success = int(sub[endpoint].sum())
            lo, hi = wilson(success, len(sub))
            summary_rows.append({
                "campaign": campaign, "endpoint": endpoint, "n_conditions": len(sub),
                "success": success, "proportion": success / len(sub),
                "ci95_low": lo, "ci95_high": hi,
                "interpretation": "bounded-box engineering probability; proper subsets outside calibration are surrogate-tier",
            })
        summary_rows.append({
            "campaign": campaign, "endpoint": "delta_H4_mean", "n_conditions": len(sub),
            "success": "", "proportion": float(sub["delta_H4"].mean()),
            "ci95_low": float(sub["delta_H4"].quantile(0.025)), "ci95_high": float(sub["delta_H4"].quantile(0.975)),
            "interpretation": "empirical 2.5/97.5 percentiles",
        })
    write_table(pd.DataFrame(summary_rows), "TX4_QMC_SUMMARY.csv")
    write_table(pd.DataFrame(summary_rows), "TX4_MC_SUMMARY.csv")

    # Replace the duplicate summary files with the campaign-specific slices.
    qmc_summary = pd.DataFrame([r for r in summary_rows if r["campaign"] == "QMC"])
    mc_summary = pd.DataFrame([r for r in summary_rows if r["campaign"] == "MC"])
    qmc_summary.to_csv(RESULTS / "TX4_QMC_SUMMARY.csv", index=False)
    mc_summary.to_csv(RESULTS / "TX4_MC_SUMMARY.csv", index=False)

    # Morris-style bounded screening from 40 preregistered bootstrap resamples.
    mc = cond[cond["campaign"] == "MC"].reset_index(drop=True)
    morris_rows = []
    rng = np.random.default_rng(20260922)
    for endpoint in ("H4_PRESENT", "EXACT_H4"):
        y = mc[endpoint].astype(float).to_numpy()
        for parameter in PARAMS:
            effects = []
            for _ in range(40):
                idx = rng.choice(len(mc), size=min(128, len(mc)), replace=False)
                x = mc.iloc[idx][parameter].to_numpy()
                median = float(np.median(x))
                effects.append(abs(float(y[idx][x >= median].mean()) - float(y[idx][x < median].mean())))
            morris_rows.append({"endpoint": endpoint, "parameter": parameter, "mu_star": float(np.mean(effects)), "sigma": float(np.std(effects, ddof=1)), "trajectories": 40, "status": "APPROX_SCREENING_FROM_MC", "deviation": "not a full Morris trajectory simulation"})
    write_table(pd.DataFrame(morris_rows), "TX4_MORRIS.csv")

    # Saltelli-form estimator driven by a deterministic calibrated surrogate.
    calibration = h4[h4["evaluator_tier"] == "EXACT_16_DAE"].dropna(subset=["alpha_all"])
    xcal = calibration[list(PARAMS)].to_numpy(float)
    sobol = qmc.Sobol(d=len(PARAMS), scramble=True, seed=20260923)
    ab = sobol.random_base2(m=11)
    a = np.asarray([BOUNDS[p][0] + ab[:, i] * (BOUNDS[p][1] - BOUNDS[p][0]) for i, p in enumerate(PARAMS)]).T
    b = np.asarray([BOUNDS[p][0] + ab[:, i + len(PARAMS)] * (BOUNDS[p][1] - BOUNDS[p][0]) for i, p in enumerate(PARAMS)]) if False else None
    # Use a second independently scrambled block for B; this keeps the design
    # deterministic without pretending the one-block campaign was Saltelli.
    bunit = qmc.Sobol(d=len(PARAMS), scramble=True, seed=20260924).random_base2(m=11)
    b = np.asarray([BOUNDS[p][0] + bunit[:, i] * (BOUNDS[p][1] - BOUNDS[p][0]) for i, p in enumerate(PARAMS)]).T
    sobol_rows = []
    for endpoint in ("H4_PRESENT", "EXACT_H4"):
        ycal = calibration[endpoint].astype(float).to_numpy()
        model = ExtraTreesRegressor(n_estimators=32, min_samples_leaf=4, random_state=20260926, n_jobs=1).fit(xcal, ycal)
        ya = model.predict(a); yb = model.predict(b); var = float(np.var(np.concatenate([ya, yb]), ddof=1)) or 1.0
        for i, parameter in enumerate(PARAMS):
            ab_i = a.copy(); ab_i[:, i] = b[:, i]
            yab = model.predict(ab_i)
            s1 = float(np.mean(yb * (yab - ya)) / var)
            st = float(np.mean((ya - yab) ** 2) / (2 * var))
            sobol_rows.append({"endpoint": endpoint, "parameter": parameter, "S1": s1, "ST": st, "N_base": 1024, "status": "SURROGATE_SALTELLI", "deviation": "indices use calibrated surrogate, not fresh exact Saltelli DAE evaluations"})
    write_table(pd.DataFrame(sobol_rows), "TX4_SOBOL.csv")

    logistic_rows = []
    mcx = mc[list(PARAMS)].to_numpy(float)
    scaler = StandardScaler().fit(mcx)
    xs = scaler.transform(mcx)
    for endpoint in ("H4_PRESENT", "EXACT_H4"):
        y = mc[endpoint].astype(int).to_numpy()
        if len(np.unique(y)) < 2:
            logistic_rows.append({"endpoint": endpoint, "status": "ONE_CLASS", "parameter": "", "coefficient": "", "odds_ratio": ""})
            continue
        model = LogisticRegression(max_iter=2000, random_state=20260927).fit(xs, y)
        for parameter, coef in zip(PARAMS, model.coef_[0]):
            logistic_rows.append({"endpoint": endpoint, "status": "FIT", "parameter": parameter, "coefficient": float(coef), "odds_ratio": float(np.exp(coef)), "n": len(y), "deviation": "predictive logistic fit; no causal interpretation"})
    write_table(pd.DataFrame(logistic_rows), "TX4_LOGISTIC_MODEL.csv")

    # g-star distribution from exact-H4 2-D phase lines and the nominal 1-D g line.
    gstar_rows = []
    def crossing(g, alpha):
        order = np.argsort(g); g = np.asarray(g)[order]; alpha = np.asarray(alpha)[order]
        for i in range(len(g) - 1):
            if alpha[i] == 0 or alpha[i] * alpha[i + 1] <= 0:
                if alpha[i + 1] == alpha[i]: return float(g[i])
                return float(g[i] + (0 - alpha[i]) * (g[i + 1] - g[i]) / (alpha[i + 1] - alpha[i]))
        return float("nan")
    one = h4[h4["campaign"] == "1D"]
    gline = one[one["condition_id"].str.startswith("1d_g_")]
    gstar_rows.append({"source": "1D_nominal", "fixed_parameter": "nominal", "fixed_value": 1.0, "g_star": crossing(gline["g"], gline["alpha_all"]), "evaluator_tier": "EXACT_H4_DAE"})
    pairs = [("g", "k"), ("g", "t"), ("g", "h"), ("g", "epsilon"), ("k", "epsilon"), ("t", "epsilon")]
    phase = h4[h4["campaign"] == "2D"].copy()
    for first, second in pairs:
        if first != "g":
            continue
        tag = f"2d_{first}x{second}_"
        sub = phase[phase["condition_id"].str.startswith(tag)]
        for fixed, group in sub.groupby(second):
            gstar_rows.append({"source": f"2D_{first}x{second}", "fixed_parameter": second, "fixed_value": float(fixed), "g_star": crossing(group["g"], group["alpha_all"]), "evaluator_tier": "EXACT_H4_DAE", "n_points": len(group)})
    write_table(pd.DataFrame(gstar_rows), "TX4_GSTAR_DISTRIBUTION.csv")

    method_rows = []
    paired = cond[cond["campaign"].isin(["QMC", "MC"])].copy()
    for campaign_name, sub in paired.groupby("campaign"):
        full = sub["stable_full"].to_numpy(bool); mode = sub["stable_mode"].to_numpy(bool)
        both = int(np.sum(full & mode)); full_only = int(np.sum(full & ~mode)); mode_only = int(np.sum(~full & mode)); neither = int(np.sum(~full & ~mode))
        table = [[both, full_only], [mode_only, neither]]
        p = float(mcnemar(table, exact=True, correction=False).pvalue)
        method_rows.append({"campaign": campaign_name, "n": len(sub), "both_stable": both, "full_only": full_only, "mode_only": mode_only, "neither": neither, "mcnemar_p": p, "status": "PAIRED_MODE_COMPARISON"})
    write_table(pd.DataFrame(method_rows), "TX4_METHOD_COMPARISON.csv")

    return_rows = []
    for campaign_name, sub in paired.groupby("campaign"):
        return_rows.append({"campaign": campaign_name, "n": len(sub), "mean_delta_H4": float(sub["delta_H4"].mean()), "mean_eta_H4": float(sub["eta_H4"].mean()), "median_q_H4": float(sub["q_H4"].median()), "min_q_H4": float(sub["q_H4"].min()), "H4_PRESENT_rate": float(sub["H4_PRESENT"].mean()), "evaluator_tier": "EXACT_H4_DAE_PLUS_SURROGATE_PROPER"})
    write_table(pd.DataFrame(return_rows), "TX4_RETURN_ROBUSTNESS.csv")
    mode_rows = []
    for campaign_name, sub in paired.groupby("campaign"):
        mode_rows.extend([
            {"campaign": campaign_name, "endpoint": "FULL_SPECTRUM", "positive": int(sub["stable_full"].sum()), "n": len(sub), "rate": float(sub["stable_full"].mean())},
            {"campaign": campaign_name, "endpoint": "EM_0.3_1.5_HZ", "positive": int(sub["stable_mode"].sum()), "n": len(sub), "rate": float(sub["stable_mode"].mean())},
        ])
    write_table(pd.DataFrame(mode_rows), "TX4_MODE_ROBUSTNESS.csv")

    # TDS endpoint is deliberately labeled a spectral trace proxy; the prior
    # PowerDynamics gate did not authorize a new same-model TDS claim.
    tds_rows = []
    for _, row in h4[h4["campaign"] == "QMC"].head(24).iterrows():
        tds_rows.append({"condition_id": row["condition_id"], "alpha_all": row["alpha_all"], "critical_frequency_hz": row["critical_frequency_hz"], "trace_status": "SPECTRAL_TRACE_PROXY", "stable_proxy": bool(row["alpha_all"] <= THRESHOLD), "evaluator_tier": row["evaluator_tier"], "note": "not a new nonlinear TDS run"})
    write_table(pd.DataFrame(tds_rows), "TX4_TDS_ROBUSTNESS.csv")

    # Cross-code file: preserve the exact nominal 16-case parity and explicitly
    # reserve the 16 random slots as not executed rather than fabricating Julia.
    cross_rows = []
    source = RESULTS / "TX4_EXACT_P4_V4_CROSSCODE.csv"
    if source.exists():
        nominal = pd.read_csv(source)
        for _, row in nominal.iterrows():
            item = row.to_dict(); item.update({"condition_id": "nominal_p4", "status": "EXACT_NOMINAL_CROSSCODE", "random_seed": 20260924})
            cross_rows.append(item)
    for index in range(16):
        cross_rows.append({"condition_id": f"random_{index:02d}", "status": "NOT_EXECUTED_JULIA_RANDOM", "random_seed": 20260924, "note": "fixed Julia script supports frozen nominal census only"})
    pd.DataFrame(cross_rows).to_csv(RESULTS / "TX4_JULIA_CROSSCODE_RANDOM.csv", index=False)

    claims = pd.DataFrame([
        {"claim_id": "R1", "status": "PASS", "evidence": "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv; frozen exact P4 closure", "scope": "nominal exact H4"},
        {"claim_id": "R2", "status": "BOUNDED_BOX_WITH_SURROGATE_PROPER", "evidence": "TX4_QMC_SUMMARY.csv; TX4_MC_SUMMARY.csv", "scope": "engineering box only"},
        {"claim_id": "R3", "status": "BOUNDED_BOX_WITH_SURROGATE_PROPER", "evidence": "TX4_QMC_SUMMARY.csv; TX4_MC_SUMMARY.csv", "scope": "engineering box only"},
        {"claim_id": "R4", "status": "NOT_CLAIMED", "evidence": "TX4_ROBUSTNESS_DEVIATIONS.md", "scope": "no physical robust radius"},
        {"claim_id": "R5", "status": "NOT_CLAIMED", "evidence": "TX4_JULIA_CROSSCODE_RANDOM.csv", "scope": "random Julia spot-check not executed"},
        {"claim_id": "R6", "status": "NOT_CLAIMED", "evidence": "TX4_ROBUSTNESS_DEVIATIONS.md", "scope": "no EMT/current-limit/hardware claim"},
        {"claim_id": "R7", "status": "NOMINAL_ONLY", "evidence": "TX4_JULIA_CROSSCODE_RANDOM.csv", "scope": "16 nominal cross-code rows; 16 random slots not executed"},
    ])
    claims.to_csv(RESULTS / "TX4_ROBUSTNESS_CLAIM_MATRIX.csv", index=False)

    headline = {
        "status": "COMPLETE_WITH_EXPLICIT_SURROGATE_TIER",
        "conditions": int(len(cond)), "portfolio_rows": int(len(df)),
        "exact_h4_conditions": int(len(h4)),
        "exact_all_portfolio_rows": int((df["evaluator_tier"] == "EXACT_16_DAE").sum()),
        "surrogate_proper_rows": int((df["evaluator_tier"] == "SURROGATE_CALIBRATED").sum()),
        "qmc_h4_present_rate": float(cond[cond.campaign == "QMC"]["H4_PRESENT"].mean()),
        "mc_h4_present_rate": float(cond[cond.campaign == "MC"]["H4_PRESENT"].mean()),
        "qmc_exact_h4_rate": float(cond[cond.campaign == "QMC"]["EXACT_H4"].mean()),
        "mc_exact_h4_rate": float(cond[cond.campaign == "MC"]["EXACT_H4"].mean()),
        "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
        "no_push": True,
    }
    (RESULTS / "TX4_ROBUSTNESS_HEADLINE.json").write_text(json.dumps(headline, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(headline, indent=2))


if __name__ == "__main__":
    main()
