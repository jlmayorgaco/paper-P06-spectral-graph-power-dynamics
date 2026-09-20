"""Build corrected blocker and continuous summaries from exact QMC/MC rows."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, norm

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_existing_robustness import condition_truth  # noqa: E402
from tx4_audit_blockers import H4_LABEL  # noqa: E402


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return (np.nan, np.nan)
    p = success / total
    denom = 1.0 + z * z / total
    centre = (p + z * z / (2 * total)) / denom
    half = z * np.sqrt((p * (1 - p) / total) + z * z / (4 * total * total)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def condition_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = [condition_truth(group) for _, group in df.groupby("condition_id", sort=False)]
    result = pd.DataFrame(rows)
    result["evaluation_source"] = "EXACT_DAE"
    result["U005_ANALYZED"] = False
    result["H005_ANALYZED"] = False
    return result


def summaries(campaign: str, truth: pd.DataFrame) -> None:
    n = len(truth)
    endpoint_rows: list[dict[str, object]] = []
    for endpoint in ("H4_UNSTABLE", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"):
        success = int(truth[endpoint].sum())
        row: dict[str, object] = {
            "campaign": campaign,
            "endpoint": endpoint,
            "n_conditions": n,
            "success": success,
            "coverage_fraction": success / n if n else np.nan,
            "interpretation": "fixed-design QMC coverage fraction" if campaign == "QMC" else "assumed independent-uniform bounded-box probability",
            "ci_method": "none" if campaign == "QMC" else "Wilson 95%",
            "ci95_low": np.nan,
            "ci95_high": np.nan,
        }
        if campaign == "MC":
            row["ci95_low"], row["ci95_high"] = wilson(success, n)
        endpoint_rows.append(row)
    kappa_counts = truth["KAPPA"].where(truth["KAPPA"].notna(), "NULL").value_counts(dropna=False).sort_index(key=lambda s: s.map(lambda x: 999 if x == "NULL" else int(x)))
    for kappa, count in kappa_counts.items():
        endpoint_rows.append(
            {
                "campaign": campaign,
                "endpoint": "KAPPA_PMF",
                "kappa": kappa,
                "n_conditions": n,
                "success": int(count),
                "coverage_fraction": int(count) / n if n else np.nan,
                "interpretation": "fixed-design QMC PMF" if campaign == "QMC" else "assumed independent-uniform bounded-box PMF",
                "ci_method": "none" if campaign == "QMC" else "Wilson 95%",
                "ci95_low": np.nan if campaign == "QMC" else wilson(int(count), n)[0],
                "ci95_high": np.nan if campaign == "QMC" else wilson(int(count), n)[1],
            }
        )
    pd.DataFrame(endpoint_rows).to_csv(RESULTS / f"TX4_{campaign}_EXACT_BLOCKER_SUMMARY.csv", index=False)


def delta_summary(campaign: str, truth: pd.DataFrame) -> None:
    values = truth.delta_H4_true_minimality.dropna().to_numpy(dtype=float)
    row: dict[str, object] = {
        "campaign": campaign,
        "n": len(values),
        "median": np.nan,
        "iqr_low": np.nan,
        "iqr_high": np.nan,
        "p05": np.nan,
        "p95": np.nan,
        "min": np.nan,
        "max": np.nan,
        "median_ci95_low": np.nan,
        "median_ci95_high": np.nan,
        "interpretation": "QMC descriptive exact true-minimality margin" if campaign == "QMC" else "MC descriptive exact true-minimality margin with BCa median interval",
    }
    if len(values):
        row.update(
            {
                "median": float(np.median(values)),
                "iqr_low": float(np.quantile(values, 0.25)),
                "iqr_high": float(np.quantile(values, 0.75)),
                "p05": float(np.quantile(values, 0.05)),
                "p95": float(np.quantile(values, 0.95)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }
        )
        if campaign == "MC":
            try:
                result = bootstrap(
                    (values,), np.median, confidence_level=0.95,
                    n_resamples=5000, method="BCa",
                    rng=np.random.default_rng(20260926),
                )
                row["median_ci95_low"] = float(result.confidence_interval.low)
                row["median_ci95_high"] = float(result.confidence_interval.high)
            except TypeError:
                result = bootstrap(
                    (values,), np.median, confidence_level=0.95,
                    n_resamples=5000, method="BCa",
                    random_state=np.random.default_rng(20260926),
                )
                row["median_ci95_low"] = float(result.confidence_interval.low)
                row["median_ci95_high"] = float(result.confidence_interval.high)
    pd.DataFrame([row]).to_csv(RESULTS / f"TX4_{campaign}_DELTA_H4_TRUE_MINIMALITY.csv", index=False)


def merge_condition_flags(campaign: str, df: pd.DataFrame, truth: pd.DataFrame) -> None:
    flag_cols = ["U0_json", "H0_json", "H4_UNSTABLE", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE", "KAPPA", "delta_H4_true_minimality", "evaluation_source"]
    endpoints = ("H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE")
    # Make the operation idempotent: a prior run may already have written
    # corrected and legacy columns into the parquet.
    base = df.drop(
        columns=[
            col for col in [*flag_cols, *(f"{endpoint}_corrected" for endpoint in endpoints), *(f"{endpoint}_legacy" for endpoint in endpoints)]
            if col in df.columns
        ],
        errors="ignore",
    ).copy()
    for endpoint in endpoints:
        if endpoint in df.columns:
            base[f"{endpoint}_legacy"] = df[endpoint].to_numpy()
    enriched = base.merge(truth[["condition_id", *flag_cols]], on="condition_id", how="left")
    enriched.to_parquet(RESULTS / f"TX4_{campaign}_ALL16_EXACT.parquet", index=False)
    enriched.to_csv(RESULTS / f"TX4_{campaign}_ALL16_EXACT.csv.gz", index=False, compression="gzip")
    truth.to_csv(RESULTS / f"TX4_{campaign}_EXACT_BLOCKER_TRUTH.csv", index=False)


def main() -> None:
    for campaign in ("QMC", "MC"):
        path = RESULTS / f"TX4_{campaign}_ALL16_EXACT.parquet"
        df = pd.read_parquet(path)
        truth = condition_table(df)
        merge_condition_flags(campaign, df, truth)
        summaries(campaign, truth)
        delta_summary(campaign, truth)
        print(json.dumps({"campaign": campaign, "rows": len(df), "conditions": len(truth)}, indent=2))
    eta_path = RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv"
    if not eta_path.exists():
        pd.DataFrame(
            [
                {
                    "campaign": campaign,
                    "analysis": "local_collective_eta_exact_stratified",
                    "n_requested_max": 128,
                    "n_clean_H4_minimal": np.nan,
                    "n_computed": 0,
                    "status": "PENDING_LOCAL_Q_RECONSTRUCTION",
                    "eta_definition": "min_i sigma_min(I+Q_{H4\\i}(s_H))",
                    "note": "Retained exact master rows do not store Q_local or s_H; no surrogate eta is promoted.",
                }
                for campaign in ("QMC", "MC")
            ]
        ).to_csv(eta_path, index=False)


if __name__ == "__main__":
    main()
