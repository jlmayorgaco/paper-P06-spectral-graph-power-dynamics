"""Audit the retained TX4 master without changing sampled coordinates."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import KFold

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tx4_audit_blockers import H4_LABEL, classify_blockers  # noqa: E402


PARAMS = ("g", "k", "t", "h", "epsilon", "damping", "inertia")
PORTFOLIOS = [
    "BASE", "30", "33", "30+33", "35", "30+35", "33+35", "30+33+35",
    "37", "30+37", "33+37", "30+33+37", "35+37", "30+35+37",
    "33+35+37", H4_LABEL,
]


def evaluation_source(tier: str) -> str:
    if tier in {"EXACT_16_DAE", "EXACT_H4_DAE"}:
        return "EXACT_DAE"
    if tier == "SURROGATE_CALIBRATED":
        return "SURROGATE_EXTRATREES"
    if tier == "DERIVED_FROM_EXACT":
        return "DERIVED_FROM_EXACT"
    return "UNKNOWN"


def read_master() -> pd.DataFrame:
    path = RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.parquet"
    df = pd.read_parquet(path)
    df["evaluation_source"] = df["evaluator_tier"].map(evaluation_source)
    df.to_parquet(path, index=False)
    df.to_csv(RESULTS / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv.gz", index=False, compression="gzip")
    provenance_cols = [
        "campaign", "condition_id", "seed", "portfolio", "evaluator_tier",
        "evaluation_source", "status", *PARAMS,
    ]
    df[provenance_cols].to_csv(RESULTS / "TX4_ROBUSTNESS_DATA_PROVENANCE.csv", index=False)
    return df


def condition_truth(group: pd.DataFrame) -> dict[str, object]:
    alpha = {
        str(row.portfolio): row.alpha_EM if row.status == "OK" else None
        for row in group.itertuples(index=False)
    }
    result = classify_blockers(alpha)
    h4 = group[group.portfolio == H4_LABEL].iloc[0]
    proper = group[group.portfolio != H4_LABEL]
    finite_proper = pd.to_numeric(proper.alpha_EM, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    max_proper = float(finite_proper.max()) if len(finite_proper) else np.nan
    h4_alpha = float(h4.alpha_EM) if pd.notna(h4.alpha_EM) else np.nan
    delta = float(min(h4_alpha, -max_proper)) if np.isfinite(h4_alpha) and np.isfinite(max_proper) else np.nan
    result.update(
        {
            "campaign": str(h4.campaign),
            "condition_id": str(h4.condition_id),
            "seed": int(h4.seed),
            **{name: float(getattr(h4, name)) for name in PARAMS},
            "n_portfolios": int(len(group)),
            "n_ok": int((group.status == "OK").sum()),
            "missing_alpha_em": int(pd.to_numeric(group.alpha_EM, errors="coerce").isna().sum()),
            "U0_json": json.dumps(result["U0"], separators=(",", ":")),
            "H0_json": json.dumps(result["H0"], separators=(",", ":")),
            "legacy_H4_UNSTABLE": bool(pd.notna(h4.alpha_EM) and float(h4.alpha_EM) > 1e-8),
            "legacy_H4_PRESENT": bool(h4.H4_PRESENT) if pd.notna(h4.H4_PRESENT) else False,
            "legacy_EXACT_H4": bool(h4.EXACT_H4) if pd.notna(h4.EXACT_H4) else False,
            "legacy_NONCOMPOSABLE": bool(h4.NONCOMPOSABLE) if pd.notna(h4.NONCOMPOSABLE) else False,
            "legacy_eta_H4": float(h4.eta_H4) if pd.notna(h4.eta_H4) else np.nan,
            "h4_alpha_EM": h4_alpha,
            "max_proper_alpha_EM": max_proper,
            "delta_H4_true_minimality": delta,
            "U005_ANALYZED": False,
            "H005_ANALYZED": False,
        }
    )
    return result


def exact_calibration_audit(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    exact = df[df.evaluator_tier == "EXACT_16_DAE"]
    truth = pd.DataFrame([condition_truth(g) for _, g in exact.groupby("condition_id", sort=True)])
    truth.to_csv(RESULTS / "TX4_EXACT_BLOCKER_CALIBRATION_TRUTH.csv", index=False)
    blocker_rows: list[dict[str, object]] = []
    for row in truth.itertuples(index=False):
        h0 = json.loads(row.H0_json)
        for rank, label in enumerate(h0, 1):
            blocker_rows.append(
                {
                    "campaign": row.campaign,
                    "condition_id": row.condition_id,
                    "blocker_rank": rank,
                    "blocker": label,
                    "blocker_size": len(label.split("+")) if label != "BASE" else 0,
                    "H4_BLOCKER": label == H4_LABEL,
                }
            )
    pd.DataFrame(blocker_rows).to_csv(RESULTS / "TX4_EXACT_BLOCKER_CALIBRATION_BLOCKERS.csv", index=False)
    return exact, truth


def h4_flag_audit(df: pd.DataFrame, truth: pd.DataFrame) -> None:
    rows: list[dict[str, object]] = []
    for scope in ("1D", "QMC", "MC"):
        subset = df[(df.campaign == scope) & (df.portfolio == H4_LABEL)]
        n = int(subset.condition_id.nunique())
        rows.append(
            {
                "scope": scope,
                "n_conditions": n,
                "h4_rows_exact_dae": int((subset.evaluation_source == "EXACT_DAE").sum()),
                "alpha_H4_ge0_n": int((pd.to_numeric(subset.alpha_EM, errors="coerce") >= 0).sum()),
                "alpha_H4_ge0_fraction": float((pd.to_numeric(subset.alpha_EM, errors="coerce") >= 0).mean()),
                "legacy_H4_PRESENT_n": int(subset.H4_PRESENT.sum()),
                "legacy_EXACT_H4_n": int(subset.EXACT_H4.sum()),
                "legacy_NONCOMPOSABLE_n": int(subset.NONCOMPOSABLE.sum()),
                "corrected_exact_truth_available": False,
                "corrected_H4_PRESENT_n": np.nan,
                "corrected_EXACT_H4_n": np.nan,
                "corrected_NONCOMPOSABLE_n": np.nan,
            }
        )
    rows.append(
        {
            "scope": "EXACT_CALIBRATION",
            "n_conditions": len(truth),
            "h4_rows_exact_dae": int(len(truth)),
            "alpha_H4_ge0_n": int(truth.H4_UNSTABLE.sum()),
            "alpha_H4_ge0_fraction": float(truth.H4_UNSTABLE.mean()),
            "legacy_H4_PRESENT_n": int(truth.legacy_H4_PRESENT.sum()),
            "legacy_EXACT_H4_n": int(truth.legacy_EXACT_H4.sum()),
            "legacy_NONCOMPOSABLE_n": int(truth.legacy_NONCOMPOSABLE.sum()),
            "corrected_exact_truth_available": True,
            "corrected_H4_PRESENT_n": int(truth.H4_PRESENT.sum()),
            "corrected_EXACT_H4_n": int(truth.EXACT_H4.sum()),
            "corrected_NONCOMPOSABLE_n": int(truth.NONCOMPOSABLE.sum()),
        }
    )
    pd.DataFrame(rows).to_csv(RESULTS / "TX4_AUDIT_H4_FLAG_COMPARISON.csv", index=False)

    confusion: list[dict[str, object]] = []
    for endpoint, legacy in (("H4_UNSTABLE", "legacy_H4_UNSTABLE"), ("H4_PRESENT", "legacy_H4_PRESENT"), ("EXACT_H4", "legacy_EXACT_H4"), ("NONCOMPOSABLE", "legacy_NONCOMPOSABLE")):
        for truth_value in (False, True):
            for legacy_value in (False, True):
                confusion.append(
                    {
                        "endpoint": endpoint,
                        "truth": truth_value,
                        "legacy": legacy_value,
                        "n_conditions": int(((truth[endpoint] == truth_value) & (truth[legacy] == legacy_value)).sum()),
                        "interpretation": "legacy flag is non-minimal or semantically mismatched" if endpoint in {"H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"} else "alpha_EM sign comparison",
                    }
                )
    pd.DataFrame(confusion).to_csv(RESULTS / "TX4_AUDIT_H4_FLAG_CONFUSION.csv", index=False)

    complement = truth[(truth.legacy_H4_PRESENT) & (~truth.H4_PRESENT)]
    pd.DataFrame(
        [
            {
                "audit": "legacy_H4_PRESENT_true_but_not_minimal_H4",
                "n": len(complement),
                "fraction_of_exact_calibration": len(complement) / len(truth),
                "conclusion": "YES" if len(complement) else "NO",
                "note": "The retained H4_PRESENT flag tests H4 instability, not membership in the H0 minimal antichain.",
            }
        ]
    ).to_csv(RESULTS / "TX4_AUDIT_COMPLEMENT_BUG.csv", index=False)


def eta_audit(df: pd.DataFrame, truth: pd.DataFrame) -> None:
    rows = truth[
        [
            "campaign", "condition_id", "h4_alpha_EM", "max_proper_alpha_EM",
            "legacy_eta_H4", "delta_H4_true_minimality", "H4_PRESENT", "EXACT_H4",
        ]
    ].copy()
    rows["legacy_formula_eta"] = -rows["max_proper_alpha_EM"]
    rows["legacy_formula_matches"] = np.isclose(rows["legacy_eta_H4"], rows["legacy_formula_eta"], equal_nan=True)
    rows["clean_exact_h4_root"] = rows["H4_PRESENT"] & rows["EXACT_H4"]
    rows["eta_corrected"] = np.nan
    rows["eta_definition"] = "min_i sigma_min(I+Q_{H4\\i}(s_H))"
    rows["eta_status"] = "NOT_COMPUTABLE_FROM_RETAINED_MASTER; Q_local_not_stored"
    rows.to_csv(RESULTS / "TX4_AUDIT_ETA_H4.csv", index=False)


def surrogate_validation(exact: pd.DataFrame, truth: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    condition_ids = sorted(truth.condition_id.tolist())
    folds = list(KFold(n_splits=5, shuffle=True, random_state=20260925).split(condition_ids))
    predictions: list[dict[str, object]] = []
    feature_by_id = exact.groupby("condition_id").first()
    for portfolio in PORTFOLIOS:
        port = exact[exact.portfolio == portfolio].set_index("condition_id").loc[condition_ids]
        y = pd.to_numeric(port.alpha_EM, errors="coerce").to_numpy(dtype=float)
        x = port[list(PARAMS)].to_numpy(dtype=float)
        pred = np.full(len(condition_ids), np.nan, dtype=float)
        for train_idx, test_idx in folds:
            train_ok = np.isfinite(y[train_idx])
            model = ExtraTreesRegressor(
                n_estimators=32,
                min_samples_leaf=2,
                random_state=20260925,
                n_jobs=1,
            )
            model.fit(x[train_idx][train_ok], y[train_idx][train_ok])
            pred[test_idx] = model.predict(x[test_idx])
        for index, condition_id in enumerate(condition_ids):
            predictions.append(
                {
                    "condition_id": condition_id,
                    "campaign": str(port.loc[condition_id, "campaign"]),
                    "portfolio": portfolio,
                    "true_alpha_EM": y[index],
                    "predicted_alpha_EM": pred[index],
                    "absolute_error": abs(y[index] - pred[index]) if np.isfinite(y[index]) else np.nan,
                    "true_unstable": bool(np.isfinite(y[index]) and y[index] >= 0),
                    "predicted_unstable": bool(np.isfinite(pred[index]) and pred[index] >= 0),
                    "prediction_source": "RECONSTRUCTED_5FOLD_EXTRA_TREES",
                }
            )
    pred_df = pd.DataFrame(predictions)
    pred_df.to_csv(RESULTS / "TX4_SURROGATE_ALPHA_VALIDATION.csv", index=False)
    pred_maps = {
        condition_id: dict(zip(g.portfolio, g.predicted_alpha_EM))
        for condition_id, g in pred_df.groupby("condition_id")
    }
    validation_rows: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []
    truth_by_id = truth.set_index("condition_id")
    for condition_id, predicted_alpha in pred_maps.items():
        pred_flags = classify_blockers(predicted_alpha)
        true_row = truth_by_id.loc[condition_id]
        h0_true = json.loads(true_row.H0_json)
        h0_pred = pred_flags["H0"]
        false_minimal = bool(pred_flags["EXACT_H4"] and not true_row.EXACT_H4)
        missed_minimal = bool(true_row.EXACT_H4 and not pred_flags["EXACT_H4"])
        validation_rows.append(
            {
                "condition_id": condition_id,
                "campaign": true_row.campaign,
                "true_H4_UNSTABLE": true_row.H4_UNSTABLE,
                "pred_H4_UNSTABLE": pred_flags["H4_UNSTABLE"],
                "true_H4_PRESENT": true_row.H4_PRESENT,
                "pred_H4_PRESENT": pred_flags["H4_PRESENT"],
                "true_EXACT_H4": true_row.EXACT_H4,
                "pred_EXACT_H4": pred_flags["EXACT_H4"],
                "true_NONCOMPOSABLE": true_row.NONCOMPOSABLE,
                "pred_NONCOMPOSABLE": pred_flags["NONCOMPOSABLE"],
                "true_KAPPA": true_row.KAPPA,
                "pred_KAPPA": pred_flags["KAPPA"],
                "true_H0_json": true_row.H0_json,
                "pred_H0_json": json.dumps(h0_pred, separators=(",", ":")),
                "h0_antichain_exact": h0_true == h0_pred,
                "false_minimal": false_minimal,
                "missed_minimal": missed_minimal,
                "prediction_source": "RECONSTRUCTED_5FOLD_EXTRA_TREES",
            }
        )
        if h0_true != h0_pred:
            errors.append(
                {
                    "condition_id": condition_id,
                    "campaign": true_row.campaign,
                    "error_type": "H0_ANTICHAIN_MISMATCH",
                    "true_H0_json": true_row.H0_json,
                    "pred_H0_json": json.dumps(h0_pred, separators=(",", ":")),
                }
            )
        if false_minimal or missed_minimal:
            errors.append(
                {
                    "condition_id": condition_id,
                    "campaign": true_row.campaign,
                    "error_type": "EXACT_H4_FALSE_MINIMAL_OR_MISSED",
                    "true_H0_json": true_row.H0_json,
                    "pred_H0_json": json.dumps(h0_pred, separators=(",", ":")),
                }
            )
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(RESULTS / "TX4_SURROGATE_MINIMALITY_VALIDATION.csv", index=False)
    pd.DataFrame(errors).to_csv(RESULTS / "TX4_SURROGATE_BLOCKER_ERRORS.csv", index=False)
    gate = {
        "false_minimal_n": int(validation.false_minimal.sum()),
        "missed_minimal_n": int(validation.missed_minimal.sum()),
        "exact_h4_flag_agreement": float((validation.true_EXACT_H4 == validation.pred_EXACT_H4).mean()),
        "h4_present_flag_agreement": float((validation.true_H4_PRESENT == validation.pred_H4_PRESENT).mean()),
        "h0_antichain_exact_fraction": float(validation.h0_antichain_exact.mean()),
        "strict_gate_pass": bool(
            validation.false_minimal.sum() == 0
            and validation.missed_minimal.sum() == 0
            and (validation.true_EXACT_H4 == validation.pred_EXACT_H4).all()
            and (validation.true_H4_PRESENT == validation.pred_H4_PRESENT).all()
            and validation.h0_antichain_exact.all()
        ),
    }
    (RESULTS / "TX4_SURROGATE_GATE.json").write_text(json.dumps(gate, indent=2) + "\n", encoding="utf-8")
    return validation, gate


def main() -> None:
    df = read_master()
    exact, truth = exact_calibration_audit(df)
    h4_flag_audit(df, truth)
    eta_audit(df, truth)
    validation, gate = surrogate_validation(exact, truth)
    print(json.dumps({"rows": len(df), "exact_conditions": len(truth), "surrogate_gate": gate}, indent=2))


if __name__ == "__main__":
    main()
