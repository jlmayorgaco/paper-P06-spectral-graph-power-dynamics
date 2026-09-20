"""Two-tier execution wrapper for the preregistered TX4 robustness design.

The H4 row is always evaluated by the exact reduced DAE.  The 16-portfolio
calibration set is also exact.  The remaining proper-subset rows are generated
by an ExtraTrees surrogate trained only on the preregistered exact calibration
set.  The tier is explicit in every row and in the final report; no surrogate
row is allowed to masquerade as an exact DAE result.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_tx4_robustness as exact  # noqa: E402

from sklearn.ensemble import ExtraTreesRegressor  # noqa: E402


OUTPUT_FIELDS = [*exact.FIELDS, "evaluator_tier"]
H4 = (30, 33, 35, 37)
SURROGATE_METRICS = ("alpha_all", "alpha_EM", "critical_frequency_hz", "ell_H4", "q_H4")


def write_rows(path: Path, rows: list[dict[str, Any]], mode: str = "a") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and mode == "a"
    with path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def calibration_id(condition_id: str) -> bool:
    if condition_id == "nominal":
        return True
    if condition_id.startswith("1d_"):
        # Exact all-portfolio calibration uses a 21-point deterministic
        # skeleton. The full 201-point 1-D design is still retained, with H4
        # evaluated exactly and proper-subset rows explicitly surrogate-tier.
        index = int(condition_id.rsplit("_", 1)[-1])
        return index % 10 == 0 or index == 200
    if condition_id.startswith("qmc_"):
        return int(condition_id.split("_")[1]) < 64
    if condition_id.startswith("mc_"):
        return int(condition_id.split("_")[1]) < 64
    if condition_id.startswith("2d_"):
        parts = condition_id.rsplit("_", 2)
        return int(parts[-2]) % 10 == 0 and int(parts[-1]) % 10 == 0
    return False


def exact_h4(item: tuple[str, str, int, dict[str, float]]) -> dict[str, Any]:
    campaign, condition_id, seed, params = item
    row = exact._one_portfolio(params, H4, campaign, condition_id, seed)
    row["evaluator_tier"] = "EXACT_H4_DAE"
    return row


def run_pool(items: list[tuple[str, str, int, dict[str, float]]], fn, workers: int) -> list[Any]:
    results: list[Any] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fn, item) for item in items]
        for index, future in enumerate(as_completed(futures), 1):
            results.append(future.result())
            if index % 100 == 0 or index == len(items):
                print(f"TX4_CAMPAIGN_PROGRESS completed={index}/{len(items)}", flush=True)
    return results


def fit_models(rows: list[dict[str, Any]]) -> dict[str, dict[str, ExtraTreesRegressor]]:
    models: dict[str, dict[str, ExtraTreesRegressor]] = {}
    x_all = np.asarray([[float(row[name]) for name in exact.PARAM_NAMES] for row in rows], dtype=float)
    for portfolio in exact.PORTFOLIO_LABELS.values():
        subset = [row for row in rows if row["portfolio"] == portfolio and row["status"] == "OK"]
        x = np.asarray([[float(row[name]) for name in exact.PARAM_NAMES] for row in subset], dtype=float)
        models[portfolio] = {}
        for metric in SURROGATE_METRICS:
            usable = [row for row in subset if row.get(metric) is not None and np.isfinite(float(row[metric]))]
            if len(usable) < 20:
                usable = subset
            xx = np.asarray([[float(row[name]) for name in exact.PARAM_NAMES] for row in usable], dtype=float)
            yy = np.asarray([float(row[metric]) if row.get(metric) is not None and np.isfinite(float(row[metric])) else 0.0 for row in usable], dtype=float)
            model = ExtraTreesRegressor(
                n_estimators=96,
                min_samples_leaf=2,
                random_state=20260925,
                n_jobs=-1,
            )
            model.fit(xx, yy)
            models[portfolio][metric] = model
    return models


def surrogate_row(
    params: dict[str, float], item: tuple[str, str, int, dict[str, float]], template: dict[str, Any],
    portfolio: str, models: dict[str, dict[str, ExtraTreesRegressor]],
) -> dict[str, Any]:
    campaign, condition_id, seed, _ = item
    row = dict(template)
    row.update(
        {
            "campaign": campaign,
            "condition_id": condition_id,
            "seed": seed,
            "portfolio": portfolio,
            **{name: float(params[name]) for name in exact.PARAM_NAMES},
            "status": "APPROX_SURROGATE",
            "error_type": "",
            "error_message": "ExtraTrees trained on exact preregistered calibration rows",
            "NONCOMPOSABLE": False,
            "evaluator_tier": "SURROGATE_CALIBRATED",
            "condition_runtime_s": 0.0,
        }
    )
    x = np.asarray([[float(params[name]) for name in exact.PARAM_NAMES]], dtype=float)
    for metric in SURROGATE_METRICS:
        row[metric] = float(models[portfolio][metric].predict(x)[0])
    row["stable_all"] = bool(row["alpha_all"] <= exact.THRESHOLD)
    row["stable_EM"] = bool(np.isfinite(row["alpha_EM"]) and row["alpha_EM"] <= exact.THRESHOLD)
    return row


def run() -> None:
    workers = max(1, min(18, (os.cpu_count() or 2) - 2))
    all_conditions = exact.make_conditions()
    calibration = [item for item in all_conditions if calibration_id(item[1])]
    remaining = [item for item in all_conditions if not calibration_id(item[1])]
    results = ROOT / "results"
    master = results / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv"
    calibration_file = results / "TX4_ROBUSTNESS_EXACT_CALIBRATION.csv"
    for path in (master, calibration_file):
        if path.exists():
            path.unlink()
    start = time.perf_counter()
    print(f"TX4_CAMPAIGN_CALIBRATION conditions={len(calibration)} workers={workers}", flush=True)
    exact_groups = run_pool(calibration, exact.evaluate_condition, workers)
    exact_rows: list[dict[str, Any]] = []
    for group in exact_groups:
        for row in group:
            row["evaluator_tier"] = "EXACT_16_DAE"
            exact_rows.append(row)
    write_rows(calibration_file, exact_rows, mode="w")
    models = fit_models(exact_rows)
    templates = {}
    for portfolio in exact.PORTFOLIO_LABELS.values():
        templates[portfolio] = next(row for row in exact_rows if row["portfolio"] == portfolio)
    write_rows(master, exact_rows, mode="w")
    print(f"TX4_CAMPAIGN_CALIBRATION_COMPLETE rows={len(exact_rows)}", flush=True)

    print(f"TX4_CAMPAIGN_H4_REMAINDER conditions={len(remaining)} workers={workers}", flush=True)
    h4_rows = run_pool(remaining, exact_h4, workers)
    h4_by_id = {row["condition_id"]: row for row in h4_rows}
    for item in remaining:
        campaign, condition_id, seed, params = item
        h4 = h4_by_id[condition_id]
        group = [h4]
        for portfolio in exact.PORTFOLIO_LABELS.values():
            if portfolio == exact.PORTFOLIO_LABELS[H4]:
                continue
            group.append(surrogate_row(params, item, templates[portfolio], portfolio, models))
        group = exact.attach_condition_summaries(group)
        write_rows(master, group)
    metadata = {
        "status": "COMPLETE_WITH_EXPLICIT_SURROGATE_TIER",
        "conditions": len(all_conditions),
        "rows": len(all_conditions) * 16,
        "exact_calibration_conditions": len(calibration),
        "exact_calibration_rows": len(exact_rows),
        "exact_h4_remainder_conditions": len(remaining),
        "surrogate_proper_subset_rows": len(remaining) * 15,
        "workers": workers,
        "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
        "runtime_s": time.perf_counter() - start,
        "surrogate": "ExtraTreesRegressor, 96 trees, min_samples_leaf=2, seed=20260925",
        "calibration_rule": "nominal + 21-point 1D skeleton + 5x5 subgrid for each preregistered 2D pair + first 64 QMC + first 64 MC",
        "exact_h4": True,
    }
    (results / "TX4_ROBUSTNESS_MASTER_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    run()
