"""Assemble checkpointed exact H4 rows and calibrated proper-subset rows."""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_tx4_robustness as exact  # noqa: E402
import run_tx4_robustness_campaign as campaign  # noqa: E402


BOOL_FIELDS = {"stable_all", "stable_EM", "H4_PRESENT", "EXACT_H4", "NONCOMPOSABLE"}
INT_FIELDS = {"seed", "replaced_count", "states", "proper_stable_count", "proper_count"}


def read_rows(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = dict(raw)
            for key in list(row):
                value = row[key]
                if value == "":
                    row[key] = None
                elif key in BOOL_FIELDS:
                    row[key] = value.strip().lower() == "true"
                elif key in INT_FIELDS:
                    row[key] = int(float(value))
                elif key not in {"campaign", "condition_id", "portfolio", "status", "error_type", "error_message", "evaluator_tier"}:
                    try:
                        row[key] = float(value)
                    except (TypeError, ValueError):
                        pass
            rows.append(row)
    return rows


def write_rows(path: Path, rows: list[dict], mode: str = "a") -> None:
    exists = path.exists() and mode == "a"
    with path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=campaign.OUTPUT_FIELDS, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def main():
    results = ROOT / "results"
    calibration_path = results / "TX4_ROBUSTNESS_EXACT_CALIBRATION.csv"
    h4_path = results / "TX4_ROBUSTNESS_H4_CHECKPOINTS.csv"
    master = results / "TX4_ROBUSTNESS_MASTER_CONDITIONS.csv"
    calibration_rows = read_rows(calibration_path)
    h4_rows = read_rows(h4_path)
    h4_by_id = {row["condition_id"]: row for row in h4_rows}
    if len(h4_by_id) != 20164:
        raise RuntimeError(f"expected 20164 unique H4 rows, found {len(h4_by_id)}")
    models = campaign.fit_models(calibration_rows)
    templates = {portfolio: next(row for row in calibration_rows if row["portfolio"] == portfolio) for portfolio in exact.PORTFOLIO_LABELS.values()}
    all_conditions = exact.make_conditions()
    calibration_ids = {row["condition_id"] for row in calibration_rows[::16]}
    if len(calibration_rows) != 6816 or len(calibration_ids) != 426:
        raise RuntimeError(f"unexpected calibration dimensions rows={len(calibration_rows)} conditions={len(calibration_ids)}")
    if master.exists():
        master.unlink()
    remaining = [item for item in all_conditions if item[1] not in calibration_ids]
    with master.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=campaign.OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(calibration_rows)
        for start in range(0, len(remaining), 500):
            batch = remaining[start : start + 500]
            x = np.asarray([[float(item[3][name]) for name in exact.PARAM_NAMES] for item in batch], dtype=float)
            predictions = {}
            for portfolio in exact.PORTFOLIO_LABELS.values():
                predictions[portfolio] = {metric: models[portfolio][metric].predict(x) for metric in campaign.SURROGATE_METRICS}
            for index, item in enumerate(batch):
                campaign_name, condition_id, seed, params = item
                group = [h4_by_id[condition_id]]
                for portfolio in exact.PORTFOLIO_LABELS.values():
                    if portfolio == exact.PORTFOLIO_LABELS[(30, 33, 35, 37)]:
                        continue
                    row = dict(templates[portfolio])
                    row.update({
                        "campaign": campaign_name,
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
                    })
                    for metric in campaign.SURROGATE_METRICS:
                        row[metric] = float(predictions[portfolio][metric][index])
                    row["stable_all"] = bool(row["alpha_all"] <= exact.THRESHOLD)
                    row["stable_EM"] = bool(np.isfinite(row["alpha_EM"]) and row["alpha_EM"] <= exact.THRESHOLD)
                    group.append(row)
                writer.writerows(exact.attach_condition_summaries(group))
    metadata = {
        "status": "COMPLETE_WITH_EXPLICIT_SURROGATE_TIER",
        "conditions": len(all_conditions),
        "rows": len(all_conditions) * 16,
        "exact_all_portfolio_conditions": len(calibration_ids),
        "exact_all_portfolio_rows": len(calibration_rows),
        "exact_h4_conditions": len(h4_by_id) + len(calibration_ids),
        "surrogate_proper_subset_rows": (len(all_conditions) - len(calibration_ids)) * 15,
        "parent": "f64db0004026ceafdb08dd13b5e2ff59d6060742",
        "assembly_runtime_s": time.time(),
        "h4_checkpoint_source": str(h4_path),
        "calibration_source": str(calibration_path),
        "surrogate": "ExtraTreesRegressor, 96 trees, min_samples_leaf=2, seed=20260925",
    }
    (results / "TX4_ROBUSTNESS_MASTER_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
