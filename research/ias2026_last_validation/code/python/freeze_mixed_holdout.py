"""Freeze and reveal the second-model mixed policy holdout.

The prediction rule is frozen before labels are copied into the reveal file:
the prior same-model V4 result predicts the proper three-target portfolio
stable and the four-target portfolio unstable, independent of second-model
policy coordinate.  If the second-model search contains no mixed class, the
script records that negative result and does not manufacture a holdout.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from campaign_root import campaign_root_from_argv

CAMPAIGN = campaign_root_from_argv()
RAW = CAMPAIGN / "raw"
SEARCH = RAW / "second_model_search" / "simplegfldc_policy_discovery.csv"
OUT = RAW / "holdout_mixed"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if not SEARCH.exists():
        raise FileNotFoundError(SEARCH)
    with SEARCH.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    passed = [row for row in rows if row.get("initialization_status") == "PASS"]
    labels = {row["stable"].strip().lower() == "true" for row in passed}
    mixed = len(labels) == 2
    status = {
        "status": "MIXED_HOLDOUT_FROZEN" if mixed else "STOPPED_NO_SECOND_MODEL_MIXED_REGION",
        "source": "raw/second_model_search/simplegfldc_policy_discovery.csv",
        "policy_search_rows": len(rows),
        "initialized_rows": len(passed),
        "classes_observed": sorted("stable" if value else "unstable" for value in labels),
        "prediction_rule": "same-model-V4 topology prior: proper three-target stable; four-target unstable, independent of second-model bandwidth coordinate",
    }
    if not mixed:
        (OUT / "holdout_status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        (OUT / "holdout_status.md").write_text(
            "# Mixed-class holdout\n\n"
            f"status: `{status['status']}`\n\n"
            "The second-model policy search did not produce both stable and unstable labels on its frozen representative grid. No artificial mixed holdout is claimed.\n",
            encoding="utf-8",
        )
        print(status["status"])
        return 0

    # Prefer rows from mixed policy coordinates, then fill to at least 12.
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in passed:
        groups.setdefault((row["axis"], row["scale"]), []).append(row)
    mixed_rows = [
        row
        for group in groups.values()
        if {item["stable"].strip().lower() == "true" for item in group} == {True, False}
        for row in group
    ]
    candidates = mixed_rows + [row for row in passed if row not in mixed_rows]
    selected: list[dict[str, str]] = []
    for row in candidates:
        if len(selected) >= 12:
            break
        selected.append(row)
    if len(selected) < 12:
        raise RuntimeError(f"mixed search produced only {len(selected)} holdout rows")

    inputs = [
        {
            "case_id": f"H{index:02d}",
            "axis": row["axis"],
            "scale": row["scale"],
            "portfolio": row["portfolio"],
            "model": row["model"],
        }
        for index, row in enumerate(selected, start=1)
    ]
    input_path = OUT / "holdout_inputs.csv"
    write_csv(input_path, list(inputs[0]), inputs)

    predictions = []
    for row in inputs:
        predictions.append(
            {
                "case_id": row["case_id"],
                "predicted_stable": row["portfolio"] != "30+33+35+37",
                "rule_frozen_before_reveal": True,
            }
        )
    prediction_path = OUT / "holdout_predictions.csv"
    write_csv(prediction_path, list(predictions[0]), predictions)
    prediction_hash = sha256(prediction_path)
    (OUT / "prediction_hash.txt").write_text(prediction_hash + "\n", encoding="utf-8")

    by_key = {(row["axis"], row["scale"], row["portfolio"]): row for row in passed}
    reveals = []
    for item in inputs:
        source = by_key[(item["axis"], item["scale"], item["portfolio"])]
        reveals.append(
            {
                "case_id": item["case_id"],
                "stable": source["stable"].strip().lower() == "true",
                "alpha_transverse": source["alpha_transverse"],
                "critical_frequency_hz": source["critical_frequency_hz"],
                "source_key": f"{item['axis']}@{item['scale']}:{item['portfolio']}",
            }
        )
    reveal_path = OUT / "holdout_reveal.csv"
    write_csv(reveal_path, list(reveals[0]), reveals)

    pred_by_id = {row["case_id"]: row["predicted_stable"] for row in predictions}
    truth = [row["stable"] for row in reveals]
    predicted = [pred_by_id[row["case_id"]] for row in reveals]
    tp = sum(p and t for p, t in zip(predicted, truth))
    tn = sum((not p) and (not t) for p, t in zip(predicted, truth))
    fp = sum(p and (not t) for p, t in zip(predicted, truth))
    fn = sum((not p) and t for p, t in zip(predicted, truth))
    accuracy = (tp + tn) / len(truth)
    tpr = tp / (tp + fn) if tp + fn else float("nan")
    tnr = tn / (tn + fp) if tn + fp else float("nan")
    metrics = {
        **status,
        "status": "MIXED_HOLDOUT_REVEALED",
        "holdout_size": len(truth),
        "stable_count": sum(truth),
        "unstable_count": len(truth) - sum(truth),
        "accuracy": accuracy,
        "balanced_accuracy": (tpr + tnr) / 2,
        "confusion": {"tp": tp, "tn": tn, "fp": fp, "fn": fn},
        "prediction_sha256": prediction_hash,
        "revealed_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "reveal_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    (OUT / "holdout_status.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    (OUT / "holdout_status.md").write_text(
        "# Mixed-class holdout\n\n"
        f"status: `{metrics['status']}`\n"
        f"size: {metrics['holdout_size']}\n"
        f"stable/unstable: {metrics['stable_count']}/{metrics['unstable_count']}\n"
        f"accuracy: {metrics['accuracy']:.6f}\n"
        f"balanced_accuracy: {metrics['balanced_accuracy']:.6f}\n"
        f"prediction_sha256: `{prediction_hash}`\n\n"
        "Inputs and predictions were written before the reveal file; the rule is deliberately a frozen topology prior, not a refit to second-model labels.\n",
        encoding="utf-8",
    )
    print(f"{metrics['status']} n={metrics['holdout_size']} accuracy={metrics['accuracy']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
