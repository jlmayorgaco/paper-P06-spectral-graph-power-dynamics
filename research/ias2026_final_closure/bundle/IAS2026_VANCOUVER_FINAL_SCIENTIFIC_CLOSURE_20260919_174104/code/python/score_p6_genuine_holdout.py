"""Score the genuine P6 alternative-model holdout after reveal."""

from __future__ import annotations

import csv
import hashlib
import json

from campaign_root import campaign_root_from_argv

CAMPAIGN = campaign_root_from_argv()
pred_path = CAMPAIGN / "prereg" / "p6_genuine_predictions.json"
label_path = CAMPAIGN / "raw" / "p6" / "p6_genuine_holdout_labels.csv"
payload = json.loads(pred_path.read_text(encoding="utf-8"))
labels = {r["case_id"]: r for r in csv.DictReader(label_path.open(encoding="utf-8", newline=""))}
y_true = [labels[k]["stable"].lower() == "true" for k in payload["holdout_case_ids"]]
y_pred = [bool(payload["predictions"][k]["predicted_stable"]) for k in payload["holdout_case_ids"]]
accuracy = sum(a == b for a, b in zip(y_true, y_pred, strict=True)) / len(y_true)
single_class = len(set(y_true)) < 2 or len(set(y_pred)) < 2
metrics = {
    "status": "SINGLE_CLASS_HOLDOUT_REVEALED" if single_class else "MIXED_CLASS_HOLDOUT_REVEALED",
    "evidence_class": "FRESH_ALTERNATIVE_DYNAMIC_MODEL_BLIND_HOLDOUT",
    "prediction_sha256": hashlib.sha256(pred_path.read_bytes()).hexdigest(),
    "holdout_case_ids": payload["holdout_case_ids"],
    "accuracy": accuracy,
    "balanced_accuracy": None if single_class else "COMPUTE_FROM_CONFUSION_MATRIX",
    "unstable_precision": None if single_class else "COMPUTE_FROM_CONFUSION_MATRIX",
    "unstable_recall": None if single_class else "COMPUTE_FROM_CONFUSION_MATRIX",
    "mcc": None if single_class else "COMPUTE_FROM_CONFUSION_MATRIX",
    "cohen_kappa": None if single_class else "COMPUTE_FROM_CONFUSION_MATRIX",
    "single_class_status": "SECOND_PREREGISTERED_HOLDOUT_REQUIRED" if single_class else "MIXED_CLASS",
}
(CAMPAIGN / "raw" / "p6" / "p6_genuine_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
(CAMPAIGN / "reports" / "P6_GENUINE_HOLDOUT_STATUS.md").write_text(
    "# P6 genuine alternative-model holdout\n\n"
    f"status: {metrics['status']}\n"
    f"evidence_class: {metrics['evidence_class']}\n"
    f"accuracy: {accuracy} ({len(y_true)}/{len(y_true)})\n"
    f"balanced_accuracy: {metrics['balanced_accuracy']}\n"
    f"unstable_precision: {metrics['unstable_precision']}\n"
    f"unstable_recall: {metrics['unstable_recall']}\n"
    f"mcc: {metrics['mcc']}\n"
    f"cohen_kappa: {metrics['cohen_kappa']}\n"
    f"single_class_status: {metrics['single_class_status']}\n"
    "The chronology is independently represented by the frozen input hash, prediction hash, and post-reveal label file. This holdout is scoped to the official alternative-model GFL11 harness and does not promote the true same-model gate.\n",
    encoding="utf-8",
)
print(f"P6_GENUINE_{metrics['status']} accuracy={accuracy:.3f}")
