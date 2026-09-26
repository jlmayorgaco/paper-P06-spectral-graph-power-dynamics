"""Reveal and score the frozen P6 holdout predictions."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from campaign_root import campaign_root_from_argv

CAMPAIGN = campaign_root_from_argv()
PREREG = CAMPAIGN / "prereg" / "blind_predictions_P6.json"
P2 = CAMPAIGN / "raw" / "p2" / "p2_powerdynamics_portfolios.csv"
RAW = CAMPAIGN / "raw" / "p6"


def antichain(portfolios: set[str], stable: dict[str, bool]) -> list[str]:
    blockers = []
    for p in portfolios:
        if stable[p]:
            continue
        ps = set(p.split("+")) if p else set()
        if not any((set(q.split("+")) if q else set()) < ps and not stable[q] for q in portfolios):
            blockers.append(p)
    return sorted(blockers)


def main() -> int:
    payload = json.loads(PREREG.read_text(encoding="utf-8"))
    digest = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    rows = {r["portfolio"]: r for r in csv.DictReader(P2.open(encoding="utf-8", newline=""))}
    holdout = payload["holdout_portfolios"]
    y_true = np.array([rows[p]["stable"].lower() == "true" for p in holdout], dtype=bool)
    y_pred = np.array([bool(payload["predictions"][p]["predicted_stable"]) for p in holdout], dtype=bool)
    tp = int(np.count_nonzero(y_true & y_pred))
    fp = int(np.count_nonzero(~y_true & y_pred))
    fn = int(np.count_nonzero(y_true & ~y_pred))
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    pbar = float(y_true.mean())
    qbar = float(y_pred.mean())
    observed = float(np.mean(y_true == y_pred))
    expected = pbar * qbar + (1.0 - pbar) * (1.0 - qbar)
    single_class = len(set(y_true.tolist())) < 2 or len(set(y_pred.tolist())) < 2
    if single_class:
        precision = None
        recall = None
    kappa = None if single_class else (observed - expected) / (1.0 - expected)
    actual_all = {p: rows[p]["stable"].lower() == "true" for p in rows}
    predicted_all = {p: bool(payload["predictions"][p]["predicted_stable"]) for p in rows}
    actual_alpha = np.array([float(rows[p]["alpha_transverse"]) for p in holdout])
    predicted_alpha = np.array([float(payload["predictions"][p]["predicted_alpha_transverse"]) for p in holdout])
    actual_freq = np.array([float(rows[p]["critical_frequency_hz"]) for p in holdout])
    predicted_freq = np.array([float(payload["predictions"][p]["predicted_critical_frequency_hz"]) for p in holdout])
    actual_blockers = antichain(set(rows), actual_all)
    predicted_blockers = antichain(set(rows), predicted_all)
    metrics = {
        "status": "POSTHOC_OR_PROCEDURAL_SINGLE_CLASS_HOLDOUT",
        "evidence_class": "POSTHOC_OR_PROCEDURAL_SINGLE_CLASS_HOLDOUT",
        "prediction_sha256": digest,
        "holdout_portfolios": holdout,
        "holdout_count": len(holdout),
        "accuracy": float(np.mean(y_true == y_pred)),
        "precision": precision,
        "recall": recall,
        "precision_recall_status": "NOT_MEANINGFUL_SINGLE_CLASS",
        "cohen_kappa": kappa,
        "cohen_kappa_status": "UNDEFINED_SINGLE_CLASS_HOLDOUT" if single_class else "DEFINED",
        "root_alpha_mae": float(np.mean(np.abs(actual_alpha - predicted_alpha))),
        "root_alpha_max_abs_error": float(np.max(np.abs(actual_alpha - predicted_alpha))),
        "root_frequency_mae_hz": float(np.mean(np.abs(actual_freq - predicted_freq))),
        "actual_antichain": actual_blockers,
        "predicted_antichain": predicted_blockers,
        "antichain_exact": None,
        "blocker_discrimination_status": "NOT_TESTED_SINGLE_CLASS_HOLDOUT",
        "repair_ranking_status": "NOT_TESTED_SINGLE_CLASS_HOLDOUT",
        "chronology_status": "NOT_INDEPENDENTLY_PROVABLE_FROM_BUNDLE",
        "known_branch_repairs_excluded": [0, 1, 13, 43],
    }
    with (RAW / "p6_revealed_holdout.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["portfolio", "split", "actual_stable", "predicted_stable", "actual_alpha_transverse", "predicted_alpha_transverse", "actual_critical_frequency_hz", "predicted_critical_frequency_hz"])
        for p in holdout:
            writer.writerow([p, "holdout", rows[p]["stable"], payload["predictions"][p]["predicted_stable"], rows[p]["alpha_transverse"], payload["predictions"][p]["predicted_alpha_transverse"], rows[p]["critical_frequency_hz"], payload["predictions"][p]["predicted_critical_frequency_hz"]])
    (RAW / "p6_reveal_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    report = CAMPAIGN / "reports" / "P6_BLIND_HOLDOUT_STATUS.md"
    report.write_text(
        "# P6 — blind holdout\n\n"
        "status: POSTHOC_OR_PROCEDURAL_SINGLE_CLASS_HOLDOUT\n"
        "evidence_class: POSTHOC_OR_PROCEDURAL_SINGLE_CLASS_HOLDOUT\n"
        f"prediction_sha256: {digest}\n"
        f"holdout_portfolios: {holdout}\n"
        f"accuracy: {metrics['accuracy']} (4/4)\n"
        "precision: NOT_MEANINGFUL_SINGLE_CLASS\n"
        "recall: NOT_MEANINGFUL_SINGLE_CLASS\n"
        "cohen_kappa: UNDEFINED_SINGLE_CLASS_HOLDOUT\n"
        f"actual_antichain: {actual_blockers}\n"
        f"predicted_antichain: {predicted_blockers}\n"
        "antichain_exact: NOT_INTERPRETABLE_SINGLE_CLASS\n"
        "blocker_discrimination: NOT_TESTED_SINGLE_CLASS_HOLDOUT\n"
        "repair_ranking: NOT_TESTED_SINGLE_CLASS_HOLDOUT\n"
        "chronology: NOT_INDEPENDENTLY_PROVABLE_FROM_BUNDLE\n"
        f"root_alpha_mae: {metrics['root_alpha_mae']}\n"
        f"root_frequency_mae_hz: {metrics['root_frequency_mae_hz']}\n"
        "known_branch_repairs_excluded: [0, 1, 13, 43]\n",
        encoding="utf-8",
    )
    print(f"P6_{metrics['status']} holdout={len(holdout)} accuracy={metrics['accuracy']:.3f} kappa=undefined")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
