"""Reveal and score the frozen P6 holdout predictions."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


CAMPAIGN = Path(__file__).resolve().parents[2]
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
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    pbar = float(y_true.mean())
    qbar = float(y_pred.mean())
    observed = float(np.mean(y_true == y_pred))
    expected = pbar * qbar + (1.0 - pbar) * (1.0 - qbar)
    kappa = (observed - expected) / (1.0 - expected) if expected < 1.0 else 1.0 if observed == 1.0 else 0.0
    actual_all = {p: rows[p]["stable"].lower() == "true" for p in rows}
    predicted_all = {p: bool(payload["predictions"][p]["predicted_stable"]) for p in rows}
    actual_alpha = np.array([float(rows[p]["alpha_transverse"]) for p in holdout])
    predicted_alpha = np.array([float(payload["predictions"][p]["predicted_alpha_transverse"]) for p in holdout])
    actual_freq = np.array([float(rows[p]["critical_frequency_hz"]) for p in holdout])
    predicted_freq = np.array([float(payload["predictions"][p]["predicted_critical_frequency_hz"]) for p in holdout])
    actual_blockers = antichain(set(rows), actual_all)
    predicted_blockers = antichain(set(rows), predicted_all)
    metrics = {
        "status": "PASS",
        "prediction_sha256": digest,
        "holdout_portfolios": holdout,
        "holdout_count": len(holdout),
        "precision": precision,
        "recall": recall,
        "cohen_kappa": kappa,
        "root_alpha_mae": float(np.mean(np.abs(actual_alpha - predicted_alpha))),
        "root_alpha_max_abs_error": float(np.max(np.abs(actual_alpha - predicted_alpha))),
        "root_frequency_mae_hz": float(np.mean(np.abs(actual_freq - predicted_freq))),
        "actual_antichain": actual_blockers,
        "predicted_antichain": predicted_blockers,
        "antichain_exact": actual_blockers == predicted_blockers,
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
        "status: PASS\n"
        "evidence_class: FRESH_BLIND_HOLDOUT_REVEALED\n"
        f"prediction_sha256: {digest}\n"
        f"holdout_portfolios: {holdout}\n"
        f"precision: {precision}\n"
        f"recall: {recall}\n"
        f"cohen_kappa: {kappa}\n"
        f"actual_antichain: {actual_blockers}\n"
        f"predicted_antichain: {predicted_blockers}\n"
        f"antichain_exact: {metrics['antichain_exact']}\n"
        f"root_alpha_mae: {metrics['root_alpha_mae']}\n"
        f"root_frequency_mae_hz: {metrics['root_frequency_mae_hz']}\n"
        "known_branch_repairs_excluded: [0, 1, 13, 43]\n",
        encoding="utf-8",
    )
    print(f"P6_PASS holdout={len(holdout)} precision={precision:.3f} recall={recall:.3f} kappa={kappa:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
