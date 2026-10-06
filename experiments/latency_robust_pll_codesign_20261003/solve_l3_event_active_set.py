"""Sparse event-aware LP around N1, using measured bus16 actuator slopes."""

from __future__ import annotations

import csv
import json
import math
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


HERE = Path(__file__).resolve().parent
BASE = "N_step1_multimode_lp"


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    probes = table(HERE / "L3_FROZEN_EVENT_GRADIENT_PROBES.csv")
    observed = {r["candidate_id"]: r for r in table(HERE / "L3_EVENT_ACTUATOR_SENSITIVITY.csv")}
    assert all(r["candidate_id"] in observed and observed[r["candidate_id"]]["complete"] == "true" for r in probes)
    base_slack = float(observed[BASE]["actuator_slack"])
    source = tomllib.loads((HERE / "designs" / f"{BASE}.toml").read_text())
    roots = table(HERE / "evaluations" / BASE / "ROOTS.csv")[:3]
    gradients = table(HERE / "evaluations" / BASE / "GRADIENTS.csv")
    taus = np.array([float(r["local_crossing_ms"]) for r in roots])
    G = np.zeros((len(roots), 20))
    for row in gradients:
        mode = int(row["mode_rank"]) - 1
        if mode >= len(roots):
            continue
        bus = int(row["bus"]) - 30
        G[mode, bus] = float(row["d_margin_ms_d_logKp"])
        G[mode, 10 + bus] = float(row["d_margin_ms_d_logKi"])
    E = np.zeros(20)
    allowed = []
    for row in probes[1:]:
        i = int(row["bus"]) - 30 + (0 if row["gain_kind"] == "Kp" else 10)
        E[i] = (float(observed[row["candidate_id"]]["actuator_slack"]) - base_slack) / float(row["delta_log_gain"])
        allowed.append(i)
    base = np.r_[source["Kp"], source["Ki"]]
    manifest = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())
    bounds = manifest["gain_bounds"]
    low = np.r_[np.repeat(bounds["Kp_min"], 10), np.repeat(bounds["Ki_min"], 10)]
    high = np.r_[np.repeat(bounds["Kp_max"], 10), np.repeat(bounds["Ki_max"], 10)]
    lb = np.maximum(-0.01, np.log(low / base))
    ub = np.minimum(0.01, np.log(high / base))
    for j in range(20):
        if j not in allowed:
            lb[j] = ub[j] = 0
    guard = 3e-6
    objective = np.r_[np.zeros(40), -1.0]
    A, b = [], []
    for m in range(len(roots)):
        a = np.zeros(41)
        a[:20] = -G[m]
        a[-1] = 1
        A.append(a)
        b.append(taus[m])
    event = np.zeros(41)
    event[:20] = -E
    A.append(event)
    b.append(base_slack - 0.002 - guard)
    for j in range(20):
        a, c = np.zeros(41), np.zeros(41)
        a[j], a[20 + j] = 1, -1
        c[j], c[20 + j] = -1, -1
        A.extend([a, c])
        b.extend([0, 0])
    budget = np.zeros(41)
    budget[20:40] = 1
    A.append(budget)
    b.append(0.05)
    result = linprog(
        objective,
        A_ub=np.asarray(A),
        b_ub=np.asarray(b),
        bounds=[(float(lb[j]), float(ub[j])) for j in range(20)] + [(0, 0.05)] * 20 + [(None, None)],
        method="highs",
    )
    assert result.success, result.message
    du = result.x[:20]
    predicted_margin = float(min(taus + G @ du))
    predicted_slack = float(base_slack + E @ du)
    name = f"{BASE}_event_active_lp"
    kp, ki = (base * np.exp(du))[:10], (base * np.exp(du))[10:]
    (HERE / "designs" / f"{name}.toml").write_text(
        "rho = [" + ", ".join(repr(float(x)) for x in source["rho"]) + "]\n"
        + "Kp = [" + ", ".join(repr(float(x)) for x in kp) + "]\n"
        + "Ki = [" + ", ".join(repr(float(x)) for x in ki) + "]\n"
        + f'parent = "{BASE}"\npredictor = "sparse_event_active_set_lp"\n',
        encoding="utf-8",
    )
    rows = [
        {
            "candidate_id": name,
            "base_actuator_slack": base_slack,
            "guard_above_limit": guard,
            "predicted_actuator_slack": predicted_slack,
            "base_crossing_ms": float(min(taus)),
            "predicted_crossing_ms": predicted_margin,
            "predicted_improvement_ms": predicted_margin - float(min(taus)),
            "allowed_gain_coordinates": len(allowed),
            "log_step_l1": float(sum(abs(du))),
            "event_constraint_linear_residual": max(0.0, 0.002 + guard - predicted_slack),
            "status": "SPARSE_ONE_SIDED_EVENT_GRADIENT_PREDICTOR_NEEDS_EXACT_EVENT_VALIDATION",
        }
    ]
    with (HERE / "L3_EVENT_ACTIVE_SET_PREDICTION.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (HERE / "L3_EVENT_ACTIVE_SET_GRADIENTS.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["gain_coordinate", "bus", "d_slack_d_log_gain", "chosen_delta_log_gain"])
        for j in allowed:
            writer.writerow(["Kp" if j < 10 else "Ki", 30 + j % 10, E[j], du[j]])
    print(name, "predicted_margin", predicted_margin, "predicted_slack", predicted_slack)


if __name__ == "__main__":
    main()
