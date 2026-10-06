"""Local four-family log-gain LP with a freshly measured actuator Jacobian.

This is a predictor only. Exact characteristic roots and all five events decide
whether its output is feasible. Unmeasured actuator directions are frozen.
"""
from __future__ import annotations

import csv
import json
import math
import sys
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    base_id = sys.argv[1] if len(sys.argv) > 1 else "pending_41ms"
    output_id = sys.argv[2] if len(sys.argv) > 2 else base_id + "_next_lp"
    actuator_minimum = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002003
    coordinate_trust = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
    l1_trust = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
    assert 0 < coordinate_trust <= 0.05 and 0 < l1_trust <= 0.25
    design = tomllib.loads((HERE / "designs" / f"{base_id}.toml").read_text())
    events = tomllib.loads((HERE / "event_validations" / base_id / "Q0_RESULT.toml").read_text())
    assert events["all_five_events_pass"], "full five-event validation required"
    actuator_slack = min(float(e["min_SG_actuator_fraction_slack"]) for e in events["events"])
    probes = rows(HERE / "L3_FROZEN_EVENT_GRADIENT_PROBES.csv")
    additional = HERE / "L3_ADDITIONAL_EVENT_GRADIENT_PROBES.csv"
    if additional.exists():
        probes.extend(rows(additional))
    measured = {r["candidate_id"]: r for r in rows(HERE / "L3_EVENT_ACTUATOR_SENSITIVITY.csv")}
    assert all(p["candidate_id"] in measured for p in probes), "fresh actuator stencil incomplete"
    base_slack = float(measured["pending_41ms"]["actuator_slack"])
    if base_id == "pending_41ms":
        assert abs(base_slack - actuator_slack) < 1e-7
    E = np.zeros(20)
    measured_coordinates = []
    for probe in probes[1:]:
        index = int(probe["bus"]) - 30 + (0 if probe["gain_kind"] == "Kp" else 10)
        E[index] = (float(measured[probe["candidate_id"]]["actuator_slack"]) - base_slack) / float(probe["delta_log_gain"])
        measured_coordinates.append(index)

    expanded = HERE / "evaluations" / base_id / "ROOTS_EXPANDED.csv"
    rootrows = rows(expanded if expanded.exists() else HERE / "evaluations" / base_id / "ROOTS.csv")
    gradients = rows(HERE / "evaluations" / base_id / "GRADIENTS.csv")
    thresholds = np.array([float(r["local_crossing_ms"]) for r in rootrows])
    G = np.zeros((len(rootrows), 20))
    for row in gradients:
        mode = int(row["mode_rank"]) - 1
        if mode >= len(rootrows):
            continue
        i = int(row["bus"]) - 30
        G[mode, i] = float(row["d_margin_ms_d_logKp"])
        G[mode, i + 10] = float(row["d_margin_ms_d_logKi"])

    nominal = np.r_[design["Kp"], design["Ki"]]
    bounds = json.loads((HERE / "PARENT_EXPERIMENT_MANIFEST.json").read_text())["gain_bounds"]
    lower = np.r_[np.full(10, bounds["Kp_min"]), np.full(10, bounds["Ki_min"])]
    upper = np.r_[np.full(10, bounds["Kp_max"]), np.full(10, bounds["Ki_max"])]
    lb = np.maximum(-coordinate_trust, np.log(lower / nominal))
    ub = np.minimum(coordinate_trust, np.log(upper / nominal))
    for index in range(20):
        if index not in measured_coordinates:
            lb[index] = ub[index] = 0.0
    A, b = [], []
    for mode, threshold in enumerate(thresholds):
        row = np.zeros(41)
        row[:20] = -G[mode]
        row[-1] = 1
        A.append(row)
        b.append(threshold)
    row = np.zeros(41)
    row[:20] = -E
    A.append(row)
    b.append(actuator_slack - actuator_minimum)
    for index in range(20):
        row = np.zeros(41)
        row[index], row[index + 20] = 1, -1
        A.append(row)
        b.append(0)
        row = np.zeros(41)
        row[index], row[index + 20] = -1, -1
        A.append(row)
        b.append(0)
    row = np.zeros(41)
    row[20:40] = 1
    A.append(row)
    b.append(l1_trust)
    objective = np.zeros(41)
    objective[-1] = -1
    result = linprog(objective, A_ub=np.asarray(A), b_ub=np.asarray(b),
                     bounds=[(float(lb[i]), float(ub[i])) for i in range(20)] + [(0, l1_trust)] * 20 + [(None, None)],
                     method="highs")
    assert result.success, result.message
    dual_ineq = -np.asarray(result.ineqlin.marginals)
    lower_dual = np.asarray(result.lower.marginals)
    upper_dual = -np.asarray(result.upper.marginals)
    residual = objective + np.asarray(A).T @ dual_ineq - lower_dual + upper_dual
    kkt_rows = []
    for index, threshold in enumerate(thresholds):
        kkt_rows.append(dict(candidate_id=output_id, constraint=f"delayed_mode_{index+1}",
                             multiplier=float(dual_ineq[index]),
                             linear_slack=float(b[index]-np.asarray(A)[index]@result.x),
                             kkt_residual_inf=float(np.max(np.abs(residual))),
                             scope="LOCAL_LP_ONLY_NOT_FULL_NONLINEAR_KKT"))
    index = len(thresholds)
    kkt_rows.append(dict(candidate_id=output_id,constraint="bus16_plus100_actuator",
                         multiplier=float(dual_ineq[index]),
                         linear_slack=float(b[index]-np.asarray(A)[index]@result.x),
                         kkt_residual_inf=float(np.max(np.abs(residual))),
                         scope="LOCAL_LP_ONLY_NOT_FULL_NONLINEAR_KKT"))
    with (HERE / f"KKT_LP_{output_id}.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(kkt_rows[0]));writer.writeheader();writer.writerows(kkt_rows)
    delta = result.x[:20]
    new = nominal * np.exp(delta)
    predicted = float(np.min(thresholds + G @ delta))
    (HERE / "designs" / f"{output_id}.toml").write_text(
        "rho = " + repr(design["rho"]) + "\n"
        + "Kp = " + repr(new[:10].tolist()) + "\n"
        + "Ki = " + repr(new[10:].tolist()) + "\n"
        + f'parent = "{base_id}"\npredictor = "fresh_actuator_four_family_lp"\n', encoding="utf-8")
    with (HERE / f"PREDICTION_{output_id}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["candidate_id", "parent", "base_tau_ms", "predicted_tau_ms", "base_actuator_slack", "predicted_actuator_slack", "actuator_target", "max_coordinate_log_step", "l1_trust", "log_gain_l1", "active_mode_count", "lp_status"])
        writer.writeheader()
        writer.writerow(dict(candidate_id=output_id, parent=base_id, base_tau_ms=float(min(thresholds)),
                             predicted_tau_ms=predicted, base_actuator_slack=actuator_slack,
                             predicted_actuator_slack=float(actuator_slack + E @ delta),
                             actuator_target=actuator_minimum,max_coordinate_log_step=coordinate_trust,
                             l1_trust=l1_trust,log_gain_l1=float(np.sum(np.abs(delta))),
                             active_mode_count=len(rootrows), lp_status=result.message))
    print(output_id, "predicted tau", predicted, "predicted slack", actuator_slack + E @ delta)


if __name__ == "__main__":
    main()
