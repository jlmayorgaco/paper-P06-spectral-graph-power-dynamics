"""Propose a bounded multi-root analytical step from an evaluated design."""

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


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_toml(path: Path, rho: list[float], kp: np.ndarray, ki: np.ndarray, parent: str, scale: str) -> None:
    text = [
        "rho = [" + ", ".join(repr(float(x)) for x in rho) + "]",
        "Kp = [" + ", ".join(repr(float(x)) for x in kp) + "]",
        "Ki = [" + ", ".join(repr(float(x)) for x in ki) + "]",
        f'parent = "{parent}"',
        f'line_search_scale = "{scale}"',
    ]
    path.write_text("\n".join(text) + "\n")


def main() -> None:
    assert len(sys.argv) == 2, "usage: propose_next.py current_design.toml"
    path = Path(sys.argv[1]).resolve()
    name = path.stem
    data = tomllib.loads(path.read_text())
    folder = HERE / "evaluations" / name
    roots = table(folder / "ROOTS.csv")[:3]
    grads = table(folder / "GRADIENTS.csv")
    assert len(roots) >= 2
    taus = np.array([float(r["local_crossing_ms"]) for r in roots])
    G = np.zeros((len(roots), 20))
    for mode in range(1, len(roots) + 1):
        relevant = [r for r in grads if int(r["mode_rank"]) == mode]
        assert len(relevant) == 10
        for row in relevant:
            i = int(row["bus"]) - 30
            G[mode - 1, i] = float(row["d_margin_ms_d_logKp"])
            G[mode - 1, i + 10] = float(row["d_margin_ms_d_logKi"])
    base = np.r_[data["Kp"], data["Ki"]]
    manifest = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())
    bounds = manifest["gain_bounds"]
    low = np.r_[np.repeat(bounds["Kp_min"], 10), np.repeat(bounds["Ki_min"], 10)]
    high = np.r_[np.repeat(bounds["Kp_max"], 10), np.repeat(bounds["Ki_max"], 10)]
    lb = np.maximum(-0.01, np.log(low / base))
    ub = np.minimum(0.01, np.log(high / base))
    # x=(du[20], absolute auxiliaries[20], common crossing t)
    objective = np.r_[np.zeros(40), -1.0]
    A, b = [], []
    for mode in range(len(roots)):
        a = np.zeros(41)
        a[:20] = -G[mode]
        a[-1] = 1
        A.append(a)
        b.append(taus[mode])
    for j in range(20):
        a, c = np.zeros(41), np.zeros(41)
        a[j], a[20 + j] = 1, -1
        c[j], c[20 + j] = -1, -1
        A += [a, c]
        b += [0, 0]
    a = np.zeros(41)
    a[20:40] = 1
    A.append(a)
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
    out = []
    for label, scale in (("full", 1.0), ("half", 0.5), ("quarter", 0.25)):
        vector = base * np.exp(scale * du)
        candidate_id = f"{name}_next_{label}"
        write_toml(HERE / "designs" / f"{candidate_id}.toml", data["rho"], vector[:10], vector[10:], name, label)
        prediction = float(np.min(taus + G @ (scale * du)))
        out.append(
            {
                "candidate_id": candidate_id,
                "parent": name,
                "line_search_scale": scale,
                "current_local_crossing_ms": float(min(taus)),
                "predicted_crossing_ms": prediction,
                "predicted_improvement_ms": prediction - float(min(taus)),
                "max_abs_log_gain_step": float(max(abs(scale * du))),
                "sum_abs_log_gain_step": float(sum(abs(scale * du))),
                "modes_in_LP": len(roots),
                "status": "ANALYTICAL_PREDICTION_NEEDS_EXACT_CORRECTOR_AND_EVENTS",
            }
        )
    with (folder / "PROPOSALS.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(out[0]))
        writer.writeheader()
        writer.writerows(out)
    for row in out:
        print(row["candidate_id"], row["predicted_improvement_ms"])


if __name__ == "__main__":
    main()
