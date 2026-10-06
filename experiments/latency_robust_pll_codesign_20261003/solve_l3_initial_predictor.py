"""Analytical ball step and two-mode log-gain LP at each frozen start."""

from __future__ import annotations

import csv
import json
import math
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


HERE = Path(__file__).resolve().parent


def data(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def save_toml(path: Path, rho: list[float], kp: np.ndarray, ki: np.ndarray, metadata: dict[str, object]) -> None:
    content = [
        "rho = [" + ", ".join(repr(float(x)) for x in rho) + "]",
        "Kp = [" + ", ".join(repr(float(x)) for x in kp) + "]",
        "Ki = [" + ", ".join(repr(float(x)) for x in ki) + "]",
    ]
    content += [f'{key} = "{value}"' for key, value in metadata.items()]
    path.write_text("\n".join(content) + "\n", encoding="utf-8")


def main() -> None:
    manifest = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())
    grads = data("L2_MULTIMODE_GRADIENTS.csv")
    scan = data("L2_FAMILY_SCAN.csv")
    output = []
    for design_id, eta, filename in [
        ("Z", 0.0, "Z_zero_delay_tuned.toml"),
        ("N", 1.0, "N_nominal.toml"),
    ]:
        design = tomllib.loads((HERE / "designs" / filename).read_text())
        active = sorted(
            (r for r in scan if abs(float(r["eta"]) - eta) < 1e-8),
            key=lambda r: int(r["rank_by_local_crossing"]),
        )[:2]
        assert len(active) == 2
        taus = np.array([float(r["local_crossing_ms"]) for r in active])
        G = np.zeros((2, 20))
        for mode in (1, 2):
            mode_rows = [r for r in grads if abs(float(r["eta"]) - eta) < 1e-8 and int(r["mode_rank"]) == mode]
            assert len(mode_rows) == 10
            for r in mode_rows:
                i = int(r["bus"]) - 30
                G[mode - 1, i] = float(r["d_margin_ms_d_logKp"])
                G[mode - 1, 10 + i] = float(r["d_margin_ms_d_logKi"])
        base = np.r_[np.asarray(design["Kp"], dtype=float), np.asarray(design["Ki"], dtype=float)]
        low = np.r_[np.repeat(manifest["gain_bounds"]["Kp_min"], 10), np.repeat(manifest["gain_bounds"]["Ki_min"], 10)]
        high = np.r_[np.repeat(manifest["gain_bounds"]["Kp_max"], 10), np.repeat(manifest["gain_bounds"]["Ki_max"], 10)]
        lb = np.maximum(-0.01, np.log(low / base))
        ub = np.minimum(0.01, np.log(high / base))
        # x = (20 signed log steps, 20 absolute-value auxiliaries, t).
        objective = np.r_[np.zeros(40), -1.0]
        A = []
        b = []
        for mode in range(2):
            row = np.zeros(41)
            row[:20] = -G[mode]
            row[-1] = 1.0
            A.append(row)
            b.append(taus[mode])
        for j in range(20):
            r1 = np.zeros(41)
            r2 = np.zeros(41)
            r1[j], r1[20 + j] = 1.0, -1.0
            r2[j], r2[20 + j] = -1.0, -1.0
            A.extend((r1, r2))
            b.extend((0.0, 0.0))
        budget = np.zeros(41)
        budget[20:40] = 1.0
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
        step = result.x[:20]
        active_g = G[0]
        ball = 0.01 * active_g / np.linalg.norm(active_g)
        ball = np.clip(ball, lb, ub)
        (HERE / "designs").mkdir(exist_ok=True)
        for method, vec in (("ball", ball), ("multimode_lp", step)):
            gains = base * np.exp(vec)
            candidate_id = f"{design_id}_step1_{method}"
            save_toml(
                HERE / "designs" / f"{candidate_id}.toml",
                design["rho"],
                gains[:10],
                gains[10:],
                {"candidate_id": candidate_id, "parent": design_id, "predictor": method},
            )
            output.append(
                {
                    "candidate_id": candidate_id,
                    "parent": design_id,
                    "method": method,
                    "base_earliest_crossing_ms": float(min(taus)),
                    "predicted_earliest_crossing_ms": float(min(taus + G @ vec)),
                    "predicted_gain_ms": float(min(taus + G @ vec) - min(taus)),
                    "log_step_l1": float(np.linalg.norm(vec, 1)),
                    "log_step_l2": float(np.linalg.norm(vec)),
                    "max_abs_log_step": float(np.max(np.abs(vec))),
                    "active_linear_modes": ";".join(str(i + 1) for i in range(2) if abs(taus[i] + G[i] @ vec - min(taus + G @ vec)) < 1e-6),
                    "step_log_vector": ";".join(repr(float(v)) for v in vec),
                    "status": "ANALYTICAL_PREDICTION_ONLY_NEEDS_EXACT_CORRECTOR",
                }
            )
    with (HERE / "L3_INITIAL_PREDICTORS.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    for row in output:
        print(row["candidate_id"], row["predicted_gain_ms"], row["active_linear_modes"])


if __name__ == "__main__":
    main()
