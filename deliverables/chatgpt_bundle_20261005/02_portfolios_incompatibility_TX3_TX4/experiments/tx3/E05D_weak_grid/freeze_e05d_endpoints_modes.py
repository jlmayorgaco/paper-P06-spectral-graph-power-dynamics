from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import andes
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
E05D = ROOT / "experiments" / "tx3" / "E05D_weak_grid"
sys.path[:0] = [str(ROOT), str(E05D)]

from weak_grid import (  # noqa: E402
    CASE,
    GeneralizedModes,
    critical_mode_index,
    damping_ratio,
    eigenvector_hash,
    evaluate_weak_grid_operator,
    generalized_modes,
    track_one_hungarian,
)


ARTIFACT = ROOT / "artifacts" / "tx3" / "E05D_weak_grid"
PREREG = ARTIFACT / "preregistration"
LIMITS = ARTIFACT / "baseline_calibration" / "E05D_BASELINE_STRESS_LIMITS.parquet"
TRACE = ARTIFACT / "baseline_calibration" / "E05D_BASELINE_CALIBRATION_TRACE.parquet"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def state_names() -> list[str]:
    system = andes.load(str(CASE), setup=True, no_output=True)
    if not system.PFlow.run() or not system.PFlow.converged:
        raise RuntimeError("state-name reference PF failed")
    system.TDS.init()
    return [str(value) for value in system.dae.x_name]


def build(seed: dict[str, Any], kappa: float, residual_max: float) -> tuple[Any, GeneralizedModes]:
    point = evaluate_weak_grid_operator(seed, kappa, {})
    if point.status != "SUCCESS":
        raise RuntimeError(f"{seed['seed_id']} empty kappa={kappa}: {point.status}")
    modes = generalized_modes(point)
    if float(np.max(modes.residuals)) > residual_max:
        raise RuntimeError(f"{seed['seed_id']} empty kappa={kappa}: EIGENPAIR_RESIDUAL_FAIL")
    return point, modes


def track_segment(
    seed: dict[str, Any], source_modes: GeneralizedModes, source_index: int,
    source_kappa: float, target_kappa: float, target_modes: GeneralizedModes,
    bmac_min: float, residual_max: float, depth: int = 0,
) -> tuple[int, float, int] | None:
    match = track_one_hungarian(source_modes, target_modes, source_index)
    if match is not None and match[1] >= bmac_min:
        return match[0], match[1], depth
    if depth >= 2:
        return None
    midpoint = 0.5 * (source_kappa + target_kappa)
    _, middle_modes = build(seed, midpoint, residual_max)
    first = track_one_hungarian(source_modes, middle_modes, source_index)
    if first is None or first[1] < bmac_min:
        return None
    second = track_segment(
        seed, middle_modes, first[0], midpoint, target_kappa, target_modes,
        bmac_min, residual_max, depth + 1,
    )
    if second is None:
        return None
    return second[0], min(first[1], second[1]), max(depth + 1, second[2])


def main() -> int:
    if (ARTIFACT / "holdout" / "e05d_connected_externalities.parquet").exists():
        raise RuntimeError("cannot freeze endpoints/modes after coalition outcomes exist")
    limits = pd.read_parquet(LIMITS).sort_values("seed_id")
    if len(limits) != 8 or not limits["status"].eq("SUCCESS").all():
        raise RuntimeError("eight successful baseline-only limits are required")
    seeds = load_json(PREREG / "E05D_HOLDOUT_FREEZE.json")["operating_points"]
    seed_map = {seed["seed_id"]: seed for seed in seeds}
    stress = load_json(PREREG / "E05D_STRESS_RULE_FREEZE.json")
    numerical = load_json(PREREG / "E05D_NUMERICAL_FREEZE.json")
    tau_grid = [float(value) for value in stress["tau_grid"]]
    bmac_min = float(numerical["stress_and_action_tracking_bmac_minimum"])
    residual_max = float(numerical["eigenpair_normalized_residual_maximum"])
    names = state_names()
    critical_rows: list[dict[str, Any]] = []
    tracking_rows: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []
    vectors: dict[str, np.ndarray] = {}
    freeze_records: list[dict[str, Any]] = []
    for limit in limits.itertuples(index=False):
        seed = seed_map[limit.seed_id]
        kappa_end = float(limit.kappa_end)
        built: dict[float, tuple[Any, GeneralizedModes]] = {}
        for tau in reversed(tau_grid):
            kappa = 1.0 + tau * (kappa_end - 1.0)
            built[tau] = build(seed, kappa, residual_max)
        endpoint_point, endpoint_modes = built[1.0]
        endpoint_index = critical_mode_index(endpoint_modes)
        endpoint_value = complex(endpoint_modes.eigenvalues[endpoint_index])
        participation = np.abs(endpoint_modes.left[:, endpoint_index].conj() * endpoint_modes.right[:, endpoint_index])
        participation /= max(float(np.sum(participation)), 1e-30)
        top = np.argsort(participation)[-8:][::-1]
        dominant = [
            {"state_index": int(index), "state_name": names[index], "normalized_participation": float(participation[index])}
            for index in top
        ]
        vectors[f"{limit.seed_id}_left"] = endpoint_modes.left[:, endpoint_index]
        vectors[f"{limit.seed_id}_right"] = endpoint_modes.right[:, endpoint_index]
        critical = {
            "seed_id": limit.seed_id,
            "tau": 1.0,
            "kappa_grid": kappa_end,
            "mode_index": endpoint_index,
            "lambda_real": float(endpoint_value.real),
            "lambda_imag": float(endpoint_value.imag),
            "frequency_hz": float(endpoint_value.imag / (2.0 * math.pi)),
            "damping_ratio": damping_ratio(endpoint_value),
            "kappa_lambda": float(endpoint_modes.condition_numbers[endpoint_index]),
            "nearest_modal_separation_per_s": float(endpoint_modes.nearest_separations[endpoint_index]),
            "eigenpair_residual": float(endpoint_modes.residuals[endpoint_index]),
            "left_right_eigenvectors_sha256": eigenvector_hash(endpoint_modes, endpoint_index),
            "dominant_device_state_participation_json": json.dumps(dominant, sort_keys=True),
            "SCR36": float(endpoint_point.metadata["SCR36"]),
            "SCR37": float(endpoint_point.metadata["SCR37"]),
            "SCR38": float(endpoint_point.metadata["SCR38"]),
            "min_SCR": float(endpoint_point.metadata["min_SCR"]),
            "status": "SUCCESS",
        }
        critical_rows.append(critical)
        current_tau = 1.0
        current_modes = endpoint_modes
        current_index = endpoint_index
        for tau in reversed(tau_grid):
            point, modes = built[tau]
            if tau == 1.0:
                index, score, depth = endpoint_index, 1.0, 0
            else:
                tracked = track_segment(
                    seed, current_modes, current_index,
                    1.0 + current_tau * (kappa_end - 1.0),
                    1.0 + tau * (kappa_end - 1.0), modes,
                    bmac_min, residual_max,
                )
                if tracked is None:
                    raise RuntimeError(f"{limit.seed_id} MODE_TRACK_FAIL at tau={tau}")
                index, score, depth = tracked
            value = complex(modes.eigenvalues[index])
            base = {
                "seed_id": limit.seed_id,
                "tau": tau,
                "kappa_grid": 1.0 + tau * (kappa_end - 1.0),
                "mode_index": int(index),
                "lambda_real": float(value.real),
                "lambda_imag": float(value.imag),
                "frequency_hz": float(value.imag / (2.0 * math.pi)),
                "damping_ratio": damping_ratio(value),
                "BMAC": float(score),
                "stress_refinement_depth": int(depth),
                "eigenpair_residual": float(modes.residuals[index]),
                "kappa_lambda": float(modes.condition_numbers[index]),
                "nearest_modal_separation_per_s": float(modes.nearest_separations[index]),
                "status": "SUCCESS",
                "SCR36": float(point.metadata["SCR36"]),
                "SCR37": float(point.metadata["SCR37"]),
                "SCR38": float(point.metadata["SCR38"]),
                "min_SCR": float(point.metadata["min_SCR"]),
            }
            tracking_rows.append(base)
            grid_rows.append({
                key: base[key] for key in (
                    "seed_id", "tau", "kappa_grid", "SCR36", "SCR37", "SCR38", "min_SCR", "status"
                )
            })
            current_tau, current_modes, current_index = tau, modes, int(index)
        freeze_records.append({**critical, "dominant_device_state_participation": dominant})

    tables = ARTIFACT / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    critical_frame = pd.DataFrame(critical_rows).sort_values("seed_id")
    tracking_frame = pd.DataFrame(tracking_rows).sort_values(["seed_id", "tau"])
    grid_frame = pd.DataFrame(grid_rows).sort_values(["seed_id", "tau"])
    critical_frame.to_parquet(tables / "e05d_critical_modes.parquet", index=False)
    tracking_frame.to_parquet(PREREG / "e05d_empty_mode_tracking_freeze.parquet", index=False)
    grid_frame.to_parquet(tables / "e05d_grid_strength.parquet", index=False)
    vector_path = PREREG / "E05D_CRITICAL_MODE_EIGENVECTORS.npz"
    np.savez_compressed(vector_path, **vectors)
    endpoint_payload = {
        "freeze_id": "TX3-E05D-STRESS-ENDPOINTS-1.0",
        "calibration_population": "EMPTY coalition only",
        "coalition_outcomes_inspected_before_freeze": False,
        "limits_path": str(LIMITS.relative_to(ROOT)).replace("\\", "/"),
        "limits_sha256": sha256(LIMITS),
        "trace_path": str(TRACE.relative_to(ROOT)).replace("\\", "/"),
        "trace_sha256": sha256(TRACE),
        "kappa_end_range": [float(limits["kappa_end"].min()), float(limits["kappa_end"].max())],
        "endpoint_reason_counts": limits["endpoint_reason"].value_counts().to_dict(),
        "limits": limits.to_dict("records"),
    }
    write_json(PREREG / "E05D_STRESS_ENDPOINT_FREEZE.json", endpoint_payload)
    mode_payload = {
        "freeze_id": "TX3-E05D-CRITICAL-MODES-1.0",
        "selection_population": "EMPTY coalition at tau=1 only",
        "coalition_outcomes_inspected_before_freeze": False,
        "selection_rule": "minimum damping-ratio valid positive-imaginary oscillatory mode in 0.1-30 Hz at each seed's empty-coalition tau=1 endpoint",
        "backward_tracking_rule": "global one-to-one Hungarian BMAC, tau=1 to tau=0, recursive half/quarter step only, no manual replacement",
        "tracking_threshold": bmac_min,
        "critical_mode_table": str((tables / "e05d_critical_modes.parquet").relative_to(ROOT)).replace("\\", "/"),
        "critical_mode_table_sha256": sha256(tables / "e05d_critical_modes.parquet"),
        "empty_tracking_table": str((PREREG / "e05d_empty_mode_tracking_freeze.parquet").relative_to(ROOT)).replace("\\", "/"),
        "empty_tracking_table_sha256": sha256(PREREG / "e05d_empty_mode_tracking_freeze.parquet"),
        "eigenvector_archive": str(vector_path.relative_to(ROOT)).replace("\\", "/"),
        "eigenvector_archive_sha256": sha256(vector_path),
        "modes": freeze_records,
    }
    write_json(PREREG / "E05D_CRITICAL_MODE_FREEZE.json", mode_payload)
    print(json.dumps({
        "status": "FROZEN", "seed_count": len(critical_frame),
        "tracking_rows": len(tracking_frame), "minimum_BMAC": float(tracking_frame["BMAC"].min()),
        "kappa_end_range": endpoint_payload["kappa_end_range"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
