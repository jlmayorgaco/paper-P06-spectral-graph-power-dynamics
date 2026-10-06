from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.linalg import eig, solve
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
sys.path[:0] = [str(ROOT), str(E03)]

from campaign_core import (  # noqa: E402
    OperatorPoint,
    action_vector,
    coalition_externality,
    evaluate_operator,
    operator_logdet,
)


@dataclass(frozen=True)
class ModeSet:
    eigenvalues: np.ndarray
    left: np.ndarray
    right: np.ndarray


def nonempty_subsets(coalition: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [
        subset
        for order in range(1, len(coalition) + 1)
        for subset in combinations(tuple(sorted(coalition)), order)
    ]


def all_subsets(coalition: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [(), *nonempty_subsets(coalition)]


def subset_key(subset: Iterable[int]) -> str:
    labels = tuple(subset)
    return "EMPTY" if not labels else "-".join(f"A{value}" for value in labels)


def reduced_state_matrix(point: OperatorPoint) -> np.ndarray:
    if point.status != "SUCCESS" or point.jacobian is None or point.descriptor is None:
        raise ValueError(f"operator point is unavailable: {point.status}")
    n = int(point.metadata["dynamic_state_count"])
    m = int(point.metadata["algebraic_variable_count"])
    jacobian = point.jacobian.toarray()
    fx = jacobian[:n, :n]
    fy = jacobian[:n, n : n + m]
    gx = jacobian[n : n + m, :n]
    gy = jacobian[n : n + m, n : n + m]
    mass = np.asarray(point.descriptor.diagonal()[:n], dtype=float)
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0.0):
        raise FloatingPointError("dynamic descriptor mass is not strictly positive")
    algebraic_response = solve(gy, gx, assume_a="gen", check_finite=True)
    schur = fx - fy @ algebraic_response
    state_matrix = schur / mass[:, None]
    if not np.isfinite(state_matrix).all():
        raise FloatingPointError("reduced state matrix is not finite")
    return state_matrix


def modal_analysis(point: OperatorPoint) -> ModeSet:
    values, left, right = eig(reduced_state_matrix(point), left=True, right=True, check_finite=True)
    if not (np.isfinite(values).all() and np.isfinite(left).all() and np.isfinite(right).all()):
        raise FloatingPointError("modal decomposition is not finite")
    return ModeSet(np.asarray(values), np.asarray(left), np.asarray(right))


def frequency_hz(value: complex) -> float:
    return float(abs(value.imag) / (2.0 * math.pi))


def damping_ratio(value: complex) -> float:
    magnitude = abs(value)
    return float(-value.real / magnitude) if magnitude > 0.0 else math.nan


def oscillatory_indices(modes: ModeSet, band_hz: tuple[float, float]) -> np.ndarray:
    values = modes.eigenvalues
    frequency = values.imag / (2.0 * math.pi)
    return np.flatnonzero(
        np.isfinite(values)
        & (frequency >= band_hz[0])
        & (frequency <= band_hz[1])
    )


def biorthogonal_mac(reference: ModeSet, target: ModeSet) -> np.ndarray:
    cross_forward = np.abs(reference.left.conj().T @ target.right)
    cross_reverse = np.abs(target.left.conj().T @ reference.right).T
    self_reference = np.abs(np.sum(reference.left.conj() * reference.right, axis=0))
    self_target = np.abs(np.sum(target.left.conj() * target.right, axis=0))
    denominator = np.maximum(self_reference[:, None] * self_target[None, :], 1e-30)
    score = np.sqrt(np.maximum(cross_forward * cross_reverse, 0.0) / denominator)
    return np.clip(np.asarray(score, dtype=float), 0.0, 1.0)


def track_modes(
    reference: ModeSet,
    target: ModeSet,
    *,
    reference_band_hz: tuple[float, float],
    target_band_hz: tuple[float, float],
) -> dict[int, tuple[int, float]]:
    reference_indices = oscillatory_indices(reference, reference_band_hz)
    target_indices = oscillatory_indices(target, target_band_hz)
    if not len(reference_indices) or not len(target_indices):
        return {}
    score = biorthogonal_mac(reference, target)[np.ix_(reference_indices, target_indices)]
    rows, columns = linear_sum_assignment(-score)
    return {
        int(reference_indices[row]): (int(target_indices[column]), float(score[row, column]))
        for row, column in zip(rows, columns, strict=True)
    }


def reference_catalog(
    modes: ModeSet,
    state_names: list[str],
    *,
    band_hz: tuple[float, float],
) -> list[dict[str, Any]]:
    indices = oscillatory_indices(modes, band_hz)
    indices = indices[np.argsort([frequency_hz(modes.eigenvalues[index]) for index in indices])]
    records = []
    for sequence, index in enumerate(indices, start=1):
        participation = np.abs(modes.left[:, index].conj() * modes.right[:, index])
        total = float(np.sum(participation))
        normalized = participation / total if total > 0.0 else participation
        top = np.argsort(normalized)[-6:][::-1]
        value = complex(modes.eigenvalues[index])
        records.append(
            {
                "reference_mode_id": f"M{sequence:03d}",
                "reference_mode_index": int(index),
                "eigenvalue_real_per_s": float(value.real),
                "eigenvalue_imag_per_s": float(value.imag),
                "frequency_hz": frequency_hz(value),
                "damping_ratio": damping_ratio(value),
                "top_state_participation": [
                    {
                        "state_index": int(position),
                        "state_name": state_names[position],
                        "normalized_participation": float(normalized[position]),
                    }
                    for position in top
                ],
            }
        )
    return records


def state_names() -> list[str]:
    import andes

    case = ROOT / "systems" / "ieee39_tx3_gfl3" / "andes_case.xlsx"
    system = andes.load(str(case), setup=True, no_output=True)
    if not system.PFlow.run() or not system.PFlow.converged:
        raise RuntimeError("canonical state-name power flow failed")
    system.TDS.init()
    return [str(value) for value in system.dae.x_name]


def build_points(
    op: dict[str, Any], candidates: dict[str, tuple[int, ...]]
) -> dict[tuple[int, ...], OperatorPoint]:
    required = {()}
    for coalition in candidates.values():
        required.update(all_subsets(coalition))
    points: dict[tuple[int, ...], OperatorPoint] = {}
    for subset in sorted(required, key=lambda item: (len(item), item)):
        alpha = action_vector({f"A{value}": 1.0 for value in subset})
        points[subset] = evaluate_operator(op, alpha)
    return points


def analyze_operating_point(
    op: dict[str, Any],
    candidates: dict[str, tuple[int, ...]],
    reference_modes: ModeSet,
    reference_records: list[dict[str, Any]],
    *,
    reference_band_hz: tuple[float, float],
    tracking_band_hz: tuple[float, float],
    density_frequencies_hz: np.ndarray,
    sigma: float,
) -> dict[str, Any]:
    points = build_points(op, candidates)
    failures = {
        subset_key(subset): point.status
        for subset, point in points.items()
        if point.status != "SUCCESS"
    }
    if failures:
        return {
            "operating_point_id": op["operating_point_id"],
            "split": op["split"],
            "status": "OPERATOR_FAIL",
            "failures": failures,
            "operator_builds": len(points),
        }
    mode_sets = {subset: modal_analysis(point) for subset, point in points.items()}
    op_baseline = mode_sets[()]
    reference_to_op = track_modes(
        reference_modes,
        op_baseline,
        reference_band_hz=reference_band_hz,
        target_band_hz=tracking_band_hz,
    )
    vertex_tracking = {
        subset: track_modes(
            op_baseline,
            modes,
            reference_band_hz=tracking_band_hz,
            target_band_hz=tracking_band_hz,
        )
        for subset, modes in mode_sets.items()
        if subset
    }
    record_by_index = {int(record["reference_mode_index"]): record for record in reference_records}
    tracking_rows: list[dict[str, Any]] = []
    consequence_rows: list[dict[str, Any]] = []
    system_rows: list[dict[str, Any]] = []
    density_rows: list[dict[str, Any]] = []

    for subset, modes in mode_sets.items():
        spectral_abscissa = float(np.max(modes.eigenvalues.real))
        system_rows.append(
            {
                "operating_point_id": op["operating_point_id"],
                "split": op["split"],
                "subset_key": subset_key(subset),
                "subset": subset,
                "spectral_abscissa_per_s": spectral_abscissa,
            }
        )

    for candidate_id, coalition in candidates.items():
        subsets = all_subsets(coalition)
        for reference_index, (op_index, cross_op_score) in reference_to_op.items():
            reference_record = record_by_index.get(reference_index)
            if reference_record is None:
                continue
            values: dict[str, dict[tuple[int, ...], float]] = {
                "damping": {},
                "real": {},
                "frequency": {},
                "log_modal_resolvent": {},
            }
            minimum_tracking = float(cross_op_score)
            complete = True
            for subset in subsets:
                if not subset:
                    target_index, vertex_score = op_index, 1.0
                else:
                    match = vertex_tracking[subset].get(op_index)
                    if match is None:
                        complete = False
                        break
                    target_index, vertex_score = match
                mode = complex(mode_sets[subset].eigenvalues[target_index])
                combined_score = min(float(cross_op_score), float(vertex_score))
                minimum_tracking = min(minimum_tracking, combined_score)
                values["damping"][subset] = damping_ratio(mode)
                values["real"][subset] = float(mode.real)
                values["frequency"][subset] = frequency_hz(mode)
                values["log_modal_resolvent"][subset] = float(-math.log(max(abs(mode.real), 1e-12)))
                tracking_rows.append(
                    {
                        "operating_point_id": op["operating_point_id"],
                        "split": op["split"],
                        "candidate_id": candidate_id,
                        "coalition_key": subset_key(coalition),
                        "reference_mode_id": reference_record["reference_mode_id"],
                        "reference_mode_index": reference_index,
                        "subset_key": subset_key(subset),
                        "subset_order": len(subset),
                        "eigenvalue_real_per_s": float(mode.real),
                        "eigenvalue_imag_per_s": float(mode.imag),
                        "frequency_hz": frequency_hz(mode),
                        "damping_ratio": damping_ratio(mode),
                        "log_modal_resolvent": values["log_modal_resolvent"][subset],
                        "cross_operating_point_bmac": float(cross_op_score),
                        "vertex_bmac": float(vertex_score),
                        "combined_bmac": combined_score,
                    }
                )
            if not complete:
                continue
            full = tuple(sorted(coalition))
            consequence_rows.append(
                {
                    "operating_point_id": op["operating_point_id"],
                    "split": op["split"],
                    "candidate_id": candidate_id,
                    "coalition_key": subset_key(coalition),
                    "order": len(coalition),
                    "reference_mode_id": reference_record["reference_mode_id"],
                    "reference_frequency_hz": reference_record["frequency_hz"],
                    "minimum_tracking_bmac": minimum_tracking,
                    "damping_externality": coalition_externality(values["damping"], coalition),
                    "real_part_externality_per_s": coalition_externality(values["real"], coalition),
                    "frequency_externality_hz": coalition_externality(values["frequency"], coalition),
                    "log_modal_resolvent_externality": coalition_externality(values["log_modal_resolvent"], coalition),
                    "joint_damping_change": values["damping"][full] - values["damping"][()],
                    "joint_real_part_change_per_s": values["real"][full] - values["real"][()],
                    "baseline_damping_ratio": values["damping"][()],
                    "joint_damping_ratio": values["damping"][full],
                    "baseline_frequency_hz": values["frequency"][()],
                    "joint_frequency_hz": values["frequency"][full],
                    "minimum_vertex_damping_ratio": min(values["damping"].values()),
                }
            )

        logdet = {}
        minimum_pivot = math.inf
        for subset in subsets:
            values, pivot = operator_logdet(points[subset], density_frequencies_hz, sigma)
            logdet[subset] = values
            minimum_pivot = min(minimum_pivot, pivot)
        density = np.zeros_like(density_frequencies_hz, dtype=float)
        for subset in subsets:
            density += (-1.0) ** (len(coalition) - len(subset)) * logdet[subset]
        log_frequency = np.log(density_frequencies_hz)
        integrated = float(np.trapezoid(density, log_frequency) / (log_frequency[-1] - log_frequency[0]))
        peak_index = int(np.argmax(np.abs(density)))
        for index, (frequency, value) in enumerate(zip(density_frequencies_hz, density, strict=True)):
            density_rows.append(
                {
                    "operating_point_id": op["operating_point_id"],
                    "split": op["split"],
                    "candidate_id": candidate_id,
                    "coalition_key": subset_key(coalition),
                    "frequency_hz": float(frequency),
                    "externality_density": float(value),
                    "integrated_externality": integrated,
                    "absolute_peak_frequency_hz": float(density_frequencies_hz[peak_index]),
                    "minimum_pivot_ratio": minimum_pivot,
                    "grid_index": index,
                }
            )

    return {
        "operating_point_id": op["operating_point_id"],
        "split": op["split"],
        "status": "SUCCESS",
        "operator_builds": len(points),
        "tracking_rows": tracking_rows,
        "consequence_rows": consequence_rows,
        "system_rows": system_rows,
        "density_rows": density_rows,
    }
