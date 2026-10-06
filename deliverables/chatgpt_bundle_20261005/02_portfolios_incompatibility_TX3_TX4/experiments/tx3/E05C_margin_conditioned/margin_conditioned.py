from __future__ import annotations

import hashlib
import math
import sys
import time
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import andes
import numpy as np
from kvxopt import matrix as kvx_matrix
from scipy.linalg import eig, solve
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csc_matrix, diags


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
sys.path[:0] = [str(ROOT), str(E03)]

from campaign_core import BASE_GFL, OperatorPoint, _uid, action_vector, coalition_externality, log_frequency_rule  # noqa: E402
from dicgrid.adapters.andes import enforce_gfl_pq_equilibrium  # noqa: E402


CASE = ROOT / "systems" / "ieee39_tx3_gfl3" / "andes_case.xlsx"
SYNCHRONOUS_PV_BUSES = (30, 31, 32, 33, 34, 35)
CANONICAL_SYNCHRONOUS_P = np.asarray([2.5, 5.7293, 6.5, 6.32, 5.08, 6.5], dtype=float)
STRESS_PARTICIPATION = CANONICAL_SYNCHRONOUS_P / np.sum(CANONICAL_SYNCHRONOUS_P)
PSEUDO_MODE_RADIUS = 1e-5


@dataclass(frozen=True)
class GeneralizedModes:
    eigenvalues: np.ndarray
    left: np.ndarray
    right: np.ndarray
    mass: np.ndarray
    residuals: np.ndarray
    condition_numbers: np.ndarray
    nearest_separations: np.ndarray


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def subset_key(subset: Iterable[int]) -> str:
    labels = tuple(subset)
    return "EMPTY" if not labels else "-".join(f"A{value}" for value in labels)


def parse_subset_key(value: str) -> tuple[int, ...]:
    if value in {"", "EMPTY"}:
        return ()
    return tuple(int(label.removeprefix("A")) for label in value.replace(",", "-").split("-"))


def mobius_vertices(coalition: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [
        tuple(coalition[index] for index in range(len(coalition)) if mask & (1 << index))
        for mask in range(1 << len(coalition))
    ]


def union_vertices(coalitions: Iterable[tuple[int, ...]]) -> list[tuple[int, ...]]:
    values: set[tuple[int, ...]] = {()}
    for coalition in coalitions:
        for order in range(1, len(coalition) + 1):
            values.update(combinations(coalition, order))
    return sorted(values, key=lambda value: (len(value), value))


def evaluate_stressed_operator(
    seed: dict[str, Any], gamma: float, alpha: dict[str, float]
) -> OperatorPoint:
    """Build a fully re-equilibrated E05C operator under frozen load stress."""

    started = time.perf_counter()
    alpha = action_vector(alpha)
    try:
        system = andes.load(str(CASE), setup=True, no_output=True)
        base_pq_p = np.asarray(system.PQ.p0.v, dtype=float).copy()
        base_pq_q = np.asarray(system.PQ.q0.v, dtype=float).copy()
        seeded_p = base_pq_p * float(seed["load_scale"])
        seeded_q = base_pq_q * float(seed["load_scale"]) * float(seed["reactive_load_scale"])
        system.PQ.p0.v[:] = seeded_p * float(gamma)
        system.PQ.q0.v[:] = seeded_q * float(gamma)

        redispatch = float(seed["redispatch_mw"]) / 100.0
        pv_buses = [int(value) for value in system.PV.bus.v]
        system.PV.p0.v[_uid(pv_buses, 30)] += redispatch
        system.PV.p0.v[_uid(pv_buses, 35)] -= redispatch
        added_demand = (float(gamma) - 1.0) * float(np.sum(seeded_p))
        for bus, participation in zip(SYNCHRONOUS_PV_BUSES, STRESS_PARTICIPATION, strict=True):
            system.PV.p0.v[_uid(pv_buses, bus)] += added_demand * participation

        targets = {bus: values.copy() for bus, values in BASE_GFL.items()}
        dispatch_scale = float(seed["gfl_dispatch_scale"])
        for target in targets.values():
            target["p_system_pu"] *= dispatch_scale
        for index in range(1, 4):
            ratio = 1.0 + 0.20 * alpha[f"A{index}"]
            system.PLL2.Kp.v[index - 1] *= ratio
            system.PLL2.Ki.v[index - 1] *= ratio * ratio
        for index in range(4, 7):
            ratio = 1.0 + 0.20 * alpha[f"A{index}"]
            system.REGCP1.Tg.v[index - 4] /= ratio
        targets[37]["p_system_pu"] += 0.20 * alpha["A7"]
        for bus, target in targets.items():
            system.PV.p0.v[_uid(pv_buses, bus)] = target["p_system_pu"]
        line45 = _uid([str(value) for value in system.Line.idx.v], "Line_45")
        system.Line.x.v[line45] *= 1.0 + 0.10 * alpha["A8"]

        enforcement = enforce_gfl_pq_equilibrium(system, targets)
        pf_return = bool(system.PFlow.run())
        if not pf_return or not system.PFlow.converged:
            return OperatorPoint(None, None, "PF_FAIL", {"elapsed_s": time.perf_counter() - started})

        voltage = np.asarray(system.Bus.v.v, dtype=float)
        pv_p = np.asarray(system.PV.p.v, dtype=float)
        pv_q = np.asarray(system.PV.q.v, dtype=float)
        slack_p = np.asarray(system.Slack.p.v, dtype=float)
        slack_q = np.asarray(system.Slack.q.v, dtype=float)
        synchronous_indices = np.asarray([_uid(pv_buses, bus) for bus in SYNCHRONOUS_PV_BUSES])
        pv_p_limit = bool(
            np.any(pv_p[synchronous_indices] > np.asarray(system.PV.pmax.v, dtype=float)[synchronous_indices] + 1e-8)
            or np.any(pv_p[synchronous_indices] < np.asarray(system.PV.pmin.v, dtype=float)[synchronous_indices] - 1e-8)
        )
        pv_q_limit = bool(
            np.any(pv_q[synchronous_indices] > np.asarray(system.PV.qmax.v, dtype=float)[synchronous_indices] + 1e-8)
            or np.any(pv_q[synchronous_indices] < np.asarray(system.PV.qmin.v, dtype=float)[synchronous_indices] - 1e-8)
        )
        slack_limit = bool(
            np.any(slack_p > np.asarray(system.Slack.pmax.v, dtype=float) + 1e-8)
            or np.any(slack_p < np.asarray(system.Slack.pmin.v, dtype=float) - 1e-8)
            or np.any(slack_q > np.asarray(system.Slack.qmax.v, dtype=float) + 1e-8)
            or np.any(slack_q < np.asarray(system.Slack.qmin.v, dtype=float) - 1e-8)
        )

        system.TDS.init()
        residual = max(
            float(np.max(np.abs(np.asarray(system.dae.f, dtype=float)))),
            float(np.max(np.abs(np.asarray(system.dae.g, dtype=float)))),
        )
        if not np.isfinite(residual) or residual > 1e-6:
            return OperatorPoint(
                None,
                None,
                "INIT_FAIL",
                {"initialization_residual": residual, "elapsed_s": time.perf_counter() - started},
            )

        device_scale = 100.0 / np.asarray(system.REGCP1.Sn.v, dtype=float)
        ip = np.asarray(system.REGCP1.Ipout.v, dtype=float) * device_scale
        iq = np.asarray(system.REGCP1.Iqout_x.v, dtype=float) * device_scale
        current = np.hypot(ip, iq)
        current_limit = np.asarray(system.REECB1.Imax.v, dtype=float)
        limiter_active = bool(np.any(current >= current_limit - 1e-8))
        saturation_active = bool(
            np.any(np.abs(np.asarray(system.REECB1.PIV_ys.v) - np.asarray(system.REECB1.PIV_y.v)) > 1e-9)
            or np.any(np.abs(np.asarray(system.REECB1.PIQ_ys.v) - np.asarray(system.REECB1.PIQ_y.v)) > 1e-9)
        )

        dense = lambda value: np.array(kvx_matrix(value), dtype=float)
        dae = system.dae
        jacobian = np.block([[dense(dae.fx), dense(dae.fy)], [dense(dae.gx), dense(dae.gy)]])
        mass = np.concatenate((np.asarray(dae.Tf, dtype=float), np.zeros(dae.m)))
        if not np.isfinite(jacobian).all() or not np.isfinite(mass).all():
            return OperatorPoint(None, None, "DAE_FAIL", {"elapsed_s": time.perf_counter() - started})
        metadata = {
            "elapsed_s": time.perf_counter() - started,
            "pf_converged": True,
            "gamma": float(gamma),
            "seeded_active_load_system_pu": float(np.sum(seeded_p)),
            "added_active_load_system_pu": added_demand,
            "scheduled_synchronous_redispatch_system_pu": added_demand,
            "slack_active_power_system_pu": float(slack_p[0]),
            "maximum_synchronous_pv_reactive_power_system_pu": float(np.max(pv_q[synchronous_indices])),
            "initialization_residual": residual,
            "voltage_min_pu": float(np.min(voltage)),
            "voltage_max_pu": float(np.max(voltage)),
            "maximum_gfl_current_device_pu": float(np.max(current)),
            "current_limit_active": limiter_active,
            "controller_saturation_active": saturation_active,
            "generator_active_limit_violation": pv_p_limit,
            "generator_reactive_limit_violation": pv_q_limit,
            "slack_limit_violation": slack_limit,
            "converted_voltage_control_disabled": enforcement["converted_voltage_control_disabled"],
            "dynamic_state_count": int(dae.n),
            "algebraic_variable_count": int(dae.m),
        }
        if limiter_active:
            status = "LIMITER_ACTIVE"
        elif saturation_active:
            status = "NONSMOOTH_REGIME"
        elif pv_p_limit or pv_q_limit or slack_limit:
            status = "LIMITER_ACTIVE"
        else:
            status = "SUCCESS"
        return OperatorPoint(csc_matrix(jacobian), diags(mass, format="csc"), status, metadata)
    except Exception as exc:
        return OperatorPoint(
            None,
            None,
            "DAE_FAIL",
            {"exception": repr(exc), "elapsed_s": time.perf_counter() - started},
        )


def generalized_modes(point: OperatorPoint) -> GeneralizedModes:
    if point.status != "SUCCESS" or point.jacobian is None or point.descriptor is None:
        raise ValueError(f"operator point unavailable: {point.status}")
    n = int(point.metadata["dynamic_state_count"])
    m = int(point.metadata["algebraic_variable_count"])
    jacobian = point.jacobian.toarray()
    fx = jacobian[:n, :n]
    fy = jacobian[:n, n : n + m]
    gx = jacobian[n : n + m, :n]
    gy = jacobian[n : n + m, n : n + m]
    mass = np.asarray(point.descriptor.diagonal()[:n], dtype=float)
    response = solve(gy, gx, assume_a="gen", check_finite=True)
    reduced = fx - fy @ response
    values, left, right = eig(reduced, np.diag(mass), left=True, right=True, check_finite=True)
    if not (np.isfinite(values).all() and np.isfinite(left).all() and np.isfinite(right).all()):
        raise FloatingPointError("non-finite generalized eigensolution")
    norm_a = np.linalg.norm(reduced, 2)
    norm_m = float(np.max(np.abs(mass)))
    residuals = np.empty(len(values), dtype=float)
    conditioning = np.empty(len(values), dtype=float)
    for index, value in enumerate(values):
        r = right[:, index]
        l = left[:, index]
        residuals[index] = np.linalg.norm(reduced @ r - value * mass * r) / (
            max((norm_a + abs(value) * norm_m) * np.linalg.norm(r), 1e-30)
        )
        denominator = abs(np.vdot(l, mass * r))
        conditioning[index] = np.linalg.norm(l) * np.linalg.norm(r) / max(denominator, 1e-30)
    separations = np.asarray(
        [float(np.min(np.abs(np.delete(values, index) - value))) for index, value in enumerate(values)],
        dtype=float,
    )
    return GeneralizedModes(values, left, right, mass, residuals, conditioning, separations)


def physical_indices(modes: GeneralizedModes) -> np.ndarray:
    return np.flatnonzero(np.abs(modes.eigenvalues) > PSEUDO_MODE_RADIUS)


def positive_oscillatory_indices(
    modes: GeneralizedModes, band_hz: tuple[float, float]
) -> np.ndarray:
    frequencies = modes.eigenvalues.imag / (2.0 * math.pi)
    return np.flatnonzero(
        (np.abs(modes.eigenvalues) > PSEUDO_MODE_RADIUS)
        & (frequencies >= band_hz[0])
        & (frequencies <= band_hz[1])
    )


def maximum_physical_real_part(modes: GeneralizedModes) -> float:
    return float(np.max(modes.eigenvalues[physical_indices(modes)].real))


def damping_ratio(value: complex) -> float:
    return float(-value.real / abs(value)) if abs(value) > 0.0 else math.nan


def biorthogonal_mac(reference: GeneralizedModes, target: GeneralizedModes) -> np.ndarray:
    cross_forward = np.abs(reference.left.conj().T @ target.right)
    cross_reverse = np.abs(target.left.conj().T @ reference.right).T
    self_reference = np.abs(np.sum(reference.left.conj() * reference.right, axis=0))
    self_target = np.abs(np.sum(target.left.conj() * target.right, axis=0))
    denominator = np.maximum(self_reference[:, None] * self_target[None, :], 1e-30)
    return np.clip(np.sqrt(np.maximum(cross_forward * cross_reverse, 0.0) / denominator), 0.0, 1.0)


def track_assignment(
    reference: GeneralizedModes,
    target: GeneralizedModes,
    reference_indices: np.ndarray,
    target_band_hz: tuple[float, float] = (0.05, 45.0),
) -> dict[int, tuple[int, float]]:
    target_indices = positive_oscillatory_indices(target, target_band_hz)
    if not len(reference_indices) or not len(target_indices):
        return {}
    score = biorthogonal_mac(reference, target)[np.ix_(reference_indices, target_indices)]
    rows, columns = linear_sum_assignment(-score)
    return {
        int(reference_indices[row]): (int(target_indices[column]), float(score[row, column]))
        for row, column in zip(rows, columns, strict=True)
    }


def track_one(
    reference: GeneralizedModes,
    target: GeneralizedModes,
    reference_index: int,
    target_band_hz: tuple[float, float] = (0.05, 45.0),
) -> tuple[int, float] | None:
    target_indices = positive_oscillatory_indices(target, target_band_hz)
    if not len(target_indices):
        return None
    scores = biorthogonal_mac(reference, target)[reference_index, target_indices]
    location = int(np.argmax(scores))
    return int(target_indices[location]), float(scores[location])


def track_one_hungarian(
    reference: GeneralizedModes,
    target: GeneralizedModes,
    reference_index: int,
    band_hz: tuple[float, float] = (0.05, 45.0),
) -> tuple[int, float] | None:
    """Track one member of a globally one-to-one positive-frequency assignment.

    The global assignment is essential when several branches have exactly tied
    biorthogonal MAC. An independent row-wise argmax can map two reference modes
    to the same target and create a spurious connected eigenvalue shift.
    """

    reference_indices = positive_oscillatory_indices(reference, band_hz)
    return track_assignment(reference, target, reference_indices, band_hz).get(int(reference_index))


def pole_factor(modes: GeneralizedModes, nodes: int = 96) -> float:
    frequencies, weights = log_frequency_rule(0.1, 30.0, nodes)
    s = 2j * math.pi * frequencies
    logdet = np.sum(np.log(np.abs(s[:, None] - modes.eigenvalues[None, :])), axis=1)
    return float(np.dot(weights, logdet))


def complex_externality(values: dict[tuple[int, ...], complex], coalition: tuple[int, ...]) -> complex:
    real = coalition_externality({key: value.real for key, value in values.items()}, coalition)
    imag = coalition_externality({key: value.imag for key, value in values.items()}, coalition)
    return complex(real, imag)


def stable_lower_order_prediction(values: dict[tuple[int, ...], complex], coalition: tuple[int, ...]) -> complex:
    full = tuple(sorted(coalition))
    return complex(values[full] - complex_externality(values, coalition))


def hard_stop_reason(point: OperatorPoint, modes: GeneralizedModes | None) -> str | None:
    if point.status != "SUCCESS":
        return point.status
    if float(point.metadata["voltage_min_pu"]) < 0.85:
        return "VOLTAGE_MIN"
    if float(point.metadata["voltage_max_pu"]) > 1.15:
        return "VOLTAGE_MAX"
    if point.metadata.get("generator_active_limit_violation") or point.metadata.get("generator_reactive_limit_violation") or point.metadata.get("slack_limit_violation"):
        return "GENERATOR_LIMIT"
    if modes is None:
        return "EIG_FAIL"
    if not np.isfinite(modes.eigenvalues).all():
        return "EIG_FAIL"
    if float(np.max(modes.residuals)) > 1e-10:
        return "EIGENPAIR_RESIDUAL_FAIL"
    if maximum_physical_real_part(modes) >= 0.0:
        return "BASELINE_UNSTABLE"
    return None


def eigenvector_hash(modes: GeneralizedModes, index: int) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(modes.left[:, index]).view(np.uint8))
    digest.update(np.ascontiguousarray(modes.right[:, index]).view(np.uint8))
    return digest.hexdigest()
