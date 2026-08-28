from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import andes
import numpy as np
from kvxopt import matrix as kvx_matrix
from scipy.sparse import csc_matrix, diags


ROOT = Path(__file__).resolve().parents[3]
E03 = ROOT / "experiments" / "tx3" / "E03_finite_externality"
E05C = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned"
sys.path[:0] = [str(ROOT), str(E03), str(E05C)]

from campaign_core import BASE_GFL, OperatorPoint, _uid, action_vector  # noqa: E402
from dicgrid.adapters.andes import enforce_gfl_pq_equilibrium  # noqa: E402
from margin_conditioned import (  # noqa: E402,F401
    GeneralizedModes,
    biorthogonal_mac,
    complex_externality,
    damping_ratio,
    eigenvector_hash,
    generalized_modes,
    maximum_physical_real_part,
    mobius_vertices,
    parse_subset_key,
    physical_indices,
    positive_oscillatory_indices,
    stable_lower_order_prediction,
    subset_key,
    track_assignment,
    track_one_hungarian,
    union_vertices,
)


CASE = ROOT / "systems" / "ieee39_tx3_gfl3" / "andes_case.xlsx"
SYNCHRONOUS_SOURCE_BUSES = (30, 31, 32, 33, 34, 35, 39)
SYNCHRONOUS_PV_BUSES = (30, 31, 32, 33, 34, 35)
GFL_POI_BUSES = (36, 37, 38)
EXPECTED_STRESSED_LINES = tuple(f"Line_{index}" for index in range(1, 35))
VALID_STATUSES = {
    "SUCCESS", "PF_FAIL", "INIT_FAIL", "EIG_FAIL", "MODE_TRACK_FAIL",
    "EIGENPAIR_RESIDUAL_FAIL", "NONSMOOTH_REGIME", "LIMITER_ACTIVE",
    "VOLTAGE_LIMIT", "GENERATOR_LIMIT", "BASELINE_UNSTABLE",
    "LOWER_ORDER_UNSTABLE", "NUMERICAL_UNRESOLVED",
}


def line_branch_records(system: Any) -> list[dict[str, Any]]:
    """Classify every branch using only frozen topology/model metadata."""

    records: list[dict[str, Any]] = []
    for index in range(system.Line.n):
        line_id = str(system.Line.idx.v[index])
        bus1 = int(system.Line.bus1.v[index])
        bus2 = int(system.Line.bus2.v[index])
        is_transformer = bool(int(system.Line.trans.v[index]))
        if line_id == "Line_45" and {bus1, bus2} == {25, 37}:
            included = False
            reason = "excluded_A8_branch_and_local_GFL37_step_up_interface"
        elif bus1 in GFL_POI_BUSES or bus2 in GFL_POI_BUSES:
            included = False
            reason = "excluded_local_GFL_terminal_step_up_interface"
        elif is_transformer:
            included = False
            reason = "excluded_fixed_tap_transformer_model"
        else:
            included = True
            reason = "included_AC_transmission_line_trans_flag_zero"
        records.append(
            {
                "line_id": line_id,
                "from_bus": bus1,
                "to_bus": bus2,
                "R_system_pu": float(system.Line.r.v[index]),
                "X_system_pu": float(system.Line.x.v[index]),
                "tap": float(system.Line.tap.v[index]),
                "transformer_flag": int(system.Line.trans.v[index]),
                "included": included,
                "reason": reason,
            }
        )
    included_ids = tuple(record["line_id"] for record in records if record["included"])
    if included_ids != EXPECTED_STRESSED_LINES:
        raise RuntimeError(f"canonical external-grid topology changed: {included_ids}")
    return records


def stressed_line_indices(system: Any) -> np.ndarray:
    records = line_branch_records(system)
    return np.asarray([index for index, record in enumerate(records) if record["included"]], dtype=int)


def network_ybus(system: Any) -> np.ndarray:
    """Passive positive-sequence network Ybus, with constant-power loads opened."""

    ybus = np.asarray(kvx_matrix(system.Line.build_ybus()), dtype=complex)
    if system.Shunt.n:
        for bus, conductance, susceptance, enabled in zip(
            system.Shunt.bus.v, system.Shunt.g.v, system.Shunt.b.v, system.Shunt.u.v, strict=True
        ):
            uid = system.Bus.idx2uid(int(bus))
            ybus[uid, uid] += float(enabled) * complex(float(conductance), float(susceptance))
    return ybus


def short_circuit_ratios_from_ybus(
    ybus: np.ndarray,
    poi_indices: Iterable[int],
    source_indices: Iterable[int],
    converter_ratings_system_pu: Iterable[float],
    poi_voltages_pu: Iterable[float] | None = None,
) -> dict[int, float]:
    """Return conventional SCR = V_POI^2 / (|Z_th| S_rated) in per unit."""

    poi_indices = tuple(int(value) for value in poi_indices)
    sources = set(int(value) for value in source_indices)
    keep = [index for index in range(ybus.shape[0]) if index not in sources]
    reduced = ybus[np.ix_(keep, keep)]
    impedance = np.linalg.inv(reduced)
    voltages = tuple(float(value) for value in (poi_voltages_pu or [1.0] * len(poi_indices)))
    ratings = tuple(float(value) for value in converter_ratings_system_pu)
    if len(voltages) != len(poi_indices) or len(ratings) != len(poi_indices):
        raise ValueError("POI, voltage, and rating vectors must have equal length")
    result: dict[int, float] = {}
    for poi, voltage, rating in zip(poi_indices, voltages, ratings, strict=True):
        local = keep.index(poi)
        zth = impedance[local, local]
        result[poi] = float(voltage * voltage / (abs(zth) * rating))
    return result


def short_circuit_ratios(system: Any, voltages: np.ndarray) -> dict[str, float]:
    poi_indices = [system.Bus.idx2uid(bus) for bus in GFL_POI_BUSES]
    source_indices = [system.Bus.idx2uid(bus) for bus in SYNCHRONOUS_SOURCE_BUSES]
    regcp_buses = [int(value) for value in system.REGCP1.bus.v]
    ratings = [
        float(system.REGCP1.Sn.v[_uid(regcp_buses, bus)]) / float(system.config.mva)
        for bus in GFL_POI_BUSES
    ]
    values = short_circuit_ratios_from_ybus(
        network_ybus(system), poi_indices, source_indices, ratings,
        [float(voltages[index]) for index in poi_indices],
    )
    result = {f"SCR{bus}": values[index] for bus, index in zip(GFL_POI_BUSES, poi_indices, strict=True)}
    result["min_SCR"] = min(result.values())
    return result


def evaluate_weak_grid_operator(
    seed: dict[str, Any], kappa_grid: float, alpha: dict[str, float]
) -> OperatorPoint:
    """Build one fully re-equilibrated E05D operator under external-grid weakening."""

    started = time.perf_counter()
    alpha = action_vector(alpha)
    try:
        if float(kappa_grid) < 1.0:
            raise ValueError("kappa_grid must be at least 1")
        system = andes.load(str(CASE), setup=True, no_output=True)
        branch_indices = stressed_line_indices(system)
        system.Line.r.v[branch_indices] *= float(kappa_grid)
        system.Line.x.v[branch_indices] *= float(kappa_grid)

        base_pq_p = np.asarray(system.PQ.p0.v, dtype=float).copy()
        base_pq_q = np.asarray(system.PQ.q0.v, dtype=float).copy()
        system.PQ.p0.v[:] = base_pq_p * float(seed["load_scale"])
        system.PQ.q0.v[:] = base_pq_q * float(seed["load_scale"]) * float(seed["reactive_load_scale"])

        pv_buses = [int(value) for value in system.PV.bus.v]
        redispatch = float(seed["redispatch_mw"]) / float(system.config.mva)
        system.PV.p0.v[_uid(pv_buses, 30)] += redispatch
        system.PV.p0.v[_uid(pv_buses, 35)] -= redispatch

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
        if not bool(system.PFlow.run()) or not system.PFlow.converged:
            return OperatorPoint(None, None, "PF_FAIL", {
                "kappa_grid": float(kappa_grid), "elapsed_s": time.perf_counter() - started,
            })

        voltage = np.asarray(system.Bus.v.v, dtype=float)
        scr = short_circuit_ratios(system, voltage)
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

        try:
            system.TDS.init()
        except Exception as exc:
            return OperatorPoint(None, None, "INIT_FAIL", {
                "kappa_grid": float(kappa_grid), **scr, "exception": repr(exc),
                "elapsed_s": time.perf_counter() - started,
            })
        residual = max(
            float(np.max(np.abs(np.asarray(system.dae.f, dtype=float)))),
            float(np.max(np.abs(np.asarray(system.dae.g, dtype=float)))),
        )
        if not np.isfinite(residual) or residual > 1e-6:
            return OperatorPoint(None, None, "INIT_FAIL", {
                "kappa_grid": float(kappa_grid), **scr, "initialization_residual": residual,
                "elapsed_s": time.perf_counter() - started,
            })

        device_scale = float(system.config.mva) / np.asarray(system.REGCP1.Sn.v, dtype=float)
        current = np.hypot(
            np.asarray(system.REGCP1.Ipout.v, dtype=float) * device_scale,
            np.asarray(system.REGCP1.Iqout_x.v, dtype=float) * device_scale,
        )
        current_limit = np.asarray(system.REECB1.Imax.v, dtype=float)
        limiter_active = bool(np.any(current >= current_limit - 1e-8))
        saturation_active = bool(
            np.any(np.abs(np.asarray(system.REECB1.PIV_ys.v) - np.asarray(system.REECB1.PIV_y.v)) > 1e-9)
            or np.any(np.abs(np.asarray(system.REECB1.PIQ_ys.v) - np.asarray(system.REECB1.PIQ_y.v)) > 1e-9)
        )

        dense = lambda value: np.asarray(kvx_matrix(value), dtype=float)
        dae = system.dae
        jacobian = np.block([[dense(dae.fx), dense(dae.fy)], [dense(dae.gx), dense(dae.gy)]])
        mass = np.concatenate((np.asarray(dae.Tf, dtype=float), np.zeros(dae.m)))
        if not np.isfinite(jacobian).all() or not np.isfinite(mass).all():
            return OperatorPoint(None, None, "NUMERICAL_UNRESOLVED", {
                "kappa_grid": float(kappa_grid), **scr, "elapsed_s": time.perf_counter() - started,
            })
        metadata = {
            "elapsed_s": time.perf_counter() - started,
            "pf_converged": True,
            "kappa_grid": float(kappa_grid),
            "stressed_branch_count": int(len(branch_indices)),
            "seeded_active_load_system_pu": float(np.sum(system.PQ.p0.v)),
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
            **scr,
        }
        if float(np.min(voltage)) < 0.85 or float(np.max(voltage)) > 1.15:
            status = "VOLTAGE_LIMIT"
        elif pv_p_limit or pv_q_limit or slack_limit:
            status = "GENERATOR_LIMIT"
        elif limiter_active:
            status = "LIMITER_ACTIVE"
        elif saturation_active:
            status = "NONSMOOTH_REGIME"
        else:
            status = "SUCCESS"
        return OperatorPoint(csc_matrix(jacobian), diags(mass, format="csc"), status, metadata)
    except Exception as exc:
        return OperatorPoint(None, None, "NUMERICAL_UNRESOLVED", {
            "kappa_grid": float(kappa_grid), "exception": repr(exc),
            "elapsed_s": time.perf_counter() - started,
        })


def critical_mode_index(modes: GeneralizedModes, band_hz: tuple[float, float] = (0.1, 30.0)) -> int:
    indices = positive_oscillatory_indices(modes, band_hz)
    if not len(indices):
        raise ValueError("no valid positive-imaginary oscillatory mode in frozen band")
    damping = np.asarray([damping_ratio(complex(modes.eigenvalues[index])) for index in indices])
    finite = np.isfinite(damping)
    if not finite.any():
        raise ValueError("no finite damping ratios in frozen band")
    eligible = indices[finite]
    values = damping[finite]
    return int(eligible[np.argmin(values)])


def baseline_boundary(point: OperatorPoint, modes: GeneralizedModes | None, damping_screen: float = 0.05) -> str | None:
    if point.status != "SUCCESS":
        return point.status
    if modes is None:
        return "EIG_FAIL"
    if float(np.max(modes.residuals)) > 1e-10:
        return "EIGENPAIR_RESIDUAL_FAIL"
    if maximum_physical_real_part(modes) >= 0.0:
        return "BASELINE_UNSTABLE"
    try:
        index = critical_mode_index(modes)
    except ValueError:
        return "EIG_FAIL"
    if damping_ratio(complex(modes.eigenvalues[index])) <= float(damping_screen):
        return "DAMPING_SCREEN"
    return None


def damping_and_frequency(value: complex) -> tuple[float, float]:
    return damping_ratio(value), float(value.imag / (2.0 * math.pi))
