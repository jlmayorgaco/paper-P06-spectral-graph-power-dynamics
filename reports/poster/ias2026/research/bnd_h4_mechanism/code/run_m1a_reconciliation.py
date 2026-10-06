"""Diagnose the DAE-to-retained-operator mismatch behind strict M1."""

from __future__ import annotations

import csv
import dataclasses
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import root as scipy_root


STENCIL = (0.0 + 0.0j, 1e-6 + 0.0j, -1e-6 + 0.0j, 0.0 + 1e-6j, 0.0 - 1e-6j)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(value))
    digest = hashlib.sha256()
    digest.update((array.dtype.str + str(array.shape)).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def normalize_hash_payload(value):
    if dataclasses.is_dataclass(value):
        return normalize_hash_payload(dataclasses.asdict(value))
    if isinstance(value, dict):
        return {str(key): normalize_hash_payload(value[key]) for key in sorted(value, key=str)}
    if isinstance(value, (list, tuple)):
        return [normalize_hash_payload(item) for item in value]
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.ndarray):
        return {"dtype": value.dtype.str, "shape": list(value.shape), "sha256": array_sha256(value)}
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def model_fingerprint(case) -> str:
    network = case.dae.network
    payload = {
        "bus_idx": normalize_hash_payload(network.bus_idx),
        "ybus_sha256": array_sha256(network.ybus),
        "loads": normalize_hash_payload(network.loads),
        "load_model": network.load_model,
        "machines": normalize_hash_payload(network.machines),
        "frequency_hz": getattr(network, "frequency_hz", None),
        "base_mva": getattr(network, "base_mva", None),
        "plan": normalize_hash_payload(case.plan),
        "slots": [
            {
                "bus": slot.bus,
                "kind": slot.kind,
                "weight": slot.weight,
                "loading": slot.loading,
                "device_type": type(slot.device).__module__ + "." + type(slot.device).__qualname__,
                "parameters": normalize_hash_payload(getattr(slot.device, "parameters", None)),
            }
            for slot in case.dae.slots
        ],
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def equilibrium_fingerprints(case) -> dict[str, str]:
    return {
        "xz_sha256": hashlib.sha256(
            (array_sha256(case.equilibrium.x) + array_sha256(case.equilibrium.z)).encode()
        ).hexdigest(),
        "reduced_A_sha256": array_sha256(case.system.A),
    }


def matched_base_port_operator(base_case, flagship_case, buses, reference_operator):
    """Linearize the base counterfactual at the matched flagship operating point.

    This changes only the comparison point used by the retained action-space
    diagnostic. It never changes either solved case, its equilibrium, dispatch,
    network, or model parameters.
    """
    from ibr_cycles.models.port_admittance import PortOperator, linearize_device, load_admittance

    base_dae, flag_dae = base_case.dae, flagship_case.dae
    base_net, flag_net = base_dae.network, flag_dae.network
    same_network = (
        base_net.bus_idx == flag_net.bus_idx
        and np.array_equal(base_net.ybus, flag_net.ybus)
        and base_net.loads == flag_net.loads
        and base_net.load_model == flag_net.load_model
        and np.array_equal(base_dae.power_flow.voltages, flag_dae.power_flow.voltages)
    )
    if not same_network:
        raise RuntimeError("matched base port construction requires the frozen identical network/PF")

    common_z = flagship_case.equilibrium.z
    common_v = flag_dae.voltages(common_z)
    flag_slots = {slot.bus: slot for slot in flag_dae.slots}
    action_buses = set(buses)
    x_base = np.empty(base_dae.n_x, dtype=float)
    for slot in base_dae.slots:
        voltage = complex(common_v[base_net.position(slot.bus)])
        if slot.bus in action_buses:
            # Keep the frozen base device parameters; use initialization only
            # to obtain its state at the shared voltage and unchanged dispatch.
            dispatch = base_dae.power_flow.injection(slot.bus, base_net.ybus)
            dispatch += base_net.loads.get(slot.bus, 0.0 + 0.0j)
            _, state = slot.device.initialize(voltage, dispatch)
        else:
            shared = flag_slots.get(slot.bus)
            if (
                shared is None
                or slot.kind != shared.kind
                or type(slot.device) is not type(shared.device)
                or slot.stop - slot.start != shared.stop - shared.start
                or getattr(slot.device, "parameters", None) != getattr(shared.device, "parameters", None)
            ):
                raise RuntimeError(f"unmatched non-action device at bus {slot.bus}")
            state = flagship_case.equilibrium.x[shared.start : shared.stop]
        state = np.asarray(state, dtype=float)
        if state.size != slot.stop - slot.start:
            raise RuntimeError(f"state-size mismatch at base device bus {slot.bus}")
        x_base[slot.start : slot.stop] = state

    ports = tuple(
        linearize_device(
            slot,
            x_base,
            complex(common_v[base_net.position(slot.bus)]),
        )
        for slot in base_dae.slots
    )
    loads = {
        bus: load_admittance(load, complex(common_v[base_net.position(bus)]))
        for bus, load in base_net.loads.items()
    }
    return PortOperator(
        ybus_real=reference_operator.ybus_real,
        ports=ports,
        load_blocks=loads,
        bus_index=reference_operator.bus_index,
        n_bus=reference_operator.n_bus,
    )


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=float), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def complex_pair(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def logdet(matrix: np.ndarray) -> complex:
    sign, logabs = np.linalg.slogdet(matrix)
    return complex(float(logabs), float(np.angle(sign)))


def multiplicative_residual(left: complex, right: complex) -> float:
    if not np.isfinite(left.real + left.imag + right.real + right.imag):
        return float("inf")
    return float(abs(np.exp(left - right) - 1.0))


def pencil(fx, fz, gx, gz, s: complex) -> np.ndarray:
    n_x = fx.shape[0]
    return np.block([[s * np.eye(n_x) - fx, -fz], [gx, gz]])


def raw_schur(fx, fz, gx, gz, s: complex) -> np.ndarray:
    return gz + gx @ np.linalg.solve(s * np.eye(fx.shape[0]) - fx, fz)


def modal_basis(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    left, singular_values, right_h = np.linalg.svd(matrix)
    l = left[:, -1]
    r = right_h.conj().T[:, -1]
    u = np.column_stack([l, null_space(l.conj().reshape(1, -1))])
    v = np.column_stack([r, null_space(r.conj().reshape(1, -1))])
    return u, v, float(singular_values[-1])


def modal_schur(space, u: np.ndarray, v: np.ndarray, s: complex, eta: float = 1.0):
    matrix = u.conj().T @ (np.eye(space.dimension, dtype=complex) + space.m(s)) @ v
    h00, h0b, hb0, hbb = matrix[0, 0], matrix[0, 1:], matrix[1:, 0], matrix[1:, 1:]
    cross = complex(h0b @ np.linalg.solve(hbb, hb0))
    h = complex(h00 - eta * cross)
    return h, cross, matrix


def solve_complex_root(function, initial: complex) -> tuple[complex, bool, float]:
    def residual(values):
        value = function(complex(values[0], values[1]))
        return [value.real, value.imag]

    result = scipy_root(residual, [initial.real, initial.imag], method="hybr", options={"xtol": 1e-11})
    candidate = complex(result.x[0], result.x[1])
    error = abs(function(candidate))
    return candidate, bool(result.success and np.isfinite(error)), float(error)


def mp_frobenius(matrix: np.ndarray) -> float:
    import mpmath as mp

    mp.mp.dps = 80
    terms = [abs(mp.mpc(float(value.real), float(value.imag))) ** 2 for value in matrix.flat]
    return float(mp.sqrt(mp.fsum(terms)))


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_m1a_reconciliation.py <repo-root> <run-root>")
    repo = Path(sys.argv[1]).resolve()
    run_root = Path(sys.argv[2]).resolve()
    if run_root.exists():
        raise FileExistsError(f"immutable M1A run root already exists: {run_root}")
    for name in ("raw/python", "raw/julia", "raw/reconciliation", "tables", "claims", "report", "environment", "logs"):
        (run_root / name).mkdir(parents=True, exist_ok=True)

    experiment = repo / "reports/poster/ias2026/research/experiments"
    source = repo / "reports/poster/ias2026/research/src"
    sys.path.insert(0, str(experiment))
    sys.path.insert(0, str(source))
    import tx4_contextual_return as tx  # type: ignore  # noqa: PLC0415
    from ibr_cycles.dynamics.linearize import central_difference_jacobians, reduce_index_one  # noqa: PLC0415
    from ibr_cycles.models.port_admittance import build_action_space, build_port_operator  # noqa: PLC0415

    case = tx.tx4_case(tx.CORE, tx.P4)
    base = tx.tx4_case((), tx.Theta(g=0.0, k=1.425, t=1.5, h=1.0))
    jac = central_difference_jacobians(case.dae, case.equilibrium.x, case.equilibrium.z, case.equilibrium.theta)
    fx, fz, gx, gz = jac.fx, jac.fz, jac.gx, jac.gz
    red = reduce_index_one(jac, case.system.labels).A
    values, vectors = np.linalg.eig(red)
    target = tx.mode_of(case).eig
    index = int(np.argmin(abs(values - target)))
    lambda_c = complex(values[index])
    vector = vectors[:, index]
    zv = -np.linalg.solve(gz, gx @ vector)
    descriptor = pencil(fx, fz, gx, gz, lambda_c)
    descriptor_vector = np.concatenate([vector, zv])
    r_a_abs = float(np.linalg.norm(red @ vector - lambda_c * vector))
    r_a_rel = r_a_abs / max(float(np.linalg.norm(red) * np.linalg.norm(vector)), 1e-300)
    r_desc_abs = float(np.linalg.norm(descriptor @ descriptor_vector))
    r_desc_rel = r_desc_abs / max(float(np.linalg.norm(descriptor) * np.linalg.norm(descriptor_vector)), 1e-300)

    raw = raw_schur(fx, fz, gx, gz, lambda_c)
    raw_u, raw_sigma, raw_vh = np.linalg.svd(raw)
    raw_q = raw_vh.conj().T[:, -1]
    raw_smin, raw_smax = float(raw_sigma[-1]), float(raw_sigma[0])
    raw_backward = float(np.linalg.norm(raw @ raw_q) / max(np.linalg.norm(raw) * np.linalg.norm(raw_q), 1e-300))
    fx_distance = float(np.min(abs(np.linalg.eigvals(fx) - lambda_c)))

    identity_rows = []
    for offset in STENCIL:
        s = lambda_c + offset
        p = pencil(fx, fz, gx, gz, s)
        t = raw_schur(fx, fz, gx, gz, s)
        lhs = logdet(p)
        right_a = logdet(gz) + logdet(s * np.eye(red.shape[0]) - red)
        right_b = logdet(s * np.eye(fx.shape[0]) - fx) + logdet(t)
        identity_rows.append({
            "offset_real": offset.real,
            "offset_imag": offset.imag,
            "is_root": bool(offset == 0.0j),
            "pencil_vs_Ared_ratio_residual": multiplicative_residual(lhs, right_a),
            "pencil_vs_raw_schur_ratio_residual": multiplicative_residual(lhs, right_b),
            "raw_schur_sigma_min": float(np.linalg.svd(t, compute_uv=False)[-1]),
        })

    candidate_operator = build_port_operator(case.dae, case.equilibrium.x, case.equilibrium.z)
    action_space = build_action_space(base, case, tx.CORE)
    matched_base = matched_base_port_operator(base, case, tx.CORE, action_space.ts)
    action_space = dataclasses.replace(action_space, t0=matched_base)
    port_rows = []
    block_rows = []
    ordered_buses = tuple(case.dae.network.bus_idx)
    for offset in STENCIL:
        s = lambda_c + offset
        t_port = candidate_operator.evaluate(s)
        t_raw = raw_schur(fx, fz, gx, gz, s)
        difference = t_port - t_raw
        port_rows.append({
            "offset_real": offset.real,
            "offset_imag": offset.imag,
            "absolute_frobenius_residual": float(np.linalg.norm(difference)),
            "relative_frobenius_residual": float(np.linalg.norm(difference) / max(np.linalg.norm(t_raw), 1e-300)),
            "candidate_port_condition": float(np.linalg.cond(t_port)),
        })
        for bus_i in ordered_buses:
            for bus_j in ordered_buses:
                pi = candidate_operator.bus_index[bus_i]
                pj = candidate_operator.bus_index[bus_j]
                rows = slice(2 * pi, 2 * pi + 2)
                cols = slice(2 * pj, 2 * pj + 2)
                block_raw = t_raw[rows, cols]
                block_delta = difference[rows, cols]
                entry_i, entry_j = np.unravel_index(np.argmax(np.abs(block_delta)), block_delta.shape)
                block_rows.append({
                    "offset_real": offset.real,
                    "offset_imag": offset.imag,
                    "bus_i": bus_i,
                    "bus_j": bus_j,
                    "block_class": "diagonal_device_load" if bus_i == bus_j else "network_offdiagonal",
                    "absolute_block_residual": float(np.linalg.norm(block_delta)),
                    "relative_block_residual": float(np.linalg.norm(block_delta) / max(np.linalg.norm(block_raw), 1e-300)),
                    "max_entry_abs": float(abs(block_delta[entry_i, entry_j])),
                    "max_entry_row_component": int(entry_i),
                    "max_entry_col_component": int(entry_j),
                    "max_entry_real": float(block_delta[entry_i, entry_j].real),
                    "max_entry_imag": float(block_delta[entry_i, entry_j].imag),
                })

    u = action_space.selector()
    delta_rows = []
    delta_block_rows = []
    for offset in STENCIL:
        s = lambda_c + offset
        t0 = action_space.t0.evaluate(s)
        th = action_space.ts.evaluate(s)
        delta = action_space.update(s)
        expected = t0 + u @ delta @ u.T
        difference = th - expected
        delta_rows.append({
            "offset_real": offset.real,
            "offset_imag": offset.imag,
            "absolute_frobenius_residual": float(np.linalg.norm(difference)),
            "relative_frobenius_residual": float(np.linalg.norm(difference) / max(np.linalg.norm(th), 1e-300)),
            "condition_base": float(np.linalg.cond(t0)),
            "condition_flagship": float(np.linalg.cond(th)),
        })
        for bus_i in ordered_buses:
            for bus_j in ordered_buses:
                pi = action_space.t0.bus_index[bus_i]
                pj = action_space.t0.bus_index[bus_j]
                block = difference[2 * pi : 2 * pi + 2, 2 * pj : 2 * pj + 2]
                entry_i, entry_j = np.unravel_index(np.argmax(np.abs(block)), block.shape)
                delta_block_rows.append({
                    "offset_real": offset.real,
                    "offset_imag": offset.imag,
                    "bus_i": bus_i,
                    "bus_j": bus_j,
                    "block_class": (
                        "action_bus_diagonal" if bus_i == bus_j and bus_i in tx.CORE
                        else "non_action_diagonal" if bus_i == bus_j
                        else "network_offdiagonal"
                    ),
                    "absolute_block_residual": float(np.linalg.norm(block)),
                    "relative_to_full_residual": float(
                        np.linalg.norm(block) / max(np.linalg.norm(difference), 1e-300)
                    ),
                    "max_entry_abs": float(abs(block[entry_i, entry_j])),
                    "max_entry_row_component": int(entry_i),
                    "max_entry_col_component": int(entry_j),
                    "max_entry_real": float(block[entry_i, entry_j].real),
                    "max_entry_imag": float(block[entry_i, entry_j].imag),
                })

    load_rows = []
    all_load_buses = sorted(set(action_space.t0.load_blocks) | set(action_space.ts.load_blocks))
    for bus in all_load_buses:
        zero = np.zeros((2, 2), dtype=complex)
        load_base = action_space.t0.load_blocks.get(bus, zero)
        load_flagship = action_space.ts.load_blocks.get(bus, zero)
        difference = load_flagship - load_base
        load_rows.append({
            "bus": int(bus),
            "absolute_residual": float(np.linalg.norm(difference)),
            "relative_residual": float(np.linalg.norm(difference) / max(np.linalg.norm(load_flagship), 1e-300)),
        })

    t0_lambda = action_space.t0.evaluate(lambda_c)
    th_lambda = action_space.ts.evaluate(lambda_c)
    device_rows = []
    for bus in ordered_buses:
        device_base = action_space.t0.bus_admittance(lambda_c, bus)
        device_flag = action_space.ts.bus_admittance(lambda_c, bus)
        load_base = action_space.t0.load_blocks.get(bus, np.zeros((2, 2), dtype=complex))
        load_flag = action_space.ts.load_blocks.get(bus, np.zeros((2, 2), dtype=complex))
        pi = action_space.t0.bus_index[bus]
        diagonal = slice(2 * pi, 2 * pi + 2)
        actual_delta = (th_lambda - t0_lambda)[diagonal, diagonal]
        device_rows.append({
            "bus": bus,
            "in_action_support": bus in tx.CORE,
            "base_kind": next((port.kind for port in action_space.t0.ports if port.bus == bus), "none"),
            "flagship_kind": next((port.kind for port in action_space.ts.ports if port.bus == bus), "none"),
            "device_admittance_delta_abs": float(np.linalg.norm(device_base - device_flag)),
            "load_block_delta_abs": float(np.linalg.norm(load_flag - load_base)),
            "operator_diagonal_delta_abs": float(np.linalg.norm(actual_delta)),
        })

    t0 = action_space.t0.evaluate(lambda_c)
    th = action_space.ts.evaluate(lambda_c)
    m = action_space.m(lambda_c)
    f_action = np.eye(action_space.dimension, dtype=complex) + m
    _, f_sigma, f_vh = np.linalg.svd(f_action)
    f_q = f_vh.conj().T[:, -1]
    action_abs = float(np.linalg.norm(f_action @ f_q))
    action_rel = action_abs / max(float(np.linalg.norm(f_action) * np.linalg.norm(f_q)), 1e-300)
    lifted = np.linalg.solve(t0, u @ f_q)
    lifted_abs = float(np.linalg.norm(th @ lifted))
    lifted_rel = lifted_abs / max(float(np.linalg.norm(th) * np.linalg.norm(lifted)), 1e-300)
    determinant_ratio_residual = float(abs(np.linalg.det(th) / np.linalg.det(t0) - np.linalg.det(f_action)))

    left_basis, right_basis, action_smin = modal_basis(f_action)
    modal_matrix = left_basis.conj().T @ f_action @ right_basis
    modal_roundtrip = left_basis @ modal_matrix @ right_basis.conj().T
    modal_roundtrip_abs = float(np.linalg.norm(f_action - modal_roundtrip))
    modal_roundtrip_rel = modal_roundtrip_abs / max(float(np.linalg.norm(f_action)), 1e-300)
    modal_orth_error = max(
        float(np.linalg.norm(left_basis.conj().T @ left_basis - np.eye(action_space.dimension))),
        float(np.linalg.norm(right_basis.conj().T @ right_basis - np.eye(action_space.dimension))),
    )
    h_value, cross, h_matrix = modal_schur(action_space, left_basis, right_basis, lambda_c)
    h_scale = max(float(np.linalg.norm(h_matrix)), 1e-300)
    schur_rows = []
    for offset in STENCIL:
        s = lambda_c + offset
        h, _, matrix = modal_schur(action_space, left_basis, right_basis, s)
        schur_rows.append({
            "offset_real": offset.real,
            "offset_imag": offset.imag,
            "h_abs": abs(h),
            "det_factorization_residual": float(abs(np.linalg.det(matrix) - np.linalg.det(matrix[1:, 1:]) * h)),
            "relative_det_factorization_residual": float(abs(np.linalg.det(matrix) - np.linalg.det(matrix[1:, 1:]) * h) / max(abs(np.linalg.det(matrix)), abs(np.linalg.det(matrix[1:, 1:]) * h), 1e-300)),
        })
    schur_root, root_success, root_residual = solve_complex_root(lambda s: modal_schur(action_space, left_basis, right_basis, s)[0], lambda_c)

    try:
        import mpmath as mp

        delta_lambda = action_space.update(lambda_c)
        expected_lambda = t0 + u @ delta_lambda @ u.T
        mp_delta_abs = mp_frobenius(th - expected_lambda)
        mp_delta_den = mp_frobenius(th)
        mp_delta_rel = mp_delta_abs / max(mp_delta_den, 1e-300)
        precision = {"backend": "mpmath", "dps": 80, "base_delta_absolute": mp_delta_abs, "base_delta_relative": mp_delta_rel}
    except Exception as error:  # pragma: no cover - environment diagnostic
        precision = {"backend": "unavailable", "error": repr(error)}

    raw_pass = raw_backward <= 1e-12
    port_pass = max(row["relative_frobenius_residual"] for row in port_rows) <= 1e-8
    delta_max = max(row["relative_frobenius_residual"] for row in delta_rows)
    delta_pass = delta_max <= 1e-12
    action_pass = action_rel <= 1e-12
    modal_pass = modal_roundtrip_rel <= 1e-12 and modal_orth_error <= 1e-12
    scalar_pass = abs(h_value) < 1e-8
    ladder = [
        {"stage": "Ared eigenpair", "matrix_dimension": red.shape[0], "absolute_residual": r_a_abs, "relative_residual": r_a_rel, "condition_number": float(np.linalg.cond(red)), "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if r_a_rel <= 1e-12 else "FAIL"},
        {"stage": "descriptor pencil", "matrix_dimension": descriptor.shape[0], "absolute_residual": r_desc_abs, "relative_residual": r_desc_rel, "condition_number": float(np.linalg.cond(descriptor)), "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if r_desc_rel <= 1e-12 else "FAIL"},
        {"stage": "raw network Schur", "matrix_dimension": raw.shape[0], "absolute_residual": raw_smin, "relative_residual": raw_backward, "condition_number": float(raw_smax / max(raw_smin, 1e-300)), "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if raw_pass else "FAIL"},
        {"stage": "candidate-port Schur", "matrix_dimension": candidate_operator.dimension, "absolute_residual": port_rows[0]["absolute_frobenius_residual"], "relative_residual": max(row["relative_frobenius_residual"] for row in port_rows), "condition_number": port_rows[0]["candidate_port_condition"], "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if port_pass else "FAIL"},
        {"stage": "base+deltaY", "matrix_dimension": th.shape[0], "absolute_residual": delta_rows[0]["absolute_frobenius_residual"], "relative_residual": delta_max, "condition_number": delta_rows[0]["condition_flagship"], "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if delta_pass else "FAIL"},
        {"stage": "action-space", "matrix_dimension": f_action.shape[0], "absolute_residual": action_abs, "relative_residual": action_rel, "condition_number": float(f_sigma[0] / max(f_sigma[-1], 1e-300)), "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if action_pass else "FAIL"},
        {"stage": "modal transform", "matrix_dimension": modal_matrix.shape[0], "absolute_residual": modal_roundtrip_abs, "relative_residual": modal_roundtrip_rel, "condition_number": 1.0, "root_real": lambda_c.real, "root_imag": lambda_c.imag, "root_distance_to_DAE": 0.0, "status": "PASS" if modal_pass else "FAIL"},
        {"stage": "scalar Schur", "matrix_dimension": 1, "absolute_residual": abs(h_value), "relative_residual": abs(h_value) / h_scale, "condition_number": float(np.linalg.cond(h_matrix[1:, 1:])), "root_real": schur_root.real, "root_imag": schur_root.imag, "root_distance_to_DAE": abs(schur_root - lambda_c), "status": "PASS" if scalar_pass else "FAIL"},
    ]
    write_csv(run_root / "tables" / "M1A_ERROR_LADDER.csv", ladder)
    write_csv(run_root / "tables" / "M1A_RAW_IDENTITY_STENCIL.csv", identity_rows)
    write_csv(run_root / "tables" / "M1A_CANDIDATE_PORT_STENCIL.csv", port_rows)
    write_csv(run_root / "tables" / "M1A_CANDIDATE_PORT_BLOCKS.csv", block_rows)
    write_csv(run_root / "tables" / "M1A_BASE_DELTAY_STENCIL.csv", delta_rows)
    write_csv(run_root / "tables" / "M1A_BASE_DELTAY_BLOCKS.csv", delta_block_rows)
    write_csv(run_root / "tables" / "M1A_LOAD_BLOCK_MISMATCH.csv", load_rows)
    write_csv(run_root / "tables" / "M1A_DEVICE_PORT_MISMATCH.csv", device_rows)
    write_csv(run_root / "tables" / "M1A_SCALAR_SCHUR_STENCIL.csv", schur_rows)
    nonroot_identity = [row for row in identity_rows if not row["is_root"]]
    write_json(run_root / "claims" / "M1A_RAW_DIAGNOSTICS.json", {"lambda_c": complex_pair(lambda_c), "fx_distance_to_lambda": fx_distance, "r_A_abs": r_a_abs, "r_A_relative": r_a_rel, "r_descriptor_abs": r_desc_abs, "r_descriptor_relative": r_desc_rel, "raw_sigma_min": raw_smin, "raw_sigma_max": raw_smax, "raw_relative_sigma": raw_smin / max(raw_smax, 1e-300), "raw_backward_error": raw_backward, "identity_stencil_max_Ared_nonroot": max(row["pencil_vs_Ared_ratio_residual"] for row in nonroot_identity), "identity_stencil_max_raw_nonroot": max(row["pencil_vs_raw_schur_ratio_residual"] for row in nonroot_identity), "identity_root_row_is_conditioned_by_zero": True})
    write_json(run_root / "claims" / "M1A_PORT_DIAGNOSTICS.json", {"candidate_port_max_relative_residual": max(row["relative_frobenius_residual"] for row in port_rows), "base_delta_max_relative_residual": delta_max, "base_delta_high_precision_relative_residual": precision.get("base_delta_relative"), "largest_load_block_mismatch": max(load_rows, key=lambda row: row["absolute_residual"]) if load_rows else None, "action_space_relative_backward_error": action_rel, "action_space_lifted_relative_backward_error": lifted_rel, "action_space_det_identity_residual": determinant_ratio_residual, "modal_roundtrip_relative_error": modal_roundtrip_rel, "modal_orthogonality_error": modal_orth_error, "scalar_h_abs": abs(h_value), "scalar_h_relative": abs(h_value) / h_scale, "scalar_schur_root": complex_pair(schur_root), "scalar_schur_root_residual": root_residual})
    write_json(run_root / "claims" / "M1A_PRECISION_AUDIT.json", precision)

    first_fail = next((row["stage"] for row in ladder if row["status"] == "FAIL"), "none")
    if first_fail == "base+deltaY":
        claim = "B_DELTAY_MISMATCH"
        next_action = "Reconcile matched load blocks between base and flagship port operators before rerunning M1."
    elif first_fail == "raw network Schur":
        claim = "G_FULL_LINEARIZATION_MISMATCH"
        next_action = "Reconcile the Jacobian blocks and raw Schur construction before any port analysis."
    elif first_fail == "candidate-port Schur":
        claim = "A_PORT_RECONSTRUCTION_BUG"
        next_action = "Reconcile the candidate-port linearization against T_raw before action compression."
    elif first_fail == "action-space":
        claim = "C_ACTION_COMPRESSION_MISMATCH"
        next_action = "Repair the action-space embedding/update identity before modal reduction."
    elif first_fail == "modal transform":
        claim = "D_MODAL_BASIS_MISMATCH"
        next_action = "Repair the modal basis transform and its metric before scalar Schur work."
    elif first_fail == "scalar Schur":
        claim = "E_SCALAR_SCHUR_MISMATCH"
        next_action = "Audit the scalar Schur partition only after all matrix identities pass."
    else:
        claim = "F_FLOATING_POINT_SCALE_LIMIT"
        next_action = "Propose a separately preregistered backward-error M1 v1.1; do not rewrite v1."

    claim_status = {
        "status": "M1A_COMPLETE",
        "first_degrading_stage": first_fail,
        "root_cause_category": claim,
        "max_absolute_residual": max(row["absolute_residual"] for row in ladder),
        "max_relative_backward_error": max(row["relative_residual"] for row in ladder),
        "recommended_m1_status": "BLOCKED_M1_STRICT",
        "safe_to_run_m2": False,
        "one_next_corrective_action": next_action,
        "ladder": ladder,
    }

    model_hashes = {"base": model_fingerprint(base), "flagship": model_fingerprint(case)}
    equilibrium_hashes = {
        "base": equilibrium_fingerprints(base),
        "flagship": equilibrium_fingerprints(case),
    }
    canonical_mode = tx.mode_of(case)
    f1_path = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/results/20260926T083621_h4_f1_f2_assets_v1/derived/F1_PORTFOLIOS.csv"
    parity_path = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/results/20260925T225832_dee6fe69_h4_modal_schur_v1/claims/CROSSCODE_GFL11_RECONCILIATION.json"
    with f1_path.open(newline="", encoding="utf-8") as handle:
        f1_rows = list(csv.DictReader(handle))
    f1_h4_rows = [row for row in f1_rows if row["portfolio"] == "30+33+35+37"]
    f1_proper_rows = [row for row in f1_rows if row["portfolio"] != "30+33+35+37"]
    parity = json.loads(parity_path.read_text(encoding="utf-8"))
    regression = {
        "current_h4": {
            "state_count": case.n_states,
            "alpha": canonical_mode.alpha,
            "frequency_hz": canonical_mode.frequency_hz,
            "model_hash_matches_pre_correction": model_hashes["flagship"] == "5bd6c2c82d876e9667ba177eff952b958902967f19b887fb1e3d7aba06e4fb5c",
            "equilibrium_hash_matches_pre_correction": equilibrium_hashes["flagship"]["xz_sha256"] == "18adf2938aa8263bbdb97165dceba885d7323a935664a1d18f16fbac7294ac34",
            "reduced_A_hash_matches_pre_correction": equilibrium_hashes["flagship"]["reduced_A_sha256"] == "da25d45096f80d872ef40caa968fa2c996ef0cc6b1bb92420ad2869c7cc92751",
        },
        "frozen_proper_subset_ledger_read_only": {
            "path": str(f1_path),
            "sha256": sha256(f1_path),
            "portfolio_rows": len(f1_rows),
            "proper_subset_rows_including_base": len(f1_proper_rows),
            "all_proper_subsets_stable": bool(f1_proper_rows) and all(row["status"] == "STABLE" for row in f1_proper_rows),
            "h4_status": f1_h4_rows[0]["status"] if len(f1_h4_rows) == 1 else "INVALID_LEDGER",
        },
        "frozen_julia_gfl11_parity_read_only": {
            "path": str(parity_path),
            "sha256": sha256(parity_path),
            "status": parity.get("status"),
            "python_states": parity["h4"]["python_state_dim"],
            "julia_states": parity["h4"]["julia_state_dim"],
            "alpha_difference_s-1": parity["h4"]["alpha_difference_s-1"],
            "frequency_difference_hz": parity["h4"]["frequency_difference_hz"],
            "pass": bool(
                parity.get("status") == "PASS_GFL11_EXACT_REPRODUCTION"
                and parity["h4"]["python_state_dim"] == 86
                and parity["h4"]["julia_state_dim"] == 86
            ),
        },
    }

    assembly_provenance = {
        "case": "frozen H4/P4 GFL11",
        "action_buses": list(tx.CORE),
        "base_counterfactual_linearization_point": "flagship solved equilibrium z; unchanged shared-device states from flagship equilibrium",
        "base_replacement_device_dispatch": "unchanged base power-flow generation at each action bus",
        "base_replacement_device_parameters": "frozen base slot parameters retained; initializer-returned retuned device discarded",
        "network_ybus_exactly_equal": bool(np.array_equal(base.dae.network.ybus, case.dae.network.ybus)),
        "bus_order_exactly_equal": bool(base.dae.network.bus_idx == case.dae.network.bus_idx),
        "bus_index_exactly_equal": bool(
            {bus: base.dae.network.position(bus) for bus in base.dae.network.bus_idx}
            == {bus: case.dae.network.position(bus) for bus in case.dae.network.bus_idx}
        ),
        "power_flow_voltage_exactly_equal": bool(
            np.array_equal(base.dae.power_flow.voltages, case.dae.power_flow.voltages)
        ),
        "load_maps_exactly_equal": bool(base.dae.network.loads == case.dae.network.loads),
        "load_model_equal": bool(base.dae.network.load_model == case.dae.network.load_model),
        "base_flag_equilibrium_voltage_distance_max_pu": float(
            np.max(abs(base.dae.voltages(base.equilibrium.z) - case.dae.voltages(case.equilibrium.z)))
        ),
        "base_equilibrium_vs_common_pf_voltage_distance_max_pu": float(
            np.max(abs(base.dae.voltages(base.equilibrium.z) - base.dae.power_flow.voltages))
        ),
        "flag_equilibrium_vs_common_pf_voltage_distance_max_pu": float(
            np.max(abs(case.dae.voltages(case.equilibrium.z) - case.dae.power_flow.voltages))
        ),
        "port_equation": "T(s)=Y_net-sum(device port admittances)-load blocks",
        "replacement_update_sign": "T_H-T_0 = E (Y_base_device-Y_flag_device) E^T",
        "port_normalization": "native rectangular current/voltage coordinates; no rescaling",
        "canonical_dae_equilibrium_unchanged": True,
    }
    write_json(run_root / "claims" / "M1A_ASSEMBLY_PROVENANCE.json", assembly_provenance)
    write_json(run_root / "claims" / "M1A_REGRESSION_CHECKS.json", regression)

    baseline_id = "20260926T100829_d0fecb32_ias26_010_reproduce"
    baseline_root = run_root.parent / baseline_id
    baseline_ladder_path = baseline_root / "tables/M1A_ERROR_LADDER.csv"
    baseline_port_path = baseline_root / "claims/M1A_PORT_DIAGNOSTICS.json"
    baseline_status_path = baseline_root / "claims/M1A_STATUS.json"
    with baseline_ladder_path.open(newline="", encoding="utf-8") as handle:
        baseline_ladder = list(csv.DictReader(handle))
    baseline_port = json.loads(baseline_port_path.read_text(encoding="utf-8"))
    baseline_status = json.loads(baseline_status_path.read_text(encoding="utf-8"))
    changed_script_rel = Path("reports/poster/ias2026/research/bnd_h4_mechanism/code/run_m1a_reconciliation.py")
    changed_diff = subprocess.check_output(
        ["git", "diff", "--no-ext-diff", "--no-color", "HEAD", "--", str(changed_script_rel)],
        cwd=repo,
        text=True,
    )
    diff_path = run_root / "raw/reconciliation/IAS26_010_CHANGED_FILES.patch"
    diff_path.write_text(changed_diff, encoding="utf-8")
    before_after = {
        "before_run_id": baseline_id,
        "after_run_id": run_root.name,
        "before_error_ladder": baseline_ladder,
        "after_error_ladder": ladder,
        "before_claim": baseline_status,
        "after_claim": claim_status,
        "before_action_space_relative_backward_error": baseline_port["action_space_relative_backward_error"],
        "after_action_space_relative_backward_error": action_rel,
        "before_base_delta_relative_residual": baseline_port["base_delta_max_relative_residual"],
        "after_base_delta_relative_residual": delta_max,
        "before_root_cause_record": "docs/IAS26-010_ROOT_CAUSE.md",
        "baseline_full_blocks": str(baseline_root / "tables/M1A_BASE_DELTAY_BLOCKS.csv"),
        "after_full_blocks": str(run_root / "tables/M1A_BASE_DELTAY_BLOCKS.csv"),
        "changed_file_diff": str(diff_path),
        "changed_file_diff_sha256": sha256(diff_path),
        "model_hashes_after": model_hashes,
        "equilibrium_and_reduced_A_hashes_after": equilibrium_hashes,
        "regression_checks": regression,
    }
    write_json(run_root / "claims" / "M1A_BEFORE_AFTER.json", before_after)

    source_files = {
        "run_m1a_reconciliation.py": repo / changed_script_rel,
        "tx4_contextual_return.py": experiment / "tx4_contextual_return.py",
        "_f7_common.py": experiment / "_f7_common.py",
        "ieee39_case.py": source / "ibr_cycles/models/ieee39_case.py",
        "ieee39_network.py": source / "ibr_cycles/models/ieee39_network.py",
        "ieee39_devices.py": source / "ibr_cycles/models/ieee39_devices.py",
        "port_admittance.py": source / "ibr_cycles/models/port_admittance.py",
        "IAS26-010_ROOT_CAUSE.md": repo / "reports/poster/ias2026/research/bnd_h4_mechanism/docs/IAS26-010_ROOT_CAUSE.md",
    }
    manifest = {
        "run_id": run_root.name,
        "run_root": str(run_root),
        "git_commit": git(repo, "rev-parse", "HEAD"),
        "git_dirty": bool(git(repo, "status", "--porcelain")),
        "python": {"version": sys.version, "implementation": platform.python_implementation(), "executable": sys.executable},
        "packages": {"numpy": np.__version__},
        "source_hashes": {name: sha256(path) for name, path in source_files.items()},
        "model_hashes": model_hashes,
        "equilibrium_and_canonical_matrix_hashes": equilibrium_hashes,
        "regression_checks": regression,
        "changed_file_diff": str(diff_path),
        "changed_file_diff_sha256": sha256(diff_path),
        "output_policy": "all diagnostics are beneath this immutable run root; no zip is produced",
    }
    write_json(run_root / "environment" / "runtime.json", manifest)
    write_json(run_root / "claims" / "M1A_STATUS.json", claim_status)
    (run_root / "report" / "RUN_SUMMARY.md").write_text("# M1A DAE-to-retained-operator reconciliation\n\n" + json.dumps(claim_status, indent=2, default=float) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
