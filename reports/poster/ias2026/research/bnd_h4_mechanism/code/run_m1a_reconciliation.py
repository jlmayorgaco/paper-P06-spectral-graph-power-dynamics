"""Diagnose the DAE-to-retained-operator mismatch behind strict M1."""

from __future__ import annotations

import csv
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
            "pencil_vs_Ared_ratio_residual": multiplicative_residual(lhs, right_a),
            "pencil_vs_raw_schur_ratio_residual": multiplicative_residual(lhs, right_b),
            "raw_schur_sigma_min": float(np.linalg.svd(t, compute_uv=False)[-1]),
        })

    candidate_operator = build_port_operator(case.dae, case.equilibrium.x, case.equilibrium.z)
    action_space = build_action_space(base, case, tx.CORE)
    port_rows = []
    block_rows = []
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
        for bus_i in tx.CORE:
            for bus_j in tx.CORE:
                pi = candidate_operator.bus_index[bus_i]
                pj = candidate_operator.bus_index[bus_j]
                rows = slice(2 * pi, 2 * pi + 2)
                cols = slice(2 * pj, 2 * pj + 2)
                block_raw = t_raw[rows, cols]
                block_delta = difference[rows, cols]
                block_rows.append({
                    "offset_real": offset.real,
                    "offset_imag": offset.imag,
                    "bus_i": bus_i,
                    "bus_j": bus_j,
                    "absolute_block_residual": float(np.linalg.norm(block_delta)),
                    "relative_block_residual": float(np.linalg.norm(block_delta) / max(np.linalg.norm(block_raw), 1e-300)),
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
        for i, bus_i in enumerate(tx.CORE):
            for j, bus_j in enumerate(tx.CORE):
                block = action_space.block(difference, i, j)
                delta_block_rows.append({
                    "offset_real": offset.real,
                    "offset_imag": offset.imag,
                    "bus_i": bus_i,
                    "bus_j": bus_j,
                    "absolute_block_residual": float(np.linalg.norm(block)),
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
    write_csv(run_root / "tables" / "M1A_SCALAR_SCHUR_STENCIL.csv", schur_rows)
    write_json(run_root / "claims" / "M1A_RAW_DIAGNOSTICS.json", {"lambda_c": complex_pair(lambda_c), "fx_distance_to_lambda": fx_distance, "r_A_abs": r_a_abs, "r_A_relative": r_a_rel, "r_descriptor_abs": r_desc_abs, "r_descriptor_relative": r_desc_rel, "raw_sigma_min": raw_smin, "raw_sigma_max": raw_smax, "raw_relative_sigma": raw_smin / max(raw_smax, 1e-300), "raw_backward_error": raw_backward, "identity_stencil_max_Ared": max(row["pencil_vs_Ared_ratio_residual"] for row in identity_rows), "identity_stencil_max_raw": max(row["pencil_vs_raw_schur_ratio_residual"] for row in identity_rows)})
    write_json(run_root / "claims" / "M1A_PORT_DIAGNOSTICS.json", {"candidate_port_max_relative_residual": max(row["relative_frobenius_residual"] for row in port_rows), "base_delta_max_relative_residual": delta_max, "base_delta_high_precision_relative_residual": precision.get("base_delta_relative"), "action_space_relative_backward_error": action_rel, "action_space_lifted_relative_backward_error": lifted_rel, "action_space_det_identity_residual": determinant_ratio_residual, "modal_roundtrip_relative_error": modal_roundtrip_rel, "modal_orthogonality_error": modal_orth_error, "scalar_h_abs": abs(h_value), "scalar_h_relative": abs(h_value) / h_scale, "scalar_schur_root": complex_pair(schur_root), "scalar_schur_root_residual": root_residual})
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

    manifest = {"run_root": str(run_root), "git_commit": git(repo, "rev-parse", "HEAD"), "git_dirty": bool(git(repo, "status", "--porcelain")), "python": {"version": sys.version, "implementation": platform.python_implementation(), "executable": sys.executable}, "packages": {"numpy": np.__version__}, "source_hashes": {"tx4_contextual_return.py": sha256(experiment / "tx4_contextual_return.py"), "ieee39_case.py": sha256(source / "ibr_cycles/models/ieee39_case.py"), "port_admittance.py": sha256(source / "ibr_cycles/models/port_admittance.py")}, "output_policy": "all diagnostics are beneath this immutable run root; no zip is produced"}
    write_json(run_root / "environment" / "runtime.json", manifest)
    claim_status = {"status": "M1A_COMPLETE", "first_degrading_stage": first_fail, "root_cause_category": claim, "max_absolute_residual": max(row["absolute_residual"] for row in ladder), "max_relative_backward_error": max(row["relative_residual"] for row in ladder), "recommended_m1_status": "BLOCKED_M1_STRICT", "safe_to_run_m2": False, "one_next_corrective_action": next_action, "ladder": ladder}
    write_json(run_root / "claims" / "M1A_STATUS.json", claim_status)
    (run_root / "report" / "RUN_SUMMARY.md").write_text("# M1A DAE-to-retained-operator reconciliation\n\n" + json.dumps(claim_status, indent=2, default=float) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
