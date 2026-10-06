"""Run the narrow H4 scalar-Schur / collective-feedback experiment.

The script reuses the frozen TX4 model and port action space. It adds only a
modal basis, scalar Schur reduction, eta continuation, derivative check, and
run-local reporting. M2/M3 are not promoted when the strict M1 gate fails.
"""

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
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mac(left: np.ndarray, right: np.ndarray) -> float:
    den = float(np.vdot(left, left).real * np.vdot(right, right).real)
    return float(abs(np.vdot(left, right)) ** 2 / den) if den else 0.0


def make_basis(frozen_operator: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    left_values, singular_values, right_h = np.linalg.svd(frozen_operator)
    left = left_values[:, -1]
    right = right_h.conj().T[:, -1]
    left_basis = np.column_stack([left, null_space(left.conj().reshape(1, -1))])
    right_basis = np.column_stack([right, null_space(right.conj().reshape(1, -1))])
    return left_basis, right_basis, float(singular_values[-1])


def transformed(space, left_basis: np.ndarray, right_basis: np.ndarray, s: complex) -> np.ndarray:
    identity = np.eye(space.dimension, dtype=complex)
    return left_basis.conj().T @ (identity + space.m(s)) @ right_basis


def scalar_schur(
    space,
    left_basis: np.ndarray,
    right_basis: np.ndarray,
    s: complex,
    eta: float,
) -> tuple[complex, complex, np.ndarray]:
    matrix = transformed(space, left_basis, right_basis, s)
    h00 = matrix[0, 0]
    h0b = matrix[0, 1:]
    hb0 = matrix[1:, 0]
    hbb = matrix[1:, 1:]
    cross = complex(h0b @ np.linalg.solve(hbb, hb0))
    return complex(h00 - eta * cross), cross, matrix


def synthetic_operator(
    space,
    left_basis: np.ndarray,
    right_basis: np.ndarray,
    s: complex,
    eta: float,
) -> np.ndarray:
    h = transformed(space, left_basis, right_basis, s)
    scale = float(np.sqrt(max(eta, 0.0)))
    h[0, 1:] *= scale
    h[1:, 0] *= scale
    return left_basis @ h @ right_basis.conj().T


def tracked_vector(
    space,
    left_basis: np.ndarray,
    right_basis: np.ndarray,
    s: complex,
    eta: float,
) -> np.ndarray:
    _, _, vh = np.linalg.svd(synthetic_operator(space, left_basis, right_basis, s, eta))
    return vh.conj().T[:, -1]


def solve_complex_root(
    function,
    initial: complex,
    *,
    tolerance: float = 1e-11,
) -> tuple[complex, bool, float, int]:
    def residual(values: np.ndarray) -> np.ndarray:
        value = function(complex(values[0], values[1]))
        return np.asarray([value.real, value.imag], dtype=float)

    result = scipy_root(residual, [initial.real, initial.imag], method="hybr", options={"xtol": tolerance})
    candidate = complex(result.x[0], result.x[1])
    error = abs(function(candidate))
    return candidate, bool(result.success and np.isfinite(error)), float(error), int(result.nfev)


def plot_outputs(run_root: Path, rows: list[dict], repo: Path) -> None:
    import matplotlib.pyplot as plt

    figures = run_root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    eta = np.asarray([float(row["eta"]) for row in rows])
    real = np.asarray([float(row["root_real_s-1"]) for row in rows])
    freq = np.asarray([float(row["root_frequency_hz"]) for row in rows])

    plt.figure(figsize=(6.2, 4.0))
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.plot(eta, real, "o-", color="#0f766e")
    plt.xlabel(r"collective feedback $\eta$")
    plt.ylabel(r"$\Re\lambda_k$ [s$^{-1}$]")
    plt.tight_layout()
    plt.savefig(figures / "F01_cross_feedback_root_locus.png", dpi=220)
    plt.close()

    plt.figure(figsize=(6.2, 4.0))
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.plot(freq, real, "o-", color="#b45309")
    plt.xlabel(r"$\Im\lambda_k/(2\pi)$ [Hz]")
    plt.ylabel(r"$\Re\lambda_k$ [s$^{-1}$]")
    plt.tight_layout()
    plt.savefig(figures / "F02_cross_feedback_complex_plane.png", dpi=220)
    plt.close()

    sweep = run_root.parent / "20260925T204017_549c3c07_h4_crossmode_v3" / "derived"
    local_path = sweep / "TX4_PHYSICAL_LOCAL_FACTORS.csv"
    collective_path = sweep / "TX4_CONTEXTUAL_RETURN_SWEEP.csv"
    if local_path.exists() and collective_path.exists():
        local_rows = list(csv.DictReader(local_path.open(encoding="utf-8")))
        local_by_g: dict[float, float] = {}
        for item in local_rows:
            g = float(item["g"])
            local_by_g[g] = min(local_by_g.get(g, float("inf")), float(item["physical_local_sigma_min"]))
        collective_rows = list(csv.DictReader(collective_path.open(encoding="utf-8")))
        gains = sorted(set(local_by_g).intersection(float(item["g"]) for item in collective_rows))
        local = [local_by_g[g] for g in gains]
        collective = [float(next(item["collective_sigma_min"] for item in collective_rows if float(item["g"]) == g)) for g in gains]
        plt.figure(figsize=(6.2, 4.0))
        plt.semilogy(gains, local, "o-", label=r"physical local $I+M_{ii}$")
        plt.semilogy(gains, collective, "o-", label=r"collective $I+Q_H$")
        plt.xlabel(r"$g$")
        plt.ylabel("smallest singular value")
        plt.legend()
        plt.tight_layout()
        plt.savefig(figures / "F03_local_collective_boundary_overlay.png", dpi=220)
        plt.close()


def write_blocked_report(run_root: Path, baseline: Path, probe: dict, manifest: dict) -> None:
    write_json(run_root / "claims" / "M1_STATUS.json", probe)
    write_json(run_root / "environment" / "runtime.json", manifest)
    (run_root / "report" / "RUN_SUMMARY.md").write_text(
        "# H4 modal-Schur run\n\n"
        "Status: **BLOCKED_M1_STRICT**\n\n"
        f"Baseline: `{baseline}`\n\n"
        "The run stopped before M2/M3 because the retained port operator did not meet the preregistered strict scalar-Schur residual. No figures were generated.\n\n"
        f"- full DAE pole: `{probe['lambda_full']}`\n"
        f"- `abs(h_k(lambda_c))`: `{probe['h_at_full_pole_abs']}`\n"
        f"- strict threshold: `{probe['m1_threshold']}`\n"
        f"- independent Schur root: `{probe['root_eta_1']}`\n"
        f"- root distance to full pole: `{probe['root_distance_to_full']}`\n"
        f"- complement condition: `{probe['complement_condition_at_full_pole']}`\n\n"
        "Interpretation: this is a retained-operator numerical gate failure, not evidence against the collective mechanism. The next safe action is to reconcile the port linearization with the full DAE before running eta continuation.\n",
        encoding="utf-8",
    )


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_modal_schur.py <repo-root> <run-root>")
    repo = Path(sys.argv[1]).resolve()
    run_root = Path(sys.argv[2]).resolve()
    for name in ("raw/python", "raw/julia", "raw/reconciliation", "tables", "figures", "claims", "report", "environment", "logs"):
        (run_root / name).mkdir(parents=True, exist_ok=True)

    experiment = repo / "reports/poster/ias2026/research/experiments"
    source = repo / "reports/poster/ias2026/research/src"
    sys.path.insert(0, str(experiment))
    sys.path.insert(0, str(source))
    import tx4_contextual_return as tx  # type: ignore  # noqa: PLC0415

    baseline = repo / "reports/poster/ias2026/research/bnd_h4_mechanism/results/20260925T204017_549c3c07_h4_crossmode_v3"
    case = tx.tx4_case(tx.CORE, tx.P4)
    base = tx.tx4_case((), tx.Theta(g=0.0, k=1.425, t=1.5, h=1.0))
    mode = tx.mode_of(case)
    space = tx.build_action_space(base, case, tx.CORE)
    frozen_operator = np.eye(space.dimension, dtype=complex) + space.m(mode.eig)
    left_basis, right_basis, operator_smin = make_basis(frozen_operator)
    h_full, _, h_matrix = scalar_schur(space, left_basis, right_basis, mode.eig, 1.0)
    root_eta_1, root_ok, root_residual, root_evals = solve_complex_root(
        lambda s: scalar_schur(space, left_basis, right_basis, s, 1.0)[0], mode.eig
    )
    complement_condition = float(np.linalg.cond(h_matrix[1:, 1:]))
    probe = {
        "lambda_full": [mode.eig.real, mode.eig.imag],
        "frequency_full_hz": mode.frequency_hz,
        "operator_dimension": space.dimension,
        "operator_singular_value_at_full_pole": operator_smin,
        "h_at_full_pole_abs": abs(h_full),
        "m1_threshold": 1e-8,
        "root_eta_1": [root_eta_1.real, root_eta_1.imag],
        "root_eta_1_residual": root_residual,
        "root_distance_to_full": abs(root_eta_1 - mode.eig),
        "root_solver_success": root_ok,
        "root_solver_evaluations": root_evals,
        "complement_condition_at_full_pole": complement_condition,
    }
    write_json(run_root / "tables" / "M1_SCALAR_SCHUR_PROBE.json", probe)
    manifest = {
        "run_root": str(run_root),
        "baseline_run": str(baseline),
        "git_commit": git(repo, "rev-parse", "HEAD"),
        "git_dirty": bool(git(repo, "status", "--porcelain")),
        "python": {"version": sys.version, "implementation": platform.python_implementation(), "executable": sys.executable},
        "packages": {"numpy": np.__version__},
        "source_hashes": {
            "tx4_contextual_return.py": sha256(experiment / "tx4_contextual_return.py"),
            "port_admittance.py": sha256(source / "ibr_cycles/models/port_admittance.py"),
        },
        "output_policy": "all files are beneath this immutable run root; no zip is produced",
    }
    strict_pass = abs(h_full) < 1e-8 and root_ok and abs(root_eta_1 - mode.eig) <= 1e-6 and np.isfinite(complement_condition) and complement_condition < 1e10
    if not strict_pass:
        probe["status"] = "BLOCKED_M1_STRICT"
        write_blocked_report(run_root, baseline, probe, manifest)
        return 2

    eta_values = np.linspace(1.0, 0.0, 21)
    continuation = []
    previous = root_eta_1
    previous_vector = tracked_vector(space, left_basis, right_basis, previous, 1.0)
    for eta in eta_values:
        if eta == 1.0:
            candidate = previous
            success, residual, evaluations = root_ok, root_residual, root_evals
        else:
            candidate, success, residual, evaluations = solve_complex_root(
                lambda s, eta=eta: scalar_schur(space, left_basis, right_basis, s, eta)[0], previous
            )
        vector = tracked_vector(space, left_basis, right_basis, candidate, float(eta))
        h_value, cross, _ = scalar_schur(space, left_basis, right_basis, candidate, float(eta))
        continuation.append(
            {
                "eta": float(eta),
                "root_real_s-1": candidate.real,
                "root_imag_rad_s-1": candidate.imag,
                "root_frequency_hz": abs(candidate.imag) / (2.0 * np.pi),
                "h_residual": abs(h_value),
                "solver_success": bool(success),
                "solver_evaluations": evaluations,
                "root_mac_previous": mac(previous_vector, vector),
                "cross_term_real": cross.real,
                "cross_term_imag": cross.imag,
            }
        )
        previous, previous_vector = candidate, vector
    continuation.sort(key=lambda row: row["eta"])
    write_csv(run_root / "tables" / "H4_MODAL_SCHUR_ETA.csv", continuation)

    derivative_rows = []
    for row in continuation[1:-1]:
        eta = float(row["eta"])
        s0 = complex(row["root_real_s-1"], row["root_imag_rad_s-1"])
        delta_s = 1e-5 * max(1.0, abs(s0))
        ds = (scalar_schur(space, left_basis, right_basis, s0 + delta_s, eta)[0] - scalar_schur(space, left_basis, right_basis, s0 - delta_s, eta)[0]) / (2.0 * delta_s)
        _, cross, _ = scalar_schur(space, left_basis, right_basis, s0, eta)
        analytic = cross / ds
        delta_eta = 1e-3
        plus, plus_ok, _, _ = solve_complex_root(lambda s: scalar_schur(space, left_basis, right_basis, s, eta + delta_eta)[0], s0)
        minus, minus_ok, _, _ = solve_complex_root(lambda s: scalar_schur(space, left_basis, right_basis, s, eta - delta_eta)[0], s0)
        finite = (plus - minus) / (2.0 * delta_eta)
        derivative_rows.append({"eta": eta, "dlam_deta_analytic_real": analytic.real, "dlam_deta_analytic_imag": analytic.imag, "dlam_deta_fd_real": finite.real, "dlam_deta_fd_imag": finite.imag, "absolute_error": abs(analytic - finite), "finite_difference_success": bool(plus_ok and minus_ok)})
    write_csv(run_root / "tables" / "H4_MODAL_SCHUR_DERIVATIVE.csv", derivative_rows)
    plot_outputs(run_root, continuation, repo)
    strong = float(continuation[0]["root_real_s-1"]) < 0.0 and float(continuation[-1]["root_real_s-1"]) > 0.0
    summary = {"status": "PASS", "strong_collective_hypothesis": "PASS" if strong else "REFUTED", "m1": probe, "eta_rows": len(continuation), "derivative_rows": len(derivative_rows), "max_derivative_error": max(row["absolute_error"] for row in derivative_rows), "min_adjacent_mac": min(row["root_mac_previous"] for row in continuation[1:])}
    write_json(run_root / "claims" / "M1_M3_SUMMARY.json", summary)
    write_json(run_root / "environment" / "runtime.json", manifest)
    (run_root / "report" / "RUN_SUMMARY.md").write_text("# H4 modal-Schur run\n\nStatus: **PASS**\n\n" + json.dumps(summary, indent=2, default=float) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
