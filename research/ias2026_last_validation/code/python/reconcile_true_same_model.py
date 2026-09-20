"""Compare frozen Python and Julia base/V4 state, spectrum, and mode MAC."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from campaign_root import campaign_root_from_argv


ROOT = campaign_root_from_argv()
OUT = ROOT / "raw" / "reconciliation"
GAUGE_TOL = 1e-4


def read_vec(path: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(path, delimiter=","))


def read_julia_eigs(label: str) -> tuple[np.ndarray, np.ndarray]:
    spectrum = np.loadtxt(OUT / f"julia_{label}_spectrum.csv", delimiter=",", skiprows=1, usecols=(1, 2))
    values = spectrum[:, 0] + 1j * spectrum[:, 1]
    real_vectors = np.loadtxt(OUT / f"julia_{label}_eigenvectors_real.csv", delimiter=",")[:, 1:]
    imag_vectors = np.loadtxt(OUT / f"julia_{label}_eigenvectors_imag.csv", delimiter=",")[:, 1:]
    return values, (real_vectors + 1j * imag_vectors).T


def mac(left: np.ndarray, right: np.ndarray) -> float:
    denominator = np.vdot(left, left).real * np.vdot(right, right).real
    return float(abs(np.vdot(left, right)) ** 2 / denominator) if denominator > 0 else 0.0


summary_rows = []
mode_rows = []
claims = {}
labels = ["base", "v4"] + ["census_" + ("none" if not buses else "+".join(str(bus) for bus in buses)) for mask in range(16) for buses in [tuple((30, 33, 35, 37)[index] for index in range(4) if (mask >> index) & 1)]]
for label in labels:
    py = np.load(OUT / f"python_{label}.npz")
    py_values = py["eigenvalues"]
    py_vectors = py["eigenvectors"]
    jl_values, jl_vectors = read_julia_eigs(label)
    row_ind, col_ind = linear_sum_assignment(abs(py_values[:, None] - jl_values[None, :]))
    state_error = float(np.max(np.abs(py["state"] - read_vec(OUT / f"julia_{label}_state.csv"))))
    algebraic_error = float(np.max(np.abs(py["algebraic"] - read_vec(OUT / f"julia_{label}_algebraic.csv"))))
    errors = []
    for i, j in zip(row_ind, col_ind, strict=True):
        delta = py_values[i] - jl_values[j]
        errors.append((i, j, delta, mac(py_vectors[:, i], jl_vectors[:, j])))
        if abs(py_values[i]) >= GAUGE_TOL and abs(py_values[i].imag) / (2.0 * np.pi) >= 0.3 and abs(py_values[i].imag) / (2.0 * np.pi) <= 1.5:
            mode_rows.append([label, i, j, py_values[i].real, jl_values[j].real, py_values[i].imag / (2.0 * np.pi), jl_values[j].imag / (2.0 * np.pi), abs(delta), mac(py_vectors[:, i], jl_vectors[:, j])])
    py_transverse = py_values[abs(py_values) >= GAUGE_TOL]
    jl_transverse = jl_values[abs(jl_values) >= GAUGE_TOL]
    py_critical = py_transverse[np.argmax(py_transverse.real)]
    jl_critical = jl_transverse[np.argmax(jl_transverse.real)]
    critical_match = min(errors, key=lambda item: abs(py_values[item[0]] - jl_critical))
    critical_mac = critical_match[3]
    claims[label] = {
        "state_max_abs_error": state_error,
        "algebraic_max_abs_error": algebraic_error,
        "critical_real_error": float(abs(py_critical.real - jl_critical.real)),
        "critical_frequency_error_hz": float(abs(abs(py_critical.imag) - abs(jl_critical.imag)) / (2.0 * np.pi)),
        "critical_mode_mac": critical_mac,
        "state_pass": bool(state_error < 1e-6),
        "algebraic_pass": bool(algebraic_error < 1e-6),
        "critical_real_pass": bool(abs(py_critical.real - jl_critical.real) < 1e-4),
        "critical_frequency_pass": bool(abs(abs(py_critical.imag) - abs(jl_critical.imag)) / (2.0 * np.pi) < 1e-4),
        "critical_mac_pass": bool(critical_mac > 0.95),
        "python_alpha_transverse": float(py_transverse.real.max()),
        "julia_alpha_transverse": float(jl_transverse.real.max()),
    }
    summary_rows.append([label, len(py["state"]), len(read_vec(OUT / f"julia_{label}_state.csv")), state_error, algebraic_error, py_critical.real, jl_critical.real, abs(py_critical.imag) / (2.0 * np.pi), abs(jl_critical.imag) / (2.0 * np.pi), critical_mac])
    with (OUT / f"{label}_spectrum_reconciliation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["python_index", "julia_index", "python_real", "julia_real", "python_imag", "julia_imag", "abs_eigenvalue_error", "mac"])
        for i, j, delta, value_mac in errors:
            writer.writerow([i, j, py_values[i].real, jl_values[j].real, py_values[i].imag, jl_values[j].imag, abs(delta), value_mac])

with (OUT / "same_model_reconciliation.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["case", "python_states", "julia_states", "state_max_abs_error", "algebraic_max_abs_error", "python_critical_real", "julia_critical_real", "python_critical_frequency_hz", "julia_critical_frequency_hz", "critical_mac"])
    writer.writerows(summary_rows)
with (OUT / "same_model_mode_family.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["case", "python_index", "julia_index", "python_real", "julia_real", "python_frequency_hz", "julia_frequency_hz", "abs_eigenvalue_error", "mac"])
    writer.writerows(mode_rows)
(OUT / "same_model_reconciliation.json").write_text(json.dumps(claims, indent=2) + "\n", encoding="utf-8")
for label, claim in claims.items():
    print(f"RECONCILIATION_{label.upper()} state={claim['state_max_abs_error']:.3e} algebraic={claim['algebraic_max_abs_error']:.3e} real={claim['critical_real_error']:.3e} freq={claim['critical_frequency_error_hz']:.3e} mac={claim['critical_mode_mac']:.6f}")
