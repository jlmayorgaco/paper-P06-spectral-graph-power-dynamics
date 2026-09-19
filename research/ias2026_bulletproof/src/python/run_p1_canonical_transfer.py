"""Generate the canonical Python 2x2 terminal transfer map for P1."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np


CAMPAIGN = Path(__file__).resolve().parents[2]
RAW = CAMPAIGN / "raw" / "gfl11"
SOURCE = RAW / "canonical_source" / "ieee39_devices.py"


def load_canonical():
    spec = importlib.util.spec_from_file_location("canonical_ieee39_devices_transfer", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(SOURCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    canonical = load_canonical()
    base = canonical.ConverterParameters(
        voltage_control=True, voltage_gain=0.03625, voltage_leak=0.05
    )
    device, xeq = base_device = canonical.GridFollowingConverter(
        bus=39, parameters=base, weight=0.75
    ).initialize(1.0 + 0.0j, 0.35 - 0.08j)
    xeq = np.asarray(xeq, dtype=float)
    ueq = np.array([1.0, 0.0], dtype=float)
    state_step = 1e-7
    input_step = 1e-7

    def response(x: np.ndarray, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        terminal = complex(float(u[0]), float(u[1]))
        f = np.asarray(device.derivatives(x, terminal), dtype=float)
        current = complex(device.injection(x, terminal))
        return f, np.array([current.real, current.imag], dtype=float)

    n = len(xeq)
    a_mat = np.zeros((n, n), dtype=float)
    b_mat = np.zeros((n, 2), dtype=float)
    c_mat = np.zeros((2, n), dtype=float)
    d_mat = np.zeros((2, 2), dtype=float)
    for j in range(n):
        h = state_step * max(1.0, abs(xeq[j]))
        xp, xm = xeq.copy(), xeq.copy()
        xp[j] += h
        xm[j] -= h
        fp, yp = response(xp, ueq)
        fm, ym = response(xm, ueq)
        a_mat[:, j] = (fp - fm) / (2.0 * h)
        c_mat[:, j] = (yp - ym) / (2.0 * h)
    for j in range(2):
        h = input_step * max(1.0, abs(ueq[j]))
        up, um = ueq.copy(), ueq.copy()
        up[j] += h
        um[j] -= h
        fp, yp = response(xeq, up)
        fm, ym = response(xeq, um)
        b_mat[:, j] = (fp - fm) / (2.0 * h)
        d_mat[:, j] = (yp - ym) / (2.0 * h)

    freqs = np.unique(np.concatenate((10.0 ** np.linspace(-2, 2, 401), np.geomspace(0.2, 2.0, 201))))
    out = RAW / "p1_transfer_canonical_python.csv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["frequency_hz", "sigma_min", "condition", "h11_re", "h11_im", "h12_re", "h12_im", "h21_re", "h21_im", "h22_re", "h22_im"]
        )
        eye = np.eye(n, dtype=complex)
        for hz in freqs:
            s = 2j * np.pi * float(hz)
            state_resolvent = s * eye - a_mat
            h_mat = c_mat @ np.linalg.solve(state_resolvent, b_mat) + d_mat
            singular_values = np.linalg.svd(state_resolvent, compute_uv=False)
            writer.writerow(
                [
                    float(hz),
                    float(singular_values[-1]),
                    float(singular_values[0] / singular_values[-1]),
                    h_mat[0, 0].real,
                    h_mat[0, 0].imag,
                    h_mat[0, 1].real,
                    h_mat[0, 1].imag,
                    h_mat[1, 0].real,
                    h_mat[1, 0].imag,
                    h_mat[1, 1].real,
                    h_mat[1, 1].imag,
                ]
            )
    print(f"P1_CANONICAL_TRANSFER_PASS points={len(freqs)} path={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
