# ruff: noqa: E501  -- test tables kept on one line
"""EMT02 diagnostic D3 (tx3-analysis): phasor reference WITH the line's electromagnetic dynamics.

Same machine test as EMT02_03_phasor_refs.py, but the T--INF line current is a dynamic phasor
(synchronous frame, R = 0): (X/w0) dI/dt = V_T - V_INF - jX I, with V_T = E' - I/Y_sg (no shunt at T)
and V_INF = E + z_s I. This is the positive-sequence content of the EMT line; the quasi-static
reference sets dI/dt = 0. Diagnostic only (not a preregistered gate).
Writes results/EMT02/dynphasor_ref.npz.
"""

from __future__ import annotations

import os
import sys
from dataclasses import replace
from pathlib import Path

# single-threaded BLAS (Radau LU): bit-reproducible reference
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import numpy as np  # noqa: E402
from EMT02_03_phasor_refs import (  # noqa: E402
    BUS,
    E0,
    S_DEV,
    Z_L,
    Z_S,
    integrate,
    operating_point,
)

import _bootstrap  # noqa: E402,F401
from ibr_cycles.models.ieee39_case import (  # noqa: E402
    _controller_payload,
    _machine_parameters,
)
from ibr_cycles.models.ieee39_devices import SynchronousMachine  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402

W0 = 2 * np.pi * 60
RESEARCH = HERE.parents[1]


def main() -> int:
    vt, vinf, i0 = operating_point()
    net = load_network()
    p = _machine_parameters(net, BUS, _controller_payload(net.config_path or None))
    p = replace(p, ka=p.ka * 1.425, ta=p.ta * 1.5)
    w = net.machines[BUS]["Sn"] / 100.0
    dev, x0 = SynchronousMachine(bus=BUS, parameters=p, weight=w).initialize(vt, S_DEV)
    pm0 = dev.parameters.pm
    y_sg = w / complex(p.ra, p.xd1)
    X = Z_L.imag

    def f(y, pm):
        x, il = y[:7], complex(y[7], y[8])
        epr = complex(x[2], -x[3]) * np.exp(1j * x[0])
        v_t = epr - il / y_sg
        v_inf = E0 + Z_S * il
        d = replace(dev, parameters=replace(dev.parameters, pm=pm))
        dil = (W0 / X) * (v_t - v_inf - 1j * X * il)
        return np.concatenate([d.derivatives(x, v_t), [dil.real, dil.imag]])

    y0 = np.concatenate([x0, [i0.real * w / w, i0.imag]])
    # steady-state line current = device current on the system base
    epr0 = complex(x0[2], -x0[3]) * np.exp(1j * x0[0])
    il0 = y_sg * (epr0 - vt)
    y0[7], y0[8] = il0.real, il0.imag
    t, y = integrate(f, y0, [(0.0, 1.0, pm0), (1.0, 1.2, pm0 * 1.02), (1.2, 10.0, pm0)])
    np.savez_compressed(RESEARCH / "results" / "EMT02" / "dynphasor_ref.npz", t=t, x=y[:, :7])
    print("dynamic-phasor reference written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
