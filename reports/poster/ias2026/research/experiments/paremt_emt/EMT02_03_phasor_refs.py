# ruff: noqa: E501  -- test tables kept on one line
"""EMT02/EMT03 phasor references (tx3-analysis; canonical device classes; prereg section 7).

Test system (identical on both sides): device at bus T -- line R = 0, X = 0.05 pu (system base) --
bus INF; the infinite bus is a Thevenin source E = 1.0 at 0 rad behind z_s = j1e-4 pu. Quasi-static
phasor network: V_T = E + (z_s + j0.05) I_dev,sys. Operating point: device output S = 4.0 + j1.0 pu at
T; bus-30 device data (Sn 1040 MVA, w = 10.4); P4 machine settings (k = 1.425, t = 1.5).
Integrator: scipy Radau, rtol 1e-10, atol 1e-12, piecewise across events, output every 1 ms.

EMT02 SG: Pm +2 % for 0.2 s at t = 1 s; 10 s.
EMT03 GFL, g in {0, 0.03625, 0.25}: (a) p_ref +1 % at 1 s; (b) |E| -1 % at 1 s; (c) E phase +0.02 rad
at 1 s; 5 s.
Writes results/EMT02/phasor_ref.npz, results/EMT03/phasor_ref.npz and results/EMT0x/op.json.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import replace
from pathlib import Path

# single-threaded BLAS (Radau LU): bit-reproducible references
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import _bootstrap  # noqa: E402,F401
from _f7_common import LEAK  # noqa: E402
from ibr_cycles.models.ieee39_case import (  # noqa: E402
    _controller_payload,
    _machine_parameters,
)
from ibr_cycles.models.ieee39_devices import (  # noqa: E402
    ConverterParameters,
    GridFollowingConverter,
    SynchronousMachine,
)
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402

RESEARCH = HERE.parents[1]
S_DEV = 4.0 + 1.0j
Z_S = 1e-4j
Z_L = 0.05j
E0 = 1.0 + 0j
BUS = 30
DT_OUT = 1e-3


def operating_point():
    z = Z_S + Z_L
    vt = E0
    for _ in range(200):
        i = np.conj(S_DEV / vt)
        vt = E0 + z * i
    i = np.conj(S_DEV / vt)
    return vt, E0 + Z_S * i, i


def integrate(f, x0, segments):
    ts, xs = [], []
    x = np.asarray(x0, float)
    for (t0, t1, param) in segments:
        tev = np.arange(np.ceil(t0 / DT_OUT - 1e-9), np.floor(t1 / DT_OUT + 1e-9) + 1) * DT_OUT
        sol = solve_ivp(lambda t, y, p=param: f(y, p), (t0, t1), x, method="Radau", rtol=1e-10, atol=1e-12, t_eval=tev)
        assert sol.success, sol.message
        keep = slice(1, None) if ts else slice(None)
        ts.append(sol.t[keep])
        xs.append(sol.y[:, keep])
        x = sol.y[:, -1]
    return np.concatenate(ts), np.concatenate(xs, axis=1).T


def sg_ref(vt):
    net = load_network()
    payload = _controller_payload(net.config_path or None)
    p = _machine_parameters(net, BUS, payload)
    p = replace(p, ka=p.ka * 1.425, ta=p.ta * 1.5)
    w = net.machines[BUS]["Sn"] / 100.0
    dev, x0 = SynchronousMachine(bus=BUS, parameters=p, weight=w).initialize(vt, S_DEV)
    pm0 = dev.parameters.pm
    y_sg = w / complex(p.ra, p.xd1)
    z = Z_S + Z_L

    def vterm(x, e):
        epr = complex(x[2], -x[3]) * np.exp(1j * x[0])
        return (e + z * y_sg * epr) / (1 + z * y_sg)

    def f(x, pm):
        d = replace(dev, parameters=replace(dev.parameters, pm=pm))
        return d.derivatives(x, vterm(x, E0))

    t, x = integrate(f, x0, [(0.0, 1.0, pm0), (1.0, 1.2, pm0 * 1.02), (1.2, 10.0, pm0)])
    vts = np.array([vterm(xi, E0) for xi in x])
    s = np.array([vi * np.conj(dev.injection(xi, vi)) for xi, vi in zip(x, vts, strict=True)])
    return t, x, vts, s, {"pm": pm0, "vref": dev.parameters.vref, "w": w}


def gfl_ref(vt, g, case):
    w = load_network().machines[BUS]["Sn"] / 100.0
    par = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK)
    dev, x0 = GridFollowingConverter(bus=BUS, parameters=par, weight=w).initialize(vt, S_DEV)
    pref0 = dev.parameters.p_ref
    z = Z_S + Z_L

    def e_of(ev):
        if case == "b" and ev:
            return E0 * 0.99
        if case == "c" and ev:
            return E0 * np.exp(0.02j)
        return E0

    def vterm(x, ev):
        i = dev.injection(x, 0j)
        return e_of(ev) + z * i

    def f(x, ev):
        pr = pref0 * 1.01 if (case == "a" and ev) else pref0
        d = replace(dev, parameters=replace(dev.parameters, p_ref=pr))
        return d.derivatives(x, vterm(x, ev))

    t, x = integrate(f, x0, [(0.0, 1.0, False), (1.0, 5.0, True)])
    vts = np.array([vterm(xi, ti >= 1.0) for xi, ti in zip(x, t, strict=True)])
    s = np.array([vi * np.conj(dev.injection(xi, vi)) for xi, vi in zip(x, vts, strict=True)])
    return t, x, vts, s, {"p_ref": pref0, "q_ref": dev.parameters.q_ref, "v_ref": dev.parameters.v_ref, "w": w}


def main() -> int:
    vt, vinf, i = operating_point()
    op = {"V_T": [vt.real, vt.imag], "V_INF": [vinf.real, vinf.imag], "S_dev": [S_DEV.real, S_DEV.imag],
          "z_s": [0.0, 1e-4], "z_line": [0.0, 0.05], "E": [1.0, 0.0], "bus_data": BUS}
    out2 = RESEARCH / "results" / "EMT02"
    out3 = RESEARCH / "results" / "EMT03"
    out2.mkdir(parents=True, exist_ok=True)
    out3.mkdir(parents=True, exist_ok=True)
    t, x, vts, s, meta = sg_ref(vt)
    np.savez_compressed(out2 / "phasor_ref.npz", t=t, x=x, vt=vts, s=s)
    (out2 / "op.json").write_text(json.dumps({**op, **meta}, indent=1), encoding="utf-8")
    refs = {}
    for g in (0.0, 0.03625, 0.25):
        for case in ("a", "b", "c"):
            t, x, vts, s, meta = gfl_ref(vt, g, case)
            key = f"g{g:g}_{case}"
            refs[f"{key}_t"], refs[f"{key}_x"], refs[f"{key}_vt"], refs[f"{key}_s"] = t, x, vts, s
            op[f"meta_{key}"] = meta
    np.savez_compressed(out3 / "phasor_ref.npz", **refs)
    (out3 / "op.json").write_text(json.dumps(op, indent=1), encoding="utf-8")
    print("references written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
