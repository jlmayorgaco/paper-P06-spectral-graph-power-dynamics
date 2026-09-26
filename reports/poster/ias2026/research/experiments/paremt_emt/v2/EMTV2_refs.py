# ruff: noqa: E501  -- test tables kept on one line
"""V2 phasor references (tx3-analysis; canonical device classes; prereg V2 section 3).

Unit-test system: device at T -- line X = 0.05 (R = 0) -- INF; INF = E behind z_s = j1e-4.
  quasi-static network : V_T = V_INF + jX I                         (the frozen phasor model)
  dynamic-line network : V_T = V_INF + (X/w0) dI/dt + jX I          (positive-sequence EMT line)
SG  : I is the line state (as v1 D3), V_T = E' - I/Y_sg.
GFL : I = w i e^{j theta} (device state), V_T solved by fixed-point iteration to 1e-14 per RHS.
Radau rtol 1e-10, atol 1e-12, piecewise across events.
Writes results/EMTV2/refs/*.npz (1-ms grids, committed) and external/paremt_runs/v2_refs/*.npz
(G4b dynamic-line references on the 50-us grid, hashed in results/EMTV2/refs/refs_manifest.json).
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from dataclasses import replace
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

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

RESEARCH = HERE.parents[2]
REPO = HERE.parents[6]
OUT = RESEARCH / "results" / "EMTV2" / "refs"
BIG = REPO / "external" / "paremt_runs" / "v2_refs"
W0 = 2 * np.pi * 60
Z_S, X_L, E0 = 1e-4j, 0.05, 1.0 + 0j
SG_CASES = {"R-SG30": (30, 4.0 + 1.0j), "B-SG36": (36, 5.79999998 + 0.60708066j)}
GFL_CASES = {"R-GFL30": (30, 4.0 + 1.0j), "B1-GFL35": (35, 2.0 + 0.5j), "B2-GFL37": (37, 3.21521338 - 0.27617116j)}
GS = (0.0, 0.03625, 0.25)


def op_point(s):
    vt = E0
    for _ in range(200):
        vt = E0 + (Z_S + 1j * X_L) * np.conj(s / vt)
    return vt, E0 + Z_S * np.conj(s / vt)


def integrate(f, x0, segments, dt_out):
    ts, xs = [], []
    x = np.asarray(x0, float)
    for (t0, t1, param) in segments:
        tev = np.arange(np.ceil(t0 / dt_out - 1e-9), np.floor(t1 / dt_out + 1e-9) + 1) * dt_out
        sol = solve_ivp(lambda t, y, p=param: f(y, p), (t0, t1), x, method="Radau", rtol=1e-10, atol=1e-12, t_eval=tev)
        assert sol.success, sol.message
        keep = slice(1, None) if ts else slice(None)
        ts.append(sol.t[keep])
        xs.append(sol.y[:, keep])
        x = sol.y[:, -1]
    return np.concatenate(ts), np.concatenate(xs, axis=1).T


def sg_refs(bus, s):
    net = load_network()
    p = _machine_parameters(net, bus, _controller_payload(net.config_path or None))
    p = replace(p, ka=p.ka * 1.425, ta=p.ta * 1.5)
    w = net.machines[bus]["Sn"] / 100.0
    vt, _ = op_point(s)
    dev, x0 = SynchronousMachine(bus=bus, parameters=p, weight=w).initialize(vt, s)
    pm0 = dev.parameters.pm
    y_sg = w / complex(p.ra, p.xd1)
    z = Z_S + 1j * X_L
    segs = [(0.0, 1.0, pm0), (1.0, 1.2, pm0 * 1.02), (1.2, 10.0, pm0)]

    def dmach(pm):
        return replace(dev, parameters=replace(dev.parameters, pm=pm))

    def vterm_qs(x):
        epr = complex(x[2], -x[3]) * np.exp(1j * x[0])
        return (E0 + z * y_sg * epr) / (1 + z * y_sg)

    t, x = integrate(lambda x, pm: dmach(pm).derivatives(x, vterm_qs(x)), x0, segs, 1e-3)
    vts = np.array([vterm_qs(xi) for xi in x])
    qs = {"t": t, "x": x, "vt": vts, "s": np.array([vi * np.conj(dev.injection(xi, vi)) for xi, vi in zip(x, vts, strict=True)])}

    def fdyn(y, pm):
        xm, il = y[:7], complex(y[7], y[8])
        epr = complex(xm[2], -xm[3]) * np.exp(1j * xm[0])
        v_t = epr - il / y_sg
        dil = (W0 / X_L) * (v_t - (E0 + Z_S * il) - 1j * X_L * il)
        return np.concatenate([dmach(pm).derivatives(xm, v_t), [dil.real, dil.imag]])

    epr0 = complex(x0[2], -x0[3]) * np.exp(1j * x0[0])
    il0 = y_sg * (epr0 - vt)
    t2, y2 = integrate(fdyn, np.concatenate([x0, [il0.real, il0.imag]]), segs, 1e-3)
    il = y2[:, 7] + 1j * y2[:, 8]
    vt2 = np.array([complex(r[2], -r[3]) * np.exp(1j * r[0]) for r in y2]) - il / y_sg
    dyn = {"t": t2, "x": y2[:, :7], "vt": vt2, "s": vt2 * np.conj(il)}
    meta = {"bus": bus, "w": w, "S": [s.real, s.imag], "V_T": [vt.real, vt.imag], "pm": pm0, "vref": dev.parameters.vref}
    return qs, dyn, meta


def gfl_refs(bus, s, g, case, dyn_dt):
    w = load_network().machines[bus]["Sn"] / 100.0
    vt, _ = op_point(s)
    par = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK)
    dev, x0 = GridFollowingConverter(bus=bus, parameters=par, weight=w).initialize(vt, s)
    pref0 = dev.parameters.p_ref

    def e_of(ev):
        if case == "b" and ev:
            return E0 * 0.99
        if case == "c" and ev:
            return E0 * np.exp(0.02j)
        return E0

    def dconv(ev):
        pr = pref0 * 1.01 if (case == "a" and ev) else pref0
        return replace(dev, parameters=replace(dev.parameters, p_ref=pr))

    def vterm_qs(x, ev):
        return e_of(ev) + (Z_S + 1j * X_L) * dev.injection(x, 0j)

    segs = [(0.0, 1.0, False), (1.0, 5.0, True)]
    t, x = integrate(lambda x, ev: dconv(ev).derivatives(x, vterm_qs(x, ev)), x0, segs, 1e-3)
    vts = np.array([vterm_qs(xi, ti >= 1.0) for xi, ti in zip(x, t, strict=True)])
    qs = {"t": t, "x": x, "vt": vts, "s": np.array([vi * np.conj(dev.injection(xi, vi)) for xi, vi in zip(x, vts, strict=True)])}

    def vterm_dyn(x, ev):
        d = dconv(ev)
        i_sys = dev.injection(x, 0j)
        idq = complex(x[6], x[7])
        rot = np.exp(1j * x[0])
        v = e_of(ev) + (Z_S + 1j * X_L) * i_sys
        for it in range(200):
            dx = d.derivatives(x, v)
            di_sys = w * rot * (complex(dx[6], dx[7]) + 1j * dx[0] * idq)
            vn = e_of(ev) + Z_S * i_sys + (X_L / W0) * di_sys + 1j * X_L * i_sys
            if abs(vn - v) <= 1e-14:
                return vn, dx, it
            v = vn
        raise RuntimeError("dynamic-line fixed point did not converge")

    def fdyn(x, ev):
        _, dx, _ = vterm_dyn(x, ev)
        return dx

    t2, x2 = integrate(fdyn, x0, segs, dyn_dt)
    vt2 = np.array([vterm_dyn(xi, ti >= 1.0)[0] for xi, ti in zip(x2, t2, strict=True)])
    s2 = vt2 * np.conj(w * (x2[:, 6] + 1j * x2[:, 7]) * np.exp(1j * x2[:, 0]))
    dyn = {"t": t2, "x": x2, "vt": vt2, "s": s2}
    meta = {"bus": bus, "w": w, "S": [s.real, s.imag], "V_T": [vt.real, vt.imag], "p_ref": pref0, "q_ref": dev.parameters.q_ref,
            "v_ref": dev.parameters.v_ref}
    return qs, dyn, meta


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    BIG.mkdir(parents=True, exist_ok=True)
    man = {"files": {}, "meta": {}}
    for cid, (bus, s) in SG_CASES.items():
        qs, dyn, meta = sg_refs(bus, s)
        for kind, d in (("qs", qs), ("dyn", dyn)):
            p = OUT / f"{cid}_{kind}.npz"
            np.savez_compressed(p, **d)
            man["files"][str(p.relative_to(RESEARCH)).replace("\\", "/")] = sha(p)
        man["meta"][cid] = meta
        print(cid, "done")
    for cid, (bus, s) in GFL_CASES.items():
        qs_all, dyn_big, dyn_ms = {}, {}, {}
        for g in GS:
            for case in "abc":
                key = f"g{g:g}_{case}"
                qs, dyn, meta = gfl_refs(bus, s, g, case, 50e-6)
                for k, v in qs.items():
                    qs_all[f"{key}_{k}"] = v
                for k, v in dyn.items():
                    dyn_big[f"{key}_{k}"] = v
                    dyn_ms[f"{key}_{k}"] = v[::20]
                man["meta"][f"{cid}_{key}"] = meta
                print(cid, key, "done", flush=True)
        for name, d, base in ((f"{cid}_qs.npz", qs_all, OUT), (f"{cid}_dyn_1ms.npz", dyn_ms, OUT), (f"{cid}_dyn_50us.npz", dyn_big, BIG)):
            p = base / name
            np.savez_compressed(p, **d)
            rel = str(p.relative_to(REPO if base == BIG else RESEARCH)).replace("\\", "/")
            man["files"][rel] = sha(p)
    (OUT / "refs_manifest.json").write_text(json.dumps(man, indent=1, default=float), encoding="utf-8")
    print("references written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
