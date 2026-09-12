# ruff: noqa: E501  -- test tables kept on one line
"""Gate SW (prereg V2 section 4; derivation section 9): software tests of the GFL voltage interface.

T1 identity: physical filter (P) driven by e_cmd of (C) reproduces f_i,TX4; e_cmd = e_TX4 + j(xf/w0) theta' i.
T2 abc-level: per-phase physical filter driven by (C) (solve_ivp, rtol 1e-11) vs dq ODE di/dt = f_i,TX4.
T3 companion: (Req, icf, Gv1) equal ParaEMT's numba_InitNet coefficients for the same R-L branch.
T4 static bases/signs at t = 0 for the three G4 cases.
Writes results/EMTV2/SW/sw_tests.json. Run with .venv/xtool-paremt.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import _emt  # noqa: E402
import numpy as np  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402
import tx4_emt_v2 as K2  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

OUT = _emt.RESULTS / "EMTV2" / "SW"
W0 = K.W0
A = np.exp(2j * np.pi / 3)
Z_S, X_L, E0 = 1e-4j, 0.05, 1.0 + 0j
CASES = {"R-GFL30": (10.4, 4.0 + 1.0j), "B1-GFL35": (10.857, 2.0 + 0.5j), "B2-GFL37": (9.702, 3.21521338 - 0.27617116j)}


def op_point(s):
    vt = E0
    for _ in range(200):
        vt = E0 + (Z_S + 1j * X_L) * np.conj(s / vt)
    return vt


def par_row(w, g, vt, s):
    par = dict(C.GFL_DEFAULTS)
    x, p_ref, q_ref, v_ref = C.gfl_init(par, w, vt, s)
    row = np.zeros(len(K.GFL_COLS))
    for n, v in dict(par, bus=0, w=w, g=g, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref).items():
        row[K.GF[n]] = v
    return row, np.array(x)


def e_tx4(x, vv, p, pref):
    """Spec section 4.2 converter voltage, written independently of the kernel."""

    th = x[0]
    vdq = vv * np.exp(-1j * th)
    v = abs(vv)
    idq = complex(x[6], x[7])
    err = p[K.GF["g"]] * (p[K.GF["v_ref"]] - v)
    q_cmd = p[K.GF["kp_v"]] * err + x[10]
    id_ref = p[K.GF["kp_p"]] * (pref - x[2]) + x[4]
    iq_ref = -(p[K.GF["kp_q"]] * (q_cmd - x[3]) + x[5])
    xi = complex(x[8], x[9])
    return vdq + p[K.GF["kp_i"]] * (complex(id_ref, iq_ref) - idq) + xi + 1j * p[K.GF["xf"]] * idq


def t1(rng):
    worst_p, worst_e = 0.0, 0.0
    for n in range(2000):
        name = list(CASES)[n % 3]
        w, s = CASES[name]
        g = (0.0, 0.03625, 0.25)[n % 3 if n < 3 else rng.integers(3)]
        vt = op_point(s)
        p, x0 = par_row(w, g, vt, s)
        x = x0 * (1.0 + 0.05 * rng.standard_normal(11)) + 0.01 * rng.standard_normal(11)
        x[0] = rng.uniform(-np.pi, np.pi)
        vv = rng.uniform(0.9, 1.1) * np.exp(1j * rng.uniform(-np.pi, np.pi))
        vq = (vv * np.exp(-1j * x[0])).imag
        x[1] = rng.uniform(-5, 5) - p[K.GF["kp_pll"]] * vq  # theta' = kp_pll v_q + x_pll in [-5, 5] rad/s
        idq = complex(x[6], x[7])
        pref = p[K.GF_PREF] * rng.uniform(0.95, 1.05)
        e, f = K2.e_command(x, vv, idq, p, pref)
        xf, rf = p[K.GF_XF], p[K.GF_RF]
        vdq = vv * np.exp(-1j * x[0])
        thd = f[0]
        phys = (W0 / xf) * (e - vdq - rf * idq) - 1j * (W0 + thd) * idq
        fi = complex(f[6], f[7])
        scale = (W0 / xf) * (abs(e) + abs(vdq) + rf * abs(idq)) + (W0 + abs(thd)) * abs(idq)
        worst_p = max(worst_p, abs(phys - fi) / scale)
        e_ref = e_tx4(x, vv, p, pref) + 1j * (xf / W0) * thd * idq
        worst_e = max(worst_e, abs(e - e_ref) / abs(e_ref))
    return {"max_rel_err_physical_vs_fi": worst_p, "max_rel_err_ecmd_vs_etx4_plus_correction": worst_e,
            "pass": bool(worst_p <= 1e-12 and worst_e <= 1e-12)}


def t2():
    """Smooth prescribed controller states and terminal voltage; theta and currents are integrated."""

    w, s = CASES["B2-GFL37"]
    vt = op_point(s)
    p, x0 = par_row(w, 0.25, vt, s)
    rf, xf = p[K.GF_RF], p[K.GF_XF]
    rf_s, lf_s = rf / w, xf / (w * W0)

    def xs(t):
        x = x0.copy()
        x[1] = 0.8 * np.sin(2 * np.pi * 1.3 * t)
        x[2:6] = x0[2:6] * (1 + 0.02 * np.sin(2 * np.pi * np.array([0.7, 1.1, 2.3, 3.1]) * t))
        x[8:11] = x0[8:11] * (1 + 0.03 * np.cos(2 * np.pi * np.array([1.7, 2.9, 0.4]) * t))
        return x

    def vv(t):
        return vt * (1 + 0.01 * np.sin(2 * np.pi * 2.0 * t)) * np.exp(1j * 0.02 * np.sin(2 * np.pi * 0.9 * t))

    def rhs_abc(t, y):
        th, ia, ib, ic = y
        x = xs(t)
        x[0] = th
        isys = (2 / 3) * (ia + A * ib + A * A * ic) * np.exp(-1j * W0 * t)
        idq = isys * np.exp(-1j * th) / w
        e, f = K2.e_command(x, vv(t), idq, p, p[K.GF_PREF])
        esys = e * np.exp(1j * th)
        v = vv(t)
        out = [f[0]]
        for k in range(3):
            ang = W0 * t - 2 * np.pi * k / 3
            ek = np.real(esys * np.exp(1j * ang))
            vk = np.real(v * np.exp(1j * ang))
            ik = (ia, ib, ic)[k]
            out.append((ek - vk - rf_s * ik) / lf_s)
        return out

    def rhs_dq(t, y):
        th, idr, iqr = y
        x = xs(t)
        x[0] = th
        x[6], x[7] = idr, iqr
        f, *_ = K._gfl_eval(x, vv(t), p, p[K.GF_PREF])
        return [f[0], f[6], f[7]]

    i0 = complex(x0[6], x0[7])
    isys0 = w * i0 * np.exp(1j * x0[0])
    y0 = [x0[0]] + [np.real(isys0 * np.exp(-1j * 2 * np.pi * k / 3)) for k in range(3)]
    tt = np.linspace(0, 0.2, 401)
    sa = solve_ivp(rhs_abc, (0, 0.2), y0, method="DOP853", rtol=1e-11, atol=1e-13, t_eval=tt)
    sq = solve_ivp(rhs_dq, (0, 0.2), [x0[0], i0.real, i0.imag], method="DOP853", rtol=1e-11, atol=1e-13, t_eval=tt)
    isys = (2 / 3) * (sa.y[1] + A * sa.y[2] + A * A * sa.y[3]) * np.exp(-1j * W0 * tt)
    idq_abc = isys * np.exp(-1j * sa.y[0]) / w
    err = float(np.max(np.abs(idq_abc - (sq.y[1] + 1j * sq.y[2]))))
    excursion = float(np.max(np.abs((sq.y[1] + 1j * sq.y[2]) - i0)))
    return {"max_abs_err_idq": err, "idq_excursion": excursion, "theta_err": float(np.max(np.abs(sa.y[0] - sq.y[0]))),
            "ok": bool(sa.success and sq.success), "pass": bool(sa.success and sq.success and err <= 1e-7)}


def t3():
    rows = []
    for name, (w, _s) in CASES.items():
        for dt in (50e-6, 25e-6):
            rf_s, lf_s, req, icf, gv1 = K2.filter_coeffs(w, 0.01, 0.15, dt)
            net = {"buses": [{"idx": 1}, {"idx": 2}], "lines": [{"bus1": 1, "bus2": 2, "r": rf_s, "x": lf_s * W0, "b": 0.0, "g": 0.0, "tap": 1.0,
                                                                   "phi": 0.0, "u": 1.0, "trans": 0.0}], "shunts": []}
            nt = K.build_network(net, {1: 1.0 + 0j, 2: 1.0 + 0j}, dt, net_damping=0.0)
            c = np.real(nt["coe0"][0])
            rows.append({"case": name, "dt": dt, "Req": [req, float(c[2])], "icf": [icf, float(c[3])], "Gv1": [gv1, float(c[4])],
                         "max_rel": float(max(abs(req - c[2]) / abs(c[2]), abs(icf - c[3]) / abs(c[3]), abs(gv1 - c[4]) / abs(c[4])))})
    return {"rows": rows, "pass": bool(max(r["max_rel"] for r in rows) <= 1e-12)}


def t4():
    rows = []
    for name, (w, s) in CASES.items():
        vt = op_point(s)
        p, x0 = par_row(w, 0.03625, vt, s)
        filt = K2.filter_coeffs(w, 0.01, 0.15, 50e-6)
        i0 = complex(x0[6], x0[7])
        ibr, ihis, esys = K2.filter_init(p, vt, i0, x0[0], filt)
        isys = K2._space_vector(ibr[0], ibr[1], ibr[2], 0.0)
        idq = isys * np.exp(-1j * x0[0]) / w
        e, f = K2.e_command(x0, vt, idq, p, p[K.GF_PREF])
        e_ref = vt * np.exp(-1j * x0[0]) + complex(0.01, 0.15) * i0
        rows.append({"case": name, "idq_err": abs(idq - i0), "S_err": abs(vt * np.conj(isys) - s) / abs(s), "e0_err": abs(e - e_ref),
                     "esys_vs_command": abs(esys - e * np.exp(1j * x0[0]))})
    worst = max(max(r["idq_err"], r["S_err"], r["e0_err"], r["esys_vs_command"]) for r in rows)
    return {"rows": rows, "pass": bool(worst <= 1e-12)}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260930)
    res = {"T1": t1(rng), "T2": t2(), "T3": t3(), "T4": t4(), "provenance": _emt.provenance()}
    res["GATE_SW"] = "PASS" if all(res[k]["pass"] for k in ("T1", "T2", "T3", "T4")) else "FAIL"
    (OUT / "sw_tests.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: (v if k in ("GATE_SW",) else {kk: vv for kk, vv in v.items() if kk != "rows"}) for k, v in res.items() if k != "provenance"}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
