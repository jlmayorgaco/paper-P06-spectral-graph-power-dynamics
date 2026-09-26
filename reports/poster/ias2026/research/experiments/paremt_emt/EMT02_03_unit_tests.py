# ruff: noqa: E501  -- test tables kept on one line
"""EMT02 (two-axis SG) and EMT03 (11-state GFL) unit tests in ParaEMT (prereg section 7, gates G3/G4).

Same test system as EMT02_03_phasor_refs.py, realized in EMT: bus T -- ParaEMT R-L line (R = 0,
X = 0.05 pu) -- bus INF; infinite bus = Thevenin source 1.0 at 0 behind j1e-4 (synchronous-frame
Norton). Device realized exactly as in the IEEE-39 runs (tx4_emt kernel). dt = 50 us, 1 ms output.

Pass rules (frozen):
  equilibrium |dP|, |dQ| <= 1e-4 pu (t = 0.9 s);
  every state: max|x_EMT - x_ph| <= 0.02 * max|x_ph - x_ph(0)| over [1 s, T];
  SG ringdown (matrix pencil on omega): |d alpha| <= 0.005 s^-1, |d f| <= 0.005 Hz;
  GFL: PLL steady state |theta - angle(V_T)| <= 1e-3 rad at T; g = 0: |q_f - q_ref| <= 1e-3 at T (case b).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _emt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402
from emt_estimator import matrix_pencil  # noqa: E402

DT = 50e-6
DS = 20
SG_NAMES = ("delta", "omega", "eq1", "ed1", "efd", "pss_w", "pss_l")


def mini_network(vt, vinf):
    net = {"buses": [{"idx": 1}, {"idx": 2}],
           "lines": [{"bus1": 1, "bus2": 2, "r": 0.0, "x": 0.05, "b": 0.0, "g": 0.0, "tap": 1.0, "phi": 0.0,
                      "u": 1.0, "trans": 0.0}], "shunts": []}
    return K.build_network(net, {1: vt, 2: vinf}, DT, net_damping=0.0)


def common(op):
    vt, vinf = complex(*op["V_T"]), complex(*op["V_INF"])
    s = complex(*op["S_dev"])
    net = mini_network(vt, vinf)
    g = net["g0"].copy()
    ys = 1.0 / complex(*op["z_s"])
    K.stamp(g, 2, 1, ys)
    # amendment A1: solidly grounded infinite bus = zero-sequence-only conductance (g_z/3) J at INF
    gz = abs(ys)
    nodes = [1, 3, 5]
    for a in nodes:
        for b in nodes:
            g[a, b] += gz / 3.0
    return vt, s, net, g, ys


def run_kernel(net, g, sg_par, sg_y, sg_x, gf_par, gf_x, ys, events, tlen):
    e_i, e_c = np.zeros(0, np.int64), np.zeros(0, np.complex128)
    return K.simulate(np.linalg.inv(g), net["coe0"], net["vsol0"], net["brch_ihis"], net["node_ihis"], 2,
                      sg_par, sg_y, sg_x, gf_par, gf_x, e_i, e_c, e_c, 1e-3,
                      np.array([1], np.int64), np.array([ys], np.complex128), np.array([1.0 + 0j]),
                      events, DT, int(round(tlen / DT)), DS, np.arange(2, dtype=np.int64))


def compare(name, t_e, x_e, t_p, x_p, labels, t_from=1.0):
    assert np.allclose(t_e, t_p[: len(t_e)], atol=1e-9)
    if len(t_e) < len(t_p):  # EMT run stopped non-finite: every state check fails (no comparison exists)
        return [{"test": name, "state": lab, "max_abs_err": float("inf"), "max_excursion": float("nan"),
                 "ratio": float("inf"), "pass": False, "emt_nonfinite_at_s": float(t_e[-1])} for lab in labels]
    m = t_e >= t_from - 1e-12
    rows = []
    for j, lab in enumerate(labels):
        err = float(np.abs(x_e[m, j] - x_p[m, j]).max())
        exc = float(np.abs(x_p[m, j] - x_p[0, j]).max())
        rows.append({"test": name, "state": lab, "max_abs_err": err, "max_excursion": exc,
                     "ratio": err / exc if exc > 0 else float("inf"), "pass": bool(err <= 0.02 * exc)})
    return rows


def main() -> int:
    data = C.Data(_emt.RESEARCH)
    out2, out3 = _emt.RESULTS / "EMT02", _emt.RESULTS / "EMT03"
    # ---------------- EMT02 ----------------
    op = json.loads((out2 / "op.json").read_text())
    ref = np.load(out2 / "phasor_ref.npz")
    vt, s, net, g, ys = common(op)
    p = C.machine_params(data, 30, (0.03625, 1.425, 1.5, 1.0), {}, {})
    w = op["w"]
    ini = C.sg_init(p, w, vt, s)
    row = np.zeros(len(K.SG_COLS))
    vals = {"bus": 0, "w": w, "ra": p["ra"], "xd": p["xd"], "xq": p["xq"], "x1": p["xd1"], "td10": p["td10"],
            "tq10": p["tq10"], "m": p["m"], "d": p["d"], "ka": p["ka"], "ta": p["ta"], "ks": p["pss_gain"],
            "t4": p["pss_lag"], "t5": p["pss_washout"], "t6": p["pss_wash_lag"], "pm": ini["pm"], "vref": ini["vref"],
            "flux_blend": 1.0, "avr_blend": 1.0, "has_gov": 0.0, "gov_r": 1.0, "gov_t1": 1.0, "gov_t2": 1.0,
            "gov_t3": 1.0, "gov_dt": 0.0, "gov_pref": ini["pm"]}
    for n, v in vals.items():
        row[K.SG[n]] = v
    ysg = w / complex(p["ra"], p["xd1"])
    K.stamp(g, 2, 0, ysg)
    x0 = np.array([[ini["delta"], 1.0, ini["eq1"], ini["ed1"], ini["efd"], ini["pm"], 0.0, ini["pm"], ini["pm"]]])
    events = np.array([[K.EV_PM, 0, 1.0, 1.2, 0.02]])
    t, sg, gf, v, fin = run_kernel(net, g, row[None, :], np.array([ysg]), x0, np.zeros((0, 20)), np.zeros((0, 11)),
                                   ys, events, 10.0)
    rows = compare("EMT02_SG_pm_pulse", t, sg[:, 0, :7], ref["t"], ref["x"], SG_NAMES)
    i09 = int(np.argmin(np.abs(t - 0.9)))
    eq2 = {"dP": float(sg[i09, 0, 10] - s.real), "dQ": float(sg[i09, 0, 11] - s.imag)}
    me = matrix_pencil(t, (sg[:, 0, 1] - 1.0)[None, :], 1e-3)
    mp = matrix_pencil(ref["t"], (ref["x"][:, 1] - 1.0)[None, :], 1e-3)
    ring = {"alpha_emt": me["alpha"], "f_emt": me["freq"], "alpha_ph": mp["alpha"], "f_ph": mp["freq"],
            "d_alpha": me["alpha"] - mp["alpha"], "d_f": me["freq"] - mp["freq"]}
    pass2 = (all(r["pass"] for r in rows) and abs(eq2["dP"]) <= 1e-4 and abs(eq2["dQ"]) <= 1e-4
             and abs(ring["d_alpha"]) <= 0.005 and abs(ring["d_f"]) <= 0.005 and fin)
    df2 = pd.DataFrame(rows)
    df2.to_csv(out2 / "EMT02_trajectory_comparison.csv", index=False)
    np.savez_compressed(out2 / "emt_traj.npz", t=t, x=sg[:, 0, :7], s=sg[:, 0, 10] + 1j * sg[:, 0, 11])
    s2 = {"equilibrium": eq2, "ringdown": ring, "finite": bool(fin), "GATE_G3": "PASS" if pass2 else "FAIL"}
    (out2 / "EMT02_summary.json").write_text(json.dumps(s2, indent=1), encoding="utf-8")
    print("EMT02", json.dumps(s2, indent=1))
    print(df2.to_string(index=False))

    # ---------------- EMT03 ----------------
    op3 = json.loads((out3 / "op.json").read_text())
    ref3 = np.load(out3 / "phasor_ref.npz")
    rows3, extra, trajs = [], [], {}
    for gg in (0.0, 0.03625, 0.25):
        for case in ("a", "b", "c"):
            key = f"g{gg:g}_{case}"
            vt, s, net, g, ys = common(op3)
            par = dict(C.GFL_DEFAULTS)
            x, p_ref, q_ref, v_ref = C.gfl_init(par, op3[f"meta_{key}"]["w"], vt, s)
            row = np.zeros(len(K.GFL_COLS))
            vals = dict(par, bus=0, w=op3[f"meta_{key}"]["w"], g=gg, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref)
            for n, v in vals.items():
                row[K.GF[n]] = v
            ev = {"a": [K.EV_PREF, 0, 1.0, 1e9, 0.01], "b": [K.EV_EMAG, 0, 1.0, 1e9, -0.01],
                  "c": [K.EV_EPHASE, 0, 1.0, 1e9, 0.02]}[case]
            t, sg, gf, v, fin = run_kernel(net, g, np.zeros((0, 27)), np.zeros(0, np.complex128), np.zeros((0, 9)),
                                           row[None, :], np.array([x]), ys, np.array([ev]), 5.0)
            rows3 += compare(f"EMT03_{key}", t, gf[:, 0, :11], ref3[f"{key}_t"], ref3[f"{key}_x"], K.GFL_STATES)
            done = bool(fin) and t[-1] >= 5.0 - 1e-9
            i09 = int(np.argmin(np.abs(t - 0.9)))
            vT = v[-1, 0, 0] + 1j * v[-1, 0, 1]
            pll_err = float(abs(np.angle(np.exp(1j * (gf[-1, 0, 0] - np.angle(vT)))))) if done else float("inf")
            qf_err = (float(abs(gf[-1, 0, 3] - q_ref)) if done else float("inf")) if (gg == 0.0 and case == "b") else float("nan")
            eq_ok = t[-1] >= 0.9 and np.isfinite(gf[i09, 0, 11:13]).all()
            extra.append({"test": key, "dP_eq": float(gf[i09, 0, 11] - s.real) if eq_ok else float("inf"),
                          "dQ_eq": float(gf[i09, 0, 12] - s.imag) if eq_ok else float("inf"),
                          "pll_steady_err_rad": pll_err, "g0_qf_minus_qref": qf_err, "finite": done,
                          "emt_last_finite_t_s": float(t[-1])})
            trajs[key] = gf[:, 0, :11]
    df3 = pd.DataFrame(rows3)
    ex3 = pd.DataFrame(extra)
    df3.to_csv(out3 / "EMT03_trajectory_comparison.csv", index=False)
    ex3.to_csv(out3 / "EMT03_checks.csv", index=False)
    np.savez_compressed(out3 / "emt_traj.npz", **trajs)
    pass3 = (bool(df3["pass"].all()) and bool((ex3.dP_eq.abs() <= 1e-4).all()) and bool((ex3.dQ_eq.abs() <= 1e-4).all())
             and bool((ex3.pll_steady_err_rad <= 1e-3).all()) and bool((ex3.g0_qf_minus_qref.dropna() <= 1e-3).all())
             and bool(ex3.finite.all()))
    s3 = {"n_state_checks": len(df3), "n_pass": int(df3["pass"].sum()), "max_ratio": float(df3.ratio.max()),
          "n_cases_nonfinite": int((~ex3.finite).sum()), "earliest_nonfinite_t_s": float(ex3.emt_last_finite_t_s.min()),
          "worst": df3.sort_values("ratio").tail(5).to_dict(orient="records"), "GATE_G4": "PASS" if pass3 else "FAIL"}
    (out3 / "EMT03_summary.json").write_text(json.dumps(s3, indent=1), encoding="utf-8")
    print("EMT03", json.dumps(s3, indent=1))
    print(ex3.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
