# ruff: noqa: E501  -- diagnostic tables kept on one line
"""Diagnostics of the preregistered unit-test failures (G3, G4). DIAGNOSTIC ONLY, not gates.

The preregistered verdicts stay as recorded in results/EMT02/EMT02_summary.json (G3) and
results/EMT03 (G4). This script localizes their cause without changing any threshold:
  D1  quasi-static-line variant: the T--INF line realized as a synchronous-frame algebraic Norton
      (no electromagnetic dynamics; a zero-sequence-only coupling (g_z/3)J keeps T referenced), so the
      EMT network is algebraic: isolates the device transcription (SG and all nine GFL cases).
  D2  SG with the preregistered EMT line at dt = 25 us: numerical vs physical attribution.
  D3  SG, preregistered EMT line (50 us) vs a phasor reference WITH the line's dI/dt term.
  D4  GFL, preregistered EMT line: size and growth of the alternating (Nyquist) component of |V_T|.
Writes results/EMT02/EMT02_diagnostics.json, results/EMT03/EMT03_diagnostics.json (+ csv).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _emt  # noqa: E402
import EMT02_03_unit_tests as U  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402

OUT2, OUT3 = _emt.RESULTS / "EMT02", _emt.RESULTS / "EMT03"


def quasi_static_g(op):
    """G with the EMT line replaced by an algebraic synchronous-frame line (+ zero-seq coupling)."""

    vt, s, net, g, ys = U.common(op)
    g = g - net["g0"]
    yl = 1.0 / complex(*op["z_line"])
    for a, b, sg in ((0, 0, 1), (1, 1, 1), (0, 1, -1), (1, 0, -1)):
        blk = K.norton_block(sg * yl)
        for k in range(3):
            for m in range(3):
                g[a + 2 * k, b + 2 * m] += blk[k, m]
    gz = abs(yl)
    for a, b, sg in ((0, 0, 1), (1, 1, 1), (0, 1, -1), (1, 0, -1)):
        for k in range(3):
            for m in range(3):
                g[a + 2 * k, b + 2 * m] += sg * gz / 3.0
    net0 = dict(net)
    net0["coe0"] = net["coe0"][:0]  # no EMT branches
    net0["brch_ihis"] = net["brch_ihis"][:0]
    net0["node_ihis"] = np.zeros_like(net["node_ihis"])
    return vt, s, net0, g, ys


def sg_setup(op, gmode, dt):
    data = C.Data(_emt.RESEARCH)
    U.DT = dt
    U.DS = int(round(1e-3 / dt))
    vt, s, net, g, ys = quasi_static_g(op) if gmode == "qs" else U.common(op)
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
    t, sg, gf, v, fin = U.run_kernel(net, g, row[None, :], np.array([ysg]), x0, np.zeros((0, 20)), np.zeros((0, 11)),
                                     ys, np.array([[K.EV_PM, 0, 1.0, 1.2, 0.02]]), 10.0)
    return t, sg[:, 0, :7], fin


def summarize(rows):
    df = pd.DataFrame(rows)
    return {"max_ratio": float(df.ratio.replace(np.inf, np.nan).max()), "all_pass_2pct": bool(df["pass"].all()),
            "per_state": df[["state", "ratio"]].to_dict(orient="records")}


def main() -> int:
    op = json.loads((OUT2 / "op.json").read_text())
    ref = np.load(OUT2 / "phasor_ref.npz")
    dyn = np.load(OUT2 / "dynphasor_ref.npz")
    d2 = {}
    t, x, fin = sg_setup(op, "qs", 50e-6)
    d2["D1_quasi_static_line_50us_vs_quasi_static_phasor"] = summarize(U.compare("D1", t, x, ref["t"], ref["x"], U.SG_NAMES)) | {"finite": fin}
    t, x, fin = sg_setup(op, "emt", 25e-6)
    d2["D2_emt_line_25us_vs_quasi_static_phasor"] = summarize(U.compare("D2", t, x, ref["t"], ref["x"], U.SG_NAMES)) | {"finite": fin}
    t, x, fin = sg_setup(op, "emt", 50e-6)
    d2["D3_emt_line_50us_vs_dynamic_phasor"] = summarize(U.compare("D3", t, x, dyn["t"], dyn["x"], U.SG_NAMES)) | {"finite": fin}
    d2["expected_dynamic_phasor_effect_f_over_f0"] = 1.177 / 60.0
    (OUT2 / "EMT02_diagnostics.json").write_text(json.dumps(d2, indent=1), encoding="utf-8")
    print("EMT02", json.dumps({k: (v["max_ratio"] if isinstance(v, dict) else v) for k, v in d2.items()}, indent=1))

    # ---------------- GFL ----------------
    op3 = json.loads((OUT3 / "op.json").read_text())
    ref3 = np.load(OUT3 / "phasor_ref.npz")
    U.DT, U.DS = 50e-6, 20
    rows, d4 = [], []
    for gg in (0.0, 0.03625, 0.25):
        for case in ("a", "b", "c"):
            key = f"g{gg:g}_{case}"
            vt, s, net, g, ys = quasi_static_g(op3)
            par = dict(C.GFL_DEFAULTS)
            w = op3[f"meta_{key}"]["w"]
            x, p_ref, q_ref, v_ref = C.gfl_init(par, w, vt, s)
            row = np.zeros(len(K.GFL_COLS))
            vals = dict(par, bus=0, w=w, g=gg, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref)
            for n, v in vals.items():
                row[K.GF[n]] = v
            ev = {"a": [K.EV_PREF, 0, 1.0, 1e9, 0.01], "b": [K.EV_EMAG, 0, 1.0, 1e9, -0.01],
                  "c": [K.EV_EPHASE, 0, 1.0, 1e9, 0.02]}[case]
            t, sg, gf, vv, fin = U.run_kernel(net, g, np.zeros((0, 27)), np.zeros(0, np.complex128), np.zeros((0, 9)),
                                              row[None, :], np.array([x]), ys, np.array([ev]), 5.0)
            rr = U.compare(f"D1_{key}", t, gf[:, 0, :11], ref3[f"{key}_t"], ref3[f"{key}_x"], K.GFL_STATES)
            for r in rr:
                r["finite"] = bool(fin)
            rows += rr
    # D4 chatter in the preregistered configuration (g = 0.25, no event), every step for 0.5 s
    vt, s, net, g, ys = U.common(op3)
    par = dict(C.GFL_DEFAULTS)
    w = op3["meta_g0.25_a"]["w"]
    x, p_ref, q_ref, v_ref = C.gfl_init(par, w, vt, s)
    row = np.zeros(len(K.GFL_COLS))
    for n, v in dict(par, bus=0, w=w, g=0.25, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref).items():
        row[K.GF[n]] = v
    U.DS = 1
    t, sg, gf, vv, fin = U.run_kernel(net, g, np.zeros((0, 27)), np.zeros(0, np.complex128), np.zeros((0, 9)),
                                      row[None, :], np.array([x]), ys, np.zeros((0, 5)), 0.7)
    vm = np.abs(vv[:, 0, 0] + 1j * vv[:, 0, 1])
    alt = np.abs(vm[1:] - vm[:-1])
    idx = [int(round(tt / 50e-6)) for tt in (0.01, 0.1, 0.3, 0.5, min(0.69, t[-2]))]
    d4 = {"alternating_amplitude_pu": {f"t={tt}": float(alt[min(i, alt.size - 1)]) for tt, i in zip((0.01, 0.1, 0.3, 0.5, 0.69), idx, strict=True)},
          "finite_0_to_0p7s": bool(fin), "note": "step-to-step alternation of |V_T| (Nyquist mode) at the GFL terminal"}
    df = pd.DataFrame(rows)
    df.to_csv(OUT3 / "EMT03_D1_quasi_static_line.csv", index=False)
    d3 = {"D1_quasi_static_line_all_9_cases": {"max_ratio": float(df.ratio.replace(np.inf, np.nan).max()),
                                               "all_pass_2pct": bool(df["pass"].all()), "all_finite": bool(df.finite.all())},
          "D4_chatter_preregistered_line": d4}
    (OUT3 / "EMT03_diagnostics.json").write_text(json.dumps(d3, indent=1), encoding="utf-8")
    print("EMT03", json.dumps(d3, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
