# ruff: noqa: E501  -- table rows kept on one line
"""EMT01 - network and operating-point equivalence (prereg section 2).

1. Official ParaEMT IEEE-39 (cases/pfd_39_1_1.json, upstream d79d735a) vs the frozen TX4 IEEE-39
   (configs/ias2026/ieee39_network.json + canonical power flow): automated diff.
2. TX4-specific ParaEMT network: 60-Hz positive-sequence Ybus assembled from the element values
   actually stored in ParaEMT's branch table (coe0: R, X or L, C, tap) vs the canonical Ybus
   (gate max|dY|/max|Y| <= 1e-8); companion coefficients checked to be pure trapezoidal;
   trapezoidal 60-Hz warping reported separately (not gated).
3. Equilibrium: base portfolio at P4, 2 s no-event run at 50 us; positive-sequence phasors averaged
   over the last 0.5 s vs the canonical power flow (gates 1e-4 pu, 1e-3 rad); device P/Q residuals.
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
from tx4_case import Data  # noqa: E402

OUT = _emt.RESULTS / "EMT01"
W0 = _emt.W0


def official_vs_tx4(data: Data) -> tuple[pd.DataFrame, pd.DataFrame]:
    off = json.loads((_emt.UP / "cases" / "pfd_39_1_1.json").read_text())
    net = data.net
    rows = []

    def add(item, key, official, tx4, note=""):
        same = None
        try:
            same = bool(
                np.isclose(complex(official), complex(tx4), rtol=1e-6, atol=1e-9)
            )
        except (TypeError, ValueError):
            same = str(official) == str(tx4)
        rows.append(
            {
                "item": item,
                "key": key,
                "official_paremt": official,
                "frozen_tx4": tx4,
                "identical": same,
                "note": note,
            }
        )

    add("system", "base_MVA", off["basemva"], 100.0)
    add("system", "frequency_hz", off["ws"] / (2 * np.pi), 60.0)
    add("buses", "count", len(off["bus_num"]), len(net["buses"]))
    add(
        "buses",
        "numbering",
        str(sorted(off["bus_num"])),
        str(sorted(int(b["idx"]) for b in net["buses"])),
    )
    add(
        "buses",
        "base_kV",
        str(sorted(set(off["bus_basekV"]))),
        "per unit (Vn 345/22 kV in source; pu model)",
        "ParaEMT case is all per unit (basekV = 1)",
    )
    # lines
    tx4_lines = {}
    for ln in net["lines"]:
        f, t = int(ln["bus1"]), int(ln["bus2"])
        is_x = float(ln["tap"]) != 1.0 or float(ln["trans"]) == 1.0
        tx4_lines[(min(f, t), max(f, t), is_x)] = (
            f,
            t,
            complex(ln["r"], ln["x"]),
            float(ln["b"]),
            float(ln["tap"]),
        )
    off_items = {}
    for f, t, rx, b in zip(
        off["line_from"], off["line_to"], off["line_RX"], off["line_chg"], strict=True
    ):
        off_items[(min(f, t), max(f, t), False)] = (f, t, complex(rx), float(b), 1.0)
    for f, t, rx, k in zip(
        off["xfmr_from"], off["xfmr_to"], off["xfmr_RX"], off["xfmr_k"], strict=True
    ):
        off_items[(min(f, t), max(f, t), True)] = (
            f,
            t,
            complex(rx),
            0.0,
            float(np.real(complex(k))),
        )
    for key in sorted(set(tx4_lines) | set(off_items)):
        name = f"{'xfmr' if key[2] else 'line'} {key[0]}-{key[1]}"
        o, x = off_items.get(key), tx4_lines.get(key)
        if o is None or x is None:
            add(name, "present", o is not None, x is not None, "endpoint set differs")
            continue
        add(name, "R", o[2].real, x[2].real)
        add(name, "X", o[2].imag, x[2].imag)
        add(name, "B", o[3], x[3])
        if key[2]:
            add(
                name,
                "tap (tx4: on bus1; official: xfmr_k, PF closes with ratio on the to-bus)",
                o[4],
                x[4],
                f"official from->to {o[0]}->{o[1]}; tx4 bus1->bus2 {x[0]}->{x[1]}",
            )
    # shunts
    for b, gb in zip(off["shnt_bus"], off["shnt_gb"], strict=True):
        tx = [s for s in net["shunts"] if int(s["bus"]) == b]
        add(f"shunt {b}", "B_pu", complex(gb).imag / 100.0, tx[0]["b"] if tx else None)
    # loads
    tx_load = {int(b): complex(*v) for b, v in data.op["loads"].items()}
    off_load = {
        int(b): complex(p, q) / 100.0
        for b, p, q in zip(
            off["load_bus"], off["load_MW"], off["load_Mvar"], strict=True
        )
    }
    for b in sorted(set(tx_load) | set(off_load)):
        add(f"load {b}", "P_pu", off_load.get(b, 0j).real, tx_load.get(b, 0j).real)
        add(f"load {b}", "Q_pu", off_load.get(b, 0j).imag, tx_load.get(b, 0j).imag)
    net_df = pd.DataFrame(rows)
    # operating point
    op = []
    gen_tx = {int(b): complex(*v) for b, v in data.op["base"]["generation"].items()}
    for b, p, q, mva in zip(
        off["gen_bus"], off["gen_MW"], off["gen_Mvar"], off["gen_MVA_base"], strict=True
    ):
        op.append(
            {
                "item": f"gen {b}",
                "key": "P_pu",
                "official_paremt": p / 100.0,
                "frozen_tx4": gen_tx[b].real,
            }
        )
        op.append(
            {
                "item": f"gen {b}",
                "key": "Q_pu",
                "official_paremt": q / 100.0,
                "frozen_tx4": gen_tx[b].imag,
            }
        )
        op.append(
            {
                "item": f"gen {b}",
                "key": "rating_MVA",
                "official_paremt": mva,
                "frozen_tx4": data.mach[b]["Sn"],
            }
        )
    vt = {int(b): complex(*v) for b, v in data.op["base"]["voltages"].items()}
    for b, vm, va in zip(off["bus_num"], off["bus_Vm"], off["bus_Va"], strict=True):
        op.append(
            {
                "item": f"bus {b}",
                "key": "Vm_pu",
                "official_paremt": vm,
                "frozen_tx4": abs(vt[b]),
            }
        )
        op.append(
            {
                "item": f"bus {b}",
                "key": "Va_rad",
                "official_paremt": va,
                "frozen_tx4": float(np.angle(vt[b])),
            }
        )
    op_df = pd.DataFrame(op)
    op_df["difference"] = op_df.official_paremt - op_df.frozen_tx4
    return net_df, op_df


def emt_ybus(built: dict, nbus: int, dt: float):
    coe0 = built["net"]["coe0"]
    nl = len(built["net"]["lines"][0])
    nx = len(built["net"]["xfmrs"][0])
    ns = len(built["net"]["shunts"][0])
    yc = np.zeros((nbus, nbus), complex)
    yd = np.zeros((nbus, nbus), complex)
    sd = 1j * (2.0 / dt) * np.tan(W0 * dt / 2.0)  # trapezoidal s at 60 Hz
    checks = []
    for i in range(nl):
        r = coe0[9 * i]  # phase A series row: [F, T, Req, icf, Gv1, R, X, 0, I]
        f, t, R, X = int(r[0].real), int(r[1].real), r[5].real, r[6].real
        L = X / W0
        checks.append(abs(r[2].real - (R + 2 * L / dt)) / (R + 2 * L / dt))
        y = 1 / complex(R, X)
        ydisc = 1 / (R + sd * L)
        for m, yy in ((yc, y), (yd, ydisc)):
            m[f, f] += yy
            m[t, t] += yy
            m[f, t] -= yy
            m[t, f] -= yy
        for rowi in (9 * i + 3, 9 * i + 6):
            rc = coe0[rowi]
            node, C = int(rc[0].real), rc[7].real
            if C > 0:
                checks.append(abs(rc[2].real - dt / (2 * C)) / (dt / (2 * C)))
                yc[node, node] += 1j * W0 * C
                yd[node, node] += sd * C
    base = 9 * nl
    for i in range(nx):
        r = coe0[base + 3 * i]
        f, t, R, L, tap = (
            int(r[0].real),
            int(r[1].real),
            r[5].real,
            r[6].real,
            r[9].real,
        )
        checks.append(abs(r[2].real - (R + 2 * L / dt)) / (R + 2 * L / dt))
        for m, yy in ((yc, 1 / complex(R, W0 * L)), (yd, 1 / (R + sd * L))):
            m[f, f] += yy / tap**2
            m[t, t] += yy
            m[f, t] -= yy / tap
            m[t, f] -= yy / tap
    base2 = base + 3 * nx
    for i in range(ns):
        r = coe0[base2 + 3 * i]
        node, C = int(r[0].real), r[7].real
        checks.append(abs(r[2].real - dt / (2 * C)) / (dt / (2 * C)))
        yc[node, node] += 1j * W0 * C
        yd[node, node] += sd * C
    return yc, yd, float(max(checks))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    data = Data(_emt.RESEARCH)
    net_df, op_df = official_vs_tx4(data)
    net_df.to_csv(OUT / "network_diff.csv", index=False)
    op_df.to_csv(OUT / "operating_point_diff.csv", index=False)
    import tx4_case

    case = _emt.CASE["EMT05_P4_BASE"]
    order = data.order
    y_can = np.array(data.op["ybus_real"]) + 1j * np.array(data.op["ybus_imag"])
    ybus = {}
    for dt in (25e-6, 50e-6, 100e-6):
        built = tx4_case.build_case(data, case, dt, amplitude=0.0)
        yc, yd, trap_check = emt_ybus(built, len(order), dt)
        ybus[dt] = {
            "rel_residual_continuous": float(
                np.abs(yc - y_can).max() / np.abs(y_can).max()
            ),
            "rel_warp_discrete_vs_continuous": float(
                np.abs(yd - yc).max() / np.abs(yc).max()
            ),
            "companion_is_pure_trapezoidal_max_rel_dev": trap_check,
        }
    # equilibrium: 2 s no event at 50 us
    out = _emt.run_spec(case, dt=50e-6, tlen=2.0, amplitude=0.0, save_raw=False)
    t, v = out["t"], out["v"]
    vend = v[t >= 1.5].mean(axis=0)
    v0 = np.array([complex(*data.op["base"]["voltages"][str(b)]) for b in order])
    i39 = order.index(39)
    dv = np.abs(np.abs(vend) - np.abs(v0))
    da = np.abs(
        np.angle(vend) - np.angle(vend[i39]) - (np.angle(v0) - np.angle(v0[i39]))
    )
    sg = out["sg"]
    gen = {int(b): complex(*x) for b, x in data.op["base"]["generation"].items()}
    pq = [
        {
            "device": f"SG {b}",
            "P_emt": float(sg[t >= 1.5, j, 10].mean()),
            "Q_emt": float(sg[t >= 1.5, j, 11].mean()),
            "P_pf": gen[b].real,
            "Q_pf": gen[b].imag,
        }
        for j, b in enumerate(out["labels_sg"])
    ]
    pq_df = pd.DataFrame(pq)
    pq_df["dP"] = pq_df.P_emt - pq_df.P_pf
    pq_df["dQ"] = pq_df.Q_emt - pq_df.Q_pf
    pq_df.to_csv(OUT / "equilibrium_device_pq.csv", index=False)
    eq = pd.DataFrame(
        {
            "bus": order,
            "Vm_pf": np.abs(v0),
            "Vm_emt": np.abs(vend),
            "dVm": np.abs(vend) - np.abs(v0),
            "dAngle_rel39": np.angle(vend)
            - np.angle(vend[i39])
            - (np.angle(v0) - np.angle(v0[i39])),
        }
    )
    eq.to_csv(OUT / "equilibrium_bus_voltages.csv", index=False)
    summary = {
        "ybus": {f"{int(k * 1e6)}us": v for k, v in ybus.items()},
        "ybus_gate_rel_1e-8": all(
            v["rel_residual_continuous"] <= 1e-8 for v in ybus.values()
        ),
        "equilibrium": {
            "max_dV_pu": float(dv.max()),
            "max_dAngle_rad": float(da.max()),
            "max_dP_pu": float(pq_df.dP.abs().max()),
            "max_dQ_pu": float(pq_df.dQ.abs().max()),
            "max_speed_deviation_pu": float(np.abs(sg[:, :, 1] - 1).max()),
        },
        "equilibrium_gate": bool(dv.max() <= 1e-4 and da.max() <= 1e-3),
        "official_vs_tx4": {
            "network_rows": len(net_df),
            "network_rows_identical": int(net_df.identical.sum()),
            "max_abs_op_difference": float(op_df.difference.abs().max()),
        },
        "manifest": out["manifest"],
    }
    summary["GATE_G1"] = (
        "PASS"
        if summary["ybus_gate_rel_1e-8"] and summary["equilibrium_gate"]
        else "FAIL"
    )
    (OUT / "EMT01_summary.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "manifest"}, indent=1))
    pd.set_option("display.width", 220)
    print(net_df[~net_df.identical.astype(bool)].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
