# ruff: noqa: E501  -- parameter tables and kernel argument lists kept on one line
"""TX4 EMT case builder (overlay; copied into external/ParaEMT_tx4).

Turns a preregistered case (results/EMT_PRED/emt_cases.json) into the arrays of tx4_emt.simulate:
device parameters from the frozen configs, device states from the phasor `initialize` formulas at
the canonical operating point, the network through ParaEMT, Norton stamps and the disturbance.
Semantics of the policy / draw / condenser / governor options follow the frozen phasor code
(ieee39_case.build_dae, F8 condenser services, governed.govern, PCV05 case_kwargs).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import tx4_emt as K

MACHINE_GROUPS = {
    "M": ("m",),
    "XP": ("xd1", "xq1"),
    "KA": ("ka",),
    "TE": ("ta",),
    "PSS_K": ("pss_gain",),
    "PSS_T": ("pss_washout", "pss_wash_lag", "pss_lag"),
}
CONV_GROUPS = {
    "PLL": ("kp_pll", "ki_pll"),
    "OUTER": ("kp_p", "ki_p", "kp_q", "ki_q"),
    "CURRENT": ("kp_i", "ki_i"),
    "TAU_P": ("tau_p",),
    "XF": ("xf",),
}
GFL_DEFAULTS = {
    "kp_pll": 53.0,
    "ki_pll": 1400.0,
    "tau_p": 0.03,
    "kp_p": 0.20,
    "ki_p": 8.0,
    "kp_q": 0.20,
    "ki_q": 8.0,
    "kp_i": 0.25,
    "ki_i": 6.0,
    "xf": 0.15,
    "rf": 0.01,
    "kp_v": 2.0,
    "ki_v": 20.0,
}
LEAK = 0.05
CONDENSER_R0010000 = {
    "inertia_scale": 0.04,
    "damping": 2.0,
    "flux_blend": 0.0,
    "avr_blend": 0.0,
    "pss_scale": 0.0,
    "q_share": 0.0,
}
TAU_M = 1e-3
PULSE = {"bus": 20, "fraction": 0.02, "t0": 1.0, "t1": 1.2}


def sha(obj) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, default=str).encode()
    ).hexdigest()[:16]


class Data:
    def __init__(self, research: Path):
        self.net = json.loads(
            (research / "configs/ias2026/ieee39_network.json").read_text()
        )
        self.gov = {
            int(g["bus"]): g
            for g in json.loads(
                (
                    research / "configs/ias2026/ieee39_governed_documented_v1.json"
                ).read_text()
            )["governors"]
        }
        self.op = json.loads(
            (research / "results/EMT_PRED/tx4_operating_points.json").read_text()
        )
        order = [int(r["bus"]) for r in self.net["machines"]]
        self.mach = {
            b: dict(
                self.net["machines"][i],
                **{"avr": self.net["avr"][i], "pss": self.net["pss"][i]},
            )
            for i, b in enumerate(order)
        }
        self.order = [int(b["idx"]) for b in self.net["buses"]]
        self.idx = {b: i for i, b in enumerate(self.order)}


def point(data: Data, variant: dict):
    if variant.get("kind") == "line":
        p = data.op["lines_x1p5"][str(variant["line_index"])]
    else:
        p = data.op["base"]
    volt = {int(b): complex(*v) for b, v in p["voltages"].items()}
    gen = {int(b): complex(*v) for b, v in p["generation"].items()}
    return volt, gen


def machine_params(
    data: Data, bus: int, theta, draw_fleet: dict, draw_unit: dict
) -> dict:
    raw = data.mach[bus]
    g, k, t, h = theta
    assert h == 1.0
    assert raw["xd1"] == raw["xq1"], bus
    p = {
        "ra": raw["ra"],
        "xd": raw["xd"],
        "xq": raw["xq"],
        "xd1": raw["xd1"],
        "xq1": raw["xq1"],
        "td10": raw["Td10"],
        "tq10": raw["Tq10"],
        "m": raw["M"],
        "d": raw["D"],
        "ka": raw["avr"]["KA"] * k,
        "ta": raw["avr"]["TE"] * t,
        "pss_gain": raw["pss"]["KS"],
        "pss_washout": raw["pss"]["T5"],
        "pss_wash_lag": raw["pss"]["T6"],
        "pss_lag": raw["pss"]["T4"],
    }
    for grp, f in draw_fleet.items():
        for fld in MACHINE_GROUPS[grp]:
            p[fld] *= f
    for grp, f in draw_unit.get(bus, {}).items():
        for fld in MACHINE_GROUPS[grp]:
            p[fld] *= f
    assert abs(p["xd1"] - p["xq1"]) < 1e-15
    return p


def sg_init(p: dict, w: float, v: complex, s: complex):
    """SynchronousMachine.initialize (phasor formulas)."""

    cur = np.conj(s / w / v)
    internal = v + complex(p["ra"], p["xq"]) * cur
    delta = float(np.angle(internal))

    def dq(z):
        return (
            z.real * np.sin(delta) - z.imag * np.cos(delta),
            z.real * np.cos(delta) + z.imag * np.sin(delta),
        )

    vd, vq = dq(v)
    id_, iq = dq(cur)
    eq1 = vq + p["ra"] * iq + p["xd1"] * id_
    ed1 = vd + p["ra"] * id_ - p["xq1"] * iq
    efd = eq1 + (p["xd"] - p["xd1"]) * id_
    power = vd * id_ + vq * iq + p["ra"] * (id_ * id_ + iq * iq)
    vref = abs(v) + efd / p["ka"]
    return {
        "delta": delta,
        "eq1": eq1,
        "ed1": ed1,
        "efd": efd,
        "pm": power,
        "vref": vref,
    }


def gfl_init(par: dict, w: float, v: complex, s: complex):
    cur = np.conj(s / w / v)
    theta = float(np.angle(v))
    dqc = cur * np.exp(-1j * theta)
    i_d, i_q = float(dqc.real), float(dqc.imag)
    v_d = float(abs(v))
    p_ref, q_ref = v_d * i_d, -v_d * i_q
    x = [
        theta,
        0.0,
        p_ref,
        q_ref,
        i_d,
        -i_q,
        i_d,
        i_q,
        par["rf"] * i_d,
        par["rf"] * i_q,
        q_ref,
    ]
    return x, p_ref, q_ref, v_d


def build_case(
    data: Data,
    case: dict,
    dt: float,
    amplitude: float = PULSE["fraction"],
    tau_m: float = TAU_M,
    net_damping: float = 0.0,
):
    members = tuple(case["members"])
    theta = tuple(case["theta"])
    variant = case.get("variant", {"kind": "none"})
    volt, gen = point(data, variant)
    fleet, unit, conv = {}, {}, {}
    if variant.get("kind") == "draw":
        fac = variant["factors"]
        fleet = dict(fac["fleet"])
        unit = {int(b): v for b, v in fac["unit"].items()}
        conv = {int(b): v for b, v in fac["conv"].items()}
    line_scale = (
        (variant["line_index"], variant["gamma"])
        if variant.get("kind") == "line"
        else None
    )
    net = K.build_network(
        data.net, volt, dt, line_scale=line_scale, net_damping=net_damping
    )
    nbus = net["nbus"]
    gmat = net["g0"].copy()
    sg_rows, sg_y, sg_x, gf_rows, gf_x, labels_sg, labels_gf = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    cond = variant.get("kind") == "condenser"
    governed = variant.get("kind") == "governed"
    for bus in sorted(data.mach):
        rating = data.mach[bus]["Sn"] / 100.0
        v = volt[bus]
        s = gen[bus]
        replaced = bus in members
        sg_specs = []
        if not replaced:
            sg_specs.append((rating, s, None))
        elif cond:
            sg_specs.append((variant["fraction"] * rating, 0j, CONDENSER_R0010000))
        for w, share, svc in sg_specs:
            p = machine_params(data, bus, theta, fleet, unit)
            if svc is not None:
                p["m"] *= svc["inertia_scale"]
                p["pss_gain"] *= svc["pss_scale"]
                p["d"] = svc["damping"]
            ini = sg_init(p, w, v, share)
            has_gov = governed and abs(ini["pm"]) > 1e-9
            gv = data.gov.get(bus, {})
            row = np.zeros(len(K.SG_COLS))
            vals = {
                "bus": data.idx[bus],
                "w": w,
                "ra": p["ra"],
                "xd": p["xd"],
                "xq": p["xq"],
                "x1": p["xd1"],
                "td10": p["td10"],
                "tq10": p["tq10"],
                "m": p["m"],
                "d": p["d"],
                "ka": p["ka"],
                "ta": p["ta"],
                "ks": p["pss_gain"],
                "t4": p["pss_lag"],
                "t5": p["pss_washout"],
                "t6": p["pss_wash_lag"],
                "pm": ini["pm"],
                "vref": ini["vref"],
                "flux_blend": svc["flux_blend"] if svc else 1.0,
                "avr_blend": svc["avr_blend"] if svc else 1.0,
                "has_gov": 1.0 if has_gov else 0.0,
                "gov_r": gv.get("R", 1.0),
                "gov_t1": gv.get("T1", 1.0),
                "gov_t2": gv.get("T2", 1.0),
                "gov_t3": gv.get("T3", 1.0),
                "gov_dt": gv.get("Dt", 0.0),
                "gov_pref": ini["pm"],
            }
            if has_gov:
                assert gv["VMIN"] <= ini["pm"] <= gv["VMAX"]
            for n, val in vals.items():
                row[K.SG[n]] = val
            y = w / complex(p["ra"], p["xd1"])
            K.stamp(gmat, nbus, data.idx[bus], y)
            sg_rows.append(row)
            sg_y.append(y)
            sg_x.append(
                [
                    ini["delta"],
                    1.0,
                    ini["eq1"],
                    ini["ed1"],
                    ini["efd"],
                    ini["pm"],
                    0.0,
                    ini["pm"],
                    ini["pm"],
                ]
            )
            labels_sg.append(bus)
        if replaced:
            par = dict(GFL_DEFAULTS)
            for grp, f in conv.get(bus, {}).items():
                for fld in CONV_GROUPS[grp]:
                    par[fld] *= f
            share = (
                complex(s.real, (1.0 - CONDENSER_R0010000["q_share"]) * s.imag)
                if cond
                else s
            )
            x, p_ref, q_ref, v_ref = gfl_init(par, rating, v, share)
            row = np.zeros(len(K.GFL_COLS))
            vals = dict(
                par,
                bus=data.idx[bus],
                w=rating,
                g=theta[0],
                leak=LEAK,
                p_ref=p_ref,
                q_ref=q_ref,
                v_ref=v_ref,
            )
            for n, val in vals.items():
                row[K.GF[n]] = val
            gf_rows.append(row)
            gf_x.append(x)
            labels_gf.append(bus)
    ld_bus, ld_s0, ld_y0 = [], [], []
    for b, sv in sorted(data.op["loads"].items(), key=lambda kv: int(kv[0])):
        b = int(b)
        s0 = complex(*sv)
        v0 = volt[b]
        y0 = np.conj(s0) / abs(v0) ** 2
        K.stamp(gmat, nbus, data.idx[b], y0)
        ld_bus.append(data.idx[b])
        ld_s0.append(s0)
        ld_y0.append(y0)
    load_list = sorted(int(b) for b in data.op["loads"])
    events = np.zeros((0, 5))
    if amplitude:
        li = load_list.index(PULSE["bus"])
        events = np.array(
            [[K.EV_LOAD, li, PULSE["t0"], PULSE["t1"], amplitude * ld_s0[li].real]]
        )
    return {
        "ginv": np.linalg.inv(gmat),
        "gmat": gmat,
        "net": net,
        "sg_par": np.array(sg_rows),
        "sg_y": np.array(sg_y, dtype=np.complex128),
        "sg_x0": np.array(sg_x),
        "gfl_par": np.array(gf_rows).reshape(len(gf_rows), len(K.GFL_COLS)),
        "gfl_x0": np.array(gf_x, dtype=np.float64).reshape(len(gf_x), 11),
        "ld_bus": np.array(ld_bus, dtype=np.int64),
        "ld_s0": np.array(ld_s0, dtype=np.complex128),
        "ld_y0": np.array(ld_y0, dtype=np.complex128),
        "tau_m": tau_m,
        "src_bus": np.zeros(0, np.int64),
        "src_y": np.zeros(0, np.complex128),
        "src_e0": np.zeros(0, np.complex128),
        "events": events,
        "labels_sg": labels_sg,
        "labels_gf": labels_gf,
        "order": data.order,
        "hashes": {
            "case": sha(case),
            "params_sg": sha(np.round(np.array(sg_rows), 12).tolist()),
            "params_gfl": sha(
                np.round(np.array(gf_rows), 12).tolist() if gf_rows else []
            ),
            "operating_point": sha({str(k): [v.real, v.imag] for k, v in volt.items()}),
        },
    }


def run(built: dict, dt: float, tlen: float, ds: int):
    nsteps = int(round(tlen / dt))
    rec = np.arange(len(built["order"]), dtype=np.int64)
    return K.simulate(
        built["ginv"],
        built["net"]["coe0"],
        built["net"]["vsol0"],
        built["net"]["brch_ihis"],
        built["net"]["node_ihis"],
        built["net"]["nbus"],
        built["sg_par"],
        built["sg_y"],
        built["sg_x0"],
        built["gfl_par"],
        built["gfl_x0"],
        built["ld_bus"],
        built["ld_s0"],
        built["ld_y0"],
        built["tau_m"],
        built["src_bus"],
        built["src_y"],
        built["src_e0"],
        built["events"],
        dt,
        nsteps,
        ds,
        rec,
    )
