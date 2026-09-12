# ruff: noqa: E501  -- kernel argument lists kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - V3 stopped at gate E0 (commit 229613de). Any use needs a new preregistration.
"""IEEE-39 case builder for the V3 scientific blocks (prereg V3 section 7).

The v1 case builder (overlay/tx4_case.build_case, frozen) is used unchanged; the only change is the
GFL realization: every GFL is the V2 average-value voltage source behind its own Rf-Lf filter
(overlay/tx4_emt_v2.py), stamped as an exact ideal-source companion at its bus and initialized from
the canonical power flow. Runs use tx4_emt_v2.simulate_v2 (identical to the v1 loop otherwise).
"""

from __future__ import annotations

import numpy as np
import tx4_case as C
import tx4_emt as K
import tx4_emt_v2 as K2


def build_v3(data, case, dt, amplitude=C.PULSE["fraction"], tau_m=C.TAU_M):
    b = C.build_case(data, case, dt, amplitude=amplitude, tau_m=tau_m)
    gmat = b["gmat"].copy()
    nbus = b["net"]["nbus"]
    volt, _ = C.point(data, case.get("variant", {"kind": "none"}))
    filt, ibr, ihis = [], [], []
    for i, row in enumerate(b["gfl_par"]):
        f = np.array(K2.filter_coeffs(row[K.GF_W], row[K.GF_RF], row[K.GF_XF], dt))
        K2.stamp_filter(gmat, nbus, int(row[K.GF_BUS]), f[2])
        x = b["gfl_x0"][i]
        ib, ih, _ = K2.filter_init(row, volt[b["labels_gf"][i]], complex(x[6], x[7]), x[0], f)
        filt.append(f)
        ibr.append(ib)
        ihis.append(ih)
    b["gmat"] = gmat
    b["ginv"] = np.linalg.inv(gmat)
    b["gf_filt"] = np.array(filt).reshape(len(filt), 5)
    b["gf_ibr0"] = np.array(ibr).reshape(len(ibr), 3)
    b["gf_ihis0"] = np.array(ihis).reshape(len(ihis), 3)
    return b


def build_native(data, case, dt, amplitude=C.PULSE["fraction"], tau_m=C.TAU_M):
    """Descriptive EMT04 holdout: every machine realized with an EMT-native stator (tx4_emt_v3)."""

    b = build_v3(data, case, dt, amplitude, tau_m)
    gmat = b["gmat"].copy()
    nbus = b["net"]["nbus"]
    volt, _ = C.point(data, case.get("variant", {"kind": "none"}))
    filt, ibr, ihis = [], [], []
    for i, row in enumerate(b["sg_par"]):
        bus_idx = int(row[K.SG_BUS])
        y = b["sg_y"][i]
        K.stamp(gmat, nbus, bus_idx, -y)
        w, ra, x1 = row[K.SG_W], row[K.SG_RA], row[K.SG_X1]
        r_s, l_s = ra / w, x1 / (w * K.W0)
        req = r_s + 2.0 * l_s / dt
        f = np.array([r_s, l_s, req, (2.0 * l_s / dt - r_s) / req, 1.0 / req])
        K2.stamp_filter(gmat, nbus, bus_idx, req)
        x = b["sg_x0"][i]
        epr = complex(x[2], -x[3]) * np.exp(1j * x[0])
        v0 = volt[b["labels_sg"][i]]
        isys = y * (epr - v0)
        ib = np.array([np.real(isys * np.exp(-2j * np.pi * k / 3)) for k in range(3)])
        u = np.array([np.real((epr - v0) * np.exp(-2j * np.pi * k / 3)) for k in range(3)])
        filt.append(f)
        ibr.append(ib)
        ihis.append(f[3] * ib + f[4] * u)
    b["gmat_native"] = gmat
    b["ginv_native"] = np.linalg.inv(gmat)
    b["sg_filt"] = np.array(filt).reshape(len(filt), 5)
    b["sg_ibr0"] = np.array(ibr).reshape(len(ibr), 3)
    b["sg_ihis0"] = np.array(ihis).reshape(len(ihis), 3)
    return b


def run_native(b, dt, tlen, ds):
    import tx4_emt_v3 as K3

    rec = np.arange(len(b["order"]), dtype=np.int64)
    t, sg, gf, v, fin = K3.simulate_v3(b["ginv_native"], b["net"]["coe0"], b["net"]["vsol0"], b["net"]["brch_ihis"], b["net"]["node_ihis"], b["net"]["nbus"],
                                       b["sg_par"], b["sg_y"], b["sg_x0"], np.ones(len(b["sg_par"])), b["sg_filt"], b["sg_ibr0"], b["sg_ihis0"],
                                       b["gfl_par"], b["gfl_x0"], b["gf_filt"], b["gf_ibr0"], b["gf_ihis0"], b["ld_bus"], b["ld_s0"], b["ld_y0"],
                                       b["tau_m"], b["src_bus"], b["src_y"], b["src_e0"], b["events"], dt, int(round(tlen / dt)), ds, rec)
    return t, sg, gf, v, np.zeros((len(t), len(b["gfl_par"]))), fin


def run_v3(b, dt, tlen, ds):
    rec = np.arange(len(b["order"]), dtype=np.int64)
    return K2.simulate_v2(b["ginv"], b["net"]["coe0"], b["net"]["vsol0"], b["net"]["brch_ihis"], b["net"]["node_ihis"], b["net"]["nbus"],
                          b["sg_par"], b["sg_y"], b["sg_x0"], b["gfl_par"], b["gfl_x0"], b["gf_filt"], b["gf_ibr0"], b["gf_ihis0"],
                          b["ld_bus"], b["ld_s0"], b["ld_y0"], b["tau_m"], b["src_bus"], b["src_y"], b["src_e0"], b["events"], dt,
                          int(round(tlen / dt)), ds, rec)
