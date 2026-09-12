# ruff: noqa: E501  -- parameter tables and kernel argument lists kept on one line
"""TX4 EMT adapter for ParaEMT (overlay file, copied into external/ParaEMT_tx4 by EMT_setup_tx4.py).

Preregistration: docs/20260911_PAREMT_EMT_PREREG_V1.md (commit c2947bd8), sections 2-4.

Network: ParaEMT's own companion models (numba_InitNet, TX4-patched for from-side taps and a
numerical-damping switch) and its history update (numba_updateIhis), trapezoidal rule, three
phases in phase-block node order, amplitude per unit on 100 MVA.

Devices (the frozen phasor equations, synchronous frame at w0 = 2 pi 60):
  - two-axis machine + first-order AVR + washout-lag PSS (+ optional TGOV1N, condenser blends);
    stator = synchronous-frame algebraic Norton (exact because x'd = x'q);
  - 11-state grid-following converter: state-driven average-value current source;
  - constant-power load: constant-Z synchronous-frame Norton + compensation from a 1 ms filtered
    voltage phasor;
  - fixed Thevenin sources (unit-test infinite bus).
Synchronous-frame algebraic Norton with admittance Y at a bus: real 3x3 conductance
G_km = (2/3) Re(a^(m-k) Y) plus source current Re(Y E exp(j(w0 t - 2 pi k/3))).
Time loop: predictor (Euler) - network solve - corrector (trapezoidal/Heun); one solve per step.
"""

from __future__ import annotations

import numba
import numpy as np
from lib_numba import numba_InitNet, numba_updateIhis

W0 = 2.0 * np.pi * 60.0
A1 = np.exp(2j * np.pi / 3.0)

# column layouts -------------------------------------------------------------
SG_COLS = (
    "bus",
    "w",
    "ra",
    "xd",
    "xq",
    "x1",
    "td10",
    "tq10",
    "m",
    "d",
    "ka",
    "ta",
    "ks",
    "t4",
    "t5",
    "t6",
    "pm",
    "vref",
    "flux_blend",
    "avr_blend",
    "has_gov",
    "gov_r",
    "gov_t1",
    "gov_t2",
    "gov_t3",
    "gov_dt",
    "gov_pref",
)
SG = {n: i for i, n in enumerate(SG_COLS)}
(
    SG_BUS,
    SG_W,
    SG_RA,
    SG_XD,
    SG_XQ,
    SG_X1,
    SG_TD10,
    SG_TQ10,
    SG_M,
    SG_D,
    SG_KA,
    SG_TA,
    SG_KS,
    SG_T4,
    SG_T5,
    SG_T6,
    SG_PM,
    SG_VREF,
    SG_FLUX,
    SG_AVRB,
    SG_HASGOV,
    SG_GR,
    SG_GT1,
    SG_GT2,
    SG_GT3,
    SG_GDT,
    SG_GPREF,
) = range(27)
SG_STATES = (
    "delta",
    "omega",
    "eq1",
    "ed1",
    "efd",
    "pss_w",
    "pss_l",
    "gov_x1",
    "gov_x2",
)
GFL_COLS = (
    "bus",
    "w",
    "kp_pll",
    "ki_pll",
    "tau_p",
    "kp_p",
    "ki_p",
    "kp_q",
    "ki_q",
    "kp_i",
    "ki_i",
    "xf",
    "rf",
    "kp_v",
    "ki_v",
    "g",
    "leak",
    "p_ref",
    "q_ref",
    "v_ref",
)
GF = {n: i for i, n in enumerate(GFL_COLS)}
(
    GF_BUS,
    GF_W,
    GF_KPPLL,
    GF_KIPLL,
    GF_TAUP,
    GF_KPP,
    GF_KIP,
    GF_KPQ,
    GF_KIQ,
    GF_KPI,
    GF_KII,
    GF_XF,
    GF_RF,
    GF_KPV,
    GF_KIV,
    GF_G,
    GF_LEAK,
    GF_PREF,
    GF_QREF,
    GF_VREF,
) = range(20)
GFL_STATES = (
    "theta",
    "x_pll",
    "p_f",
    "q_f",
    "x_p",
    "x_q",
    "i_d",
    "i_q",
    "x_id",
    "x_iq",
    "x_v",
)
# event kinds: 1 load dP (abs pu) on [t0,t1); 2 SG Pm x(1+v) on [t0,t1); 3 GFL p_ref x(1+v) t>=t0;
#              4 source |E| x(1+v) t>=t0; 5 source phase +v t>=t0
EV_LOAD, EV_PM, EV_PREF, EV_EMAG, EV_EPHASE = 1, 2, 3, 4, 5


# ------------------------------------------------------------------ network --
def build_network(
    net: dict, voltages: dict, ts: float, line_scale=None, net_damping: float = 0.0
):
    """TX4 network through ParaEMT's numba_InitNet (TX4 patch). voltages: bus -> complex."""

    order = [int(b["idx"]) for b in net["buses"]]
    bus_num = np.array(order, dtype=np.int64)
    v = np.array([voltages[b] for b in order], dtype=np.complex128)
    lf, lt, lrx, lchg, xf_, xt, xrx, xk = [], [], [], [], [], [], [], []
    for ell, ln in enumerate(net["lines"]):
        assert (
            float(ln["u"]) == 1.0 and float(ln["phi"]) == 0.0 and float(ln["g"]) == 0.0
        )
        s = line_scale[1] if (line_scale is not None and ell == line_scale[0]) else 1.0
        rx = complex(float(ln["r"]), float(ln["x"])) / s
        tap = float(ln["tap"])
        if tap != 1.0 or float(ln["trans"]) == 1.0:
            assert float(ln["b"]) == 0.0
            xf_.append(int(ln["bus1"]))
            xt.append(int(ln["bus2"]))
            xrx.append(rx)
            xk.append(tap)
        else:
            lf.append(int(ln["bus1"]))
            lt.append(int(ln["bus2"]))
            lrx.append(rx)
            lchg.append(float(ln["b"]) * s)
    shb = np.array([int(x["bus"]) for x in net["shunts"]], dtype=np.int64)
    shg = np.array(
        [100.0 * complex(float(x["g"]), float(x["b"])) for x in net["shunts"]],
        dtype=np.complex128,
    )
    for x in net["shunts"]:
        assert float(x["g"]) == 0.0
    e_i, e_f, e_c = (
        np.zeros(0, np.int64),
        np.zeros(0, np.float64),
        np.zeros(0, np.complex128),
    )
    out = numba_InitNet(
        100.0,
        W0,
        bus_num,
        np.ones(len(order)),
        np.abs(v),
        np.angle(v),
        e_i,
        e_f,
        e_f,
        np.array(lf, np.int64),
        np.array(lt, np.int64),
        np.array(lrx, np.complex128),
        np.array(lchg, np.float64),
        np.array(xf_, np.int64),
        np.array(xt, np.int64),
        np.array(xrx, np.complex128),
        e_i,
        e_f,
        e_f,
        shb,
        shg,
        e_i,
        e_c,
        ts,
        2,
        np.array(xk, np.float64),
        float(net_damping),
    )
    (
        _,
        _,
        _,
        _,
        _,
        vt,
        _,
        n_nodes,
        coe0,
        vnet,
        _,
        brch_ipre,
        node_ihis,
        brch_ihis,
        rows,
        cols,
        data,
    ) = out
    g0 = np.zeros((n_nodes, n_nodes))
    np.add.at(g0, (rows, cols), data)
    return {
        "order": order,
        "nbus": len(order),
        "coe0": coe0,
        "g0": g0,
        "vsol0": np.real(vt).copy(),
        "brch_ihis": brch_ihis.copy(),
        "node_ihis": node_ihis.copy(),
        "brch_ipre": brch_ipre.copy(),
        "lines": (lf, lt, lrx, lchg),
        "xfmrs": (xf_, xt, xrx, xk),
        "shunts": (shb.tolist(), shg.tolist()),
    }


def norton_block(y: complex) -> np.ndarray:
    g = np.zeros((3, 3))
    for k in range(3):
        for m in range(3):
            g[k, m] = (2.0 / 3.0) * np.real(A1 ** (m - k) * y)
    return g


def stamp(gmat, nbus, bus_idx, y):
    blk = norton_block(y)
    nodes = [bus_idx, bus_idx + nbus, bus_idx + 2 * nbus]
    for k in range(3):
        for m in range(3):
            gmat[nodes[k], nodes[m]] += blk[k, m]


# ------------------------------------------------------------------ devices --
@numba.njit(cache=False)
def _sg_eval(x, vv, p, pm_eff):
    """Two-axis machine: returns (derivative[9], Pe, P_sys, Q_sys, I_sys)."""

    delta = x[0]
    sd, cd = np.sin(delta), np.cos(delta)
    ra, x1, w = p[SG_RA], p[SG_X1], p[SG_W]
    epr = complex(x[2], -x[3]) * complex(
        np.cos(delta), np.sin(delta)
    )  # (e'q - j e'd) e^{j delta}
    idev = (epr - vv) / complex(ra, x1)
    id_ = idev.real * sd - idev.imag * cd
    iq = idev.real * cd + idev.imag * sd
    vd = vv.real * sd - vv.imag * cd
    vq = vv.real * cd + vv.imag * sd
    pe = vd * id_ + vq * iq + ra * (id_ * id_ + iq * iq)
    term = abs(vv)
    f = np.zeros(9)
    omega = x[1]
    pmech = pm_eff
    if p[SG_HASGOV] > 0.5:
        y = (p[SG_GT2] / p[SG_GT3]) * (x[7] - x[8]) + x[8]
        pmech = y - p[SG_GDT] * (omega - 1.0)
        pd = p[SG_GPREF] - (omega - 1.0) / p[SG_GR]
        f[7] = (pd - x[7]) / p[SG_GT1]
        f[8] = (x[7] - x[8]) / p[SG_GT3]
    f[0] = W0 * (omega - 1.0)
    f[1] = (pmech - pe - p[SG_D] * (omega - 1.0)) / p[SG_M]
    f[2] = p[SG_FLUX] * (x[4] - x[2] - (p[SG_XD] - x1) * id_) / p[SG_TD10]
    f[3] = p[SG_FLUX] * (-x[3] + (p[SG_XQ] - x1) * iq) / p[SG_TQ10]
    f[4] = (
        p[SG_AVRB]
        * (p[SG_KA] * (p[SG_VREF] + p[SG_KS] * x[6] - term) - x[4])
        / p[SG_TA]
    )
    f[5] = (pe - x[5]) / p[SG_T6]
    f[6] = (p[SG_T5] * (pe - x[5]) / p[SG_T6] - x[6]) / p[SG_T4]
    isys = w * idev
    s = vv * np.conj(isys)
    return f, pe, s.real, s.imag, epr


@numba.njit(cache=False)
def _gfl_eval(x, vv, p, pref_eff):
    """11-state GFL: returns (derivative[11], P_sys, Q_sys, I_sys, theta_dot)."""

    th = x[0]
    rot = vv * complex(np.cos(th), -np.sin(th))
    v_d, v_q = rot.real, rot.imag
    i_d, i_q = x[6], x[7]
    power = v_d * i_d + v_q * i_q
    reactive = v_q * i_d - v_d * i_q
    err = p[GF_G] * (p[GF_VREF] - abs(vv))
    q_cmd = p[GF_KPV] * err + x[10]
    id_ref = p[GF_KPP] * (pref_eff - x[2]) + x[4]
    iq_ref = -(p[GF_KPQ] * (q_cmd - x[3]) + x[5])
    xf, rf, kp_i = p[GF_XF], p[GF_RF], p[GF_KPI]
    e_d = v_d + kp_i * (id_ref - i_d) + x[8] - xf * i_q
    e_q = v_q + kp_i * (iq_ref - i_q) + x[9] + xf * i_d
    f = np.zeros(11)
    f[0] = p[GF_KPPLL] * v_q + x[1]
    f[1] = p[GF_KIPLL] * v_q
    f[2] = (power - x[2]) / p[GF_TAUP]
    f[3] = (reactive - x[3]) / p[GF_TAUP]
    f[4] = p[GF_KIP] * (pref_eff - x[2])
    f[5] = p[GF_KIQ] * (q_cmd - x[3])
    f[6] = (W0 / xf) * (e_d - v_d - rf * i_d + xf * i_q)
    f[7] = (W0 / xf) * (e_q - v_q - rf * i_q - xf * i_d)
    f[8] = p[GF_KII] * (id_ref - i_d)
    f[9] = p[GF_KII] * (iq_ref - i_q)
    f[10] = p[GF_KIV] * err - p[GF_LEAK] * (x[10] - p[GF_QREF])
    w = p[GF_W]
    isys = w * complex(i_d, i_q) * complex(np.cos(th), np.sin(th))
    return f, w * power, w * reactive, isys, f[0]


@numba.njit(cache=False)
def _phasors(vsol, nbus, t):
    out = np.empty(nbus, dtype=np.complex128)
    a = complex(np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3))
    rot = complex(np.cos(W0 * t), -np.sin(W0 * t))
    for b in range(nbus):
        s = (2.0 / 3.0) * (vsol[b] + a * vsol[b + nbus] + a * a * vsol[b + 2 * nbus])
        out[b] = s * rot
    return out


@numba.njit(cache=False)
def _inject(irhs, nbus, bus, cur, t):
    """Add Re(cur exp(j(w0 t - 2 pi k/3))) to the three phase nodes of `bus`."""

    for k in range(3):
        ang = W0 * t - 2.0 * np.pi * k / 3.0
        irhs[bus + k * nbus] += cur.real * np.cos(ang) - cur.imag * np.sin(ang)


@numba.njit(cache=False)
def _event_value(ev, kind, idx, t):
    out = 0.0
    hit = False
    for e in range(ev.shape[0]):
        if int(ev[e, 0]) == kind and int(ev[e, 1]) == idx:
            if kind == 1 or kind == 2:
                if ev[e, 2] <= t < ev[e, 3]:
                    out += ev[e, 4]
                    hit = True
            else:
                if t >= ev[e, 2]:
                    out += ev[e, 4]
                    hit = True
    return out, hit


@numba.njit(cache=False)
def simulate(
    ginv,
    coe0,
    vsol0,
    brch_ihis0,
    node_ihis0,
    nbus,
    sg_par,
    sg_y,
    sg_x0,
    gfl_par,
    gfl_x0,
    ld_bus,
    ld_s0,
    ld_y0,
    tau_m,
    src_bus,
    src_y,
    src_e0,
    events,
    dt,
    nsteps,
    ds,
    rec_buses,
):
    nsg, ngf, nld, nsrc = (
        sg_par.shape[0],
        gfl_par.shape[0],
        ld_bus.shape[0],
        src_bus.shape[0],
    )
    nn = 3 * nbus
    vsol = vsol0.copy()
    brch_ihis = brch_ihis0.copy()
    node_ihis = node_ihis0.copy()
    sgx = sg_x0.copy()
    gfx = gfl_x0.copy()
    vph = _phasors(vsol, nbus, 0.0)
    vf = np.empty(nld, dtype=np.complex128)
    for i in range(nld):
        vf[i] = vph[ld_bus[i]]
    nout = nsteps // ds + 1
    t_out = np.zeros(nout)
    sg_out = np.zeros((nout, nsg, 12))
    gf_out = np.zeros((nout, ngf, 14))
    v_out = np.zeros((nout, rec_buses.shape[0], 2))
    sgf0 = np.zeros((nsg, 9))
    gff0 = np.zeros((ngf, 11))
    sgxp = np.zeros((nsg, 9))
    gfxp = np.zeros((ngf, 11))
    irhs = np.zeros(nn)
    k_out = 0
    finite = True
    for step in range(nsteps + 1):
        t_prev = step * dt
        # ---- derivatives at (x_{n-1}, V_{n-1}), and recording ----
        for i in range(nsg):
            dpm, hit = _event_value(events, 2, i, t_prev)
            pm_eff = sg_par[i, SG_PM] * (1.0 + dpm)
            f, pe, ps, qs, epr = _sg_eval(
                sgx[i], vph[int(sg_par[i, SG_BUS])], sg_par[i], pm_eff
            )
            sgf0[i] = f
            if step % ds == 0:
                sg_out[k_out, i, :9] = sgx[i]
                sg_out[k_out, i, 9] = pe
                sg_out[k_out, i, 10] = ps
                sg_out[k_out, i, 11] = qs
        for i in range(ngf):
            dpr, hit = _event_value(events, 3, i, t_prev)
            f, ps, qs, isys, thd = _gfl_eval(
                gfx[i],
                vph[int(gfl_par[i, GF_BUS])],
                gfl_par[i],
                gfl_par[i, GF_PREF] * (1.0 + dpr),
            )
            gff0[i] = f
            if step % ds == 0:
                gf_out[k_out, i, :11] = gfx[i]
                gf_out[k_out, i, 11] = ps
                gf_out[k_out, i, 12] = qs
                gf_out[k_out, i, 13] = thd
        if step % ds == 0:
            t_out[k_out] = t_prev
            for r in range(rec_buses.shape[0]):
                v_out[k_out, r, 0] = vph[rec_buses[r]].real
                v_out[k_out, r, 1] = vph[rec_buses[r]].imag
            k_out += 1
        if step == nsteps:
            break
        t = (step + 1) * dt
        # ---- predictor ----
        for i in range(nsg):
            sgxp[i] = sgx[i] + dt * sgf0[i]
        for i in range(ngf):
            gfxp[i] = gfx[i] + dt * gff0[i]
        vfp = np.empty(nld, dtype=np.complex128)
        for i in range(nld):
            vfp[i] = vf[i] + (dt / tau_m) * (vph[ld_bus[i]] - vf[i])
        # ---- source currents at t ----
        irhs[:] = node_ihis
        for i in range(nsg):
            d = sgxp[i, 0]
            epr = complex(sgxp[i, 2], -sgxp[i, 3]) * complex(np.cos(d), np.sin(d))
            _inject(irhs, nbus, int(sg_par[i, SG_BUS]), sg_y[i] * epr, t)
        for i in range(ngf):
            th = gfxp[i, 0]
            isys = (
                gfl_par[i, GF_W]
                * complex(gfxp[i, 6], gfxp[i, 7])
                * complex(np.cos(th), np.sin(th))
            )
            _inject(irhs, nbus, int(gfl_par[i, GF_BUS]), isys, t)
        for i in range(nld):
            dp, hit = _event_value(events, 1, i, t)
            s = ld_s0[i] + dp
            icomp = -(np.conj(s / vfp[i]) - ld_y0[i] * vfp[i])
            _inject(irhs, nbus, ld_bus[i], icomp, t)
        for i in range(nsrc):
            dm, h1 = _event_value(events, 4, i, t)
            dph, h2 = _event_value(events, 5, i, t)
            e = src_e0[i] * (1.0 + dm) * complex(np.cos(dph), np.sin(dph))
            _inject(irhs, nbus, src_bus[i], src_y[i] * e, t)
        # ---- network solve ----
        vsol = ginv @ irhs
        vnew = _phasors(vsol, nbus, t)
        # ---- corrector ----
        for i in range(nsg):
            dpm, hit = _event_value(events, 2, i, t)
            pm_eff = sg_par[i, SG_PM] * (1.0 + dpm)
            f1, pe, ps, qs, epr = _sg_eval(
                sgxp[i], vnew[int(sg_par[i, SG_BUS])], sg_par[i], pm_eff
            )
            sgx[i] = sgx[i] + 0.5 * dt * (sgf0[i] + f1)
        for i in range(ngf):
            dpr, hit = _event_value(events, 3, i, t)
            f1, ps, qs, isys, thd = _gfl_eval(
                gfxp[i],
                vnew[int(gfl_par[i, GF_BUS])],
                gfl_par[i],
                gfl_par[i, GF_PREF] * (1.0 + dpr),
            )
            gfx[i] = gfx[i] + 0.5 * dt * (gff0[i] + f1)
        for i in range(nld):
            vf[i] = vf[i] + 0.5 * (dt / tau_m) * (
                (vph[ld_bus[i]] - vf[i]) + (vnew[ld_bus[i]] - vfp[i])
            )
        vph = vnew
        # ---- network history (ParaEMT) ----
        brch_ipre, node_ihis = numba_updateIhis(brch_ihis, vsol, coe0, nn)
        if step % ds == 0 and not np.all(np.isfinite(vsol)):
            finite = False
            return t_out[:k_out], sg_out[:k_out], gf_out[:k_out], v_out[:k_out], finite
    return t_out[:k_out], sg_out[:k_out], gf_out[:k_out], v_out[:k_out], finite
