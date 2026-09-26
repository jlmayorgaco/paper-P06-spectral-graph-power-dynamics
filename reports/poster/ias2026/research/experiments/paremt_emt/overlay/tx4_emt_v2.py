# ruff: noqa: E501  -- kernel argument lists kept on one line
"""TX4 EMT kernel V2 (overlay; copied into external/ParaEMT_tx4 by EMT_setup_tx4.py).

Preregistration V2: docs/20260911_PAREMT_EMT_PREREG_V2.md (commit 0ab8efd4).
Interface: docs/20260911_PAREMT_GFL_VOLTAGE_INTERFACE_DERIVATION.md.

Identical to the v1 kernel (tx4_emt.simulate, left unchanged) except the grid-following converter:
an average-value controlled voltage source behind the physical Rf-Lf filter branch.
  - filter branch per phase: trapezoidal R-L companion (ParaEMT formulas), node E eliminated
    exactly (ideal source); the caller stamps 1/Req on the three phase diagonals of the GFL bus;
  - i_d, i_q are the MEASURED filter currents (never integrated); the controller integrates its
    other 9 states with the v1 Heun scheme;
  - command at t_n (derivation (C), section 7): e = v + rf i + (xf/w0)[f_i,TX4 + j(w0 + theta') i]
    evaluated at the predicted controller states and at linearly extrapolated terminal voltage
    and filter current phasors.
"""

from __future__ import annotations

import numba
import numpy as np
from lib_numba import numba_updateIhis
from tx4_emt import (
    GF_BUS,
    GF_PREF,
    GF_RF,
    GF_W,
    GF_XF,
    SG_BUS,
    SG_PM,
    W0,
    _event_value,
    _gfl_eval,
    _inject,
    _phasors,
    _sg_eval,
)


def filter_coeffs(w: float, rf: float, xf: float, dt: float):
    """System-base filter and its trapezoidal companion: (Rf, Lf, Req, icf, Gv1)."""

    r = rf / w
    ell = xf / (w * W0)
    req = r + 2.0 * ell / dt
    return r, ell, req, (2.0 * ell / dt - r) / req, 1.0 / req


@numba.njit(cache=False)
def _space_vector(x0, x1, x2, t):
    a = complex(np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3))
    return (2.0 / 3.0) * (x0 + a * x1 + a * a * x2) * complex(np.cos(W0 * t), -np.sin(W0 * t))


@numba.njit(cache=False)
def e_command(x, vv, idq, p, pref_eff):
    """Derivation (C). x: controller vector (11 entries, x[0] = theta); vv: terminal voltage phasor
    (network synchronous frame); idq: filter current, device base, PLL frame. Returns (e_dq, f)."""

    xx = x.copy()
    xx[6] = idq.real
    xx[7] = idq.imag
    f, ps, qs, isys, thd = _gfl_eval(xx, vv, p, pref_eff)
    th = xx[0]
    vdq = vv * complex(np.cos(th), -np.sin(th))
    fi = complex(f[6], f[7])
    e = vdq + p[GF_RF] * idq + (p[GF_XF] / W0) * (fi + 1j * (W0 + thd) * idq)
    return e, f


@numba.njit(cache=False)
def simulate_v2(
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
    gf_filt,
    gf_ibr0,
    gf_ihis0,
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
    """gf_filt[i] = (Rf, Lf, Req, icf, Gv1) system base; gf_ibr0/gf_ihis0: (ngf, 3) initial branch
    currents and histories. Returns t, sg_out, gf_out[n, ngf, 14] (11 quantities with measured
    i_d, i_q; P, Q system base from the branch-current phasor; theta'), v_out, nyq_out[n, ngf]
    (max ||V_T,n| - |V_T,n-1|| over the steps since the previous record), finite."""

    nsg, ngf, nld, nsrc = sg_par.shape[0], gfl_par.shape[0], ld_bus.shape[0], src_bus.shape[0]
    nn = 3 * nbus
    vsol = vsol0.copy()
    brch_ihis = brch_ihis0.copy()
    node_ihis = node_ihis0.copy()
    sgx = sg_x0.copy()
    gfx = gfl_x0.copy()
    ibr = gf_ibr0.copy()
    ihis = gf_ihis0.copy()
    vph = _phasors(vsol, nbus, 0.0)
    vf = np.empty(nld, dtype=np.complex128)
    for i in range(nld):
        vf[i] = vph[ld_bus[i]]
    nout = nsteps // ds + 1
    t_out = np.zeros(nout)
    sg_out = np.zeros((nout, nsg, 12))
    gf_out = np.zeros((nout, ngf, 14))
    v_out = np.zeros((nout, rec_buses.shape[0], 2))
    nyq_out = np.zeros((nout, ngf))
    sgf0 = np.zeros((nsg, 9))
    gff0 = np.zeros((ngf, 11))
    sgxp = np.zeros((nsg, 9))
    gfxp = np.zeros((ngf, 11))
    irhs = np.zeros(nn)
    ecur = np.zeros((ngf, 3))
    i_prev = np.zeros(ngf, dtype=np.complex128)
    i_prev2 = np.zeros(ngf, dtype=np.complex128)
    v_prev2 = np.zeros(ngf, dtype=np.complex128)
    vmag_last = np.zeros(ngf)
    nyq_blk = np.zeros(ngf)
    for i in range(ngf):
        b = int(gfl_par[i, GF_BUS])
        cur = _space_vector(ibr[i, 0], ibr[i, 1], ibr[i, 2], 0.0)
        i_prev[i] = cur
        i_prev2[i] = cur
        v_prev2[i] = vph[b]
        vmag_last[i] = abs(vph[b])
        th = gfx[i, 0]
        idq = cur * complex(np.cos(th), -np.sin(th)) / gfl_par[i, GF_W]
        gfx[i, 6] = idq.real
        gfx[i, 7] = idq.imag
    k_out = 0
    finite = True
    for step in range(nsteps + 1):
        t_prev = step * dt
        # ---- derivatives at (x_{n-1}, V_{n-1}, i_{n-1}), and recording ----
        for i in range(nsg):
            dpm, hit = _event_value(events, 2, i, t_prev)
            pm_eff = sg_par[i, SG_PM] * (1.0 + dpm)
            f, pe, ps, qs, epr = _sg_eval(sgx[i], vph[int(sg_par[i, SG_BUS])], sg_par[i], pm_eff)
            sgf0[i] = f
            if step % ds == 0:
                sg_out[k_out, i, :9] = sgx[i]
                sg_out[k_out, i, 9] = pe
                sg_out[k_out, i, 10] = ps
                sg_out[k_out, i, 11] = qs
        for i in range(ngf):
            b = int(gfl_par[i, GF_BUS])
            dpr, hit = _event_value(events, 3, i, t_prev)
            f, ps, qs, isys, thd = _gfl_eval(gfx[i], vph[b], gfl_par[i], gfl_par[i, GF_PREF] * (1.0 + dpr))
            gff0[i] = f
            if step % ds == 0:
                s = vph[b] * np.conj(i_prev[i])
                gf_out[k_out, i, :11] = gfx[i]
                gf_out[k_out, i, 11] = s.real
                gf_out[k_out, i, 12] = s.imag
                gf_out[k_out, i, 13] = thd
                nyq_out[k_out, i] = nyq_blk[i]
                nyq_blk[i] = 0.0
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
            b = int(gfl_par[i, GF_BUS])
            thp = gfxp[i, 0]
            rot = complex(np.cos(thp), np.sin(thp))
            v_hat = 2.0 * vph[b] - v_prev2[i]
            i_hat = (2.0 * i_prev[i] - i_prev2[i]) * np.conj(rot) / gfl_par[i, GF_W]
            dpr, hit = _event_value(events, 3, i, t)
            edq, fcmd = e_command(gfxp[i], v_hat, i_hat, gfl_par[i], gfl_par[i, GF_PREF] * (1.0 + dpr))
            esys = edq * rot
            for k in range(3):
                ang = W0 * t - 2.0 * np.pi * k / 3.0
                ek = esys.real * np.cos(ang) - esys.imag * np.sin(ang)
                ecur[i, k] = ek
                irhs[b + k * nbus] += ek / gf_filt[i, 2] + ihis[i, k]
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
        # ---- filter branches (companion update) and measured currents ----
        i_new = np.zeros(ngf, dtype=np.complex128)
        for i in range(ngf):
            b = int(gfl_par[i, GF_BUS])
            for k in range(3):
                u = ecur[i, k] - vsol[b + k * nbus]
                ik = u / gf_filt[i, 2] + ihis[i, k]
                ibr[i, k] = ik
                ihis[i, k] = gf_filt[i, 3] * ik + gf_filt[i, 4] * u
            i_new[i] = _space_vector(ibr[i, 0], ibr[i, 1], ibr[i, 2], t)
            dv = abs(abs(vnew[b]) - vmag_last[i])
            if dv > nyq_blk[i]:
                nyq_blk[i] = dv
            vmag_last[i] = abs(vnew[b])
        # ---- corrector ----
        for i in range(nsg):
            dpm, hit = _event_value(events, 2, i, t)
            pm_eff = sg_par[i, SG_PM] * (1.0 + dpm)
            f1, pe, ps, qs, epr = _sg_eval(sgxp[i], vnew[int(sg_par[i, SG_BUS])], sg_par[i], pm_eff)
            sgx[i] = sgx[i] + 0.5 * dt * (sgf0[i] + f1)
        for i in range(ngf):
            b = int(gfl_par[i, GF_BUS])
            thp = gfxp[i, 0]
            xc = gfxp[i].copy()
            idq_p = i_new[i] * complex(np.cos(thp), -np.sin(thp)) / gfl_par[i, GF_W]
            xc[6] = idq_p.real
            xc[7] = idq_p.imag
            dpr, hit = _event_value(events, 3, i, t)
            f1, ps, qs, isys, thd = _gfl_eval(xc, vnew[b], gfl_par[i], gfl_par[i, GF_PREF] * (1.0 + dpr))
            for j in range(11):
                if j != 6 and j != 7:
                    gfx[i, j] = gfx[i, j] + 0.5 * dt * (gff0[i, j] + f1[j])
            th = gfx[i, 0]
            idq = i_new[i] * complex(np.cos(th), -np.sin(th)) / gfl_par[i, GF_W]
            gfx[i, 6] = idq.real
            gfx[i, 7] = idq.imag
            v_prev2[i] = vph[b]
            i_prev2[i] = i_prev[i]
            i_prev[i] = i_new[i]
        for i in range(nld):
            vf[i] = vf[i] + 0.5 * (dt / tau_m) * ((vph[ld_bus[i]] - vf[i]) + (vnew[ld_bus[i]] - vfp[i]))
        vph = vnew
        # ---- network history (ParaEMT) ----
        brch_ipre, node_ihis = numba_updateIhis(brch_ihis, vsol, coe0, nn)
        if step % ds == 0 and not np.all(np.isfinite(vsol)):
            finite = False
            return t_out[:k_out], sg_out[:k_out], gf_out[:k_out], v_out[:k_out], nyq_out[:k_out], finite
    return t_out[:k_out], sg_out[:k_out], gf_out[:k_out], v_out[:k_out], nyq_out[:k_out], finite


def filter_init(par_row: np.ndarray, vv: complex, idq0: complex, theta0: float, filt, t0: float = 0.0):
    """Section 8: branch currents and histories from the continuous phasor steady state."""

    w = par_row[GF_W]
    isys = w * idq0 * np.exp(1j * theta0)
    esys = vv + complex(filt[0], W0 * filt[1]) * isys
    ibr = np.array([np.real(isys * np.exp(1j * (W0 * t0 - 2 * np.pi * k / 3))) for k in range(3)])
    u = np.array([np.real((esys - vv) * np.exp(1j * (W0 * t0 - 2 * np.pi * k / 3))) for k in range(3)])
    ihis = filt[3] * ibr + filt[4] * u
    return ibr, ihis, esys


def stamp_filter(gmat: np.ndarray, nbus: int, bus_idx: int, req: float):
    for k in range(3):
        gmat[bus_idx + k * nbus, bus_idx + k * nbus] += 1.0 / req
