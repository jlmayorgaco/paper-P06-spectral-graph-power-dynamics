# ruff: noqa: E501  -- kernel argument lists kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - V3 stopped at gate E0 (commit 229613de). Any use needs a new preregistration.
"""TX4 EMT kernel V3 (overlay; copied into external/ParaEMT_tx4 by EMT_setup_tx4.py).

Used ONLY for the descriptive EMT04 holdout "EMT-native stator" (prereg v1 section 3.1 and prereg V3
section 7): selected synchronous machines are realized as E' (predicted states) behind an EMT R-L
stator branch R = ra/w, L = x'/(w w0), with the same exact ideal-source trapezoidal companion as the
V2 GFL filter; the machine equations then use the MEASURED stator current instead of the algebraic
(E' - V)/(ra + j x'). Everything else is identical to tx4_emt_v2.simulate_v2 (which, with no native
machine, it reproduces). The primary V3 scientific runs use tx4_emt_v2.simulate_v2, not this file.
"""

from __future__ import annotations

import numba
import numpy as np
from lib_numba import numba_updateIhis
from tx4_emt import (
    GF_BUS,
    GF_PREF,
    GF_W,
    SG_AVRB,
    SG_BUS,
    SG_D,
    SG_FLUX,
    SG_GDT,
    SG_GPREF,
    SG_GR,
    SG_GT1,
    SG_GT2,
    SG_GT3,
    SG_HASGOV,
    SG_KA,
    SG_KS,
    SG_M,
    SG_PM,
    SG_RA,
    SG_T4,
    SG_T5,
    SG_T6,
    SG_TA,
    SG_TD10,
    SG_TQ10,
    SG_VREF,
    SG_W,
    SG_X1,
    SG_XD,
    SG_XQ,
    W0,
    _event_value,
    _gfl_eval,
    _inject,
    _phasors,
    _sg_eval,
)
from tx4_emt_v2 import _space_vector, e_command


@numba.njit(cache=False)
def _sg_eval_meas(x, vv, idev, p, pm_eff):
    """Two-axis machine with the stator current given (device base). Same equations as _sg_eval."""

    delta = x[0]
    sd, cd = np.sin(delta), np.cos(delta)
    ra = p[SG_RA]
    x1 = p[SG_X1]
    epr = complex(x[2], -x[3]) * complex(np.cos(delta), np.sin(delta))
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
    f[4] = p[SG_AVRB] * (p[SG_KA] * (p[SG_VREF] + p[SG_KS] * x[6] - term) - x[4]) / p[SG_TA]
    f[5] = (pe - x[5]) / p[SG_T6]
    f[6] = (p[SG_T5] * (pe - x[5]) / p[SG_T6] - x[6]) / p[SG_T4]
    s = vv * np.conj(p[SG_W] * idev)
    return f, pe, s.real, s.imag, epr


@numba.njit(cache=False)
def simulate_v3(ginv, coe0, vsol0, brch_ihis0, node_ihis0, nbus, sg_par, sg_y, sg_x0, sg_native, sg_filt, sg_ibr0, sg_ihis0,
                gfl_par, gfl_x0, gf_filt, gf_ibr0, gf_ihis0, ld_bus, ld_s0, ld_y0, tau_m, src_bus, src_y, src_e0, events, dt, nsteps, ds, rec_buses):
    """simulate_v2 plus native-stator machines (sg_native[i] = 1: E' behind the R-L companion sg_filt[i];
    the caller must NOT stamp that machine's Norton admittance but stamp 1/Req on its phase diagonals)."""

    nsg, ngf, nld, nsrc = sg_par.shape[0], gfl_par.shape[0], ld_bus.shape[0], src_bus.shape[0]
    nn = 3 * nbus
    vsol = vsol0.copy()
    brch_ihis = brch_ihis0.copy()
    node_ihis = node_ihis0.copy()
    sgx = sg_x0.copy()
    gfx = gfl_x0.copy()
    ibr = gf_ibr0.copy()
    ihis = gf_ihis0.copy()
    sbr = sg_ibr0.copy()
    shis = sg_ihis0.copy()
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
    ecur = np.zeros((ngf, 3))
    escur = np.zeros((nsg, 3))
    i_prev = np.zeros(ngf, dtype=np.complex128)
    i_prev2 = np.zeros(ngf, dtype=np.complex128)
    v_prev2 = np.zeros(ngf, dtype=np.complex128)
    is_meas = np.zeros(nsg, dtype=np.complex128)
    for i in range(nsg):
        if sg_native[i] > 0.5:
            is_meas[i] = _space_vector(sbr[i, 0], sbr[i, 1], sbr[i, 2], 0.0) / sg_par[i, SG_W]
    for i in range(ngf):
        b = int(gfl_par[i, GF_BUS])
        cur = _space_vector(ibr[i, 0], ibr[i, 1], ibr[i, 2], 0.0)
        i_prev[i] = cur
        i_prev2[i] = cur
        v_prev2[i] = vph[b]
        th = gfx[i, 0]
        idq = cur * complex(np.cos(th), -np.sin(th)) / gfl_par[i, GF_W]
        gfx[i, 6] = idq.real
        gfx[i, 7] = idq.imag
    k_out = 0
    finite = True
    for step in range(nsteps + 1):
        t_prev = step * dt
        for i in range(nsg):
            dpm, hit = _event_value(events, 2, i, t_prev)
            pm_eff = sg_par[i, SG_PM] * (1.0 + dpm)
            b = int(sg_par[i, SG_BUS])
            if sg_native[i] > 0.5:
                f, pe, ps, qs, epr = _sg_eval_meas(sgx[i], vph[b], is_meas[i], sg_par[i], pm_eff)
            else:
                f, pe, ps, qs, epr = _sg_eval(sgx[i], vph[b], sg_par[i], pm_eff)
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
        if step % ds == 0:
            t_out[k_out] = t_prev
            for r in range(rec_buses.shape[0]):
                v_out[k_out, r, 0] = vph[rec_buses[r]].real
                v_out[k_out, r, 1] = vph[rec_buses[r]].imag
            k_out += 1
        if step == nsteps:
            break
        t = (step + 1) * dt
        for i in range(nsg):
            sgxp[i] = sgx[i] + dt * sgf0[i]
        for i in range(ngf):
            gfxp[i] = gfx[i] + dt * gff0[i]
        vfp = np.empty(nld, dtype=np.complex128)
        for i in range(nld):
            vfp[i] = vf[i] + (dt / tau_m) * (vph[ld_bus[i]] - vf[i])
        irhs[:] = node_ihis
        for i in range(nsg):
            d = sgxp[i, 0]
            epr = complex(sgxp[i, 2], -sgxp[i, 3]) * complex(np.cos(d), np.sin(d))
            b = int(sg_par[i, SG_BUS])
            if sg_native[i] > 0.5:
                for k in range(3):
                    ang = W0 * t - 2.0 * np.pi * k / 3.0
                    ek = epr.real * np.cos(ang) - epr.imag * np.sin(ang)
                    escur[i, k] = ek
                    irhs[b + k * nbus] += ek / sg_filt[i, 2] + shis[i, k]
            else:
                _inject(irhs, nbus, b, sg_y[i] * epr, t)
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
        vsol = ginv @ irhs
        vnew = _phasors(vsol, nbus, t)
        for i in range(nsg):
            if sg_native[i] > 0.5:
                b = int(sg_par[i, SG_BUS])
                for k in range(3):
                    u = escur[i, k] - vsol[b + k * nbus]
                    ik = u / sg_filt[i, 2] + shis[i, k]
                    sbr[i, k] = ik
                    shis[i, k] = sg_filt[i, 3] * ik + sg_filt[i, 4] * u
                is_meas[i] = _space_vector(sbr[i, 0], sbr[i, 1], sbr[i, 2], t) / sg_par[i, SG_W]
        i_new = np.zeros(ngf, dtype=np.complex128)
        for i in range(ngf):
            b = int(gfl_par[i, GF_BUS])
            for k in range(3):
                u = ecur[i, k] - vsol[b + k * nbus]
                ik = u / gf_filt[i, 2] + ihis[i, k]
                ibr[i, k] = ik
                ihis[i, k] = gf_filt[i, 3] * ik + gf_filt[i, 4] * u
            i_new[i] = _space_vector(ibr[i, 0], ibr[i, 1], ibr[i, 2], t)
        for i in range(nsg):
            dpm, hit = _event_value(events, 2, i, t)
            pm_eff = sg_par[i, SG_PM] * (1.0 + dpm)
            b = int(sg_par[i, SG_BUS])
            if sg_native[i] > 0.5:
                f1, pe, ps, qs, epr = _sg_eval_meas(sgxp[i], vnew[b], is_meas[i], sg_par[i], pm_eff)
            else:
                f1, pe, ps, qs, epr = _sg_eval(sgxp[i], vnew[b], sg_par[i], pm_eff)
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
        brch_ipre, node_ihis = numba_updateIhis(brch_ihis, vsol, coe0, nn)
        if step % ds == 0 and not np.all(np.isfinite(vsol)):
            finite = False
            return t_out[:k_out], sg_out[:k_out], gf_out[:k_out], v_out[:k_out], finite
    return t_out[:k_out], sg_out[:k_out], gf_out[:k_out], v_out[:k_out], finite
