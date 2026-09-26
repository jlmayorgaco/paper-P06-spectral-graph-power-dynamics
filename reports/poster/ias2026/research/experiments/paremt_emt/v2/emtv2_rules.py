# ruff: noqa: E501  -- rule tables kept on one line
"""Comparison operators of preregistration V2 (docs/20260911_PAREMT_EMT_PREREG_V2.md section 4), frozen.

v1-style rule (G3a, G3b, G4a): raw 1-ms samples, max|x - x_ref| <= 0.02 * max|x_ref - x_ref(0)|.
G4b: LP = zero-phase 4th-order Butterworth low-pass, 15 Hz, at the native rate, sampled on the 1-ms
grid; window W = [1.0, 4.8] s; tolerances 0.02*exc + 1e-9 (trajectories), 0.01*exc + 1e-9 (50 vs 25 us).
Nyquist indicator A(k) = max over the steps in 1-ms block k of ||V_T,n| - |V_T,n-1||.
Estimator R: multi-channel matrix pencil on LP'd channels decimated to 100 Hz over [1.05, 4.8] s,
band 0.1-15 Hz, energy >= 1 % of in-band total, resolved iff residual <= 0.05.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

LP_FC, LP_ORDER = 15.0, 4
W_LO, W_HI = 1.0, 4.8
FLOOR = 1e-9
TOL_TRAJ, TOL_CONV = 0.02, 0.01
R_T0, R_T1, R_DT = 1.05, 4.8, 0.01
R_BAND = (0.1, 15.0)
SV_REL, MAX_ORDER, ENERGY_FRAC, RESID_MAX = 1e-4, 30, 0.01, 0.05


def v1_rule(t_e, x_e, t_p, x_p, labels, t_from, name):
    """v1 EMT02/EMT03 rule on raw 1-ms samples (no floor)."""

    n = min(len(t_e), len(t_p))
    if len(t_e) < len(t_p):
        return [{"test": name, "state": lab, "max_abs_err": float("inf"), "max_excursion": float("nan"), "ratio": float("inf"),
                 "pass": False, "emt_nonfinite_at_s": float(t_e[-1])} for lab in labels]
    assert np.allclose(t_e[:n], t_p[:n], atol=1e-9)
    m = t_e[:n] >= t_from - 1e-12
    rows = []
    for j, lab in enumerate(labels):
        err = float(np.abs(x_e[:n][m, j] - x_p[:n][m, j]).max())
        exc = float(np.abs(x_p[:n][m, j] - x_p[0, j]).max())
        rows.append({"test": name, "state": lab, "max_abs_err": err, "max_excursion": exc,
                     "ratio": err / exc if exc > 0 else (0.0 if err == 0 else float("inf")), "pass": bool(err <= TOL_TRAJ * exc)})
    return rows


def lp(x, dt):
    sos = signal.butter(LP_ORDER, LP_FC, fs=1.0 / dt, output="sos")
    return signal.sosfiltfilt(sos, x, axis=0)


def on_ms_grid(x, dt, tmax):
    step = int(round(1e-3 / dt))
    n = int(round(tmax / 1e-3)) + 1
    return x[: (n - 1) * step + 1 : step][:n]


def window_mask(t):
    return (t >= W_LO - 1e-12) & (t <= W_HI + 1e-12)


def lp_compare(t_ms, x_lp, xref_lp, x0_ref, labels, tol, name):
    """x_lp, xref_lp: LP'd signals on the 1-ms grid (n x k); x0_ref: reference values at t = 0."""

    m = window_mask(t_ms)
    rows = []
    for j, lab in enumerate(labels):
        err = float(np.abs(x_lp[m, j] - xref_lp[m, j]).max())
        exc = float(np.abs(xref_lp[m, j] - x0_ref[j]).max())
        rows.append({"test": name, "quantity": lab, "max_abs_err": err, "excursion": exc, "tolerance": tol * exc + FLOOR,
                     "ratio": err / exc if exc > 0 else float("nan"), "pass": bool(err <= tol * exc + FLOOR)})
    return rows


def nyquist_blocks(a_step, dt):
    """a_step: per-step indicator (length n_steps + 1, entry 0 = 0). Returns block maxima per ms."""

    step = int(round(1e-3 / dt))
    a = np.asarray(a_step[1:])
    nb = len(a) // step
    return np.arange(1, nb + 1) * 1e-3, a[: nb * step].reshape(nb, step).max(axis=1)


def nyquist_event_check(tb, ab):
    late = ab[(tb > 4.5) & (tb <= 5.0 + 1e-12)].max()
    early = ab[(tb > 1.0) & (tb <= 1.5 + 1e-12)].max()
    m = (tb > 1.5) & (ab > 0)
    slope = float(np.polyfit(tb[m], np.log(ab[m]), 1)[0]) if m.sum() > 10 else float("nan")
    return {"A_max_4.5_5.0": float(late), "A_max_1.0_1.5": float(early), "log_slope_1.5_5.0": slope,
            "pass": bool(late <= 1e-5 and late <= early)}


def nyquist_noevent_check(tb, ab):
    tot = ab[tb > 0.1].max()
    late = ab[(tb > 4.5) & (tb <= 5.0 + 1e-12)].max()
    early = ab[(tb > 0.5) & (tb <= 1.0 + 1e-12)].max()
    return {"A_max_0.1_5.0": float(tot), "A_max_4.5_5.0": float(late), "A_max_0.5_1.0": float(early),
            "pass": bool(tot <= 1e-5 and late <= max(early, 1e-9))}


def estimator_r(t_ms, chans_lp):
    """Estimator R on LP'd channels sampled at 1 ms (n x k): decimate to 100 Hz by sampling."""

    stride = int(round(R_DT / 1e-3))
    tt = t_ms[::stride]
    x = chans_lp[::stride]
    m = (tt >= R_T0 - 1e-12) & (tt <= R_T1 + 1e-12)
    ys = []
    for j in range(x.shape[1]):
        y = x[m, j] - x[m, j].mean()
        r = float(np.sqrt(np.mean(y**2)))
        if r >= 1e-12:
            ys.append(y / r)
    if not ys:
        return {"resolved": False, "alpha": float("nan"), "freq": float("nan"), "resid": float("nan"), "order": 0, "modes": []}
    y = np.array(ys)
    n = y.shape[1]
    L = n // 3
    H = np.vstack([np.array([yc[i : i + L + 1] for i in range(n - L)]) for yc in y])
    _, sv, vh = np.linalg.svd(H, full_matrices=False)
    order = int(np.clip(np.sum(sv / sv[0] >= SV_REL), 2, MAX_ORDER))
    v = vh[:order].T
    z = np.linalg.eigvals(np.linalg.pinv(v[:-1]) @ v[1:])
    s = np.log(z.astype(complex)) / R_DT
    vand = z[None, :] ** np.arange(n)[:, None]
    res, *_ = np.linalg.lstsq(vand, y.T, rcond=None)
    resid = float(np.linalg.norm(y - (vand @ res).real.T) / np.linalg.norm(y))
    T = n * R_DT
    integ = np.where(np.abs(2 * s.real) > 1e-12, (np.exp(2 * s.real * T) - 1) / (2 * s.real), T)
    energy = np.sum(np.abs(res) ** 2, axis=1) * integ
    f = np.abs(s.imag) / (2 * np.pi)
    band = (f >= R_BAND[0]) & (f <= R_BAND[1]) & (s.imag > 0)
    if not band.any():
        return {"resolved": False, "alpha": float("nan"), "freq": float("nan"), "resid": resid, "order": order, "modes": []}
    tot = energy[band].sum()
    sel = band & (energy >= ENERGY_FRAC * tot)
    idx = np.flatnonzero(sel)
    best = idx[np.argmax(s[idx].real)]
    modes = sorted(((float(s[i].real), float(f[i]), float(energy[i] / tot)) for i in np.flatnonzero(band)), key=lambda q: -q[2])
    return {"resolved": bool(sel.any() and resid <= RESID_MAX), "alpha": float(s[best].real), "freq": float(f[best]), "resid": resid,
            "order": order, "modes": modes[:6]}


def ringdown_check(r_emt, r_ref):
    if not r_ref["resolved"]:
        return {"applicable": False, "pass": True}
    da, dfreq = abs(r_emt["alpha"] - r_ref["alpha"]), abs(r_emt["freq"] - r_ref["freq"])
    ok = bool(r_emt["resolved"] and da <= 0.005 + 0.02 * abs(r_ref["alpha"]) and dfreq <= 0.005 + 0.02 * r_ref["freq"])
    return {"applicable": True, "d_alpha": da, "d_f": dfreq, "pass": ok}
