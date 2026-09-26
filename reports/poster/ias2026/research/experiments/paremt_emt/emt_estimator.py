# ruff: noqa: E501  -- parameter tables and kernel argument lists kept on one line
"""Preregistered EMT mode estimators (prereg section 5; frozen at commit c2947bd8).

Input: relative channels sampled every 1 ms (from _emt.estimator_channels) and their time axis.
Window [2.2 s, T_end]; each channel detrended (least-squares line) and scaled to unit RMS;
zero-phase FIR anti-alias decimation to 50 Hz (estimator B) and 10 Hz (estimator A).

A  multi-channel matrix pencil (10 Hz): stacked Hankel, L = floor(N/3), order = #(sigma_i/sigma_1 >= 1e-4)
   clipped to [2, 30]; modes s = ln(z)/dt; residues by least squares; energy = sum_ch |res|^2 * int|e^{st}|^2;
   in-band (0.2-1.2 Hz) modes with energy >= 1 % of the in-band total; alpha_A = max Re, f_A its frequency;
   resolved iff >= 1 selected mode and relative reconstruction residual <= 0.05.
B  analytic-signal envelope (50 Hz): channel with the largest in-band Welch power; f_pk from the Hann-windowed
   zero-padded FFT in band; 4th-order Butterworth band-pass [f_pk-0.15, f_pk+0.15] n [0.2, 1.2], zero phase;
   alpha_B = slope of ln|analytic| on the window trimmed 3 s at each end, f_B = phase slope / 2 pi;
   resolved iff R^2 >= 0.90 and >= 5 cycles.
Classification (eps = 0.002): STABLE iff both resolved and both < -eps; UNSTABLE iff both > +eps; else UNRESOLVED.
Agreement flag: |alpha_A - alpha_B| <= 0.02 and |f_A - f_B| <= 0.03.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

BAND = (0.2, 1.2)
EPS = 0.002
T_START = 2.2
SV_REL = 1e-4
MAX_ORDER = 30
ENERGY_FRAC = 0.01
RESID_MAX = 0.05
R2_MIN = 0.90
MIN_CYCLES = 5
TRIM = 3.0


def _prep(t, chans):
    m = t >= T_START - 1e-12
    tt = t[m]
    x = chans[:, m].astype(np.float64)
    a = np.vstack([np.ones_like(tt), tt - tt[0]]).T
    coef, *_ = np.linalg.lstsq(a, x.T, rcond=None)
    x = x - (a @ coef).T
    rms = np.sqrt(np.mean(x**2, axis=1))
    keep = rms > 1e-14
    x = x[keep] / rms[keep, None]
    return tt, x, keep


def _decimate(x, factor):
    out = x
    for q in _factors(factor):
        out = signal.decimate(out, q, ftype="fir", axis=-1, zero_phase=True)
    return out


def _factors(n):
    fs = []
    for q in (5, 4, 2):
        while n % q == 0 and n > 1:
            fs.append(q)
            n //= q
    assert n == 1
    return fs


def matrix_pencil(t, chans, dt_in):
    tt, x, _ = _prep(t, chans)
    fac = int(round(0.1 / dt_in))
    y = _decimate(x, fac)
    dt = dt_in * fac
    n = y.shape[1]
    L = n // 3
    blocks = [np.array([yc[i : i + L + 1] for i in range(n - L)]) for yc in y]
    H = np.vstack(blocks)
    _, sv, vh = np.linalg.svd(H, full_matrices=False)
    order = int(np.clip(np.sum(sv / sv[0] >= SV_REL), 2, MAX_ORDER))
    v = vh[:order].T
    v1, v2 = v[:-1], v[1:]
    z = np.linalg.eigvals(np.linalg.pinv(v1) @ v2)
    s = np.log(z.astype(complex)) / dt
    k = np.arange(n)
    vand = z[None, :] ** k[:, None]
    res, *_ = np.linalg.lstsq(vand, y.T, rcond=None)
    recon = (vand @ res).real.T
    resid = float(np.linalg.norm(y - recon) / np.linalg.norm(y))
    T = n * dt
    integ = np.where(
        np.abs(2 * s.real) > 1e-12, (np.exp(2 * s.real * T) - 1) / (2 * s.real), T
    )
    energy = np.sum(np.abs(res) ** 2, axis=1) * integ
    f = np.abs(s.imag) / (2 * np.pi)
    inband = (f >= BAND[0]) & (f <= BAND[1]) & (s.imag >= 0)
    if not inband.any():
        return {
            "resolved": False,
            "alpha": float("nan"),
            "freq": float("nan"),
            "order": order,
            "resid": resid,
            "modes": [],
        }
    tot = energy[inband].sum()
    sel = inband & (energy >= ENERGY_FRAC * tot)
    idx = np.flatnonzero(sel)
    best = idx[np.argmax(s[idx].real)]
    modes = sorted(
        (
            (float(s[i].real), float(f[i]), float(energy[i] / tot))
            for i in np.flatnonzero(inband)
        ),
        key=lambda m: -m[2],
    )
    return {
        "resolved": bool(sel.any() and resid <= RESID_MAX),
        "alpha": float(s[best].real),
        "freq": float(f[best]),
        "order": order,
        "resid": resid,
        "modes": modes[:8],
    }


def hilbert_envelope(t, chans, dt_in):
    tt, x, _ = _prep(t, chans)
    fac = int(round(0.02 / dt_in))
    y = _decimate(x, fac)
    dt = dt_in * fac
    fs = 1.0 / dt
    fw, pw = signal.welch(y, fs=fs, nperseg=min(y.shape[1], 1024), axis=-1)
    band = (fw >= BAND[0]) & (fw <= BAND[1])
    ch = int(np.argmax(pw[:, band].sum(axis=1)))
    yc = y[ch]
    nfft = 16 * int(2 ** np.ceil(np.log2(yc.size)))
    spec = np.abs(np.fft.rfft(yc * np.hanning(yc.size), n=nfft))
    ff = np.fft.rfftfreq(nfft, dt)
    bm = (ff >= BAND[0]) & (ff <= BAND[1])
    fpk = float(ff[bm][np.argmax(spec[bm])])
    lo, hi = max(BAND[0], fpk - 0.15), min(BAND[1], fpk + 0.15)
    sos = signal.butter(4, [lo, hi], btype="bandpass", fs=fs, output="sos")
    yb = signal.sosfiltfilt(sos, yc)
    an = signal.hilbert(yb)
    tloc = np.arange(yc.size) * dt
    m = (tloc >= TRIM) & (tloc <= tloc[-1] - TRIM)
    le = np.log(np.abs(an[m]) + 1e-300)
    ph = np.unwrap(np.angle(an[m]))
    a = np.vstack([np.ones(m.sum()), tloc[m]]).T
    cl, *_ = np.linalg.lstsq(a, le, rcond=None)
    fit = a @ cl
    r2 = 1.0 - np.sum((le - fit) ** 2) / max(np.sum((le - le.mean()) ** 2), 1e-300)
    cp, *_ = np.linalg.lstsq(a, ph, rcond=None)
    fb = float(cp[1] / (2 * np.pi))
    span = tloc[m][-1] - tloc[m][0]
    cycles = span * abs(fb)
    return {
        "resolved": bool(r2 >= R2_MIN and cycles >= MIN_CYCLES),
        "alpha": float(cl[1]),
        "freq": abs(fb),
        "r2": float(r2),
        "cycles": float(cycles),
        "channel": ch,
        "f_peak": fpk,
    }


def classify(t, chans, dt_in):
    a = matrix_pencil(t, chans, dt_in)
    b = hilbert_envelope(t, chans, dt_in)
    both = a["resolved"] and b["resolved"]
    if both and a["alpha"] < -EPS and b["alpha"] < -EPS:
        verdict = "STABLE"
    elif both and a["alpha"] > EPS and b["alpha"] > EPS:
        verdict = "UNSTABLE"
    else:
        verdict = "UNRESOLVED"
    agree = bool(
        abs(a["alpha"] - b["alpha"]) <= 0.02 and abs(a["freq"] - b["freq"]) <= 0.03
    )
    return {
        "verdict": verdict,
        "alpha_A": a["alpha"],
        "f_A": a["freq"],
        "resolved_A": a["resolved"],
        "order_A": a["order"],
        "resid_A": a["resid"],
        "alpha_B": b["alpha"],
        "f_B": b["freq"],
        "resolved_B": b["resolved"],
        "r2_B": b["r2"],
        "cycles_B": b["cycles"],
        "agree": agree,
        "modes_A": a["modes"],
    }
