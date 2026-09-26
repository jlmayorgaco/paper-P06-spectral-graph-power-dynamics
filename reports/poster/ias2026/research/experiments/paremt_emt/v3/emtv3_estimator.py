# ruff: noqa: E501  -- parameter tables kept on one line
"""V3 instrument: multi-channel bounded variable-projection damped-mode estimator ("VP").

Specification: docs/20260912_PAREMT_EMT_PREREG_V3.md sections 2-3 (commit 6d339006). New code; it
does not import the frozen v1 emt_estimator.py.

Model, channel c, centred time tau:  y_c = e^{a tau}[A_c cos 2 pi f tau + B_c sin 2 pi f tau] + C_c + D_c tau.
For fixed (a, f) the 4 linear coefficients of every channel come from one Householder QR of the
common column-scaled design matrix; only (a, f) are optimized, inside explicit bounds (grid, then
bounded trf refinement). No full-record z^k Vandermonde, no exp(2 a T) energy weighting.
Uncertainty: channel and time-block jackknife with floors; resolved flag R1-R4; CI verdict (EPS).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.optimize import least_squares

EPS = 0.002
BLOCK = 0.1
ALPHA_BOUNDS = (-0.8, 0.8)
GRID_DA, GRID_DF = 0.02, 0.005
N_STARTS = 3
K_TIME = 5
SE_FLOOR_A, SE_FLOOR_F = 5e-4, 5e-4
R3_MAX = 0.9
R4_A, R4_F = 0.02, 0.02
INTERIOR_A, INTERIOR_F = 0.01, 0.002
RANK_TOL = 1e-12
T0_SCI = 2.2
DURATIONS = (30.0, 60.0, 90.0)
NEIGH = 0.15
BAND = (0.2, 1.2)
SENT_EXPLAINED, SENT_DF = 0.10, 0.05
Z95 = 1.96


# ------------------------------------------------------------------ preprocessing --
def block_mean(t, y, t0, t1):
    """Non-overlapping 0.1-s block means inside [t0, t1); returns block centres and (n_ch, M) data."""

    t = np.asarray(t, float)
    y = np.atleast_2d(np.asarray(y, float))
    h = float(t[1] - t[0])
    if not np.allclose(np.diff(t), h, rtol=0, atol=1e-9 * max(1.0, abs(t[-1]))):
        raise ValueError("non-uniform sampling")
    spb = int(round(BLOCK / h))
    if spb < 1 or abs(spb * h - BLOCK) > 1e-9:
        raise ValueError("0.1 s is not an integer multiple of the sample step")
    i0 = int(np.searchsorted(t, t0 - 1e-9))
    if abs(t[i0] - t0) > 1e-9:
        raise ValueError("window start not on the sample grid")
    nb = int(np.floor((t1 - t0) / BLOCK + 1e-9))
    nb = min(nb, (len(t) - i0) // spb)
    seg = y[:, i0 : i0 + nb * spb].reshape(y.shape[0], nb, spb).mean(axis=2)
    tb = t0 + (np.arange(nb) + 0.5) * BLOCK
    return tb, seg


def normalize(tau, y):
    a = np.column_stack([np.ones_like(tau), tau])
    coef, *_ = np.linalg.lstsq(a, y.T, rcond=None)
    r = np.sqrt(np.mean((y - (a @ coef).T) ** 2, axis=1))
    keep = (r >= 1e-12) & (r >= 1e-8 * r.max()) if r.size else np.zeros(0, bool)
    return y[keep] / r[keep, None], keep, r


# ------------------------------------------------------------------ core VP ---------
# Amendment V3-AM1 (docs/20260912_PAREMT_EMT_PREREG_V3_AMENDMENTS.md): envelope-equalizing weight
# W(tau) = exp(-alpha_w tau), alpha_w = clip(first-pass alpha, -GAMMA_W/L, +GAMMA_W/L), L = window length.
GAMMA_W = 3.0


def _q(alpha, f, tau, wt=None):
    e = np.exp(alpha * tau)
    w = 2 * np.pi * f * tau
    b = np.column_stack([e * np.cos(w), e * np.sin(w), np.ones_like(tau), tau])
    if wt is not None:
        b = b * wt[:, None]
    b = b / np.linalg.norm(b, axis=0)
    q, r = np.linalg.qr(b)
    d = np.abs(np.diag(r))
    if not np.all(np.isfinite(d)) or d.min() < RANK_TOL * d.max():
        return None
    return q


def objective(alpha, f, tau, y, wt=None):
    q = _q(alpha, f, tau, wt)
    if q is None:
        return np.inf
    yw = y if wt is None else y * wt
    yq = yw @ q
    return float(np.sum(yw * yw) - np.sum(yq * yq))


def residual(alpha, f, tau, y, wt=None):
    """Weighted projected residual W (y - model) (unweighted when wt is None)."""

    q = _q(alpha, f, tau, wt)
    if q is None:
        return np.full(y.size, 1e6)
    yw = y if wt is None else y * wt
    return (yw - (yw @ q) @ q.T).ravel()


def _refine(x0, tau, y, lo, hi, max_nfev=None, wt=None):
    x0 = np.clip(np.asarray(x0, float), np.array(lo) + [1e-6, 1e-6], np.array(hi) - [1e-6, 1e-6])
    r = least_squares(lambda p: residual(p[0], p[1], tau, y, wt), x0, bounds=(lo, hi), method="trf", xtol=1e-12, ftol=1e-12,
                      gtol=1e-12, x_scale=np.array([GRID_DA, GRID_DF]), max_nfev=max_nfev)
    return r


def _grid_starts(tau, y, lo, hi, wt=None):
    ag = np.arange(lo[0], hi[0] + 1e-12, GRID_DA)
    fg = np.arange(lo[1], hi[1] + 1e-12, GRID_DF)
    J = np.array([[objective(a, f, tau, y, wt) for f in fg] for a in ag])
    starts = []
    order = np.argsort(J, axis=None)
    for flat in order:
        i, j = np.unravel_index(flat, J.shape)
        if not np.isfinite(J[i, j]):
            break
        nb = J[max(i - 1, 0) : i + 2, max(j - 1, 0) : j + 2]
        if J[i, j] <= nb.min() and all(abs(i - a) > 2 or abs(j - b) > 2 for a, b in starts):
            starts.append((i, j))
        if len(starts) == N_STARTS:
            break
    return [(ag[i], fg[j]) for i, j in starts]


def _best(tau, y, lo, hi, wt=None):
    best = None
    for x0 in _grid_starts(tau, y, lo, hi, wt):
        r = _refine(x0, tau, y, lo, hi, wt=wt)
        if np.all(np.isfinite(r.x)) and (best is None or r.cost < best.cost):
            best = r
    return best


def fit(tau, y, f_bounds, a_bounds=ALPHA_BOUNDS, jackknife=True, weighted=True):
    """VP fit on prepared data (tau centred, y normalized (n_ch, M)). Returns the result dictionary.
    weighted=True applies amendment V3-AM1 (first pass unweighted, final pass envelope-equalized)."""

    lo, hi = [a_bounds[0], f_bounds[0]], [a_bounds[1], f_bounds[1]]
    out = {"n_ch": int(y.shape[0]), "n_samples": int(y.shape[1]), "f_bounds": list(f_bounds), "a_bounds": list(a_bounds)}
    if y.shape[0] == 0:
        return {**out, "converged": False, "resolved": False, "reasons": ["no channel"], "alpha": np.nan, "freq": np.nan}
    best = _best(tau, y, lo, hi)
    if best is None:
        return {**out, "converged": False, "resolved": False, "reasons": ["no feasible start"], "alpha": np.nan, "freq": np.nan}
    wt = None
    out["alpha_first_pass"] = float(best.x[0])
    if weighted:
        span = float(tau[-1] - tau[0] + BLOCK)
        aw = float(np.clip(best.x[0], -GAMMA_W / span, GAMMA_W / span))
        wt = np.exp(-aw * tau)
        b2 = _best(tau, y, lo, hi, wt)
        if b2 is None:
            return {**out, "converged": False, "resolved": False, "reasons": ["no feasible weighted start"], "alpha": np.nan, "freq": np.nan}
        best = b2
        out["alpha_weight"] = aw
    a, f = float(best.x[0]), float(best.x[1])
    conv = bool(best.status > 0 and np.all(np.isfinite(best.fun)))
    interior = bool(a - lo[0] > INTERIOR_A and hi[0] - a > INTERIOR_A and f - lo[1] > INTERIOR_F and hi[1] - f > INTERIOR_F)
    m = y.shape[1]
    wv = np.ones(m) if wt is None else wt
    res_w = residual(a, f, tau, y, wt).reshape(y.shape)
    qa = np.linalg.qr(np.column_stack([wv, wv * tau]))[0]
    yw = y * wv
    ya_w = yw - (yw @ qa) @ qa.T
    per_ch = np.linalg.norm(res_w, axis=1) / np.maximum(np.linalg.norm(ya_w, axis=1), 1e-300)
    out.update({"alpha": a, "freq": f, "converged": conv, "interior": interior, "cost": float(best.cost),
                "resid_median": float(np.median(per_ch)), "resid_global": float(np.linalg.norm(res_w) / max(np.linalg.norm(ya_w), 1e-300))})
    out["_residual"] = res_w / wv
    se_ch = se_t = 0.0
    if jackknife and conv:
        th = []
        if y.shape[0] >= 3:
            for c in range(y.shape[0]):
                rr = _refine([a, f], tau, np.delete(y, c, axis=0), lo, hi, max_nfev=200, wt=wt)
                th.append(rr.x)
            th = np.array(th)
            n = len(th)
            se_ch = np.sqrt((n - 1) / n * np.sum((th - th.mean(axis=0)) ** 2, axis=0))
        tk = []
        edges = np.linspace(0, m, K_TIME + 1).astype(int)
        for k in range(K_TIME):
            keep = np.ones(m, bool)
            keep[edges[k] : edges[k + 1]] = False
            rr = _refine([a, f], tau[keep], y[:, keep], lo, hi, max_nfev=200, wt=None if wt is None else wt[keep])
            tk.append(rr.x)
        tk = np.array(tk)
        se_t = np.sqrt((K_TIME - 1) / K_TIME * np.sum((tk - tk.mean(axis=0)) ** 2, axis=0))
        se_ch = se_ch if np.ndim(se_ch) else np.zeros(2)
    else:
        se_ch, se_t = np.zeros(2), np.zeros(2)
    se_a = float(max(se_ch[0], se_t[0], SE_FLOOR_A))
    se_f = float(max(se_ch[1], se_t[1], SE_FLOOR_F))
    out.update({"se_alpha": se_a, "se_freq": se_f, "se_alpha_channel": float(se_ch[0]), "se_alpha_time": float(se_t[0]),
                "ci_alpha": [a - Z95 * se_a, a + Z95 * se_a], "ci_freq": [f - Z95 * se_f, f + Z95 * se_f]})
    reasons = []
    if not conv:
        reasons.append("R1 not converged")
    if not interior:
        reasons.append("R2 at bound")
    if out["resid_median"] > R3_MAX:
        reasons.append("R3 residual")
    if Z95 * se_a > R4_A or Z95 * se_f > R4_F:
        reasons.append("R4 CI width")
    out["reasons"] = reasons
    out["resolved"] = bool(conv and not reasons)
    out["valid_for_extension"] = bool(conv and interior and out["resid_median"] <= R3_MAX)
    return out


def prepare(t, y, t0, t1):
    tb, yb = block_mean(t, y, t0, t1)
    tau = tb - 0.5 * (t0 + t1)
    yn, keep, r = normalize(tau, yb)
    return tau, yn, keep


def verdict(res, flagged=False):
    if flagged:
        return "UNEXPECTED_MODE"
    if not res.get("resolved"):
        return "UNRESOLVED"
    lo, hi = res["ci_alpha"]
    if hi < -EPS:
        return "STABLE"
    if lo > EPS:
        return "UNSTABLE"
    return "UNRESOLVED"


def sentinel(tau, resid, f_target):
    s = fit(tau, resid, BAND, ALPHA_BOUNDS, jackknife=False, weighted=False)
    s.pop("_residual", None)
    if not (s.get("converged") and s.get("interior")):
        return {**s, "flag": False}
    r2 = residual(s["alpha"], s["freq"], tau, resid)
    expl = 1.0 - float(np.sum(r2**2) / max(np.sum(resid**2), 1e-300))
    s["explained"] = expl
    flag = False
    if s["alpha"] > EPS and expl >= SENT_EXPLAINED and abs(s["freq"] - f_target) > SENT_DF:
        sj = fit(tau, resid, BAND, ALPHA_BOUNDS, jackknife=True, weighted=False)
        lo_ci = sj["alpha"] - Z95 * sj["se_alpha"]
        s["ci_alpha_lower"] = lo_ci
        flag = bool(sj.get("converged") and sj.get("interior") and lo_ci > EPS)
    return {**s, "flag": flag}


def track(t, y, f_pred, t0=T0_SCI, t1=30.0, use_sentinel=True):
    fb = (max(BAND[0], f_pred - NEIGH), min(BAND[1], f_pred + NEIGH))
    tau, yn, keep = prepare(t, y, t0, t1)
    res = fit(tau, yn, fb)
    resid = res.pop("_residual", None)
    sent = sentinel(tau, resid, res["freq"]) if (use_sentinel and resid is not None and res.get("converged")) else {"flag": False}
    res["sentinel"] = {k: v for k, v in sent.items() if not k.startswith("_")}
    res["verdict"] = verdict(res, sent.get("flag", False))
    res["window"] = [t0, t1]
    res["channels_used"] = int(keep.sum())
    return res


def classify_adaptive(get_record: Callable[[float], tuple], f_pred, t0=T0_SCI, durations=DURATIONS):
    """get_record(T) -> (t, y) covering [0, T] (a shorter record is used as available). Adaptive rule
    of prereg V3 section 2.6. Returns the final result plus the history."""

    hist = []
    last_end = -np.inf
    for T in durations:
        t, y = get_record(T)
        t_end = float(min(T, t[-1] + (t[1] - t[0])))
        if t_end <= last_end + 1e-9:
            break
        last_end = t_end
        r = track(t, y, f_pred, t0, t_end)
        r["duration"] = t_end
        hist.append(r)
        extend = r["verdict"] == "UNRESOLVED" and r.get("valid_for_extension", False) and not r["sentinel"].get("flag", False)
        if not extend:
            break
    final = dict(hist[-1])
    final["history"] = [{k: h[k] for k in ("duration", "verdict", "alpha", "freq", "ci_alpha", "resolved", "reasons")} for h in hist]
    return final


def unit_ringdown(t, y, t0=1.5, t1=10.0, f_bounds=(0.5, 2.0)):
    """Unit-test settings (prereg V3 section 6): single or multi channel, time-block jackknife applies."""

    tau, yn, keep = prepare(t, np.atleast_2d(y), t0, t1)
    res = fit(tau, yn, f_bounds)
    res.pop("_residual", None)
    res["verdict"] = verdict(res)
    return res
