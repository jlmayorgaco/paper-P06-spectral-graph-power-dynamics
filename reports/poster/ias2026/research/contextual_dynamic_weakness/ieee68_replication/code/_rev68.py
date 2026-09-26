# ruff: noqa: E501
"""R7/R8 reversal machinery on a 64-portfolio census (docs/CDW68_PREREG_V1.md sections 8-9).

A census is {cond: {'alpha','status','hz','gap2','em_re','em_hz','modes'}} with arrays indexed by bitmask over
V68 = (3, 4, 6, 9, 11, 12) (bit b <-> unit V68[b]). Level definitions follow the hardening campaign exactly;
level T (EM-tracked) is the preregistered secondary.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import _r68 as R

V = R.V68
NB = len(V)
N = 1 << NB
RESOLVED = ("STABLE", "UNSTABLE")
TAUS = (0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05)


def compact(m):
    ph = np.asarray(m["phi"], np.float64)
    n = ph.size // 2
    return {"re": float(m["re"]), "hz": float(m["hz"]), "crit": bool(m.get("crit")), "em_top": bool(m.get("em_top")), "phi": ph[:n] + 1j * ph[n:]}


def mac_c(a, b) -> float:
    num = abs(np.vdot(a, b)) ** 2
    den = float(np.vdot(a, a).real * np.vdot(b, b).real) or 1e-300
    return float(num / den)


def match(mode, cands, band_only=False):
    best, bm = None, 0.0
    for c in cands:
        if band_only and not R.in_band(c["hz"]):
            continue
        if abs(c["hz"] - mode["hz"]) > R.DF_MAX:
            continue
        m = mac_c(mode["phi"], c["phi"])
        if m > bm:
            best, bm = c, m
    if bm < R.MAC_MIN:
        return None, bm
    return best, bm


def pick(ms, flag):
    if not ms:
        return None
    for m in ms:
        if m[flag]:
            return m
    return None


def census_from_records(records) -> dict:
    """records: iterable of (cond, S_label, record)."""
    cen = {}
    for cond, S, r in records:
        d = cen.setdefault(cond, {"alpha": np.full(N, np.nan), "status": np.array(["MISSING"] * N, dtype=object), "hz": np.full(N, np.nan),
                                  "gap2": np.full(N, np.nan), "em_re": np.full(N, np.nan), "em_hz": np.full(N, np.nan), "modes": [None] * N})
        m = R.mask_of(() if S == "BASE" else tuple(int(b) for b in S.split("+")))
        st = r.get("status")
        d["status"][m] = st
        a = r.get("alpha")
        d["alpha"][m] = float(a) if a is not None and st not in ("INFEASIBLE", "ERROR", "PF_FAIL") else np.nan
        d["hz"][m] = r.get("lam_hz", np.nan) if r.get("lam_hz") is not None else np.nan
        d["gap2"][m] = r.get("gap2", np.nan) if r.get("gap2") is not None else np.nan
        d["em_re"][m] = r.get("em_top_re", np.nan) if r.get("em_top_re") is not None else np.nan
        d["em_hz"][m] = r.get("em_top_hz", np.nan) if r.get("em_top_hz") is not None else np.nan
        d["modes"][m] = [compact(x) for x in r.get("modes", [])] if r.get("modes") else None
    return cen


def marginals(cond, d, tau=R.TAU_MAT) -> pd.DataFrame:
    rows = []
    alpha, status, hz, modes, em_re = d["alpha"], d["status"], d["hz"], d["modes"], d["em_re"]
    for m in range(N):
        cm = pick(modes[m], "crit")
        et = pick(modes[m], "em_top")
        for b in range(NB):
            if m >> b & 1:
                continue
            mi = m | (1 << b)
            delta = alpha[mi] - alpha[m]
            same, mac = False, np.nan
            if cm is not None and modes[mi]:
                mm, mac = match(cm, modes[mi])
                same = bool(mm is not None and mm["crit"])
            t_ok, t_delta, t_mac = False, np.nan, np.nan
            if et is not None and modes[mi]:
                mt, t_mac = match(et, modes[mi], band_only=True)
                if mt is not None and mt["em_top"]:
                    t_ok, t_delta = True, mt["re"] - et["re"]
            rows.append((cond, m, V[b], alpha[m], alpha[mi], delta, status[m], status[mi], same, mac, hz[m], hz[mi], t_ok, t_delta, t_mac, em_re[m], em_re[mi]))
    df = pd.DataFrame(rows, columns=["cond", "mask", "i", "alpha_S", "alpha_Si", "delta", "status_S", "status_Si", "same_mode", "mac", "hz_S", "hz_Si",
                                     "t_ok", "t_delta", "t_mac", "em_re_S", "em_re_Si"])
    df["em"] = df.hz_S.apply(R.in_band) & df.hz_Si.apply(R.in_band)
    df["lvA"] = (df.status_S == "STABLE") & df.status_Si.isin(RESOLVED) & np.isfinite(df.delta)
    df["lvB"] = df.lvA & df.same_mode
    df["lvC"] = df.lvB & df.em
    df["lvT"] = (df.status_S == "STABLE") & df.status_Si.isin(RESOLVED) & df.t_ok & np.isfinite(df.t_delta)
    df["cls"] = np.where(df.delta <= -tau, -1, np.where(df.delta >= tau, 1, 0))
    df["cls_T"] = np.where(df.t_delta <= -tau, -1, np.where(df.t_delta >= tau, 1, 0))
    df.loc[~np.isfinite(df.delta), "cls"] = 0
    df.loc[~np.isfinite(df.t_delta), "cls_T"] = 0
    return df


def _nested(lo, hi):
    if not len(lo) or not len(hi):
        return np.empty((0, 2), int)
    a = np.asarray(lo)[:, None]
    b = np.asarray(hi)[None, :]
    ok = ((a & b) == a) & (a != b)
    ia, ib = np.nonzero(ok)
    return np.stack([np.asarray(lo)[ia], np.asarray(hi)[ib]], axis=1)


def chain(f, S1, S2, bi):
    """Canonical chain S1 -> S2 (added units in increasing bus order): identity residual and anatomy (R8)."""
    add = [b for b in range(NB) if (S2 >> b & 1) and not (S1 >> b & 1)]
    Sr, terms = S1, []
    for b in add:
        i_, j_, ij = Sr | (1 << bi), Sr | (1 << b), Sr | (1 << bi) | (1 << b)
        terms.append(f[ij] - f[i_] - f[j_] + f[Sr])
        Sr = j_
    d1 = f[S1 | (1 << bi)] - f[S1]
    d2 = f[S2 | (1 << bi)] - f[S2]
    resid = d2 - (d1 + sum(terms))
    partial = d1
    first_cross = None
    for k, t in enumerate(terms, 1):
        new = partial + t
        if first_cross is None and np.sign(new) != np.sign(partial) and new != 0:
            first_cross = k
        partial = new
    return {"m": len(add), "resid": float(resid) if np.isfinite(resid) else np.nan, "cum": float(sum(terms)), "pos": float(sum(t for t in terms if t > 0)),
            "neg": float(sum(t for t in terms if t < 0)), "first_cross": first_cross,
            "wit": float(max(terms) if d2 > d1 else min(terms)) if terms else np.nan}


def pairs_for(cond, d, mg, tau=R.TAU_MAT):
    alpha, modes, gap2, em_re = d["alpha"], d["modes"], d["gap2"], d["em_re"]
    out = []
    for lvl in ("A", "B", "C", "T"):
        col, dcol, ccol = f"lv{lvl}", ("t_delta" if lvl == "T" else "delta"), ("cls_T" if lvl == "T" else "cls")
        cls = np.where(mg[dcol] <= -tau, -1, np.where(mg[dcol] >= tau, 1, 0))
        e = mg[mg[col].to_numpy() & np.isfinite(mg[dcol].to_numpy())].assign(_c=cls[mg[col].to_numpy() & np.isfinite(mg[dcol].to_numpy())])
        f = em_re if lvl == "T" else alpha
        for i in V:
            ei = e[e.i == i]
            dl = dict(zip(ei["mask"], ei[dcol], strict=True))
            bi = V.index(i)
            for direction, lo, hi in (("s2d", ei[ei._c == -1]["mask"].to_numpy(), ei[ei._c == 1]["mask"].to_numpy()),
                                      ("d2s", ei[ei._c == 1]["mask"].to_numpy(), ei[ei._c == -1]["mask"].to_numpy())):
                for S1, S2 in _nested(lo, hi):
                    S1, S2 = int(S1), int(S2)
                    ch = chain(f, S1, S2, bi)
                    four = (S1, S1 | (1 << bi), S2, S2 | (1 << bi))
                    rec = {"cond": cond, "level": lvl, "i": i, "S1": S1, "S2": S2, "dir": direction, "d1": dl[S1], "d2": dl[S2],
                           "mag": min(abs(dl[S1]), abs(dl[S2])), "size1": bin(S1).count("1"), "size2": bin(S2).count("1"),
                           "gap_min": float(np.nanmin([gap2[x] for x in four])), **{f"chain_{k}": v for k, v in ch.items()}}
                    if lvl in ("C", "T"):
                        flag = "crit" if lvl == "C" else "em_top"
                        c1, c2 = pick(modes[S1], flag), pick(modes[S2], flag)
                        ok = c1 is not None and c2 is not None and abs(c1["hz"] - c2["hz"]) <= R.DF_MAX
                        rec["mac12"] = mac_c(c1["phi"], c2["phi"]) if ok else np.nan
                        rec["same12"] = bool(ok and rec["mac12"] >= R.MAC_MIN)
                        rec["hz1"], rec["hz2"] = (c1["hz"] if c1 else np.nan), (c2["hz"] if c2 else np.nan)
                    out.append(rec)
    df = pd.DataFrame(out)
    if len(df):
        df["lvD"] = (df.level == "C") & df.get("same12", False).fillna(False).astype(bool)
        df["lvT2"] = (df.level == "T") & df.get("same12", False).fillna(False).astype(bool)
    return df


def hypergraph(d) -> tuple[str, float]:
    unst = [m for m in range(N) if d["status"][m] == "UNSTABLE"]
    minimal = [m for m in unst if not any((u & m) == u and u != m for u in unst)]
    if not minimal:
        return "EMPTY", np.nan
    lab = ";".join(sorted(R.label(R.members_of(m)) for m in minimal))
    return lab, float(min(bin(m).count("1") for m in minimal))


def analyse(cen: dict, tau=R.TAU_MAT):
    mgs, prs, rows = [], [], []
    for cond in sorted(cen):
        d = cen[cond]
        mg = marginals(cond, d, tau)
        pr = pairs_for(cond, d, mg, tau)
        mgs.append(mg)
        if len(pr):
            prs.append(pr)
        hg, kap = hypergraph(d)
        rec = {"cond": cond, "base_status": d["status"][0], "base_stable": d["status"][0] == "STABLE", "base_alpha": d["alpha"][0],
               "base_hz": d["hz"][0], "base_em_re": d["em_re"][0], "n_stable": int((d["status"] == "STABLE").sum()),
               "n_unstable": int((d["status"] == "UNSTABLE").sum()), "n_unresolved": int(np.isin(d["status"], ["BOUNDARY_OR_UNRESOLVED"]).sum()),
               "n_infeasible": int(np.isin(d["status"], ["INFEASIBLE", "ERROR", "PF_FAIL", "LIMIT_ACTIVE", "MISSING"]).sum()),
               "frac_crit_em": float(np.mean([R.in_band(h) for h, s in zip(d["hz"], d["status"], strict=True) if s == "STABLE"])) if (d["status"] == "STABLE").any() else np.nan,
               "H_RHP": hg, "kappa_RHP": kap}
        for lvl, sel in (("A", lambda p: p.level == "A"), ("B", lambda p: p.level == "B"), ("C", lambda p: p.level == "C"),
                         ("D", lambda p: p.lvD), ("T", lambda p: p.lvT2)):
            q = pr[sel(pr)] if len(pr) else pr
            rec[f"{lvl}_has"] = bool(len(q) > 0)
            rec[f"{lvl}_n_pairs"] = int(len(q))
            rec[f"{lvl}_n_s2d"] = int((q.dir == "s2d").sum()) if len(q) else 0
            rec[f"{lvl}_n_d2s"] = int((q.dir == "d2s").sum()) if len(q) else 0
            rec[f"{lvl}_both_dirs"] = bool(rec[f"{lvl}_n_s2d"] > 0 and rec[f"{lvl}_n_d2s"] > 0)
            rec[f"{lvl}_max_mag"] = float(q.mag.max()) if len(q) else np.nan
            rec[f"{lvl}_med_mag"] = float(q.mag.median()) if len(q) else np.nan
        rows.append(rec)
    mg = pd.concat(mgs, ignore_index=True) if mgs else pd.DataFrame()
    pr = pd.concat(prs, ignore_index=True) if prs else pd.DataFrame()
    return mg, pr, pd.DataFrame(rows)


def coverage_curve(cen, conds, taus=TAUS):
    out = []
    for tau in taus:
        _, pr, pol = analyse({c: cen[c] for c in conds}, tau)
        el = pol[pol.base_stable]
        out.append({"tau": tau, "n": int(len(el)), "C": int(el.C_has.sum()), "D": int(el.D_has.sum()), "T": int(el.T_has.sum()),
                     "D_both": int(el.D_both_dirs.sum()), "T_both": int(el.T_both_dirs.sum())})
    return pd.DataFrame(out)
