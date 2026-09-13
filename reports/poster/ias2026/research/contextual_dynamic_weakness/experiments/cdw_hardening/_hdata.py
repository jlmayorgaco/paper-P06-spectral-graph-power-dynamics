# ruff: noqa: E501
"""Census loaders and per-marginal features for the hardening analyses (no model evaluation).

A census is the 512 V9 portfolios at one policy, indexed by bitmask over V9 = (30..38) in
order (bit b <-> unit V9[b]). Old census: raw/E01 (D01-D15, H01-H24); new: raw/H_H01.
"""

from __future__ import annotations

import _hinfra as HI  # noqa: I001  (first: sets sys.path for the CDW modules)

import numpy as np
import pandas as pd

import _cdw as C
import _infra as I
import E01_census as E1

V9 = C.V9
BIT = {u: b for b, u in enumerate(V9)}
EMLO, EMHI = HI.EM_BAND


def mask_of(label: str) -> int:
    if label == "BASE":
        return 0
    m = 0
    for u in label.split("+"):
        m |= 1 << BIT[int(u)]
    return m


def label_of(mask: int) -> str:
    return C.label(tuple(V9[b] for b in range(9) if mask >> b & 1))


def load_census(store: str, device="gfl") -> dict:
    """{pid: {'alpha','status','hz','modes'}} arrays of length 512 (bitmask index)."""

    out = {}
    for rec in I.Store(store).all():
        if not rec.get("ok"):
            raise RuntimeError(f"failed task in {store}: {rec.get('error')}")
        t = rec["task"]
        if t.get("device", "gfl") != device:
            continue
        d = out.setdefault(t["pid"], {"alpha": np.full(512, np.nan), "status": np.array(["MISSING"] * 512, dtype=object),
                                      "hz": np.full(512, np.nan), "modes": [None] * 512})
        for r in rec["records"]:
            m = mask_of(r["S"])
            d["alpha"][m] = r.get("alpha", np.nan) if r.get("alpha") is not None else np.nan
            d["status"][m] = r.get("status")
            d["hz"][m] = r.get("lam_hz", np.nan) if r.get("lam_hz") is not None else np.nan
            d["modes"][m] = [E1.compact(x) for x in r.get("modes", [])]
    return out


def load_all_census() -> dict:
    cen = load_census("E01")
    cen.update(load_census("H_H01"))
    return cen


def crit(modes):
    return E1.crit_mode(modes) if modes else None


def in_em(hz) -> bool:
    return bool(np.isfinite(hz) and EMLO <= hz <= EMHI)


def marginals(pid: str, d: dict) -> pd.DataFrame:
    rows = []
    alpha, status, hz, modes = d["alpha"], d["status"], d["hz"], d["modes"]
    for m in range(512):
        cm = crit(modes[m])
        for b in range(9):
            if m >> b & 1:
                continue
            mi = m | (1 << b)
            delta = alpha[mi] - alpha[m]
            same, mac = False, np.nan
            if cm is not None and modes[mi]:
                mm, mac = E1.match(cm, modes[mi])
                same = bool(mm is not None and mm["crit"])
            rows.append((pid, m, V9[b], alpha[m], alpha[mi], delta, status[m], status[mi], status[m] == "STABLE",
                         same, mac, hz[m], hz[mi]))
    df = pd.DataFrame(rows, columns=["pid", "mask", "i", "alpha_S", "alpha_Si", "delta", "status_S", "status_Si", "stable_S",
                                     "same_mode", "mac", "hz_S", "hz_Si"])
    df["em"] = df.hz_S.between(EMLO, EMHI) & df.hz_Si.between(EMLO, EMHI)
    df["cls"] = np.where(df.delta <= -C.TAU_MAT, -1, np.where(df.delta >= C.TAU_MAT, 1, 0))
    df.loc[~np.isfinite(df.delta), "cls"] = 0
    df["lvA"] = df.stable_S & np.isfinite(df.delta)
    df["lvB"] = df.lvA & df.same_mode
    df["lvC"] = df.lvB & df.em
    return df


def all_marginals(cen: dict) -> pd.DataFrame:
    return pd.concat([marginals(pid, cen[pid]) for pid in sorted(cen)], ignore_index=True)
