# ruff: noqa: E501
"""CDW hardening infrastructure: paths, policies, draws, status. Reuses experiments/cdw/_infra.

Prereg: docs/CDW_HARDENING_PREREG_V1.md (commit 05b507e3). Writes only under
contextual_dynamic_weakness/{raw/H_*, results/hardening, logs/hardening, figures/hardening}.
"""

from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
CDW = HERE.parent / "cdw"
for p in (str(HERE), str(CDW)):
    if p not in sys.path:
        sys.path.insert(0, p)

import _infra as I  # noqa: E402  (pins BLAS before numpy)

PREREG_COMMIT = "05b507e3"
RESULTS = I.RESULTS / "hardening"
LOGS = I.LOGS / "hardening"
FIGS = I.FIGURES / "hardening"
INPUTS = RESULTS / "prereg_inputs"
STATUS = RESULTS / "CDWH_RUN_STATUS.json"
for p in (RESULTS, LOGS, FIGS):
    p.mkdir(parents=True, exist_ok=True)

HPOL = [f"HARDENING_H{i:02d}" for i in range(1, 25)]
OLD_HOLD = [f"H{i:02d}" for i in range(1, 25)]
DISC = [f"D{i:02d}" for i in range(1, 16)]
ENVS = ("EM-f", "EM-u", "EC", "EMC")
OMEGA_REF_HZ = 0.6222796695779355
SEED_BOOT = 20260932
SEED_PERM = 20260933
EM_BAND = (0.1, 2.0)


@lru_cache(maxsize=1)
def hardening_policies() -> dict:
    rows = json.loads((INPUTS / "hardening_policies.json").read_text(encoding="utf-8"))
    return {r["id"]: (r["g"], r["k"], r["t"], r["h"]) for r in rows}


@lru_cache(maxsize=1)
def fresh_draws() -> list:
    return json.loads((INPUTS / "hardening_envelope_draws.json").read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def null_groups() -> dict:
    return json.loads((INPUTS / "corridor_null_groups.json").read_text(encoding="utf-8"))


def theta_of(pid: str) -> tuple:
    import _cdw as C

    if pid.startswith("HARDENING_"):
        return hardening_policies()[pid]
    return C.policy(pid)


def split_of(pid: str) -> str:
    if pid.startswith("HARDENING_"):
        return "new"
    if pid.startswith("H"):
        return "old"
    return "discovery"


def draw_of(source: str | None, env: str | None, k: int | None):
    """source: None/'NOMINAL', 'CDW' (old holdout draws) or 'FRESH' (hardening draws)."""

    import _cdw as C

    if source in (None, "NOMINAL") or env is None:
        return None
    pool = C.cdw_draws() if source == "CDW" else fresh_draws()
    for d in pool:
        if d["envelope"] == env and d["draw"] == k:
            return d
    raise KeyError((source, env, k))


def in_band(hz) -> bool:
    return hz is not None and EM_BAND[0] <= hz <= EM_BAND[1]


# ------------------------------------------------------------------ status --
def load_status() -> dict:
    if STATUS.exists():
        return json.loads(STATUS.read_text(encoding="utf-8"))
    return {"phases": {}}


def set_status(phase: str, **fields) -> None:
    import time

    st = load_status()
    rec = st["phases"].get(phase, {})
    rec.update(fields)
    rec["updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    st["phases"][phase] = rec
    I.atomic_write_json(STATUS, st)


def write_json(name: str, payload) -> None:
    I.atomic_write_json(RESULTS / name, payload)
