# ruff: noqa: E501, E402  -- long argument lists; BLAS threads must be pinned before numpy is imported
"""Shared helpers for the TX4 ParaEMT campaign scripts (run with .venv/xtool-paremt).

Puts the working copy external/ParaEMT_tx4 on sys.path (after EMT_setup_tx4.py has refreshed
it), and provides the run manifest, raw-output saving and the estimator-signal extraction.
"""

from __future__ import annotations

import os

# single-threaded BLAS: the dense solve / inverse is then bit-reproducible (checked on EMT01)
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import hashlib  # noqa: E402
import json  # noqa: E402
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
REPO = HERE.parents[5]
TX4 = REPO / "external" / "ParaEMT_tx4"
UP = REPO / "external" / "ParaEMT_upstream"
RAW = REPO / "external" / "paremt_runs" / "raw"
RESULTS = RESEARCH / "results"
if str(TX4) not in sys.path:
    sys.path.insert(0, str(TX4))
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import tx4_case  # noqa: E402
import tx4_emt  # noqa: E402

W0 = tx4_emt.W0
GEN_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
CASES = json.loads((RESULTS / "EMT_PRED" / "emt_cases.json").read_text())["cases"]
CASE = {c["case_id"]: c for c in CASES}


def _git(path, *args):
    return subprocess.run(
        ["git", "-C", str(path), *args], capture_output=True, text=True
    ).stdout.strip()


def provenance() -> dict:
    diff = _git(TX4, "diff", "HEAD")
    untracked = sorted(
        (TX4 / f).read_bytes()
        for f in _git(TX4, "ls-files", "--others", "--exclude-standard").split()
        if f
    )
    h = hashlib.sha256(diff.encode())
    for b in untracked:
        h.update(b)
    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True
    ).stdout
    return {
        "git_head": _git(RESEARCH, "rev-parse", "HEAD"),
        "paremt_upstream_sha": _git(UP, "rev-parse", "HEAD"),
        "paremt_working_copy_diff_sha": h.hexdigest(),
        "environment_hash": hashlib.sha256(freeze.encode()).hexdigest(),
    }


def run_spec(
    case: dict,
    dt: float = 50e-6,
    tlen: float = 30.0,
    amplitude: float = 0.02,
    tau_m: float = 1e-3,
    store_ms: float = 1e-3,
    save_raw: bool = True,
    tag: str = "",
):
    data = tx4_case.Data(RESEARCH)
    built = tx4_case.build_case(data, case, dt, amplitude=amplitude, tau_m=tau_m)
    ds = int(round(store_ms / dt))
    t0 = time.time()
    t, sg, gf, v, finite = tx4_case.run(built, dt, tlen, ds)
    wall = time.time() - t0
    out = {
        "t": t,
        "sg": sg,
        "gf": gf,
        "v": v[..., 0] + 1j * v[..., 1],
        "finite": bool(finite),
        "labels_sg": built["labels_sg"],
        "labels_gf": built["labels_gf"],
        "order": built["order"],
    }
    run_id = f"{case['case_id']}{tag}_dt{int(round(dt * 1e6))}us_a{amplitude:g}"
    manifest = {
        "run_id": run_id,
        **provenance(),
        **{f"{k}_hash": v for k, v in built["hashes"].items()},
        "dt_s": dt,
        "duration_s": tlen,
        "disturbance": {
            "bus": tx4_case.PULSE["bus"],
            "fraction": amplitude,
            "t0": tx4_case.PULSE["t0"],
            "t1": tx4_case.PULSE["t1"],
        },
        "tau_m_s": tau_m,
        "stored_interval_s": store_ms,
        "seed": None,
        "wall_s": round(wall, 2),
        "finite": bool(finite),
    }
    if save_raw:
        RAW.mkdir(parents=True, exist_ok=True)
        path = RAW / f"{run_id}.npz"
        np.savez_compressed(
            path,
            t=t,
            sg=sg.astype(np.float32),
            gf=gf.astype(np.float32),
            v=out["v"].astype(np.complex64),
            labels_sg=np.array(built["labels_sg"]),
            labels_gf=np.array(built["labels_gf"]),
        )
        manifest["output_path"] = str(path.relative_to(REPO)).replace("\\", "/")
        manifest["output_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    out["manifest"] = manifest
    return out


def estimator_channels(out: dict) -> tuple[np.ndarray, list[str]]:
    """Relative signals (prereg section 5): rel speeds, GFL freq deviations, angle diffs, |V|."""

    sg, gf, v, order = out["sg"], out["gf"], out["v"], out["order"]
    labels_sg, labels_gf = list(out["labels_sg"]), list(out["labels_gf"])
    i39 = labels_sg.index(39)
    w39 = sg[:, i39, 1]
    chans, names = [], []
    for j, b in enumerate(labels_sg):
        if b != 39:
            chans.append(sg[:, j, 1] - w39)
            names.append(f"relspeed_{b}")
    for j, b in enumerate(labels_gf):
        chans.append(gf[:, j, 13] / W0 - (w39 - 1.0))
        names.append(f"gflfreq_{b}")
    ib39 = order.index(39)
    for b in GEN_BUSES:
        if b != 39:
            chans.append(
                np.unwrap(np.angle(v[:, order.index(b)]) - np.angle(v[:, ib39]))
            )
            names.append(f"angdiff_{b}")
    for b in GEN_BUSES:
        chans.append(np.abs(v[:, order.index(b)]))
        names.append(f"vmag_{b}")
    return np.array(chans), names
