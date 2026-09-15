"""Blind reduced-port prediction for the frozen TX4 V4/V9 candidate sets.

The module intentionally contains no import of historical full-order result
tables. It builds one all-SG baseline port kernel and one single-replacement
local Delta-Y model per candidate, then predicts portfolio closure roots.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import root

HERE = Path(__file__).resolve()
EXP = HERE.parent
SRC = EXP.parent / "src"
ROOT = HERE.parents[5]
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(SRC))

from _f7_common import LEAK, Theta  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.port_admittance import build_port_operator  # noqa: E402

OUT = ROOT / "results"
BAND_HZ = (0.3, 1.5)
THETA = Theta(g=0.03625, k=1.425, t=1.5, h=1.0)
V4 = (30, 33, 35, 37)
V9 = tuple(range(30, 39))
SIGMA_GRID = np.array([0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4])
FREQ_GRID = np.linspace(BAND_HZ[0], BAND_HZ[1], 25)
FREQ_CLOSURE_GRID = np.linspace(BAND_HZ[0], BAND_HZ[1], 121)
LOCAL_REGULARITY_FLOOR = 1e-7


def tx4_case(members: tuple[int, ...]):
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        converter=ConverterParameters(
            voltage_control=True, voltage_gain=THETA.g, voltage_leak=LEAK
        ),
        machine_scaling={"ka": THETA.k, "ta": THETA.t},
    )


@dataclass
class Kernel:
    candidates: tuple[int, ...]
    base_case: object
    local_cases: dict[int, object]
    base_port: object
    local_ports: dict[int, object]

    def __post_init__(self):
        self._cache: dict[tuple[float, float], tuple[np.ndarray, dict[int, np.ndarray]]] = {}
        self._positions = {b: self.base_port.bus_index[b] for b in self.candidates}
        selector = np.zeros((self.base_port.dimension, 2 * len(self.candidates)), dtype=complex)
        for k, b in enumerate(self.candidates):
            p = self._positions[b]
            selector[2 * p : 2 * p + 2, 2 * k : 2 * k + 2] = np.eye(2)
        self._selector = selector

    def at(self, s: complex):
        key = (round(float(s.real), 12), round(float(s.imag), 12))
        if key in self._cache:
            return self._cache[key]
        t0 = self.base_port.evaluate(s)
        kfull = self._selector.T @ np.linalg.solve(t0, self._selector)
        deltas = {
            b: self.base_port.bus_admittance(s, b) - self.local_ports[b].bus_admittance(s, b)
            for b in self.candidates
        }
        self._cache[key] = (kfull, deltas)
        return self._cache[key]

    def qdata(self, members: tuple[int, ...], s: complex) -> dict:
        if not members:
            return {
                "collective_det": 1.0 + 0j,
                "collective_sigma_min": 1.0,
                "local_sigma_min": 1.0,
                "local_det_min_abs": 1.0,
                "q": np.zeros((0, 0), dtype=complex),
            }
        kfull, deltas = self.at(s)
        pos = [self.candidates.index(b) for b in members]
        idx = sum(([2 * p, 2 * p + 1] for p in pos), [])
        kss = kfull[np.ix_(idx, idx)]
        d = np.zeros((2 * len(members), 2 * len(members)), dtype=complex)
        for j, b in enumerate(members):
            d[2 * j : 2 * j + 2, 2 * j : 2 * j + 2] = deltas[b]
        m = d @ kss
        total = np.eye(2 * len(members), dtype=complex) + m
        local = np.zeros_like(total)
        local_smin, local_det = [], []
        for j in range(len(members)):
            block = total[2 * j : 2 * j + 2, 2 * j : 2 * j + 2]
            local[2 * j : 2 * j + 2, 2 * j : 2 * j + 2] = block
            local_smin.append(float(np.linalg.svd(block, compute_uv=False)[-1]))
            local_det.append(complex(np.linalg.det(block)))
        q = np.linalg.solve(local, total) - np.eye(2 * len(members), dtype=complex)
        c = np.eye(2 * len(members), dtype=complex) + q
        return {
            "collective_det": complex(np.linalg.det(c)),
            "collective_sigma_min": float(np.linalg.svd(c, compute_uv=False)[-1]),
            "local_sigma_min": min(local_smin),
            "local_det_min_abs": min(abs(x) for x in local_det),
            "q": q,
        }


def make_kernel(candidates: tuple[int, ...]) -> Kernel:
    base = tx4_case(())
    locals_ = {b: tx4_case((b,)) for b in candidates}
    bp = build_port_operator(base.dae, base.equilibrium.x, base.equilibrium.z)
    lp = {b: build_port_operator(c.dae, c.equilibrium.x, c.equilibrium.z) for b, c in locals_.items()}
    return Kernel(candidates, base, locals_, bp, lp)


def root_for(kernel: Kernel, members: tuple[int, ...]):
    if not members:
        return None

    def f(x):
        s = complex(float(x[0]), 2.0 * math.pi * float(x[1]))
        z = kernel.qdata(members, s)["collective_det"]
        scale = max(1.0, abs(z))
        return [z.real / scale, z.imag / scale]

    seeds = []
    for sigma in SIGMA_GRID:
        for freq in FREQ_GRID:
            z = kernel.qdata(members, complex(float(sigma), 2 * math.pi * float(freq)))
            seeds.append((z["collective_sigma_min"], sigma, freq))
    seeds.sort(key=lambda x: x[0])
    roots = []
    for _, sigma, freq in seeds[:8]:
        sol = root(f, [sigma, freq], method="hybr", options={"xtol": 1e-9})
        if not sol.success or not np.all(np.isfinite(sol.x)):
            continue
        sig, fr = float(sol.x[0]), float(sol.x[1])
        if sig < -1e-6 or fr < BAND_HZ[0] - 1e-6 or fr > BAND_HZ[1] + 1e-6:
            continue
        z = kernel.qdata(members, complex(sig, 2 * math.pi * fr))
        residual = abs(z["collective_det"])
        if residual > 1e-6 or z["local_sigma_min"] < LOCAL_REGULARITY_FLOOR:
            continue
        if not any(abs(sig - a) < 1e-5 and abs(fr - b) < 1e-5 for a, b in roots):
            roots.append((sig, fr))
    if not roots:
        return None
    # The least-damped RHP/boundary root is the relevant portfolio root.
    return min(roots, key=lambda x: x[0])


def predict_one(kernel: Kernel, members: tuple[int, ...]) -> dict:
    started = time.perf_counter()
    root_hit = root_for(kernel, members)
    closure = []
    for freq in FREQ_CLOSURE_GRID:
        s = complex(0.0, 2 * math.pi * float(freq))
        closure.append(kernel.qdata(members, s))
    cmin = min(closure, key=lambda x: x["collective_sigma_min"])
    if root_hit is None:
        verdict, sigma, freq = "STABLE_PREDICTED", float("nan"), float("nan")
    else:
        sigma, freq = root_hit
        verdict = "UNSTABLE_PREDICTED" if sigma > 1e-6 else "BOUNDARY_PREDICTED"
    return {
        "portfolio": "+".join(map(str, members)) or "BASE",
        "members": "+".join(map(str, members)) or "BASE",
        "cardinality": len(members),
        "predicted_verdict": verdict,
        "predicted_root_real": sigma,
        "predicted_root_frequency_hz": freq,
        "predicted_closure_min_sigma": cmin["collective_sigma_min"],
        "predicted_closure_frequency_hz": FREQ_CLOSURE_GRID[int(np.argmin([x["collective_sigma_min"] for x in closure]))],
        "predicted_local_factor_min_sigma": cmin["local_sigma_min"],
        "predicted_local_factor_min_det_abs": cmin["local_det_min_abs"],
        "reduced_runtime_s": time.perf_counter() - started,
    }


def all_subsets(candidates):
    return [tuple(s) for n in range(len(candidates) + 1) for s in combinations(candidates, n)]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def run_set(label: str, candidates: tuple[int, ...]):
    started = time.perf_counter()
    kernel = make_kernel(candidates)
    setup = time.perf_counter() - started
    rows = [predict_one(kernel, s) for s in all_subsets(candidates)]
    path = OUT / f"TX4_{label}_BLIND_PREDICTIONS.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return {
        "label": label,
        "candidates": list(candidates),
        "n_portfolios": len(rows),
        "setup_s": setup,
        "total_reduced_s": time.perf_counter() - started,
        "median_per_portfolio_s": float(np.median([r["reduced_runtime_s"] for r in rows])),
        "prediction_file": str(path.relative_to(ROOT)),
    }


def main():
    OUT.mkdir(exist_ok=True)
    runtimes = [run_set("V4", V4), run_set("V9", V9)]
    pd.DataFrame(runtimes).to_csv(OUT / "TX4_BLIND_PREDICTION_RUNTIME.csv", index=False)
    hashes = []
    for name in ("TX4_V4_BLIND_PREDICTIONS.csv", "TX4_V9_BLIND_PREDICTIONS.csv", "TX4_BLIND_PREDICTION_RUNTIME.csv"):
        p = OUT / name
        hashes.append(f"{sha256(p)}  {p.relative_to(ROOT)}")
    (OUT / "TX4_BLIND_PREDICTION_HASHES.txt").write_text("\n".join(hashes) + "\n", encoding="utf-8")
    (OUT / "TX4_BLIND_PREDICTION_SETUP.json").write_text(
        json.dumps({"candidates": {"V4": list(V4), "V9": list(V9)}, "theta": THETA.as_dict(), "band_hz": BAND_HZ, "grid_sigma": SIGMA_GRID.tolist(), "grid_frequency_hz": FREQ_GRID.tolist(), "closure_grid_frequency_hz": FREQ_CLOSURE_GRID.tolist(), "answer_key_loaded": False}, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(runtimes, indent=2))


if __name__ == "__main__":
    main()
