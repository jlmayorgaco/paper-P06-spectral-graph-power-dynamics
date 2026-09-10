"""F11 - computational scaling of three evaluation methods.

M1  full DAE build + equilibrium + central-difference Jacobian + eigensolve,
    per portfolio and per parameter point (the reference)
M2  localized A_red assembly: three M1 solves per portfolio at initialization,
    then per point only the converter rows (affine in g) and each machine's efd
    row (K/T, 1/T) move, followed by one eigensolve per portfolio
M3  reduced port core: one base and one all-replaced operator per point, every
    portfolio a principal minor of I + M(s), counted by the argument principle

Lattices: the 16 core subsets and the full 2^9 = 512 portfolios over the nine
candidate buses. Points: the six F8 points and random frozen-F7A points.
Times are CPU seconds per task (process_time) so that parallel execution does
not flatter any method. Accuracy is always against M1.

Regimes. FROZEN equilibrium: the matched dispatch, where M2 and M3 are exact
reformulations. RE-EQUILIBRATED: unity power factor, where every portfolio has
its own operating point; M2 and M3 are then applied naively (their references
are the matched ones) to measure how wrong reuse is. No speed-up is claimed
there.
"""

from __future__ import annotations

import json
import math
import sys
import time
import tracemalloc
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import LEAK, Theta, in_gamma
from _overnight import choose_workers, pin_blas_threads
from _postfreeze import PostFreezeExperiment
from _v2c_common import BAND_HZ, CORE
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    incompatibility_hypergraph,
)
from ibr_cycles.models.ieee39_case import ReplacementPlan, build_dae, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_core import PortCore, all_subsets

NAME = "F11_scaling"
OUT = RESULTS / "F11"
LO, HI = 2 * math.pi * BAND_HZ[0], 2 * math.pi * BAND_HZ[1]
CANDIDATES = tuple(load_network().replacement_candidates)
SEED = 20260913
N_RANDOM = 18


def _converter(theta, policy):
    if policy == "unity_pf":
        return ConverterParameters()
    return ConverterParameters(
        voltage_control=True, voltage_gain=theta.g, voltage_leak=LEAK
    )


def _plan(members, policy):
    return ReplacementPlan.of(
        {b: 1.0 for b in members},
        q_policy="unity_pf" if policy == "unity_pf" else "matched",
    )


def _scaling(theta):
    return {"ka": theta.k, "ta": theta.t}


def count(values) -> int:
    return int(np.count_nonzero(in_gamma(values)))


# -------------------------------------------------------------------- M1 ----


def m1_task(task):
    members, theta_d, policy = task
    theta = Theta(**theta_d)
    t0 = time.process_time()
    case = solve_case(
        _plan(members, policy),
        converter=_converter(theta, policy),
        machine_scaling=_scaling(theta),
    )
    values = np.linalg.eigvals(case.system.A)
    return {
        "members": members,
        "N": count(values),
        "cpu": time.process_time() - t0,
        "n_states": case.n_states,
        "band": _band(values),
    }


def _band(values):
    f = np.abs(values.imag) / (2 * math.pi)
    keep = (
        (values.imag >= 0)
        & (np.abs(values) > 1e-3)
        & (f >= BAND_HZ[0])
        & (f <= BAND_HZ[1])
    )
    return sorted(values[keep].tolist(), key=lambda v: (v.imag, v.real))


# -------------------------------------------------------------------- M2 ----


def m2_init_task(task):
    """Reference parts of one portfolio: the F7 localized assembly, generic lattice."""

    members, policy = task
    t0 = time.process_time()
    ka = {int(b): float(v["KA"]) for b, v in _payload()["avr_by_bus"].items()}
    te = {int(b): float(v["TE"]) for b, v in _payload()["avr_by_bus"].items()}

    def solve(g, k):
        th = Theta(g=g, k=k, t=1.0)
        return solve_case(
            _plan(members, policy),
            converter=_converter(th, "matched"),
            machine_scaling=_scaling(th),
        )

    ref, g1, k2 = solve(0.0, 1.0), solve(1.0, 1.0), solve(0.0, 2.0)
    a0 = ref.system.A.copy()
    rows = np.array(["_gfl" in n for n in ref.system.labels])
    dg = np.zeros_like(a0)
    dg[rows] = g1.system.A[rows] - a0[rows]
    efd = []
    for i, name in enumerate(ref.system.labels):
        if name.startswith("efd_sg"):
            bus = int(name.removeprefix("efd_sg"))
            r1, r2 = a0[i], k2.system.A[i]
            efd.append((i, bus, (r2 - r1) * te[bus] / ka[bus], (2 * r1 - r2) * te[bus]))
    return {
        "members": members,
        "parts": (a0, dg, efd),
        "cpu": time.process_time() - t0,
        "bytes": a0.nbytes * 2 + sum(u.nbytes + w.nbytes for _, _, u, w in efd),
    }


def _payload():
    from ibr_cycles.models.ieee39_case import _controller_payload

    return _controller_payload()


def m2_point(parts: dict, theta: Theta) -> dict:
    ka = {int(b): float(v["KA"]) for b, v in _payload()["avr_by_bus"].items()}
    te = {int(b): float(v["TE"]) for b, v in _payload()["avr_by_bus"].items()}
    out = {}
    for members, (a0, dg, efd) in parts.items():
        a = a0 + theta.g * dg
        for i, bus, u, w in efd:
            kk, tt = theta.k * ka[bus], theta.t * te[bus]
            a[i] = (kk / tt) * u + (1.0 / tt) * w
        values = np.linalg.eigvals(a)
        out[members] = (count(values), _band(values))
    return out


# -------------------------------------------------------------------- M3 ----


def m3_point(theta: Theta, buses, subsets, policy="matched"):
    conv = _converter(theta, policy)
    base = build_dae(_plan((), policy), converter=conv, machine_scaling=_scaling(theta))
    repl = build_dae(
        _plan(buses, policy), converter=conv, machine_scaling=_scaling(theta)
    )
    core = PortCore.from_daes(*base, *repl, buses)
    counts, info = core.counts(subsets, omega_lo=LO, omega_hi=HI)
    return {s: counts[s]["delta_N"] for s in subsets}, info, core


def label_of(counts: dict) -> str:
    return hypergraph_label(
        incompatibility_hypergraph({frozenset(s): n for s, n in counts.items()})
    )


# ------------------------------------------------------------------ driver ----


def sample_points():
    pts = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    chosen = {n: {k: p[k] for k in ("g", "k", "t", "h")} for n, p in pts.items()}
    spots = pd.read_csv(RESULTS / "F7" / "F7A_spotcheck.csv")
    rng = np.random.default_rng(SEED)
    for i, r in enumerate(
        spots.iloc[rng.choice(len(spots), N_RANDOM, replace=False)].itertuples()
    ):
        chosen[f"R{i:02d}"] = {"g": float(r.x), "k": float(r.y), "t": 1.5, "h": 1.0}
    return chosen


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    points = sample_points()
    workers = choose_workers()
    exp = PostFreezeExperiment(
        name=NAME,
        question="What does each method cost, and how accurate is it?",
        config={"points": points, "candidates": CANDIDATES, "seed": SEED},
        workers=workers,
    )
    started = time.time()
    rows, summary = [], {}
    scenarios = [
        ("core16", CORE, list(points)),
        ("lattice512", CANDIDATES, list(points)[:6]),
    ]
    with Pool(workers, initializer=pin_blas_threads) as pool:
        for scen, buses, pnames in scenarios:
            subsets = all_subsets(buses)
            # M2 initialisation (once per lattice)
            w0 = time.time()
            init = pool.map(
                m2_init_task, [(s, "matched") for s in subsets], chunksize=2
            )
            m2_init_wall = time.time() - w0
            parts = {r["members"]: r["parts"] for r in init}
            m2_init_cpu = sum(r["cpu"] for r in init)
            m2_bytes = sum(r["bytes"] for r in init)
            summary[scen] = {
                "subsets": len(subsets),
                "m2_init_cpu_s": m2_init_cpu,
                "m2_init_wall_s": m2_init_wall,
                "m2_reference_MB": m2_bytes / 2**20,
            }
            for pname in pnames:
                theta = Theta(**points[pname])
                w0 = time.time()
                m1 = pool.map(
                    m1_task,
                    [(s, points[pname], "matched") for s in subsets],
                    chunksize=2,
                )
                m1_wall = time.time() - w0
                truth = {r["members"]: r["N"] for r in m1}
                t0 = time.process_time()
                m2 = m2_point(parts, theta)
                m2_cpu = time.process_time() - t0
                tracemalloc.start()
                t0 = time.process_time()
                m3, info, core = m3_point(theta, buses, subsets)
                m3_cpu = time.process_time() - t0
                m3_peak = tracemalloc.get_traced_memory()[1] / 2**20
                tracemalloc.stop()
                band_err = (
                    max(
                        (
                            min(abs(complex(a) - complex(b)) for b in m2[s][1])
                            if m2[s][1]
                            else 0.0
                        )
                        for r in m1
                        for s in [r["members"]]
                        for a in r["band"]
                    )
                    if any(r["band"] for r in m1)
                    else 0.0
                )
                rows.append(
                    {
                        "scenario": scen,
                        "point": pname,
                        **points[pname],
                        "subsets": len(subsets),
                        "m1_cpu_s": sum(r["cpu"] for r in m1),
                        "m1_wall_s": m1_wall,
                        "m1_states_min": min(r["n_states"] for r in m1),
                        "m1_states_max": max(r["n_states"] for r in m1),
                        "m2_cpu_s": m2_cpu,
                        "m3_cpu_s": m3_cpu,
                        "m3_peak_MB": m3_peak,
                        "m3_contour_points": info["contour_points"],
                        "m3_operator_dim": core.base.dimension,
                        "m3_action_dim": 2 * len(buses),
                        "m2_band_eig_max_err": band_err,
                        "m2_count_errors": sum(m2[s][0] != truth[s] for s in subsets),
                        "m3_count_errors": sum(m3[s] != truth[s] for s in subsets),
                        "H_m1": label_of(truth),
                        "H_m2": label_of({s: m2[s][0] for s in subsets}),
                        "H_m3": label_of(m3),
                    }
                )
                exp.note(
                    f"{scen} {pname}: M1 {rows[-1]['m1_cpu_s']:.1f}s M2 {m2_cpu:.2f}s "
                    f"M3 {m3_cpu:.2f}s errors M2 {rows[-1]['m2_count_errors']} "
                    f"M3 {rows[-1]['m3_count_errors']}"
                )
        # re-equilibrated regime: unity power factor, core lattice, naive reuse
        subsets = all_subsets(CORE)
        reeq = []
        for pname in list(points)[:6]:
            theta = Theta(**points[pname])
            m1 = pool.map(
                m1_task, [(s, points[pname], "unity_pf") for s in subsets], chunksize=2
            )
            truth = {r["members"]: r["N"] for r in m1}
            m3, _, _ = m3_point(theta, CORE, subsets, policy="unity_pf")
            reeq.append(
                {
                    "point": pname,
                    "H_m1_unity": label_of(truth),
                    "H_m3_naive": label_of(m3),
                    "m3_naive_count_errors": sum(m3[s] != truth[s] for s in subsets),
                    "m1_cpu_s": sum(r["cpu"] for r in m1),
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "F11_benchmark.csv", index=False)
    pd.DataFrame(reeq).to_csv(OUT / "F11_reequilibrated.csv", index=False)
    for scen in summary:
        f = frame[frame.scenario == scen]
        summary[scen].update(
            {
                "points": int(len(f)),
                "m1_cpu_per_point_s": float(f.m1_cpu_s.mean()),
                "m2_cpu_per_point_s": float(f.m2_cpu_s.mean()),
                "m3_cpu_per_point_s": float(f.m3_cpu_s.mean()),
                "speedup_m1_over_m2": float(f.m1_cpu_s.mean() / f.m2_cpu_s.mean()),
                "speedup_m1_over_m3": float(f.m1_cpu_s.mean() / f.m3_cpu_s.mean()),
                "m2_count_errors": int(f.m2_count_errors.sum()),
                "m3_count_errors": int(f.m3_count_errors.sum()),
                "H_errors_m2": int((f.H_m1 != f.H_m2).sum()),
                "H_errors_m3": int((f.H_m1 != f.H_m3).sum()),
                "m2_band_eig_max_err": float(f.m2_band_eig_max_err.max()),
                "m3_peak_MB": float(f.m3_peak_MB.max()),
            }
        )
    summary["reequilibrated"] = pd.DataFrame(reeq).to_dict("records")
    (OUT / "F11_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    exp.finish("COMPUTED", summary=summary, elapsed_s=time.time() - started)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
