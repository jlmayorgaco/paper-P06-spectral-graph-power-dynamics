"""F8C - mechanism-conditioned port closure at service-induced boundaries.

For every pair of F8 configurations that differ in ONE service and have
different H_Gamma (both ends feasible, base stable), the service is moved
continuously between its two levels, each subset whose indicator N > 0 flips is
bisected to its physical boundary on the direct path, the crossing is
classified (Safeguard B), and the port closure of that witness subset is
evaluated at the crossing frequency against the base of the SAME configuration.

Also: the electromagnetic-presence path (condenser rating 1e-3 -> 0.25 at fixed
per-unit services) at every F8 point, the synthetic-inertia path, the AVR blend
path on the surviving fleet, and the AVR-necessity sweep along the policy gain.

Continuous service levels: inertia (kinetic fraction or scale) geometric,
everything else linear. Every configuration keeps the operating point.
"""

from __future__ import annotations

import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import SUBSETS, Theta, in_gamma
from _overnight import choose_workers, pin_blas_threads
from _postfreeze import PostFreezeExperiment
from _v2c_common import BAND_HZ
from F8_service_attribution import A_FACTORS, R_FACTORS, solve_config
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    incompatibility_hypergraph,
)
from ibr_cycles.models.ieee39_case import InfeasibleReplacement
from ibr_cycles.models.port_admittance import (
    build_action_space,
    closure_at,
    refine_closure,
)

NAME = "F8C_service_closure"
F8 = RESULTS / "F8"
OUT = RESULTS / "F8C"
ITER = 26
GEOMETRIC = {"S1_inertia", "A1_inertia"}


def _levels(config: dict) -> dict:
    return dict(config.get("condenser") or config.get("services") or {})


def interpolate(cfg_a: dict, cfg_b: dict, factor: str, s: float) -> dict:
    out = json.loads(json.dumps(cfg_a))
    if factor == "virtual_h":
        out["virtual_h"] = (1 - s) * cfg_a.get("virtual_h", 0.0) + s * cfg_b.get(
            "virtual_h", 0.0
        )
        return out
    if factor == "EM_presence":
        out["condenser"] = dict(cfg_b["condenser"])
        out["condenser"]["EM_rating"] = 1e-3 + s * (
            cfg_b["condenser"]["EM_rating"] - 1e-3
        )
        out["em_fixed_pu"] = True
        return out
    key = "condenser" if cfg_a["class"] == "R" else "services"
    a, b = cfg_a[key][factor], cfg_b[key][factor]
    value = (
        a * (b / a) ** s
        if factor in GEOMETRIC and a > 0 and b > 0
        else (1 - s) * a + s * b
    )
    out[key][factor] = value
    return out


def _solve(members, theta, config):
    if config.get("em_fixed_pu") and members:
        # EM-presence path: per-unit services fixed at the 0.25-rating values
        cfg = json.loads(json.dumps(config))
        cfg["condenser"]["S1_inertia"] = (
            cfg["condenser"]["S1_inertia"] * cfg["condenser"]["EM_rating"] / 0.25
        )
        return solve_config(members, theta, cfg)
    return solve_config(members, theta, config)


def _count(members, theta, config):
    case = _solve(members, theta, config)
    values = np.linalg.eigvals(case.system.A)
    return int(np.count_nonzero(in_gamma(values))), values


def _violation(v: complex) -> str:
    f = abs(v.imag) / (2 * math.pi)
    parts = []
    if v.real <= 0:
        parts.append("IMAGINARY_AXIS")
    if f < BAND_HZ[0]:
        parts.append("LOWER_BAND_EDGE")
    if f > BAND_HZ[1]:
        parts.append("UPPER_BAND_EDGE")
    return "+".join(parts) or "NONE"


def locate(task) -> list[dict]:
    """All subset boundaries along one single-service path, with port closure."""

    info, theta_d, cfg_a, cfg_b, factor = task
    theta = Theta(**theta_d)
    out = []
    try:
        ends = []
        for s in (0.0, 1.0):
            cfg = interpolate(cfg_a, cfg_b, factor, s)
            ends.append({m: _count(m, theta, cfg) for m in SUBSETS})
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
        return [{**info, "status": f"INFEASIBLE: {str(error)[:60]}"}]
    for members in SUBSETS[1:]:
        na, nb = ends[0][members][0], ends[1][members][0]
        if (na > 0) == (nb > 0):
            continue
        lo, hi, v_lo, v_hi = 0.0, 1.0, ends[0][members][1], ends[1][members][1]
        unsafe_lo = na > 0
        for _ in range(ITER):
            mid = 0.5 * (lo + hi)
            n, v = _count(members, theta, interpolate(cfg_a, cfg_b, factor, mid))
            if (n > 0) == unsafe_lo:
                lo, v_lo = mid, v
            else:
                hi, v_hi = mid, v
        unsafe, safe = (v_lo, v_hi) if unsafe_lo else (v_hi, v_lo)
        inside = unsafe[in_gamma(unsafe)]
        upper = safe[safe.imag >= 0]
        best = None
        for lam in inside:
            partner = upper[np.argmin(np.abs(upper - lam))]
            kind = _violation(partner)
            if kind != "NONE":
                best = (lam, partner, kind)
                break
        row = {
            **info,
            "status": "LOCATED",
            "factor": factor,
            "witness": "+".join(map(str, members)),
            "s_star": 0.5 * (lo + hi),
            "direction": "LEAVES" if unsafe_lo else "ENTERS",
            "N_a": na,
            "N_b": nb,
        }
        if best is None:
            row["boundary_type"] = "UNRESOLVED"
            out.append(row)
            continue
        lam, partner, kind = best
        star = 0.5 * (lam + partner)
        row.update(
            boundary_type=kind,
            lambda_re=float(star.real),
            lambda_im=float(star.imag),
            freq_hz=float(abs(star.imag) / (2 * math.pi)),
            rhp_unsafe=int(
                np.count_nonzero(
                    (unsafe.real > 0) & (np.abs(unsafe) > 1e-3) & (unsafe.imag >= 0)
                )
            ),
            rhp_safe=int(
                np.count_nonzero(
                    (safe.real > 0) & (np.abs(safe) > 1e-3) & (safe.imag >= 0)
                )
            ),
        )
        if kind == "IMAGINARY_AXIS":
            cfg = interpolate(cfg_a, cfg_b, factor, row["s_star"])
            try:
                base = _solve((), theta, cfg)
                case = _solve(members, theta, cfg)
                space = build_action_space(base, case, members)
                omega = abs(star.imag)
                split = space.split(complex(0.0, omega))
                ref = np.median(
                    [
                        abs(complex(space.split(complex(0.0, w))["full"]))
                        for w in np.linspace(
                            2 * math.pi * BAND_HZ[0], 2 * math.pi * BAND_HZ[1], 25
                        )
                    ]
                )
                total = complex(split["full"])
                product = complex(split["individual"]) * complex(split["collective"])
                row.update(
                    port_det_ratio=abs(total) / ref,
                    port_visible=bool(abs(total) / ref < 1e-3),
                    individual_min=float(min(abs(d) for d in split["diagonals"])),
                    closure_distance=float(
                        np.min(np.abs(np.asarray(split["q_eigenvalues"]) + 1.0))
                    )
                    if len(members) > 1
                    else float("nan"),
                    factorization_residual=abs(total - product)
                    / max(abs(complex(split["individual"])), 1e-300),
                    solve_residual=closure_at(space, omega).solve_residual,
                )
                if len(members) > 1:
                    refined = refine_closure(
                        space, omega, half_width=0.05 * omega, tolerance=1e-9
                    )
                    row.update(
                        closure_refined=refined.margin,
                        freq_port_hz=refined.frequency_hz,
                    )
            except (np.linalg.LinAlgError, ValueError, InfeasibleReplacement) as error:
                row["port_error"] = str(error)[:80]
        out.append(row)
    if not out:
        out.append({**info, "status": "NO_SUBSET_FLIP"})
    return out


def _configs():
    from F8_service_attribution import a_configs, r_configs

    return {c["name"]: c for c in r_configs() + a_configs()}


def build_tasks() -> list:
    frame = pd.read_csv(F8 / "F8_interventions.csv", low_memory=False)
    points = json.loads((F8 / "F8_points.json").read_text(encoding="utf-8"))
    configs = _configs()
    tasks = []
    for fname, factors in (
        ("F8_flips_R.csv", R_FACTORS),
        ("F8_flips_A.csv", A_FACTORS),
    ):
        flips = pd.read_csv(F8 / fname)
        flips = flips[
            (flips.status_from == "OK") & (flips.status_to == "OK") & flips.H_changes
        ]
        prefix = "R_" if fname.endswith("R.csv") else "A_"
        for _, f in flips.iterrows():
            i = list(factors).index(f.factor)
            code = [int(c) for c in str(f.context).zfill(len(factors))]
            other = list(code)
            other[i] += 1
            a, b = (
                configs[prefix + "".join(map(str, code))],
                configs[prefix + "".join(map(str, other))],
            )
            theta = {k: points[f.point][k] for k in ("g", "k", "t", "h")}
            tasks.append(
                (
                    {
                        "point": f.point,
                        "pair": f"{a['name']}->{b['name']}",
                        "H_a": f.H_from,
                        "H_b": f.H_to,
                        "path": "factorial_flip",
                    },
                    theta,
                    a,
                    b,
                    f.factor,
                )
            )
    # electromagnetic presence: tiny condenser -> 0.25 rating, minimal services
    minimal = configs["R_0000000"]
    for name, p in points.items():
        theta = {k: p[k] for k in ("g", "k", "t", "h")}
        tasks.append(
            (
                {
                    "point": name,
                    "pair": "EM 1e-3 -> 0.25 (R_0000000)",
                    "path": "EM_presence",
                    "H_a": "",
                    "H_b": frame[
                        (frame.point == name) & (frame.config == "R_0000000")
                    ].label.iloc[0],
                },
                theta,
                minimal,
                minimal,
                "EM_presence",
            )
        )
        vi = configs["R_virtual_inertia"]
        tasks.append(
            (
                {
                    "point": name,
                    "pair": "virtual inertia 0 -> H",
                    "path": "virtual_inertia",
                    "H_a": frame[
                        (frame.point == name) & (frame.config == "R_none")
                    ].label.iloc[0],
                    "H_b": frame[
                        (frame.point == name) & (frame.config == "R_virtual_inertia")
                    ].label.iloc[0],
                },
                theta,
                {**vi, "virtual_h": 1e-6},
                vi,
                "virtual_h",
            )
        )
    return tasks


def avr_necessity(pool) -> pd.DataFrame:
    """H along the policy gain, surviving fleet on AVR or on manual excitation."""

    points = json.loads((F8 / "F8_points.json").read_text(encoding="utf-8"))
    rows = []
    jobs = []
    for pname in ("P4", "P1"):
        k = points[pname]["k"]
        for damping in (0.0, 2.0):
            for avr in (0.0, 1.0):
                cfg = {
                    "class": "A",
                    "name": f"A_nec_D{damping}_avr{avr}",
                    "services": {
                        "A1_inertia": 1.0,
                        "A2_damping": damping,
                        "A3_flux": 1.0,
                        "A4_avr": avr,
                        "A5_pss": 1.0,
                    },
                }
                for g in np.round(np.linspace(0, 1, 41), 3):
                    jobs.append((pname, k, damping, avr, float(g), cfg))
    results = pool.map(_necessity_job, jobs, chunksize=2)
    for job, res in zip(jobs, results, strict=True):
        rows.append(
            {
                "point_k": job[0],
                "k": job[1],
                "D": job[2],
                "avr": job[3],
                "g": job[4],
                **res,
            }
        )
    return pd.DataFrame(rows)


def _necessity_job(job):
    pname, k, damping, avr, g, cfg = job
    theta = Theta(g=g, k=k, t=1.5)
    counts = {}
    for members in SUBSETS:
        try:
            n, v = _count(members, theta, cfg)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
            return {"label": "INFEASIBLE"}
        if not members:
            dyn = v[np.abs(v) > 1e-3]
            if dyn.real.max() >= 0:
                return {"label": "BASE_UNSTABLE"}
        counts[frozenset(members)] = n
    return {"label": hypergraph_label(incompatibility_hypergraph(counts))}


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = build_tasks()
    workers = choose_workers()
    exp = PostFreezeExperiment(
        name=NAME,
        question="Does port closure mark service-induced boundaries?",
        config={"transitions": len(tasks), "bisection": ITER},
        workers=workers,
    )
    started = time.time()
    exp.note(f"{len(tasks)} single-service paths")
    with Pool(workers, initializer=pin_blas_threads) as pool:
        results = pool.map(locate, tasks, chunksize=1)
        necessity = avr_necessity(pool)
    events = pd.DataFrame([r for group in results for r in group])
    events.to_csv(OUT / "F8C_events.csv", index=False)
    necessity.to_csv(OUT / "F8C_avr_necessity.csv", index=False)
    loc = events[events.status == "LOCATED"]
    imag = loc[loc.boundary_type == "IMAGINARY_AXIS"]
    multi = imag[imag.witness.str.contains(r"\+")]
    summary = {
        "paths": len(tasks),
        "located": int(len(loc)),
        "boundary_types": loc.boundary_type.value_counts().to_dict(),
        "port_visible": int(imag.port_visible.fillna(False).astype(bool).sum()),
        "imaginary_axis": int(len(imag)),
        "multi_port": int(len(multi)),
        "multi_closure_lt_1e-4": int((multi.closure_distance < 1e-4).sum()),
        "closure_max": float(multi.closure_distance.max()) if len(multi) else None,
        "single_port": int(len(imag) - len(multi)),
        "single_individual_max": float(
            imag[~imag.witness.str.contains(r"\+")].individual_min.max()
        )
        if len(imag) > len(multi)
        else None,
        "port_errors": int(
            events.get("port_error", pd.Series(dtype=str)).notna().sum()
        ),
        "by_factor": loc.groupby("factor").size().to_dict(),
    }
    (OUT / "F8C_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    exp.finish("COMPUTED", **summary, elapsed_s=time.time() - started)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
