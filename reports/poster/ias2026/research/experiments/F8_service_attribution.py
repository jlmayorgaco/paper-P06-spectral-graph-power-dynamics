"""F8 - synchronous-service causal attribution of the incompatibility hypergraph.

Representative points are taken from the FROZEN F7A map (tag
IAS2026_TRACKA_F7_POLICY_HYPERGRAPH_FREEZE), read only. At each point every
intervention is a physical model configuration applied consistently to all 16
subsets of the core, and the full hypergraph is recomputed on the direct path.

Class R - restore services of the RETIRED machines at their buses. Each retired
bus of every subset carries a synchronous condenser (no active power; the
converter keeps every megawatt). Full 2^7 factorial over

    EM   condenser rating                 0.25 / 1.00 of the retired machine
    S1   inertia (kinetic energy)          1 % / 100 % of the retired machine
    S2   mechanical damping D              0 / 2 pu on the condenser base
    S3   electromagnetic transients        frozen EMFs (classical) / dynamic
    S4   AVR                               manual (field held) / native AVR
    S5   PSS                               off / native
    S6   reactive support                  condenser carries 0 / all bus Q

plus two arms without electromagnetic presence: the plain replacement, and the
converters with synthetic inertia equal to the retired machine's H.

Class A - ablate or augment services of the SURVIVING synchronous fleet
(inertia x0.5/1/2, D 0/2, EMF transients, AVR, PSS). Full 3 x 2^4 factorial.

S7 governor: not represented in the harmonized model (TGOV1N dropped, reheat
time constant outside the band, Dt = 0), so it is excluded, not ablated.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import ndimage

from _bootstrap import RESULTS
from _f7_common import (
    LABELS,
    LEAK,
    SUBSETS,
    Theta,
    in_gamma,
    native_te,
    te_geometric_mean,
)
from _overnight import choose_workers, pin_blas_threads
from _postfreeze import PostFreezeExperiment as Experiment
from _v2c_common import (
    BAND_HZ,
    CORE,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from F7_report import kappa_of, points_file, raster
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_order,
    incompatibility_hypergraph,
)
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space, closure_profile

NAME = "F8_service_attribution"
OUT = RESULTS / "F8"
BAND_SAMPLES = 61
ZERO = 1e-3

R_FACTORS = {
    "EM_rating": (0.25, 1.00),
    "S1_inertia": (0.01, 1.00),
    "S2_damping": (0.0, 2.0),
    "S3_flux": (0.0, 1.0),
    "S4_avr": (0.0, 1.0),
    "S5_pss": (0.0, 1.0),
    "S6_qshare": (0.0, 1.0),
}
A_FACTORS = {
    "A1_inertia": (0.5, 1.0, 2.0),
    "A2_damping": (0.0, 2.0),
    "A3_flux": (0.0, 1.0),
    "A4_avr": (0.0, 1.0),
    "A5_pss": (0.0, 1.0),
}

_CONTEXT: dict = {}


# ---------------------------------------------------------------- points ----


def choose_points() -> dict[str, dict]:
    """Interior points of the frozen F7A map, one per kappa class, plus P_fold."""

    points = pd.read_csv(points_file("F7A"), low_memory=False)
    leaves = pd.read_csv(RESULTS / "F7" / "F7A_leaves.csv")
    grid = raster(points, leaves)
    xs = (
        points.groupby("I")
        .x.first()
        .reindex(range(grid.shape[0]))
        .interpolate()
        .to_numpy()
    )
    ys = (
        points.groupby("J")
        .y.first()
        .reindex(range(grid.shape[1]))
        .interpolate()
        .to_numpy()
    )
    chosen = {}
    for name, kappa in (("P_inf", -1), ("P4", 4), ("P3", 3), ("P2", 2), ("P1", 1)):
        labels = [s for s in set(grid.ravel()) if s and kappa_of(s) == kappa]
        label = max(labels, key=lambda s: (grid == s).sum())
        mask = grid == label
        comp, n = ndimage.label(mask, structure=np.ones((3, 3)))
        biggest = 1 + int(np.argmax(ndimage.sum(mask, comp, range(1, n + 1))))
        depth = ndimage.distance_transform_edt(comp == biggest)
        i, j = np.unravel_index(np.argmax(depth), depth.shape)
        chosen[name] = {
            "g": float(xs[i]),
            "k": float(ys[j]),
            "t": 1.5,
            "h": 1.0,
            "frozen_label": label,
            "depth_nodes": float(depth[i, j]),
        }
    tongue = json.loads((RESULTS / "F7" / "F7_tongue.json").read_text(encoding="utf-8"))
    tip = tongue["tip_fold"]
    chosen["P_fold"] = {
        "g": round(tip["g"], 5),
        "k": round(tip["k"] + 0.01, 5),
        "t": 1.5,
        "h": 1.0,
        "frozen_label": "30+33+35+37 (inside tongue, 0.01 above tip)",
        "depth_nodes": 0.0,
    }
    return chosen


# --------------------------------------------------------- configurations ----


def r_configs() -> list[dict]:
    configs = [{"class": "R", "name": "R_none", "condenser": None, "virtual_h": 0.0}]
    configs.append(
        {"class": "R", "name": "R_virtual_inertia", "condenser": None, "virtual_h": 1.0}
    )
    for levels in itertools.product(*[(0, 1)] * len(R_FACTORS)):
        spec = {f: R_FACTORS[f][lv] for f, lv in zip(R_FACTORS, levels, strict=True)}
        code = "".join(str(lv) for lv in levels)
        configs.append(
            {
                "class": "R",
                "name": f"R_{code}",
                "condenser": spec,
                "virtual_h": 0.0,
                "code": code,
            }
        )
    return configs


def a_configs() -> list[dict]:
    out = []
    for levels in itertools.product(*[range(len(v)) for v in A_FACTORS.values()]):
        spec = {f: A_FACTORS[f][lv] for f, lv in zip(A_FACTORS, levels, strict=True)}
        code = "".join(str(lv) for lv in levels)
        out.append({"class": "A", "name": f"A_{code}", "services": spec, "code": code})
    return out


def _machine_m(bus: int) -> float:
    return float(load_network().machines[bus]["M"])


def solve_config(members, theta: Theta, config: dict):
    """One subset under one intervention, at one F7 point, direct path."""

    mean = te_geometric_mean()
    bus_scaling = (
        {b: {"ta": (mean / te) ** (1.0 - theta.h)} for b, te in native_te().items()}
        if theta.h != 1.0
        else None
    )
    converter = ConverterParameters(
        voltage_control=True,
        voltage_gain=theta.g,
        voltage_leak=LEAK,
        inertia_emulation=config.get("virtual_h", 0.0) > 0.0,
    )
    converters = None
    if config.get("virtual_h", 0.0) > 0.0 and members:
        # synthetic inertia equal to each retired machine's H = M / 2, own base
        converters = {
            b: ConverterParameters(
                voltage_control=True,
                voltage_gain=theta.g,
                voltage_leak=LEAK,
                inertia_emulation=True,
                h_virtual=config["virtual_h"] * _machine_m(b) / 2.0,
            )
            for b in members
        }
    plan_kwargs = {}
    machine_services = None
    if config["class"] == "R" and config.get("condenser") and members:
        spec = config["condenser"]
        rating = spec["EM_rating"]
        plan_kwargs = {
            "condenser": {b: rating for b in members},
            "condenser_services": {
                # kinetic-energy fraction -> per-unit inertia on the condenser base
                "inertia": spec["S1_inertia"] / rating,
                "damping": spec["S2_damping"],
                "flux_blend": spec["S3_flux"],
                "avr_blend": spec["S4_avr"],
                "pss": spec["S5_pss"],
                "q_share": spec["S6_qshare"],
            },
        }
    if config["class"] == "A":
        spec = config["services"]
        machine_services = {
            "inertia": spec["A1_inertia"],
            "damping": spec["A2_damping"],
            "flux_blend": spec["A3_flux"],
            "avr_blend": spec["A4_avr"],
            "pss": spec["A5_pss"],
        }
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}, **plan_kwargs),
        converter=converter,
        converters=converters,
        machine_scaling={"ka": theta.k, "ta": theta.t},
        machine_bus_scaling=bus_scaling,
        machine_services=machine_services,
    )


def evaluate(task) -> dict:
    point_name, theta_dict, config = task
    theta = Theta(**theta_dict)
    row = {
        "point": point_name,
        **theta_dict,
        "config": config["name"],
        "class": config["class"],
        **(config.get("condenser") or {}),
        **config.get("services", {}),
        "virtual_h": config.get("virtual_h", 0.0),
    }
    counts, cases = {}, {}
    for members, label in zip(SUBSETS, LABELS, strict=True):
        try:
            case = solve_config(members, theta, config)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row.update(
                status="INFEASIBLE",
                label="INFEASIBLE",
                reason=f"{label}: {str(error)[:70]}",
            )
            return row
        values = np.linalg.eigvals(case.system.A)
        dynamic = values[np.abs(values) > ZERO]
        row[f"N_{label}"] = int(np.count_nonzero(in_gamma(values)))
        f = np.abs(values.imag) / (2 * math.pi)
        band = (
            (values.imag >= 0)
            & (np.abs(values) > ZERO)
            & (f >= BAND_HZ[0])
            & (f <= BAND_HZ[1])
        )
        row[f"bandmax_{label}"] = (
            float(values.real[band].max()) if band.any() else float("-inf")
        )
        row[f"maxload_{label}"] = case.max_loading
        counts[frozenset(members)] = row[f"N_{label}"]
        cases[members] = (case, values)
        if not members:
            row["base_abscissa"] = float(dynamic.real.max())
            if row["base_abscissa"] >= 0.0:
                row.update(status="BASE_UNSTABLE", label="BASE_UNSTABLE")
                return row
    edges = incompatibility_hypergraph(counts)
    row.update(
        label=hypergraph_label(edges), kappa=hypergraph_order(edges), status="OK"
    )
    base, base_values = cases[()]
    flag, flag_values = cases[CORE]
    band = [
        (v, abs(v.imag) / (2 * math.pi))
        for v in flag_values
        if v.imag >= 0
        and abs(v) > ZERO
        and BAND_HZ[0] <= abs(v.imag) / (2 * math.pi) <= BAND_HZ[1]
    ]
    worst = max(band, key=lambda p: p[0].real)
    row["band_alpha"] = float(worst[0].real)
    row["band_f_hz"] = float(worst[1])
    from ibr_cycles.dynamics.modes import eigen_analysis

    base_spec, flag_spec = eigen_analysis(base.system.A), eigen_analysis(flag.system.A)
    anchor = base_anchor(base_spec, base, _CONTEXT["nominal"])
    reading = (
        read_family(anchor, base, flag, flag_spec, common_labels=_CONTEXT["pinned"])
        if anchor
        else None
    )
    row["alpha_IA"] = reading.alpha if reading else float("nan")
    row["f_IA_hz"] = reading.frequency_worst_hz if reading else float("nan")
    try:
        profile = closure_profile(
            build_action_space(base, flag, CORE), band_hz=BAND_HZ, samples=BAND_SAMPLES
        )
        best = min(profile, key=lambda s: s.margin)
        row["m_cl"], row["f_port_hz"] = best.margin, best.frequency_hz
    except (np.linalg.LinAlgError, ValueError):
        row["m_cl"], row["f_port_hz"] = float("nan"), float("nan")
    return row


def _init():
    pin_blas_threads()
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    _CONTEXT.update(
        pinned=sorted(machine_labels(flagship.system.labels)),
        nominal=nominal_reference(),
    )


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    points = choose_points()
    (OUT / "F8_points.json").write_text(json.dumps(points, indent=2), encoding="utf-8")
    configs = r_configs() + a_configs()
    tasks = [
        (name, {k: p[k] for k in ("g", "k", "t", "h")}, config)
        for name, p in points.items()
        for config in configs
    ]
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Which synchronous services reshape the incompatibility hypergraph?",
        config={
            "points": points,
            "R_factors": R_FACTORS,
            "A_factors": A_FACTORS,
            "n_configs": len(configs),
            "source_map": "frozen F7A",
            "S7": "governor not represented; excluded",
        },
        workers=workers,
    )
    started = time.time()
    experiment.note(f"{len(tasks)} (point, intervention) cases on {workers} workers")
    with Pool(workers, initializer=_init) as pool:
        rows = pool.map(evaluate, tasks, chunksize=2)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "F8_interventions.csv", index=False)
    experiment.finish("COMPUTED", cases=len(frame), elapsed_s=time.time() - started)
    print(frame.groupby(["point", "class"]).label.nunique())
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
