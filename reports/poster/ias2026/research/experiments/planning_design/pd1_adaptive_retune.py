# ruff: noqa: E402
"""PD1 - adaptive per-point converter retune at the 240 unstable E35 samples.

Preregistration: docs/20260913_PLANNING_DESIGN_PREREG.md (rules frozen in
configs/ias2026/planning_design_prereg_v1.yaml before any campaign run).

For every E35 sample whose four-unit portfolio is unstable (240 of 1004 accepted):

1. regenerate the operating point exactly as E35 did (same seed rule, same draws),
   and check the frozen E35 flagship label (alpha_IA reproduction gate);
2. label the flagship, the frozen RC retune (E21/E32) and the 25 % condenser with
   (a) the E35 band criterion (no inter-area-band RHP eigenvalue) and
   (b) the series metric: the transverse whole-RHP classification of the FC01
   direct path (certification.physical.physical_report, status_rest);
3. re-optimize the E34 M1 converter-only retune at this operating point: minimum
   ||v||^2 over the five E34 physical log-coordinates, bounds +-ln 8, SLSQP, start
   -0.15, maxiter 90, ftol 1e-7 (E34 verbatim), subject to the transverse spectral
   abscissa <= TARGET (the E34 disc rule |lambda| > 1e-3 is replaced by the exact
   transverse quotient, which is the only change to the E34 optimizer);
4. label the adaptive solution with (a) and (b).
"""

from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

import time

import _pd_common as P
import numpy as np
from scipy.optimize import minimize

from _v2c_common import CORE, base_anchor, machine_labels, nominal_reference, read_family, BAND_HZ
from E32_tds_validation import rc_converter
from E34_mitigation_frontier import PHYSICAL, converter_from_physical
from E35_mc_operating import DISPATCH_JITTER, LOAD_SCATTER, SEED
from ibr_cycles.certification.physical import physical_report
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator
from ibr_cycles.certification.transverse import transverse_operator
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import InfeasibleReplacement, ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.uncertainty.sampling import sample_operating_point

FLAGSHIP = ReplacementPlan.of({b: 1.0 for b in CORE})
CONDENSER = ReplacementPlan.of({b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE})
BOUND = float(np.log(8.0))
START = -0.15
MAXITER = 90
FTOL = 1e-7
TARGET = -0.02  # E34 first declared margin; frozen in the config
FAILED = (InfeasibleReplacement, ValueError, np.linalg.LinAlgError, RuntimeError)

_CTX: dict = {}


def initialise():
    if _CTX:
        return
    flagship = solve_case(FLAGSHIP)
    _CTX.update(
        network=load_network(),
        pinned=sorted(machine_labels(flagship.system.labels)),
        nominal=nominal_reference(),
        rc=rc_converter(),
    )


def operating_point(row) -> dict:
    initialise()
    rng = np.random.default_rng(SEED + 100003 * int(row["sample"]))
    point = sample_operating_point(
        _CTX["network"],
        active_load=float(row["load"]),
        reactive_load=float(row["reactive_load"]),
        availability=float(row["availability"]),
        availability_buses=CORE,
        rng=rng,
        load_scatter=LOAD_SCATTER,
        dispatch_jitter=DISPATCH_JITTER,
    )
    if not point.dispatchable:
        raise RuntimeError("regenerated point not dispatchable")
    base = solve_case(ReplacementPlan.of({}), network=point.network)
    anchor = base_anchor(eigen_analysis(base.system.A), base, _CTX["nominal"])
    if anchor is None:
        raise RuntimeError("no anchor")
    return {"network": point.network, "base": base, "anchor": anchor}


def _band(ctx, case) -> dict:
    spectrum = eigen_analysis(case.system.A)
    reading = read_family(ctx["anchor"], ctx["base"], case, spectrum, common_labels=_CTX["pinned"])
    rhp = sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0.0)
    return {"alpha_IA": reading.family.alpha if reading else float("nan"), "band_rhp": int(rhp)}


def _perp(case) -> dict:
    rep = physical_report(case)
    return {"status_perp": rep.status_rest, "alpha_perp": float(rep.abscissa_rest),
            "n_dead": rep.n_dead, "partner": rep.partner}


def labels_of(ctx, case) -> dict:
    return {**_band(ctx, case), **_perp(case)}


def flagship_labels(ctx) -> dict:
    return labels_of(ctx, solve_case(FLAGSHIP, network=ctx["network"]))


def abscissa_perp(ctx, vector) -> float:
    """Fast transverse abscissa used inside the optimizer (no classification)."""

    case = solve_case(FLAGSHIP, network=ctx["network"], converter=converter_from_physical(vector))
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    return float(np.linalg.eigvals(transverse_operator(case.system.A, r_x, w).a_perp).real.max())


def adaptive(ctx, target=TARGET) -> dict:
    cache: dict = {}

    def value(v):
        key = np.round(v, 14).tobytes()
        if key not in cache:
            try:
                cache[key] = abscissa_perp(ctx, v)
            except FAILED:
                cache[key] = None
        return cache[key]

    def constraint(v):
        a = value(v)
        return -1.0 if a is None else target - a

    t0 = time.time()
    sol = minimize(
        lambda v: float(v @ v),
        np.full(len(PHYSICAL), START),
        method="SLSQP",
        bounds=[(-BOUND, BOUND)] * len(PHYSICAL),
        constraints=[{"type": "ineq", "fun": constraint}],
        options={"maxiter": MAXITER, "ftol": FTOL},
    )
    out = {
        "slsqp_success": bool(sol.success),
        "slsqp_message": str(sol.message)[:80],
        "slsqp_nit": int(sol.nit),
        "n_evaluations": len(cache),
        "wall_s": time.time() - t0,
        "change_norm": float(np.linalg.norm(sol.x)),
        "bound_active": bool(np.any(np.abs(sol.x) >= BOUND - 1e-4)),
        "x": [float(v) for v in sol.x],
    }
    try:
        case = solve_case(FLAGSHIP, network=ctx["network"], converter=converter_from_physical(sol.x))
        out.update({f"ad_{k}": v for k, v in labels_of(ctx, case).items()})
    except FAILED as e:
        out.update(ad_status_perp="INFEASIBLE", ad_alpha_perp=float("nan"), ad_error=str(e)[:80])
    return out


def run_sample(row: dict) -> dict:
    rec = {"sample": int(row["sample"]), "alpha_IA_e35": float(row["alpha_IA_exact"]),
           "rc_success_e35": bool(row["rc_success"]), "sc_success_e35": bool(row["sc_success"])}
    t0 = time.time()
    try:
        ctx = operating_point(row)
        rec.update({f"flag_{k}": v for k, v in flagship_labels(ctx).items()})
        for tag, kw in (("rc", {"converter": _CTX["rc"]}), ("sc", {})):
            plan = FLAGSHIP if tag == "rc" else CONDENSER
            try:
                rec.update({f"{tag}_{k}": v for k, v in labels_of(ctx, solve_case(plan, network=ctx["network"], **kw)).items()})
            except FAILED as e:
                rec[f"{tag}_status_perp"] = "INFEASIBLE"
                rec[f"{tag}_error"] = str(e)[:80]
        rec.update(adaptive(ctx))
        rec["ok"] = True
    except FAILED as e:
        rec["ok"] = False
        rec["error"] = str(e)[:200]
    rec["sample_wall_s"] = time.time() - t0
    return rec
