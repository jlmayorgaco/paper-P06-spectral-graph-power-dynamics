"""F2/F2B - numerical corroboration of the composability-region theorem.

The theorem predicts that ``kappa_Gamma`` is constant on connected components of
the regular parameter set, and that it can change only where some subset has an
eigenvalue on the boundary of Gamma. A one-dimensional continuation is the
cleanest test: sweep a single physical parameter finely, record for EVERY subset
the largest real part inside the frozen band, and check that every change of
kappa coincides with a zero crossing of one of those curves.

The sweep runs along the excitation-gain scale at a fixed excitation time-scale,
preserving each machine's own relative excitation data.
"""

from __future__ import annotations

import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _overnight import Experiment, choose_workers, pin_blas_threads
from _v2c_common import BAND_HZ, CORE
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.port_admittance import build_action_space, closure_profile

NAME = "F2_composability_continuation"
ZERO = 1e-3
BAND_SAMPLES = 61
SUBSETS = [tuple(sorted(s)) for k in range(len(CORE) + 1) for s in combinations(CORE, k)]
T_SCALE = 1.5           # fixed, chosen before the sweep: the F1B row richest in kappa values
GAINS = np.round(np.arange(0.50, 2.001, 0.025), 4)


def _initialise():
    pin_blas_threads()


def band_max_real(spectrum) -> float:
    """Largest real part inside the frozen band: the signed distance to dGamma."""

    modes = band_candidates(spectrum, BAND_HZ)
    return float(max(m.real for m in modes)) if modes else float("-inf")


def gamma_count(spectrum) -> int:
    return sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0.0)


def _evaluate(gain):
    scaling = {"ka": float(gain), "ta": float(T_SCALE)}
    row = {"gain_scale": float(gain), "t_scale": T_SCALE}
    counts, margins = {}, {}
    base = None
    flagship = None
    for members in SUBSETS:
        try:
            case = solve_case(
                ReplacementPlan.of({b: 1.0 for b in members}), machine_scaling=scaling
            )
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row.update({"status": f"INFEASIBLE {members}: {str(error)[:40]}"})
            return row
        spectrum = eigen_analysis(case.system.A)
        label = "+".join(map(str, members)) or "BASE"
        row[f"maxre_{label}"] = band_max_real(spectrum)
        row[f"count_{label}"] = gamma_count(spectrum)
        counts[members] = gamma_count(spectrum)
        margins[members] = band_max_real(spectrum)
        if not members:
            base = case
            row["base_abscissa"] = float(
                max(m.real for m in spectrum.modes if abs(m.value) > ZERO)
            )
        if len(members) == len(CORE):
            flagship = case
    row["base_stable"] = bool(row["base_abscissa"] < 0.0)
    reference = counts[()]
    changed = {
        len(S): max(
            (abs(counts[T] - reference) for T in counts if len(T) == len(S)), default=0
        )
        for S in SUBSETS
        if S
    }
    order = next((k for k in range(1, len(CORE) + 1) if changed.get(k, 0) > 0), None)
    row["kappa_gamma"] = order if order is not None else -1
    try:
        space = build_action_space(base, flagship, CORE)
        profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
        row["m_cl"] = min(s.margin for s in profile)
    except (np.linalg.LinAlgError, ValueError):
        row["m_cl"] = float("nan")
    row["status"] = "OK"
    return row


def main() -> int:
    pin_blas_threads()
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Is kappa piecewise constant, changing only at spectral boundaries?",
        config={
            "sweep": "excitation gain scale, multiplicative on the native fleet",
            "t_scale": T_SCALE,
            "gain_range": [float(GAINS.min()), float(GAINS.max())],
            "step": 0.025,
            "gamma": "right half plane intersected with 0.3-1.5 Hz",
        },
        workers=workers,
    )
    started = time.time()
    experiment.note(f"{len(GAINS)} continuation points on {workers} workers")
    with Pool(processes=workers, initializer=_initialise) as pool:
        rows = pool.map(_evaluate, list(GAINS), chunksize=2)

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "F2_policy_composability_order.csv", index=False)
    experiment.save_table(table, "F2_composability_continuation.csv")
    good = table[(table.status == "OK") & table.base_stable].sort_values("gain_scale")

    # every change of kappa must coincide with a band-max-real zero crossing
    labels = ["+".join(map(str, s)) or "BASE" for s in SUBSETS]
    events = []
    values = good.kappa_gamma.to_numpy()
    gains = good.gain_scale.to_numpy()
    for i in range(1, len(values)):
        if values[i] == values[i - 1]:
            continue
        crossers = [
            label
            for label in labels
            if np.sign(good.iloc[i][f"maxre_{label}"])
            != np.sign(good.iloc[i - 1][f"maxre_{label}"])
        ]
        events.append(
            {
                "gain_from": float(gains[i - 1]),
                "gain_to": float(gains[i]),
                "kappa_from": int(values[i - 1]),
                "kappa_to": int(values[i]),
                "subsets_crossing_dGamma": ";".join(crossers),
                "n_crossing": len(crossers),
                "explained": bool(crossers),
            }
        )
    event_frame = pd.DataFrame(events)
    experiment.save_table(event_frame, "F2_kappa_transitions.csv")

    runs = (good.kappa_gamma != good.kappa_gamma.shift()).cumsum().nunique()
    checks = {
        "points": int(len(table)),
        "base_stable_points": int(len(good)),
        "kappa_values_seen": sorted(set(int(v) for v in values)),
        "constant_runs": int(runs),
        "transitions": int(len(events)),
        "transitions_explained_by_a_boundary_crossing": int(
            sum(1 for e in events if e["explained"])
        ),
        "all_transitions_explained": bool(all(e["explained"] for e in events)),
        "m_cl_by_kappa": {
            int(k): float(v)
            for k, v in good.groupby("kappa_gamma").m_cl.median().items()
        },
    }
    print()
    if len(event_frame):
        print(event_frame.to_string(index=False))
    print()
    for key, value in checks.items():
        print(f"  {key:44s} {value}")
    experiment.finish(
        "THEOREM_CORROBORATED" if checks["all_transitions_explained"] else "MISMATCH",
        checks=checks,
        elapsed_s=time.time() - started,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
