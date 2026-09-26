"""E37 - overnight 9. Machine and synchronous-support Monte Carlo, N = 1000.

Seed 20260915, never used before. The operating point and the converter are the
frozen nominal ones; the SYNCHRONOUS machines are uncertain. Track A claims the
failure is a loss of synchronous dynamic support, so the machines' own parameters
are the uncertainty that most directly threatens the claim.

    inertia H                      +- 20 %
    transient reactances Xd', Xq'  +- 10 %
    AVR gain and time constant     +- 15 %
    PSS gain and time constants    +- 20 %

Machine damping D is ZERO throughout the source case, so a +-25 % perturbation of
it is identically no perturbation. It is reported as not applicable rather than
sampled and presented as though it had been varied.

Every draw carries the whole 16-subset lattice, and the base case is rebuilt per
draw because the machines themselves are what changed.
"""

from __future__ import annotations

import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import beta

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, choose_workers, pin_blas_threads
from _v2c_common import (
    BAND_HZ,
    CORE,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.port_admittance import build_action_space, closure_profile

NAME = "E37_MC_machine"
SEED = 20260915
N_DRAW = 1000
BAND_SAMPLES = 61
SUBSETS = [
    tuple(sorted(s)) for k in range(len(CORE) + 1) for s in combinations(CORE, k)
]
NAMES = tuple(map(str, CORE))
RANGES = {
    "inertia": (0.80, 1.20),
    "transient_reactance": (0.90, 1.10),
    "avr_gain": (0.85, 1.15),
    "avr_time_constant": (0.85, 1.15),
    "pss_gain": (0.80, 1.20),
    "pss_time_constants": (0.80, 1.20),
}
DAMPING_NOTE = "machine damping D is zero in the source case; not sampled"

_CONTEXT: dict = {}


def _initialise():
    pin_blas_threads()
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    _CONTEXT.update(
        pinned=sorted(machine_labels(flagship.system.labels)),
        nominal=nominal_reference(),
    )


def build_scaling(draw):
    return {
        "m": draw["inertia"],
        "xd1": draw["transient_reactance"],
        "xq1": draw["transient_reactance"],
        "ka": draw["avr_gain"],
        "ta": draw["avr_time_constant"],
        "pss_gain": draw["pss_gain"],
        "pss_washout": draw["pss_time_constants"],
        "pss_wash_lag": draw["pss_time_constants"],
        "pss_lag": draw["pss_time_constants"],
    }


def _evaluate(task):
    index, draw = task
    pinned = _CONTEXT["pinned"]
    scaling = build_scaling(draw)
    row = {"sample": index, **draw}
    try:
        base = solve_case(ReplacementPlan.of({}), machine_scaling=scaling)
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
        row["status"] = "REJECTED"
        row["reason"] = f"base: {str(error)[:60]}"
        return row
    anchor = base_anchor(eigen_analysis(base.system.A), base, _CONTEXT["nominal"])
    if anchor is None:
        row["status"] = "REJECTED"
        row["reason"] = "TRACKING_FAILURE: no anchor"
        return row
    row["alpha_base"] = float(anchor.real)

    lattice: dict[frozenset, float] = {}
    sizes, overlaps = [], []
    proper_stable = True
    flagship_case = None
    for members in SUBSETS:
        try:
            case = (
                base
                if not members
                else solve_case(
                    ReplacementPlan.of({b: 1.0 for b in members}),
                    machine_scaling=scaling,
                )
            )
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row["status"] = "REJECTED"
            row["reason"] = f"subset {members}: {str(error)[:55]}"
            return row
        spectrum = eigen_analysis(case.system.A)
        reading = read_family(anchor, base, case, spectrum, common_labels=pinned)
        if reading is None:
            row["status"] = "REJECTED"
            row["reason"] = f"TRACKING_FAILURE at {members}"
            return row
        lattice[frozenset(map(str, members))] = reading.family.alpha
        sizes.append(reading.family.size)
        overlaps.append(min(reading.family.overlaps))
        if len(members) == len(CORE):
            flagship_case = case
            row["alpha_IA_exact"] = reading.family.alpha
            row["freq_IA_hz"] = reading.family.frequency_worst_hz
            row["band_rhp_flagship"] = sum(
                1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0
            )
        elif members and reading.family.alpha > 0.0:
            proper_stable = False

    decomposition = decompose(lambda s: lattice[frozenset(s)], NAMES)
    truncations = {k: float(decomposition.truncated(NAMES, k)) for k in range(5)}
    row.update(
        {
            "status": "ACCEPTED",
            "alpha_IA_le_1": truncations[1],
            "alpha_IA_le_2": truncations[2],
            "alpha_IA_le_3": truncations[3],
            "mu4": float(decomposition.terms[frozenset(NAMES)]),
            "proper_subsets_stable": bool(proper_stable),
            "portfolio_unstable": bool(row["alpha_IA_exact"] > 0.0),
            "fourth_order_crossing": bool(
                truncations[1] < 0
                and truncations[2] < 0
                and truncations[3] < 0
                and truncations[4] > 0
            ),
            "family_size_min": int(min(sizes)),
            "family_size_max": int(max(sizes)),
            "overlap_min": float(min(overlaps)),
        }
    )
    try:
        space = build_action_space(base, flagship_case, CORE)
        profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
        best = min(profile, key=lambda s: s.margin)
        row["m4"] = best.margin
        row["freq_port_hz"] = best.frequency_hz
    except (np.linalg.LinAlgError, ValueError):
        row["m4"] = float("nan")

    for tag, plan in (
        (
            "sc",
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE}
            ),
        ),
    ):
        try:
            case = solve_case(plan, machine_scaling=scaling)
            spectrum = eigen_analysis(case.system.A)
            reading = read_family(anchor, base, case, spectrum, common_labels=pinned)
            row[f"alpha_{tag}"] = reading.family.alpha if reading else float("nan")
            row[f"{tag}_success"] = bool(
                sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0) == 0
            )
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
            row[f"{tag}_success"] = False
    return row


def clopper_pearson(successes, trials, level=0.95):
    if trials == 0:
        return float("nan"), 0.0, 1.0
    a = 1.0 - level
    low = beta.ppf(a / 2.0, successes, trials - successes + 1) if successes else 0.0
    high = (
        beta.ppf(1.0 - a / 2.0, successes + 1, trials - successes)
        if successes < trials
        else 1.0
    )
    return successes / trials, float(low), float(high)


def main() -> int:
    pin_blas_threads()
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Is the result robust to reasonable synchronous-machine uncertainty?",
        seed=SEED,
        config={
            "n_draw": N_DRAW,
            "ranges": RANGES,
            "damping": DAMPING_NOTE,
            "coordinates": "physical machine parameters",
        },
        workers=workers,
    )
    started = time.time()
    rng = np.random.default_rng(SEED)
    tasks = [
        (i, {name: float(rng.uniform(*bounds)) for name, bounds in RANGES.items()})
        for i in range(N_DRAW)
    ]
    experiment.note(f"{N_DRAW} machine-parameter draws on {workers} workers")

    with Pool(processes=workers, initializer=_initialise) as pool:
        rows = pool.map(_evaluate, tasks, chunksize=4)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E37_MC_machine.csv", parquet=True)
    accepted = table[table.status == "ACCEPTED"]
    n = len(accepted)

    def rate(mask, label, trials=None):
        successes = int(mask.sum())
        total = trials if trials is not None else n
        point, low, high = clopper_pearson(successes, total)
        return {
            "question": label,
            "successes": successes,
            "trials": total,
            "estimate": point,
            "ci_low": low,
            "ci_high": high,
        }

    unstable = accepted[accepted.portfolio_unstable]
    summary = [
        rate(accepted.proper_subsets_stable, "all proper subsets stable"),
        rate(accepted.portfolio_unstable, "full portfolio unstable"),
        rate(accepted.fourth_order_crossing, "genuine order-4 crossing"),
        rate(accepted.family_size_min > 0, "inter-area family retained"),
    ]
    if len(unstable):
        summary.append(
            rate(
                unstable.sc_success,
                "condenser repair succeeds (unstable draws)",
                trials=len(unstable),
            )
        )
    frame = pd.DataFrame(summary)
    experiment.save_table(frame, "E37_MC_machine_summary.csv")

    checks = {
        "drawn": int(len(table)),
        "accepted": n,
        "rejected": int(len(table) - n),
        "damping": DAMPING_NOTE,
        "alpha_exact_min": float(accepted.alpha_IA_exact.min()),
        "alpha_exact_max": float(accepted.alpha_IA_exact.max()),
        "alpha_base_range": [
            float(accepted.alpha_base.min()),
            float(accepted.alpha_base.max()),
        ],
        "m4_median": float(accepted.m4.median()),
        "overlap_min": float(accepted.overlap_min.min()),
        "family_sizes": sorted(set(accepted.family_size_max.dropna().astype(int))),
    }
    print()
    print(frame.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))
    print()
    for key, value in checks.items():
        print(f"  {key:24s} {value}")
    order4 = frame[frame.question == "genuine order-4 crossing"].iloc[0]
    verdict = (
        "PASS" if order4.ci_low > 0.5 else ("MIXED" if order4.estimate > 0.2 else "FAIL")
    )
    experiment.finish(
        verdict, summary=summary, checks=checks, elapsed_s=time.time() - started
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
