"""E36 - overnight 8. Fresh held-out controller Monte Carlo, N = 1500.

Seed 20260914, never used before. The operating point is the frozen nominal one;
only the converter is uncertain, so this separates controller uncertainty from
operating-point uncertainty rather than mixing them.

Coordinates are PHYSICAL, not raw gains. Sampling kp and ki independently would
wander outside the region where the loops mean anything; sampling a natural
frequency, a damping ratio and three bandwidths keeps every draw a controller a
person could commission.

    PLL natural frequency        nominal +- 20 %
    PLL damping ratio            0.50 to 1.00, a frozen plausible range
    outer active bandwidth       +- 25 %
    outer reactive bandwidth     +- 25 %
    measurement bandwidth        +- 25 %
    current-loop bandwidth       +- 15 %

Every draw carries the whole 16-subset lattice.
"""

from __future__ import annotations

import time
from dataclasses import replace
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
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.port_admittance import build_action_space, closure_profile

NAME = "E36_MC_controller"
SEED = 20260914
N_DRAW = 1500
BAND_SAMPLES = 61
SUBSETS = [
    tuple(sorted(s)) for k in range(len(CORE) + 1) for s in combinations(CORE, k)
]
NAMES = tuple(map(str, CORE))

_BASE = ConverterParameters()
NOMINAL_WN = float(np.sqrt(_BASE.ki_pll))
NOMINAL_ZETA = float(_BASE.kp_pll / (2.0 * np.sqrt(_BASE.ki_pll)))
RANGES = {
    "pll_wn": (0.80, 1.20),
    "pll_zeta_absolute": (0.50, 1.00),
    "p_bandwidth": (0.75, 1.25),
    "q_bandwidth": (0.75, 1.25),
    "measurement_bandwidth": (0.75, 1.25),
    "current_bandwidth": (0.85, 1.15),
}

_CONTEXT: dict = {}


def _initialise():
    pin_blas_threads()
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    base = solve_case(ReplacementPlan.of({}))
    _CONTEXT.update(
        base=base,
        pinned=sorted(machine_labels(flagship.system.labels)),
        anchor=base_anchor(eigen_analysis(base.system.A), base, nominal_reference()),
    )


def build_converter(draw):
    wn = NOMINAL_WN * draw["pll_wn"]
    pll = ConverterParameters.pll_from_design(float(wn), float(draw["pll_zeta"]))
    return replace(
        _BASE,
        kp_pll=pll["kp_pll"],
        ki_pll=pll["ki_pll"],
        kp_p=_BASE.kp_p * draw["p_bandwidth"],
        ki_p=_BASE.ki_p * draw["p_bandwidth"],
        kp_q=_BASE.kp_q * draw["q_bandwidth"],
        ki_q=_BASE.ki_q * draw["q_bandwidth"],
        tau_p=_BASE.tau_p / draw["measurement_bandwidth"],
        kp_i=_BASE.kp_i * draw["current_bandwidth"],
        ki_i=_BASE.ki_i * draw["current_bandwidth"],
    )


def _evaluate(task):
    index, draw = task
    base = _CONTEXT["base"]
    anchor = _CONTEXT["anchor"]
    pinned = _CONTEXT["pinned"]
    converter = build_converter(draw)
    row = {"sample": index, **draw}

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
                    ReplacementPlan.of({b: 1.0 for b in members}), converter=converter
                )
            )
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row["status"] = "REJECTED"
            row["reason"] = f"subset {members}: {str(error)[:60]}"
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

    # repair robustness: the condenser, which does not depend on the converter
    try:
        case = solve_case(
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE}
            ),
            converter=converter,
        )
        spectrum = eigen_analysis(case.system.A)
        reading = read_family(anchor, base, case, spectrum, common_labels=pinned)
        row["alpha_sc"] = reading.family.alpha if reading else float("nan")
        row["sc_success"] = bool(
            sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0) == 0
        )
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
        row["sc_success"] = False
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
        question="Does the order-4 behaviour survive controller uncertainty?",
        seed=SEED,
        config={"n_draw": N_DRAW, "ranges": RANGES, "coordinates": "physical"},
        workers=workers,
    )
    started = time.time()
    rng = np.random.default_rng(SEED)
    tasks = []
    for i in range(N_DRAW):
        draw = {
            "pll_wn": float(rng.uniform(*RANGES["pll_wn"])),
            "pll_zeta": float(rng.uniform(*RANGES["pll_zeta_absolute"])),
            "p_bandwidth": float(rng.uniform(*RANGES["p_bandwidth"])),
            "q_bandwidth": float(rng.uniform(*RANGES["q_bandwidth"])),
            "measurement_bandwidth": float(
                rng.uniform(*RANGES["measurement_bandwidth"])
            ),
            "current_bandwidth": float(rng.uniform(*RANGES["current_bandwidth"])),
        }
        tasks.append((i, draw))
    experiment.note(f"{N_DRAW} controller draws on {workers} workers")

    with Pool(processes=workers, initializer=_initialise) as pool:
        rows = pool.map(_evaluate, tasks, chunksize=4)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E36_MC_controller.csv", parquet=True)
    accepted = table[table.status == "ACCEPTED"]
    n = len(accepted)

    def rate(mask, label, trials=None, subset=None):
        block = accepted if subset is None else subset
        successes = int(mask.sum())
        total = trials if trials is not None else len(block)
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
                subset=unstable,
            )
        )
    frame = pd.DataFrame(summary)
    experiment.save_table(frame, "E36_MC_controller_summary.csv")

    checks = {
        "drawn": int(len(table)),
        "accepted": n,
        "rejected": int(len(table) - n),
        "alpha_exact_min": float(accepted.alpha_IA_exact.min()),
        "alpha_exact_max": float(accepted.alpha_IA_exact.max()),
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
    verdict = "PASS" if order4.ci_low > 0.5 else ("MIXED" if order4.estimate > 0.2 else "FAIL")
    experiment.finish(verdict, summary=summary, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
