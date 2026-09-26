"""E35 - overnight 7. Fresh held-out operating-point Monte Carlo, N = 2000.

Seed 20260913, never used in v1 (20260909), v2A (20260910), v2B (20260911) or
v2C (20260912). No earlier sample is reused.

Every accepted sample carries the WHOLE 16-subset lattice, so the order-4
statement is tested per operating point rather than once at the nominal point:

    P1  every proper subset stable
    P2  full portfolio unstable
    P3  reconstruction through order 3 stable while the exact portfolio is not
    P4  the same inter-area family is retained
    P5  the closure margin is smaller in unstable or near-boundary cases
    P6  the RC converter repair succeeds
    P7  the synchronous condenser succeeds

Non-dispatchable draws are rejected explicitly and counted as neither stable nor
unstable.
"""

from __future__ import annotations

import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import beta, mannwhitneyu

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
from E32_tds_validation import rc_converter
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space, closure_profile
from ibr_cycles.uncertainty.sampling import sample_operating_point

NAME = "E35_MC_operating"
SEED = 20260913
N_DRAW = 2000
BAND_SAMPLES = 61
LOAD_RANGE = (0.85, 1.15)
REACTIVE_RANGE = (0.90, 1.10)
AVAILABILITY_RANGE = (0.70, 1.00)
LOAD_SCATTER = 0.05
DISPATCH_JITTER = 0.10
NEAR_BOUNDARY = 0.05
SUBSETS = [
    tuple(sorted(s))
    for k in range(len(CORE) + 1)
    for s in __import__("itertools").combinations(CORE, k)
]
NAMES = tuple(map(str, CORE))

_CONTEXT: dict = {}


def _initialise():
    pin_blas_threads()
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    _CONTEXT.update(
        network=load_network(),
        pinned=sorted(machine_labels(flagship.system.labels)),
        nominal=nominal_reference(),
        rc=rc_converter(),
    )


def _read(anchor, base, case, pinned):
    spectrum = eigen_analysis(case.system.A)
    reading = read_family(anchor, base, case, spectrum, common_labels=pinned)
    band = band_candidates(spectrum, BAND_HZ)
    rhp = sum(1 for m in band if m.real > 0.0)
    if reading is None:
        return None, rhp
    return reading, rhp


def _evaluate(task):
    index, load, reactive, availability = task
    network = _CONTEXT["network"]
    pinned = _CONTEXT["pinned"]
    row = {
        "sample": index,
        "load": load,
        "reactive_load": reactive,
        "availability": availability,
    }
    rng = np.random.default_rng(SEED + 100003 * index)
    point = sample_operating_point(
        network,
        active_load=load,
        reactive_load=reactive,
        availability=availability,
        availability_buses=CORE,
        rng=rng,
        load_scatter=LOAD_SCATTER,
        dispatch_jitter=DISPATCH_JITTER,
    )
    row["dispatchable"] = bool(point.dispatchable)
    if not point.dispatchable:
        row["status"] = "REJECTED"
        row["reason"] = "not dispatchable"
        return row
    try:
        base = solve_case(ReplacementPlan.of({}), network=point.network)
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
        row["status"] = "REJECTED"
        row["reason"] = f"base infeasible: {str(error)[:60]}"
        return row
    anchor = base_anchor(eigen_analysis(base.system.A), base, _CONTEXT["nominal"])
    if anchor is None:
        row["status"] = "REJECTED"
        row["reason"] = "TRACKING_FAILURE: no anchor"
        return row

    lattice: dict[frozenset, float] = {}
    family_sizes = []
    overlaps = []
    proper_stable = True
    flagship_case = None
    for members in SUBSETS:
        try:
            case = (
                base
                if not members
                else solve_case(
                    ReplacementPlan.of({b: 1.0 for b in members}), network=point.network
                )
            )
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row["status"] = "REJECTED"
            row["reason"] = f"subset {members} infeasible: {str(error)[:50]}"
            return row
        reading, rhp = _read(anchor, base, case, pinned)
        if reading is None:
            row["status"] = "REJECTED"
            row["reason"] = f"TRACKING_FAILURE at subset {members}"
            return row
        lattice[frozenset(map(str, members))] = reading.family.alpha
        family_sizes.append(reading.family.size)
        overlaps.append(min(reading.family.overlaps))
        if len(members) == len(CORE):
            flagship_case = case
            row["band_rhp_flagship"] = rhp
            row["alpha_IA_exact"] = reading.family.alpha
            row["freq_IA_hz"] = reading.family.frequency_worst_hz
        elif members and reading.family.alpha > 0.0:
            proper_stable = False
        if not members:
            row["alpha_base"] = reading.family.alpha

    decomposition = decompose(lambda s: lattice[frozenset(s)], NAMES)
    truncations = {k: float(decomposition.truncated(NAMES, k)) for k in range(5)}
    row.update(
        {
            "alpha_IA_le_1": truncations[1],
            "alpha_IA_le_2": truncations[2],
            "alpha_IA_le_3": truncations[3],
            "mu4": float(decomposition.terms[frozenset(NAMES)]),
            "proper_subsets_stable": bool(proper_stable),
            "worst_proper_subset_alpha": float(
                max(v for k, v in lattice.items() if 0 < len(k) < len(CORE))
            ),
            "family_size_min": int(min(family_sizes)),
            "family_size_max": int(max(family_sizes)),
            "overlap_min": float(min(overlaps)),
        }
    )
    row["portfolio_unstable"] = bool(row["alpha_IA_exact"] > 0.0)
    row["fourth_order_crossing"] = bool(
        truncations[1] < 0
        and truncations[2] < 0
        and truncations[3] < 0
        and truncations[4] > 0
    )

    try:
        space = build_action_space(base, flagship_case, CORE)
        profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
        best = min(profile, key=lambda s: s.margin)
        row.update(
            {
                "m4": best.margin,
                "freq_port_hz": best.frequency_hz,
                "cond_T0": best.condition,
                "solve_residual": best.solve_residual,
            }
        )
    except (np.linalg.LinAlgError, ValueError):
        row["m4"] = float("nan")

    for tag, plan, converter in (
        ("rc", ReplacementPlan.of({b: 1.0 for b in CORE}), _CONTEXT["rc"]),
        (
            "sc",
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE}
            ),
            None,
        ),
    ):
        try:
            case = solve_case(plan, network=point.network, converter=converter)
            reading, rhp = _read(anchor, base, case, pinned)
            row[f"alpha_{tag}"] = reading.family.alpha if reading else float("nan")
            row[f"band_rhp_{tag}"] = rhp
            row[f"{tag}_success"] = bool(rhp == 0)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row[f"alpha_{tag}"] = float("nan")
            row[f"{tag}_success"] = False
            row[f"{tag}_error"] = str(error)[:60]

    row["status"] = "ACCEPTED"
    return row


def clopper_pearson(successes, trials, level=0.95):
    if trials == 0:
        return float("nan"), 0.0, 1.0
    alpha = 1.0 - level
    low = beta.ppf(alpha / 2.0, successes, trials - successes + 1) if successes else 0.0
    high = (
        beta.ppf(1.0 - alpha / 2.0, successes + 1, trials - successes)
        if successes < trials
        else 1.0
    )
    return successes / trials, float(low), float(high)


def main() -> int:
    pin_blas_threads()
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Do the Track-A statements hold across a fresh held-out operating envelope?",
        seed=SEED,
        config={
            "n_draw": N_DRAW,
            "load_range": LOAD_RANGE,
            "reactive_range": REACTIVE_RANGE,
            "availability_range": AVAILABILITY_RANGE,
            "load_scatter": LOAD_SCATTER,
            "dispatch_jitter": DISPATCH_JITTER,
            "subsets_per_sample": len(SUBSETS),
            "seed_never_used_before": True,
        },
        workers=workers,
    )
    started = time.time()
    rng = np.random.default_rng(SEED)
    tasks = [
        (
            i,
            float(rng.uniform(*LOAD_RANGE)),
            float(rng.uniform(*REACTIVE_RANGE)),
            float(rng.uniform(*AVAILABILITY_RANGE)),
        )
        for i in range(N_DRAW)
    ]
    experiment.note(f"{N_DRAW} draws on {workers} workers, {len(SUBSETS)} subsets each")

    with Pool(processes=workers, initializer=_initialise) as pool:
        rows = pool.map(_evaluate, tasks, chunksize=4)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E35_MC_operating_2000.csv", parquet=True)
    accepted = table[table.status == "ACCEPTED"].copy()
    n = len(accepted)

    def rate(mask, label):
        successes = int(mask.sum())
        point, low, high = clopper_pearson(successes, n)
        return {
            "probability": label,
            "successes": successes,
            "trials": n,
            "estimate": point,
            "ci_low": low,
            "ci_high": high,
        }

    unstable = accepted[accepted.portfolio_unstable]
    summary = [
        rate(accepted.proper_subsets_stable, "P1 all proper subsets stable"),
        rate(accepted.portfolio_unstable, "P2 full portfolio unstable"),
        rate(accepted.fourth_order_crossing, "P3 order<=3 stable while exact unstable"),
        rate(accepted.family_size_min > 0, "P4 inter-area family retained"),
    ]
    if len(unstable):
        successes = int(unstable.rc_success.sum())
        point, low, high = clopper_pearson(successes, len(unstable))
        summary.append(
            {
                "probability": "P6 RC repair succeeds (unstable samples)",
                "successes": successes,
                "trials": len(unstable),
                "estimate": point,
                "ci_low": low,
                "ci_high": high,
            }
        )
        successes = int(unstable.sc_success.sum())
        point, low, high = clopper_pearson(successes, len(unstable))
        summary.append(
            {
                "probability": "P7 condenser succeeds (unstable samples)",
                "successes": successes,
                "trials": len(unstable),
                "estimate": point,
                "ci_low": low,
                "ci_high": high,
            }
        )
    valid = accepted[np.isfinite(accepted.m4)]
    near = valid[valid.alpha_IA_exact.abs() <= NEAR_BOUNDARY]
    far = valid[valid.alpha_IA_exact.abs() > NEAR_BOUNDARY]
    p5 = (
        float(mannwhitneyu(near.m4, far.m4, alternative="less").pvalue)
        if len(near) >= 5 and len(far) >= 5
        else float("nan")
    )
    summary.append(
        {
            "probability": "P5 closure smaller near the boundary (Mann-Whitney p)",
            "successes": len(near),
            "trials": len(valid),
            "estimate": p5,
            "ci_low": float(near.m4.median()) if len(near) else float("nan"),
            "ci_high": float(far.m4.median()) if len(far) else float("nan"),
        }
    )
    frame = pd.DataFrame(summary)
    experiment.save_table(frame, "E35_MC_operating_summary.csv")

    checks = {
        "drawn": int(len(table)),
        "accepted": n,
        "rejected": int(len(table) - n),
        "rejection_reasons": table[table.status == "REJECTED"].reason.value_counts().to_dict(),
        "tracking_failures": int(
            table.reason.fillna("").str.startswith("TRACKING_FAILURE").sum()
        ),
        "p5_mannwhitney": p5,
        "near_median_m4": float(near.m4.median()) if len(near) else None,
        "far_median_m4": float(far.m4.median()) if len(far) else None,
        "family_sizes": sorted(set(accepted.family_size_min.dropna().astype(int)))
        + sorted(set(accepted.family_size_max.dropna().astype(int))),
        "overlap_min": float(accepted.overlap_min.min()),
    }
    print()
    print(frame.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))
    print()
    for key, value in checks.items():
        print(f"  {key:24s} {value}")

    p3 = frame[frame.probability.str.startswith("P3")].iloc[0]
    verdict = "PASS" if p3.ci_low > 0.5 else ("G3_MIXED" if p3.estimate > 0.2 else "G3_FAIL")
    experiment.finish(verdict, summary=summary, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
