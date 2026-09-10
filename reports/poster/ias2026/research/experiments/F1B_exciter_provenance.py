"""F1B - excitation provenance audit. Is the order-4 phenomenon an artifact?

F1 established that the internal simulator and ANDES agree once both use the same
first-order AVR. It did NOT establish that the inherited parameters are a
defensible benchmark. They were obtained by pairing the IEEEX1 amplifier gain
``KA = 10.1`` with the IEEEX1 exciter time constant ``TE``, from two different
blocks, so nothing yet licenses calling them representative.

This maps the phenomenon over a physically admissible excitation region and asks
whether it occupies a finite plausible area or only the inherited point.

Documented ranges, taken from the ANDES SEXS model card, which is the reference
implementation of the first-order excitation family used here:

    K    default 20,  vrange [20, 100]
    TE   default 1,   vrange [0, 0.5]
    design note in the source:  5 <= K * TA / TB <= 15

With the lead-lag collapsed (TA/TB = 1) the SEXS transfer function is exactly
``K / (1 + s TE)``, the internal AVR, and the design product ``K * TA/TB``
reduces to ``K``.

The branch-independent order used here is the one F2 formalises:

    Gamma   = right half plane intersected with the frozen 0.3-1.5 Hz band
    N(S)    = number of modes of the replaced system inside Gamma
    kappa   = min { |S| : N(S) != N(empty) },  infinity if no subset qualifies

No modal labels are involved.
"""

from __future__ import annotations

import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _overnight import Experiment, choose_workers, pin_blas_threads
from _v2c_common import BAND_HZ, CORE, base_anchor, machine_labels, nominal_reference, read_family
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.port_admittance import build_action_space, closure_profile

NAME = "F1B_exciter_provenance"
ZERO = 1e-3
BAND_SAMPLES = 61
SUBSETS = [tuple(sorted(s)) for k in range(len(CORE) + 1) for s in combinations(CORE, k)]

#: documented SEXS ranges, for admissibility labelling only
SEXS_K_RANGE = (20.0, 100.0)
SEXS_T_RANGE = (0.0, 0.5)
SEXS_DESIGN_PRODUCT = (5.0, 15.0)

K_GRID = (5.0, 8.0, 10.1, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0)
T_GRID = (0.02, 0.05, 0.10, 0.25, 0.40, 0.50, 1.00)
STABILIZER = {"pss on": None, "pss off": {"pss": 0.0}}

#: Multiplicative sweep around the case's OWN heterogeneous excitation data.
#: Flattening every machine to a common time constant is itself a model change:
#: at K = 10.1 with a uniform T = 0.25 s the flagship is STABLE, while the case's
#: heterogeneous TE in 0.25-0.50 s gives kappa = 4. A real fleet has diverse
#: exciters, so the multiplicative sweep is the physically meaningful one and the
#: uniform grid is reported beside it as the comparison with documented ranges.
K_SCALE = (0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 4.00)
T_SCALE = (0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00)

_CONTEXT: dict = {}


def _initialise():
    pin_blas_threads()
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    _CONTEXT.update(
        pinned=sorted(machine_labels(flagship.system.labels)),
        nominal=nominal_reference(),
    )


def gamma_count(spectrum) -> int:
    """Modes inside Gamma: right half plane, inside the frozen band."""

    return sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0.0)


def _evaluate(task):
    mode, gain, constant, stabilizer = task
    if mode == "absolute":
        overrides, scaling = {"ka": float(gain), "ta": float(constant)}, None
        effective_k, effective_t = gain, constant
    else:
        overrides, scaling = None, {"ka": float(gain), "ta": float(constant)}
        effective_k = 10.1 * gain          # every machine carries KA = 10.1
        effective_t = float("nan")         # heterogeneous by construction
    services = STABILIZER[stabilizer]
    row = {
        "mode": mode,
        "K_exc": gain,
        "T_exc": constant,
        "effective_K": effective_k,
        "stabilizer": stabilizer,
        "K_in_documented_range": bool(
            SEXS_K_RANGE[0] <= effective_k <= SEXS_K_RANGE[1]
        ),
        "T_in_documented_range": bool(
            SEXS_T_RANGE[0] < constant <= SEXS_T_RANGE[1]
            if mode == "absolute"
            else constant * 0.50 <= SEXS_T_RANGE[1]
        ),
        "design_product_ok": bool(
            SEXS_DESIGN_PRODUCT[0] <= effective_k <= SEXS_DESIGN_PRODUCT[1]
        ),
    }
    try:
        base = solve_case(
            ReplacementPlan.of({}),
            machine_overrides=overrides,
            machine_scaling=scaling,
            machine_services=services,
        )
    except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
        row.update({"status": "BASE_INFEASIBLE", "reason": str(error)[:60]})
        return row
    base_spectrum = eigen_analysis(base.system.A)
    base_abscissa = float(
        max(m.real for m in base_spectrum.modes if abs(m.value) > ZERO)
    )
    base_gamma = gamma_count(base_spectrum)
    row.update(
        {
            "base_abscissa": base_abscissa,
            "base_gamma_count": base_gamma,
            "base_stable": bool(base_abscissa < 0.0),
        }
    )
    if base_abscissa >= 0.0:
        # F1B item 6: reject regions where the BASE system is unstable. They are
        # recorded so the rejected area is visible, and excluded from every
        # statement about the phenomenon.
        row["status"] = "BASE_UNSTABLE"
        return row

    anchor = base_anchor(base_spectrum, base, _CONTEXT["nominal"])
    counts: dict[int, int] = {}
    flagship_case = None
    for members in SUBSETS:
        if not members:
            counts[0] = base_gamma
            continue
        try:
            case = solve_case(
                ReplacementPlan.of({b: 1.0 for b in members}),
                machine_overrides=overrides,
                machine_scaling=scaling,
                machine_services=services,
            )
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row.update({"status": "SUBSET_INFEASIBLE", "reason": str(error)[:60]})
            return row
        spectrum = eigen_analysis(case.system.A)
        count = gamma_count(spectrum)
        size = len(members)
        counts[size] = max(counts.get(size, 0), abs(count - base_gamma))
        if size == len(CORE):
            flagship_case = case
            row["flagship_gamma_count"] = count
            row["flagship_abscissa"] = float(
                max(m.real for m in spectrum.modes if abs(m.value) > ZERO)
            )
            reading = read_family(
                anchor, base, case, spectrum, common_labels=_CONTEXT["pinned"]
            )
            row["flagship_alpha_IA"] = reading.family.alpha if reading else float("nan")
            row["flagship_freq_IA"] = (
                reading.family.frequency_worst_hz if reading else float("nan")
            )

    order = next((k for k in range(1, len(CORE) + 1) if counts.get(k, 0) > 0), None)
    row["kappa_gamma"] = order if order is not None else -1  # -1 encodes infinity
    row["kappa_is_infinite"] = order is None

    if flagship_case is not None:
        try:
            space = build_action_space(base, flagship_case, CORE)
            profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
            best = min(profile, key=lambda s: s.margin)
            row["m_cl"] = best.margin
            row["freq_port_hz"] = best.frequency_hz
        except (np.linalg.LinAlgError, ValueError):
            row["m_cl"] = float("nan")
    row["status"] = "OK"
    return row


def main() -> int:
    pin_blas_threads()
    workers = choose_workers()
    experiment = Experiment(
        name=NAME,
        question="Does the order-4 phenomenon occupy a physically credible excitation region?",
        config={
            "K_grid": list(K_GRID),
            "T_grid": list(T_GRID),
            "stabilizer_settings": list(STABILIZER),
            "documented_K_range": list(SEXS_K_RANGE),
            "documented_T_range": list(SEXS_T_RANGE),
            "documented_design_product": list(SEXS_DESIGN_PRODUCT),
            "provenance_source": "ANDES SEXS model card vrange metadata",
            "inherited_point": {"K": 10.1, "T": 0.25},
            "kappa_definition": "min |S| with a change in the right-half-plane count inside the frozen band",
            "model": "configs/ieee39_harmonized_dynamic_model.yaml",
        },
        workers=workers,
    )
    started = time.time()
    tasks = [
        ("absolute", gain, constant, stabilizer)
        for stabilizer in STABILIZER
        for gain in K_GRID
        for constant in T_GRID
    ] + [
        ("scaled", gain, constant, stabilizer)
        for stabilizer in STABILIZER
        for gain in K_SCALE
        for constant in T_SCALE
    ]
    experiment.note(f"{len(tasks)} excitation points on {workers} workers")

    with Pool(processes=workers, initializer=_initialise) as pool:
        rows = pool.map(_evaluate, tasks, chunksize=2)

    table = pd.DataFrame(rows)
    target = RESULTS / "F1B_exciter_parameter_map.csv"
    table.to_csv(target, index=False)
    experiment.save_table(table, "F1B_exciter_parameter_map.csv")

    admissible = table[table.status == "OK"]
    checks = {"points": int(len(table))}
    for stabilizer in STABILIZER:
        block = table[(table.stabilizer == stabilizer) & (table["mode"] == "scaled")]
        good = block[block.status == "OK"]
        checks[stabilizer] = {
            "evaluated": int(len(block)),
            "base_unstable_rejected": int((block.status == "BASE_UNSTABLE").sum()),
            "admissible": int(len(good)),
            "kappa_distribution": {
                int(k): int(v) for k, v in good.kappa_gamma.value_counts().items()
            },
            "kappa4_points": int((good.kappa_gamma == 4).sum()),
            "kappa4_inside_documented_K": int(
                ((good.kappa_gamma == 4) & good.K_in_documented_range).sum()
            ),
            "kappa4_K_range": (
                [
                    float(good[good.kappa_gamma == 4].K_exc.min()),
                    float(good[good.kappa_gamma == 4].K_exc.max()),
                ]
                if (good.kappa_gamma == 4).any()
                else None
            ),
            "kappa4_T_range": (
                [
                    float(good[good.kappa_gamma == 4].T_exc.min()),
                    float(good[good.kappa_gamma == 4].T_exc.max()),
                ]
                if (good.kappa_gamma == 4).any()
                else None
            ),
        }

    print()
    for mode in ("scaled", "absolute"):
      for stabilizer in STABILIZER:
        good = table[
            (table.stabilizer == stabilizer)
            & (table.status == "OK")
            & (table["mode"] == mode)
        ]
        if good.empty:
            continue
        grid = good.pivot_table(
            index="K_exc", columns="T_exc", values="kappa_gamma", aggfunc="first"
        )
        print(f"=== kappa_Gamma, {mode} sweep, {stabilizer} ===")
        print("    (-1 means no destabilizing subset; blank means base unstable)")
        print(grid.to_string(float_format=lambda v: f"{v:4.0f}"))
        print()
    for key, value in checks.items():
        print(f"  {key}: {value}")
    experiment.finish("REPORTED", checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
