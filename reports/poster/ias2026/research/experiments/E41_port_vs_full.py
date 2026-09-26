"""E41 - overnight 4. Does the port reduction reproduce the full-state conclusion?

The whole theoretical apparatus of Track A lives in port space: an 8x8 action
operator standing in for a system of several hundred states. If the port
representation did not reproduce the full model's engineering conclusion over the
entire subset lattice, the theory would be an elegant description of something
other than the system.

Two questions per case:

  1. VISIBILITY. Is each inter-area family mode of the full DAE a zero of
     det T(s)? A mode confined to a device and invisible at its terminals is
     absent from the port operator by construction, and this reports which ones
     those are instead of hiding them. The test is the smallest singular value of
     T(lambda) relative to the norm of T(lambda): near zero means the full model's
     eigenvalue is a zero of the port operator.

  2. ACCURACY. Solving the port nonlinear eigenvalue problem in its own right,
     by secant iteration on the smallest eigenvalue of T(s), how far is the port
     model's own answer from the full model's?

No explicit inverse of T is formed anywhere.
"""

from __future__ import annotations

import time
from itertools import combinations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, pin_blas_threads
from _v2c_common import (
    BAND_HZ,
    CORE,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from E32_tds_validation import rc_converter
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.port_admittance import build_port_operator

NAME = "E41_port_vs_full"
VISIBLE_TOLERANCE = 1e-6
"""sigma_min(T)/||T|| below this counts the full model's mode as a port zero."""


def smallest_eigenvalue(operator, s):
    values = np.linalg.eigvals(operator.evaluate(s))
    return complex(values[int(np.argmin(np.abs(values)))])


def relative_sigma_min(operator, s):
    matrix = operator.evaluate(s)
    singular = np.linalg.svd(matrix, compute_uv=False)
    return float(singular[-1] / singular[0])


def nep_zero(operator, start, *, iterations=60, tolerance=1e-12):
    """Secant iteration on the smallest eigenvalue of T(s).

    The determinant itself is hopelessly scaled at this dimension; its smallest
    eigenvalue is not, and it vanishes at exactly the same points.
    """

    s0 = complex(start)
    s1 = s0 * (1.0 + 1e-4) + 1e-6
    f0 = smallest_eigenvalue(operator, s0)
    f1 = smallest_eigenvalue(operator, s1)
    for _ in range(iterations):
        if abs(f1 - f0) < 1e-300:
            break
        s2 = s1 - f1 * (s1 - s0) / (f1 - f0)
        if not np.isfinite(s2):
            break
        s0, f0, s1 = s1, f1, s2
        f1 = smallest_eigenvalue(operator, s1)
        if abs(s1 - s0) < tolerance * max(abs(s1), 1.0):
            break
    return s1, abs(f1)


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Does the reduced port representation reproduce the full-state conclusion?",
        config={
            "cases": "base, 15 proper subsets, flagship, RC repair, condenser 25%",
            "visible_tolerance": VISIBLE_TOLERANCE,
            "band_hz": list(BAND_HZ),
            "no_explicit_inverse": True,
        },
        workers=1,
    )
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship_plan = ReplacementPlan.of({b: 1.0 for b in CORE})
    flagship = solve_case(flagship_plan)
    pinned = sorted(machine_labels(flagship.system.labels))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)

    cases = {}
    for size in range(len(CORE) + 1):
        for members in combinations(CORE, size):
            label = "+".join(map(str, members)) or "BASE"
            cases[label] = (
                base if not members else ReplacementPlan.of({b: 1.0 for b in members})
            )
    extra = {
        "RC repair": (flagship_plan, rc_converter()),
        "condenser 25%": (
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE}
            ),
            None,
        ),
    }

    rows = []

    def evaluate(label, case):
        spectrum = eigen_analysis(case.system.A)
        reading = read_family(anchor, base, case, spectrum, common_labels=pinned)
        operator = build_port_operator(
            case.dae, case.equilibrium.x, case.equilibrium.z
        )
        band = band_candidates(spectrum, BAND_HZ)
        family_modes = reading.family.modes if reading else ()
        visible = 0
        for mode in family_modes:
            if relative_sigma_min(operator, mode.value) < VISIBLE_TOLERANCE:
                visible += 1
        worst = reading.family.worst.value if reading else complex("nan")
        zero, residual = nep_zero(operator, worst)
        band_visible = sum(
            1
            for mode in band
            if relative_sigma_min(operator, mode.value) < VISIBLE_TOLERANCE
        )
        port_rhp = 0
        for mode in band:
            if relative_sigma_min(operator, mode.value) < VISIBLE_TOLERANCE and mode.real > 0:
                port_rhp += 1
        return {
            "case": label,
            "alpha_full": float(worst.real),
            "freq_full_hz": float(abs(worst.imag) / (2 * np.pi)),
            "alpha_port": float(zero.real),
            "freq_port_hz": float(abs(zero.imag) / (2 * np.pi)),
            "eigenvalue_error": float(abs(zero - worst)),
            "alpha_error": float(abs(zero.real - worst.real)),
            "freq_error_hz": float(abs(abs(zero.imag) - abs(worst.imag)) / (2 * np.pi)),
            "nep_residual": float(residual),
            "sigma_min_at_full_mode": relative_sigma_min(operator, worst),
            "family_size": len(family_modes),
            "family_modes_visible_at_ports": visible,
            "band_modes": len(band),
            "band_modes_visible_at_ports": band_visible,
            "band_rhp_full": sum(1 for m in band if m.real > 0),
            "band_rhp_port_visible": port_rhp,
            "n_states_full": int(case.system.n),
            "n_states_port": int(operator.dimension),
        }

    for label, item in cases.items():
        try:
            case = item if not isinstance(item, ReplacementPlan) else solve_case(item)
            rows.append(evaluate(label, case))
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            rows.append({"case": label, "error": str(error)[:100]})
        print(f"  {label:14s} done", flush=True)
    for label, (plan, converter) in extra.items():
        try:
            case = solve_case(plan, converter=converter)
            rows.append(evaluate(label, case))
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            rows.append({"case": label, "error": str(error)[:100]})
        print(f"  {label:14s} done", flush=True)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E41_port_vs_full.csv")
    good = table[table.get("alpha_full").notna()] if "alpha_full" in table else table

    checks = {
        "cases": int(len(good)),
        "max_eigenvalue_error": float(good.eigenvalue_error.max()),
        "median_eigenvalue_error": float(good.eigenvalue_error.median()),
        "max_alpha_error": float(good.alpha_error.max()),
        "max_freq_error_hz": float(good.freq_error_hz.max()),
        "max_sigma_min_at_full_mode": float(good.sigma_min_at_full_mode.max()),
        "all_family_modes_visible": bool(
            (good.family_modes_visible_at_ports == good.family_size).all()
        ),
        "invisible_band_modes_total": int(
            (good.band_modes - good.band_modes_visible_at_ports).sum()
        ),
        "rhp_count_agrees_everywhere": bool(
            (good.band_rhp_full == good.band_rhp_port_visible).all()
        ),
        "state_reduction": f"{int(good.n_states_full.max())} -> {int(good.n_states_port.iloc[0])}",
    }
    verdict = (
        "PASS"
        if checks["rhp_count_agrees_everywhere"]
        and checks["all_family_modes_visible"]
        and checks["max_eigenvalue_error"] < 1e-6
        else "G6_MIXED"
    )

    print()
    print(
        good[
            [
                "case",
                "alpha_full",
                "alpha_port",
                "eigenvalue_error",
                "sigma_min_at_full_mode",
                "family_size",
                "family_modes_visible_at_ports",
                "band_rhp_full",
                "band_rhp_port_visible",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:11.4e}")
    )
    print()
    for key, value in checks.items():
        print(f"  {key:34s} {value}")
    experiment.finish(verdict, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
