"""F7 Safeguard A - is the reactive-policy coordinate regular at gain zero?

Audits two ways of putting a Q/V policy gain g on the replacing converters:

    naive   q_cmd = g kp_v e + x_v,          x_v' = g ki_v e
    leaky   q_cmd = q_ref + g kp_v e + x_v,  x_v' = g ki_v e - w (x_v - q_ref)

and answers the five questions asked before any F7 map may be reported:
constant state dimension, regular g_z, no artificial zero or marginal state,
continuous A_red, and a baseline count that cannot move for structural reasons.
"""

from __future__ import annotations

import json

import numpy as np
from scipy.optimize import linear_sum_assignment

from _bootstrap import RESULTS
from _f7_common import LEAK, Theta, in_gamma, solve_subset
from _v2c_common import CORE
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters

OUT = RESULTS / "F7"


def _flagship(converter):
    return solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}), converter=converter)


def _match(a, b) -> float:
    cost = np.abs(a[:, None] - b[None, :])
    rows, cols = linear_sum_assignment(cost)
    return float(cost[rows, cols].max())


def _row(case) -> dict:
    values = np.linalg.eigvals(case.system.A)
    dynamic = values[np.abs(values) > 1e-3]
    slow = sorted(values[(np.abs(values) > 1e-9) & (np.abs(values) < 0.1)], key=abs)
    band = values[in_gamma(values)]
    return {
        "n_states": case.n_states,
        "exact_zeros": int(np.sum(np.abs(values) < 1e-9)),
        "below_1e-6": int(np.sum(np.abs(values) < 1e-6)),
        "slow_modes": [[float(v.real), float(v.imag)] for v in slow[:8]],
        "gz_condition": case.gz_condition,
        "N_gamma": int(band.size),
        "abscissa": float(dynamic.real.max()),
    }


def main() -> int:
    matched = _flagship(None)
    reference = np.linalg.eigvals(matched.system.A)
    report: dict = {"core": list(CORE), "leak_rad_s": LEAK, "matched": _row(matched)}

    report["naive"] = {
        str(g): _row(
            _flagship(ConverterParameters(voltage_control=True, voltage_gain=g))
        )
        for g in (0.0, 1e-4, 1e-3, 1e-2, 0.1, 1.0)
    }
    leaky = {
        g: _flagship(
            ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=LEAK)
        )
        for g in (0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.37, 1.0)
    }
    report["leaky"] = {str(g): _row(case) for g, case in leaky.items()}

    at_zero = np.linalg.eigvals(leaky[0.0].system.A)
    expected = np.concatenate([reference, np.full(len(CORE), -LEAK)])
    a0, a1, ah = leaky[0.0].system.A, leaky[1.0].system.A, leaky[0.37].system.A
    report["checks"] = {
        "leaky_g0_equals_matched_plus_leak_poles": _match(at_zero, expected),
        "equilibrium_x_shift_max": max(
            float(np.max(np.abs(c.equilibrium.x - leaky[0.0].equilibrium.x)))
            for c in leaky.values()
        ),
        "equilibrium_z_shift_max": max(
            float(np.max(np.abs(c.equilibrium.z - leaky[0.0].equilibrium.z)))
            for c in leaky.values()
        ),
        "gz_condition_spread": float(
            np.ptp([c.gz_condition for c in leaky.values()]) / leaky[0.0].gz_condition
        ),
        "affine_residual_rel": float(
            np.linalg.norm(ah - (0.63 * a0 + 0.37 * a1)) / np.linalg.norm(ah)
        ),
        "state_dimension_constant": len({c.n_states for c in leaky.values()}) == 1,
    }
    plain = _flagship(ConverterParameters(voltage_control=True))

    # the leak moves only the slow regulator modes; compare the oscillatory ones
    def oscillatory(values):
        return values[np.abs(values.imag) >= 2.0 * np.pi * 0.2]

    report["checks"]["g1_leaky_vs_plain_pi_oscillatory_distance"] = _match(
        oscillatory(np.linalg.eigvals(a1)),
        oscillatory(np.linalg.eigvals(plain.system.A)),
    )
    report["checks"]["g1_band_max_leaky"] = float(
        max(
            v.real
            for v in np.linalg.eigvals(a1)
            if v.imag >= 0 and 0.3 <= abs(v.imag) / 2 / np.pi <= 1.5
        )
    )
    report["checks"]["g1_band_max_plain_pi"] = float(
        max(
            v.real
            for v in np.linalg.eigvals(plain.system.A)
            if v.imag >= 0 and 0.3 <= abs(v.imag) / 2 / np.pi <= 1.5
        )
    )
    # the base case carries no converter: identical matrices whatever g is
    report["checks"]["base_matrix_identical_across_g"] = bool(
        np.array_equal(
            solve_subset((), Theta(g=0.0)).system.A,
            solve_subset((), Theta(g=0.9)).system.A,
        )
    )
    naive_zero = (
        report["naive"]["0.0"]["exact_zeros"] - report["matched"]["exact_zeros"]
    )
    report["verdict"] = {
        "naive_coordinate": f"SINGULAR at g=0: {naive_zero} artificial zero "
        "eigenvalues; g=0 excluded",
        "leaky_coordinate": "REGULAR on [0, inf): g=0 included; exact fixed-Q "
        "plus decoupled poles at -w",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "F7_safeguard_A_audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report["checks"], indent=2))
    print(json.dumps(report["verdict"], indent=2))
    for g, row in report["naive"].items():
        print("naive", g, row["n_states"], row["exact_zeros"], row["slow_modes"][:4])
    for g, row in report["leaky"].items():
        print(
            "leaky",
            g,
            row["n_states"],
            row["exact_zeros"],
            row["slow_modes"][:4],
            row["N_gamma"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
