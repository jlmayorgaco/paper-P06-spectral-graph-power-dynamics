"""The port-core all-subset count equals the full eigenproblem count.

One base operator, one replaced operator, principal minors of I + M, winding on
the boundary of the band region plus the open-loop device-pole correction.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ibr_cycles.models.ieee39_case import ReplacementPlan, build_dae, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.port_core import PortCore, all_subsets

CORE = (30, 33, 35, 37)
LO, HI = 2 * math.pi * 0.3, 2 * math.pi * 1.5


def _in_gamma(v):
    f = np.abs(v.imag) / (2 * math.pi)
    return (v.real > 0) & (v.imag >= 0) & (np.abs(v) > 1e-3) & (f >= 0.3) & (f <= 1.5)


@pytest.mark.parametrize("gain,k", [(0.0, 1.0), (0.036, 1.425), (0.11, 1.85)])
def test_port_core_counts_match_eigenvalue_counts(gain, k):
    conv = ConverterParameters(
        voltage_control=True, voltage_gain=gain, voltage_leak=0.05
    )
    scaling = {"ka": k, "ta": 1.5}
    base = build_dae(ReplacementPlan.of({}), converter=conv, machine_scaling=scaling)
    repl = build_dae(
        ReplacementPlan.of({b: 1.0 for b in CORE}),
        converter=conv,
        machine_scaling=scaling,
    )
    core = PortCore.from_daes(*base, *repl, CORE)
    subsets = all_subsets(CORE)
    counts, info = core.counts(subsets, omega_lo=LO, omega_hi=HI)
    for s in subsets:
        case = solve_case(
            ReplacementPlan.of({b: 1.0 for b in s}),
            converter=conv,
            machine_scaling=scaling,
        )
        direct = int(np.count_nonzero(_in_gamma(np.linalg.eigvals(case.system.A))))
        assert counts[s]["delta_N"] == direct, s
        assert counts[s]["winding_residual"] < 1e-6
