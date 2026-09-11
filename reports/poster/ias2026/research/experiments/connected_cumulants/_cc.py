# ruff: noqa: E501  -- formulas in docstrings and test labels kept on one line
"""Shared definitions for the connected-cumulant extension (CC02-CC06).

Frozen inputs only:
- the ten FC18 boundaries (theta and s* read from FC18_summary.json; NOT re-located);
- the FC04 Pg-matched portfolios (10 failing, 25 stable) at the E12 census policy;
- the frozen F7B line t = 0.8515625 (points and labels from results/F7/F7B_points.csv.gz);
- the frozen 12-line holdout (holdout_lines.json).

The port machinery is a copy of FC18's (C2 common realization, T_S, M = D K,
Q = (I + D K_d)^-1 D K_o) so that nothing here imports a module that creates output
directories in a frozen run. Outputs go to results/CC/<ID>/.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

import numpy as np  # noqa: E402

import _bootstrap  # noqa: E402,F401  (puts src/ on the path)
from _f7_common import (  # noqa: E402
    LEAK,
    Theta,
    native_te,
    solve_subset,
    te_geometric_mean,
)
from _overnight import RESEARCH  # noqa: E402
from ibr_cycles.certification.binary import (  # noqa: E402
    build_common_realization,
    reduced,
)
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402

ROOT = RESEARCH
RESULTS = ROOT / "results" / "CC"
FC_RUN = ROOT / "outputs/ias2026/final_math_nonlinear_validation_20260910T231539"
FC18_SUMMARY = FC_RUN / "FC18_targeted_port_checks/FC18_summary.json"
FC04_TABLE = FC_RUN / "FC04_e14_n6_pg/FC04_n6_pg.csv"
E13_PREDICTORS = ROOT / "results/tables/E13_baseline_challenge_predictors.csv"
E12_CENSUS = ROOT / "results/tables/E12_compatibility_census_portfolios.csv"
F7B_POINTS = ROOT / "results/F7/F7B_points.csv.gz"
HOLDOUT = (
    ROOT
    / "validation_inputs/cross_tool_handoff/claude_cross_tool_handoff"
    / "dynamic_forest_ieee39_validation/dynamic_forest_ieee39_validation/holdout_lines.json"
)
PREREG = ROOT / "configs/ias2026/connected_cumulants_prereg_v1.yaml"
CORE = (30, 33, 35, 37)
PATH_T = 0.8515625
NAMES = ("g", "k", "t", "h")
STEPS = {"g": 1e-4, "k": 1e-4, "t": 1e-4, "h": 1e-4}  # FC18 steps
SCALE = np.array([1.0, 1.8, 2.5, 2.0])  # FC18 policy-box ranges (g, k, t, h)
BAND = (0.3, 1.5)


def out_dir(name: str) -> Path:
    path = RESULTS / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")


def label(members) -> str:
    return "+".join(map(str, sorted(members))) or "BASE"


# ------------------------------------------------------------ frozen machinery --
def kwargs_of(th) -> dict:
    """FC18 kwargs_of: the (g, k, t, h) policy point as solve_case keyword arguments."""

    g, k, t, h = th
    kw = {
        "converter": ConverterParameters(
            voltage_control=True, voltage_gain=g, voltage_leak=LEAK
        ),
        "machine_scaling": {"ka": k, "ta": t},
    }
    if h != 1.0:
        mean = te_geometric_mean()
        kw["machine_bus_scaling"] = {
            b: {"ta": (mean / te) ** (1.0 - h)} for b, te in native_te().items()
        }
    return kw


def port_t(jac, s):
    n = jac.fx.shape[0]
    return jac.gz + jac.gx @ np.linalg.solve(s * np.eye(n) - jac.fx, jac.fz)


class Realization:
    """FC18 Realization: C2 common realization of ``members`` at policy ``th``.

    ``extra`` passes further solve_case keywords (e.g. ``network`` for line scaling).
    ``th=None`` is the solve_case default policy (the E12 census policy).
    """

    def __init__(self, members, th=None, **extra):
        self.members = tuple(members)
        kw = {} if th is None else kwargs_of(th)
        kw.update(extra)
        self.cr = build_common_realization(self.members, **kw)
        self._jac = {}
        net = self.cr.case.dae.network
        self.ch = {
            b: [2 * net.position(b), 2 * net.position(b) + 1] for b in self.members
        }
        self.idx = [c for b in self.members for c in self.ch[b]]

    def jac(self, subset):
        """C2 vertex Jacobian of ``subset``, computed on first use.

        Identical to FC18's eager dict (same ``cr.jacobian(delta, "C2")`` call); only
        the vertices actually used are evaluated.
        """

        key = tuple(sorted(subset))
        if key not in self._jac:
            delta = {b: (1 if b in key else 0) for b in self.members}
            self._jac[key], _, _ = self.cr.jacobian(delta, "C2")
        return self._jac[key]

    def m_matrix(self, s):
        t0 = port_t(self.jac(()), s)
        kfull = np.linalg.inv(t0)
        k = kfull[np.ix_(self.idx, self.idx)]
        n = 2 * len(self.members)
        d = np.zeros((n, n), complex)
        for i, b in enumerate(self.members):
            c = self.ch[b]
            d[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = (port_t(self.jac((b,)), s) - t0)[
                np.ix_(c, c)
            ]
        return d, k

    def q_matrix(self, s):
        """Q = (I + D K_d)^-1 D K_o (zero 2x2 diagonal blocks) and the local blocks."""

        d, k = self.m_matrix(s)
        n = k.shape[0]
        kd = np.zeros_like(k)
        for i in range(0, n, 2):
            kd[i : i + 2, i : i + 2] = k[i : i + 2, i : i + 2]
        loc = np.eye(n) + d @ kd
        q = np.linalg.solve(loc, d @ (k - kd))
        return q, d, k, loc

    def cols(self, subset):
        return [2 * self.members.index(b) + c for b in subset for c in (0, 1)]

    def vertex_spectrum(self, subset):
        return np.linalg.eigvals(reduced(self.jac(subset)))


def direct_perp(members, th):
    case = solve_subset(tuple(members), Theta(*th))
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    tr = transverse_operator(case.system.A, r_x, frequency_partner(case.dae).w)
    return np.linalg.eigvals(tr.a_perp)


def band_critical(ev):
    """Rightmost transverse eigenvalue with 0.3 <= f <= 1.5 Hz (upper half plane)."""

    f = ev.imag / (2 * np.pi)
    band = ev[(f >= BAND[0]) & (f <= BAND[1])]
    return complex(band[np.argmax(band.real)])


def frozen_events() -> list[dict]:
    data = json.loads(FC18_SUMMARY.read_text(encoding="utf-8"))
    out = []
    for e in data["events"]:
        out.append(
            {
                "event": e["event"],
                "subset": tuple(int(b) for b in e["subset"].split("+")),
                "theta": tuple(e["theta"]),
                "s_star": complex(*e["s_star_realization"]),
                "freq_hz": e["freq_hz"],
                "grad_port": {k: complex(*v) for k, v in e["grad_port"].items()},
                "grad_direct": {k: complex(*v) for k, v in e["grad_direct"].items()},
            }
        )
    return out
