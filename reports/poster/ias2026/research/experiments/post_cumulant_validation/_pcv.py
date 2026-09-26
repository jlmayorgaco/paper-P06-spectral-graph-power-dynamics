# ruff: noqa: E501  -- definitions quoted from the preregistration kept on one line
"""Shared definitions for the post-cumulant validation campaign (PCV).

Preregistration: docs/20260911_POST_CUMULANT_VALIDATION_PREREG.md (commit 5d0b1986).

- B0 truth: FC01 part-3 direct path (solve_subset, physical_matrices h/2h error,
  transverse operator, four-state classifier), recomputed.
- B3: F10 B5 first-order modal sensitivity (1 % partial replacement), with the full
  theta including h (F10 ignored h; declared correction).
- B8/B9: non-oracle (band minimum over 121 frequencies) and oracle-assisted (at the
  true critical j omega*) closure distance / |chi|, on the FC18 C2 realization.
Outputs go to results/PCV/<ID>/.
"""

from __future__ import annotations

import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
for p in (EXPERIMENTS, EXPERIMENTS / "connected_cumulants"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
from _cc import (  # noqa: E402
    BAND,
    CORE,
    FC18_SUMMARY,
    HOLDOUT,
    ROOT,
    SCALE,
    Realization,
    band_critical,
    kwargs_of,
)

from _f7_common import Theta, solve_subset  # noqa: E402
from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.physical import physical_matrices  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.cycles.connected import characteristic_values, cumulants  # noqa: E402
from ibr_cycles.diagnosis.baselines import (  # noqa: E402
    generalized_scr,
    interaction_factors,
    nodal_metrics,
    short_circuit_ybus,
)
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_network import (  # noqa: E402
    load_network,
    solve_power_flow,
)

RESULTS = ROOT / "results" / "PCV"
PREREG_COMMIT = "5d0b1986"
FC_RUN = ROOT / "outputs/ias2026/final_math_nonlinear_validation_20260910T231539"
FC03_SUBSETS = FC_RUN / "FC03_governed_replication/FC03_subsets.csv.gz"
FC03_POINTS = FC_RUN / "FC03_governed_replication/FC03_points.csv"
FC10_CENSUS = FC_RUN / "FC10_monotone_and_paths/FC10_census_transverse.csv"
F10_BASELINES = ROOT / "results/F10/F10_baselines.csv"
H4 = CORE
POINTS = {
    "P4": (0.03625, 1.425, 1.5, 1.0),
    "G_S": (0.25, 1.425, 1.5, 1.0),
    "G_S2": (1.0, 1.425, 1.5, 1.0),
    "P_inf": (1.0, 0.5, 1.5, 1.0),
}
N_GRID = 121
SUBSETS4 = [tuple(sorted(s)) for r in range(5) for s in combinations(CORE, r)]


def out_dir(name: str) -> Path:
    path = RESULTS / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, payload) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")


def label(members) -> str:
    return "+".join(map(str, sorted(members))) or "BASE"


def sha256_array(a) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.round(np.asarray(a), 12)).tobytes()
    ).hexdigest()


# ------------------------------------------------------------------ B0 truth --
def transverse_of(case, *, classify: bool = True) -> dict:
    """FC01 direct path on one solved case."""

    a, d, _ = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    tr = transverse_operator(a, r_x, w)
    ev = np.linalg.eigvals(tr.a_perp)
    status = (
        classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
        if classify
        else ("UNSTABLE" if (ev.real > 0).any() else "STABLE")
    )
    f = ev.imag / (2 * np.pi)
    band = ev[(f >= BAND[0]) & (f <= BAND[1])]
    crit = complex(band[np.argmax(band.real)]) if band.size else complex("nan")
    return {
        "status": status,
        "alpha": float(ev.real.max()),
        "rhp": int((ev.real > 0).sum()),
        "crit_re": crit.real,
        "crit_hz": crit.imag / (2 * np.pi),
    }


def truth_task(task):
    """(point_id, theta, members) -> B0 record."""

    pid, theta, members = task
    case = solve_subset(tuple(members), Theta(*theta))
    return {
        "point": pid,
        "subset": label(members),
        "size": len(members),
        **transverse_of(case),
    }


def hypergraph(statuses: dict) -> dict:
    """FC01 h_of on {members_tuple: status}."""

    base = statuses[()]
    if base != "STABLE":
        return {"H": f"BASE_{base}", "kappa": float("nan"), "exact": False}
    unsafe = [s for s, v in statuses.items() if v == "UNSTABLE"]
    unres = [s for s, v in statuses.items() if v == "BOUNDARY_OR_UNRESOLVED"]
    minimal = sorted(
        (s for s in unsafe if not any(set(r) < set(s) for r in unsafe)),
        key=lambda s: (len(s), s),
    )
    return {
        "H": "|".join(label(e) for e in minimal) or "EMPTY",
        "kappa": float(min(len(e) for e in minimal)) if minimal else float("inf"),
        "exact": not unres,
    }


# ------------------------------------------------------------ lower orders --
def mobius_truncation(alpha: dict, order: int) -> dict:
    """Order-r Moebius truncation of a set function alpha (keys: sorted tuples)."""

    keys = list(alpha)
    mu = {}
    for t in keys:
        mu[t] = sum(
            (-1) ** (len(t) - len(r)) * alpha[r] for r in keys if set(r) <= set(t)
        )
    return {
        s: sum(mu[t] for t in keys if set(t) <= set(s) and len(t) <= order)
        for s in keys
    }


def predicted_hypergraph(pred_alpha: dict) -> dict:
    statuses = {s: ("UNSTABLE" if v > 0 else "STABLE") for s, v in pred_alpha.items()}
    return hypergraph(statuses)


def modal_sensitivity(theta, candidates, *, rho: float = 0.01) -> dict:
    """F10 B5 (with the full theta): base least-damped band mode and d lambda / d rho_a."""

    kw = kwargs_of(theta) if theta is not None else {}
    base = solve_case(ReplacementPlan.of({}), **kw)
    values = np.linalg.eigvals(base.system.A)
    f = values.imag / (2 * np.pi)
    band = [
        i
        for i in range(values.size)
        if values[i].imag > 0 and BAND[0] <= f[i] <= BAND[1]
    ]
    i0 = max(band, key=lambda i: values[i].real)
    lam0 = complex(values[i0])
    dlam = {}
    for b in candidates:
        case = solve_case(ReplacementPlan.of({b: rho}), **kw)
        ev = np.linalg.eigvals(case.system.A)
        dlam[b] = complex((ev[np.argmin(np.abs(ev - lam0))] - lam0) / rho)
    return {"lam0": lam0, "dlam": dlam}


def b3_prediction(ms: dict, members) -> tuple[float, bool]:
    lam = ms["lam0"] + sum(ms["dlam"][b] for b in members)
    unstable = lam.real > 0 and BAND[0] <= abs(lam.imag) / (2 * np.pi) <= BAND[1]
    return float(lam.real), bool(unstable)


def participation(theta, candidates) -> dict:
    kw = kwargs_of(theta) if theta is not None else {}
    base = solve_case(ReplacementPlan.of({}), **kw)
    values, vectors = np.linalg.eig(base.system.A)
    f = values.imag / (2 * np.pi)
    band = [
        i
        for i in range(values.size)
        if values[i].imag > 0 and BAND[0] <= f[i] <= BAND[1]
    ]
    i0 = max(band, key=lambda i: values[i].real)
    left = np.linalg.inv(vectors)[i0, :]
    part = np.abs(left * vectors[:, i0])
    part /= part.sum()
    labels = base.system.labels
    return {
        b: float(
            sum(
                part[j]
                for j, n in enumerate(labels)
                if n in (f"delta_sg{b}", f"omega_sg{b}")
            )
        )
        for b in candidates
    }


# ----------------------------------------------------------------- static --
class Static:
    """B1 / B2 / B7 (policy-independent, network only)."""

    def __init__(self):
        self.net = load_network()
        self.flow = solve_power_flow(self.net)
        self.nodal = nodal_metrics(self.net, self.flow)
        payload = json.loads((ROOT / "configs/ias2026/ieee39_network.json").read_text())
        self.p0 = {int(r["bus"]): 100.0 * float(r["p0"]) for r in payload["pv"]}
        self.sn = {
            b: float(self.net.machines[b]["Sn"]) for b in self.net.generator_buses
        }
        self.total_pg = sum(
            (self.flow.injection(b, self.net.ybus) + self.net.loads.get(b, 0j)).real
            * 100.0
            for b in self.net.generator_buses
        )

    def row(self, members) -> dict:
        m = tuple(sorted(members))
        pg = sum(self.p0[b] for b in m)
        out = {
            "pg_mw": pg,
            "sn_mva": sum(self.sn[b] for b in m),
            "penetration": pg / self.total_pg,
            "gscr": float(generalized_scr(self.net, m, self.flow)),
            "min_scr": float(min(self.nodal[b].scr for b in m)),
        }
        if len(m) >= 2:
            f = interaction_factors(self.net, m)
            out["max_miif"] = float(np.max(f[~np.eye(len(m), dtype=bool)]))
            z = np.linalg.inv(short_circuit_ybus(self.net, exclude=m))
            pos = [self.net.position(b) for b in m]
            dist = [
                abs(z[i, i] + z[j, j] - z[i, j] - z[j, i])
                for i, j in combinations(pos, 2)
            ]
            out["compactness_score"] = -float(np.mean(dist))
        else:
            out["max_miif"] = float("nan")
            out["compactness_score"] = float("nan")
        return out


# --------------------------------------------------------------- closure --
def closure_scan(r: Realization, members_all, subsets, crit_hz: dict) -> dict:
    """B8a/B9a (non-oracle, band grid) and B8b/B9b (oracle, at j omega*) per subset."""

    grid = np.linspace(2 * np.pi * BAND[0], 2 * np.pi * BAND[1], N_GRID)
    qs = [r.q_matrix(1j * w)[0] for w in grid]
    out = {}
    for s in subsets:
        if len(s) < 2:
            continue
        idx = [members_all.index(b) for b in s]
        cols = [2 * i + c for i in idx for c in (0, 1)]
        dists = [
            float(np.min(np.abs(np.linalg.eigvals(q[np.ix_(cols, cols)]) + 1.0)))
            for q in qs
        ]
        k = int(np.argmin(dists))
        chi_a = _chi_top(qs[k], cols, len(s))
        rec = {
            "closure_nonoracle": dists[k],
            "closure_freq_hz": grid[k] / (2 * np.pi),
            "abs_chi_nonoracle": chi_a,
        }
        hz = crit_hz.get(label(s))
        if hz is not None and np.isfinite(hz):
            q = r.q_matrix(1j * 2 * np.pi * hz)[0]
            rec["closure_oracle"] = float(
                np.min(np.abs(np.linalg.eigvals(q[np.ix_(cols, cols)]) + 1.0))
            )
            rec["abs_chi_oracle"] = _chi_top(q, cols, len(s))
        out[label(s)] = rec
    return out


def _chi_top(q, cols, n) -> float:
    vals = characteristic_values(q[np.ix_(cols, cols)], [2] * n)
    return float(abs(cumulants(vals, range(n))[frozenset(range(n))]))


__all__ = [
    "BAND",
    "CORE",
    "FC18_SUMMARY",
    "HOLDOUT",
    "ROOT",
    "SCALE",
    "H4",
    "POINTS",
    "SUBSETS4",
    "Realization",
    "Static",
    "band_critical",
    "b3_prediction",
    "closure_scan",
    "hypergraph",
    "kwargs_of",
    "label",
    "mobius_truncation",
    "modal_sensitivity",
    "participation",
    "predicted_hypergraph",
    "transverse_of",
    "truth_task",
]
