"""F7 shared definitions: the parameter point, its evaluation, and boundary events.

Coordinates of one parameter point ``theta = (g, k, t, h)``:

    g   reactive-policy gain. Every replacing converter runs a LEAKY voltage
        regulator  q_cmd = q_ref + g kp_v e + x_v,  x_v' = g ki_v e - w (x_v - q_ref)
        with the E27/E30 gains kp_v = 2, ki_v = 20 and the frozen leak w = LEAK.
        g = 0 is exact fixed-Q (the matched dispatch) plus one decoupled stable
        pole per converter at -w; g = 1 is the E30 voltage regulator up to the
        leak. The matched operating point is an equilibrium for every g.
    k   excitation-gain scale, multiplicative on every machine's own KA
    t   excitation time-constant scale, multiplicative on every machine's own TE
    h   excitation heterogeneity amplitude. TE_i(h) = G (TE_i / G)^h with G the
        geometric mean of the ten native TE: h = 1 is the native fleet, h = 0 a
        uniform fleet at G, h > 1 an exaggerated spread. Applied BEFORE t.

Two evaluation paths.

``direct``  solve_case at theta: equilibrium, central-difference Jacobians,
            index-1 reduction. The reference.
``fast``    exact parametric assembly. Along every coordinate the equilibrium,
            the state dimension and the algebraic block are fixed (Safeguard A),
            and the only rows of A_red that move are
                converter rows   affine in g
                efd row of machine i   (K_i / T_i) u_i + (1 / T_i) w_i
            so three direct solves per subset determine A_red everywhere. The
            direct path validates it at random points (F7 random spot-check).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from _v2c_common import BAND_HZ, CORE
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_order,
    incompatibility_hypergraph,
    upward_closure_defect,
)
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    _controller_payload,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.port_admittance import (
    build_action_space,
    closure_at,
    refine_closure,
)

LEAK = 0.05
"""Regulator leak, rad/s. Frozen before any F7 map was computed. Pole at -0.05
(0.008 Hz) is far below the band; at 0.57 Hz the leak alters the integral path
by 1.4 percent."""

ZERO = 1e-3
ARTIFICIAL_ZERO = 1e-6
SUBSETS = [
    tuple(sorted(s)) for k in range(len(CORE) + 1) for s in combinations(CORE, k)
]
LABELS = ["+".join(map(str, s)) or "BASE" for s in SUBSETS]
OMEGA_LO = 2.0 * math.pi * BAND_HZ[0]
OMEGA_HI = 2.0 * math.pi * BAND_HZ[1]


def native_te() -> dict[int, float]:
    return {
        int(b): float(r["TE"]) for b, r in _controller_payload()["avr_by_bus"].items()
    }


def native_ka() -> dict[int, float]:
    return {
        int(b): float(r["KA"]) for b, r in _controller_payload()["avr_by_bus"].items()
    }


def te_geometric_mean() -> float:
    values = list(native_te().values())
    return float(math.exp(sum(math.log(v) for v in values) / len(values)))


@dataclass(frozen=True)
class Theta:
    g: float
    k: float = 1.0
    t: float = 1.0
    h: float = 1.0

    def lerp(self, other: Theta, s: float) -> Theta:
        return Theta(
            *(
                a + s * (b - a)
                for a, b in zip(
                    (self.g, self.k, self.t, self.h),
                    (other.g, other.k, other.t, other.h),
                    strict=True,
                )
            )
        )

    def shifted(self, name: str, delta: float) -> Theta:
        values = self.as_dict()
        values[name] += delta
        return Theta(**values)

    def as_dict(self) -> dict[str, float]:
        return {"g": self.g, "k": self.k, "t": self.t, "h": self.h}


def machine_te(theta: Theta) -> dict[int, float]:
    mean = te_geometric_mean()
    return {b: theta.t * mean * (te / mean) ** theta.h for b, te in native_te().items()}


def solve_subset(members, theta: Theta):
    """The direct path."""

    bus_scaling = None
    if theta.h != 1.0:
        mean = te_geometric_mean()
        bus_scaling = {
            bus: {"ta": (mean / te) ** (1.0 - theta.h)}
            for bus, te in native_te().items()
        }
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        converter=ConverterParameters(
            voltage_control=True, voltage_gain=theta.g, voltage_leak=LEAK
        ),
        machine_scaling={"ka": theta.k, "ta": theta.t},
        machine_bus_scaling=bus_scaling,
    )


class FastModel:
    """Exact parametric assembly of A_red(S; theta) from three direct solves."""

    def __init__(self):
        self.ka, self.te = native_ka(), native_te()
        self.parts = {}
        for members in SUBSETS:
            ref = solve_subset(members, Theta(0.0))
            g1 = solve_subset(members, Theta(1.0))
            k2 = solve_subset(members, Theta(0.0, k=2.0))
            a0 = ref.system.A.copy()
            labels = ref.system.labels
            converter_rows = np.array([("_gfl" in n) for n in labels])
            delta_g = np.zeros_like(a0)
            delta_g[converter_rows] = g1.system.A[converter_rows] - a0[converter_rows]
            efd = []
            for i, name in enumerate(labels):
                if name.startswith("efd_sg"):
                    bus = int(name.removeprefix("efd_sg"))
                    r1, r2 = a0[i], k2.system.A[i]
                    u = (r2 - r1) * self.te[bus] / self.ka[bus]
                    w = (2.0 * r1 - r2) * self.te[bus]
                    efd.append((i, bus, u, w))
            self.parts[members] = (a0, delta_g, efd, ref.gz_condition, labels)

    def matrix(self, members, theta: Theta) -> np.ndarray:
        a0, delta_g, efd, _, _ = self.parts[members]
        a = a0 + theta.g * delta_g
        te = machine_te(theta)
        for i, bus, u, w in efd:
            kk, tt = theta.k * self.ka[bus], te[bus]
            a[i] = (kk / tt) * u + (1.0 / tt) * w
        return a


_FAST: FastModel | None = None


def fast_model() -> FastModel:
    global _FAST
    if _FAST is None:
        _FAST = FastModel()
    return _FAST


def eigenvalues(members, theta: Theta, *, fast: bool = True) -> np.ndarray:
    if fast:
        return np.linalg.eigvals(fast_model().matrix(members, theta))
    return np.linalg.eigvals(solve_subset(members, theta).system.A)


def in_gamma(values: np.ndarray) -> np.ndarray:
    f = np.abs(values.imag) / (2.0 * math.pi)
    return (
        (values.real > 0.0)
        & (values.imag >= 0.0)
        & (np.abs(values) > ZERO)
        & (f >= BAND_HZ[0])
        & (f <= BAND_HZ[1])
    )


def summarize(values: np.ndarray) -> dict[str, float]:
    f = np.abs(values.imag) / (2.0 * math.pi)
    upper = (values.imag >= 0.0) & (np.abs(values) > ZERO)
    band = upper & (f >= BAND_HZ[0]) & (f <= BAND_HZ[1])
    unstable = upper & (values.real > 0.0)
    edge_gap = np.minimum(np.abs(f - BAND_HZ[0]), np.abs(f - BAND_HZ[1]))
    return {
        "N": int(np.count_nonzero(in_gamma(values))),
        "band_max": float(values.real[band].max()) if band.any() else float("-inf"),
        "rhp": int(np.count_nonzero(unstable)),
        "abscissa": float(values.real[np.abs(values) > ZERO].max()),
        "artificial_zeros": int(np.count_nonzero(np.abs(values) < ARTIFICIAL_ZERO)),
        "edge_gap_hz": float(edge_gap[unstable].min())
        if unstable.any()
        else float("inf"),
    }


def evaluate(theta: Theta, *, fast: bool = True) -> dict:
    """Every subset of the core at one parameter point."""

    row: dict = {**theta.as_dict(), "path": "fast" if fast else "direct"}
    counts: dict[frozenset, int] = {}
    axis_gap, edge_gap = float("inf"), float("inf")
    for members, label in zip(SUBSETS, LABELS, strict=True):
        try:
            if fast:
                values = eigenvalues(members, theta)
                nx, gz = len(values), fast_model().parts[members][3]
            else:
                case = solve_subset(members, theta)
                values = np.linalg.eigvals(case.system.A)
                nx, gz = case.n_states, case.gz_condition
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            row.update(
                status="INFEASIBLE",
                label="INFEASIBLE",
                reason=f"{label}: {str(error)[:60]}",
            )
            return row
        summary = summarize(values)
        row[f"N_{label}"] = summary["N"]
        row[f"bandmax_{label}"] = summary["band_max"]
        row[f"rhp_{label}"] = summary["rhp"]
        row[f"nx_{label}"] = nx
        row[f"zeros_{label}"] = summary["artificial_zeros"]
        row[f"gzcond_{label}"] = gz
        counts[frozenset(members)] = summary["N"]
        axis_gap = min(axis_gap, abs(summary["band_max"]))
        edge_gap = min(edge_gap, summary["edge_gap_hz"])
        if not members:
            row["base_abscissa"] = summary["abscissa"]
            if summary["abscissa"] >= 0.0:
                row.update(status="BASE_UNSTABLE", label="BASE_UNSTABLE")
                row["signature"] = "BASE_UNSTABLE"
                row["axis_gap"] = abs(summary["abscissa"])
                row["edge_gap_hz"] = float("inf")
                return row
    edges = incompatibility_hypergraph(counts)
    row["label"] = hypergraph_label(edges)
    row["kappa"] = hypergraph_order(edges)
    row["n_edges"] = len(edges)
    row["closure_defect"] = ";".join(
        "+".join(map(str, sorted(s))) for s in upward_closure_defect(counts)
    )
    row["signature"] = (
        row["label"] + "#" + "".join(str(counts[frozenset(s)]) for s in SUBSETS)
    )
    row["axis_gap"] = axis_gap
    row["edge_gap_hz"] = edge_gap
    row["status"] = "OK"
    return row


# ------------------------------------------------------------ boundary events --


def _gamma_violation(value: complex) -> list[str]:
    """Which constraints of Gamma a (upper half plane) eigenvalue violates."""

    f = abs(value.imag) / (2.0 * math.pi)
    out = []
    if value.real <= 0.0:
        out.append("IMAGINARY_AXIS")
    if f < BAND_HZ[0]:
        out.append("LOWER_BAND_EDGE")
    if f > BAND_HZ[1]:
        out.append("UPPER_BAND_EDGE")
    return out


def _nearest(values: np.ndarray, target: complex) -> complex:
    upper = values[values.imag >= 0.0]
    return complex(upper[int(np.argmin(np.abs(upper - target)))])


def locate_event(
    members,
    a: Theta,
    b: Theta,
    *,
    scales: dict[str, float],
    axes: tuple[str, str],
    iterations: int = 40,
    tangency_edges: float = 1.0,
) -> dict:
    """Bisect the indicator ``N(S) > 0`` along a -> b on the fast path and classify.

    Safeguard B: the crossing eigenvalue is matched across the final bracket; the
    constraint of Gamma its safe-side partner violates names the segment.
    Imaginary-axis crossings are split into TRANSVERSAL and TANGENCY by whether
    Re(lambda) along the path has a stationary point within ``tangency_edges``
    edge lengths of the crossing. The in-plane gradient is reported beside it.
    """

    def count(theta):
        values = eigenvalues(members, theta)
        return int(np.count_nonzero(in_gamma(values))), values

    n_a, values_a = count(a)
    n_b, values_b = count(b)
    unsafe_a = n_a > 0
    if unsafe_a == (n_b > 0):
        return {"status": "NO_INDICATOR_CHANGE"}
    lo, hi = 0.0, 1.0
    v_lo, v_hi, n_lo, n_hi = values_a, values_b, n_a, n_b
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        n_mid, v_mid = count(a.lerp(b, mid))
        if (n_mid > 0) == unsafe_a:
            lo, v_lo, n_lo = mid, v_mid, n_mid
        else:
            hi, v_hi, n_hi = mid, v_mid, n_mid
    unsafe_values, safe_values = (v_lo, v_hi) if unsafe_a else (v_hi, v_lo)
    inside = unsafe_values[in_gamma(unsafe_values)]
    matches = [(lam, _nearest(safe_values, lam)) for lam in inside]
    crossing = [
        (lam, p, _gamma_violation(p)) for lam, p in matches if _gamma_violation(p)
    ]
    s_star = 0.5 * (lo + hi)
    theta_star = a.lerp(b, s_star)
    event = {
        "status": "LOCATED",
        "witness": "+".join(map(str, members)),
        "s_star": s_star,
        **{f"{k}_star": v for k, v in theta_star.as_dict().items()},
        "bracket": hi - lo,
        "n_crossing_modes": len(crossing),
        "direction": "ENTERS_GAMMA" if not unsafe_a else "LEAVES_GAMMA",
        "N_from": n_lo,
        "N_to": n_hi,
        "rhp_from": int(
            np.count_nonzero((v_lo.real > 0) & (np.abs(v_lo) > ZERO) & (v_lo.imag >= 0))
        ),
        "rhp_to": int(
            np.count_nonzero((v_hi.real > 0) & (np.abs(v_hi) > ZERO) & (v_hi.imag >= 0))
        ),
    }
    if not crossing:
        event["boundary_type"] = "UNRESOLVED"
        return event
    lam, partner, violated = crossing[0]
    star = 0.5 * (lam + partner)
    event.update(
        lambda_re=float(star.real),
        lambda_im=float(star.imag),
        freq_star_hz=float(abs(star.imag) / (2 * math.pi)),
        matching_gap=float(abs(lam - partner)),
    )
    if len(crossing) > 1 or len(violated) > 1:
        event["boundary_type"] = "CORNER_OR_MULTIPLE"
        return event
    segment = violated[0]
    if segment != "IMAGINARY_AXIS":
        event["boundary_type"] = segment
        return event
    # In-plane gradient of Re(lambda), physical units, on the smooth fast model.
    grad = {}
    for name in axes:
        step = 1e-6 * scales[name]
        plus = _nearest(eigenvalues(members, theta_star.shifted(name, step)), star)
        minus = _nearest(eigenvalues(members, theta_star.shifted(name, -step)), star)
        grad[name] = (plus.real - minus.real) / (2 * step)
    # Transversality along the PATH, scale-free: fit Re(lambda)(s) = r0 + a u + c u^2
    # around the crossing (u in edge lengths). The crossing is a tangency when the
    # stationary point -a / 2c of that parabola lies within one edge length, i.e.
    # the path meets a fold of the boundary at the lattice resolution.
    du = 1e-3
    values = [
        _nearest(eigenvalues(members, a.lerp(b, s_star + d)), star).real
        for d in (-du, 0.0, du)
    ]
    slope = (values[2] - values[0]) / (2 * du)
    curvature = (values[2] - 2 * values[1] + values[0]) / (2 * du * du)
    stationary = abs(slope / (2 * curvature)) if curvature != 0 else float("inf")
    event.update(
        {f"dre_d{n}": grad[n] for n in axes},
        dre_ds_per_edge=slope,
        d2re_ds2_per_edge=curvature,
        fold_distance_edges=stationary,
        boundary_type="IMAGINARY_AXIS_TANGENCY"
        if stationary <= tangency_edges
        else "IMAGINARY_AXIS_TRANSVERSAL",
    )
    return event


def port_diagnostics(members, theta: Theta, omega: float) -> dict:
    """Closure quantities of subset S at s = i omega against the base, direct path."""

    if not members:
        return {}
    base = solve_subset((), theta)
    case = solve_subset(members, theta)
    direct = np.linalg.eigvals(case.system.A)
    nearest = direct[np.argmin(np.abs(direct - complex(0.0, omega)))]
    space = build_action_space(base, case, members)
    s = complex(0.0, omega)
    split = space.split(s)
    total = complex(split["full"])
    product = complex(split["individual"]) * complex(split["collective"])
    scale = max(abs(complex(split["individual"])), 1e-300) * max(
        1.0, abs(complex(split["collective"]))
    )
    sample = closure_at(space, omega)
    reference = [
        abs(complex(space.split(complex(0.0, w))["full"]))
        for w in np.linspace(OMEGA_LO, OMEGA_HI, 25)
    ]
    out = {
        "direct_re_at_star": float(nearest.real),
        "direct_freq_at_star_hz": float(abs(nearest.imag) / (2 * math.pi)),
        "port_det_full": abs(total),
        "port_det_ratio_to_band_median": abs(total) / float(np.median(reference)),
        "port_min_individual": float(min(abs(d) for d in split["diagonals"])),
        "closure_distance": float(
            np.min(np.abs(np.asarray(split["q_eigenvalues"]) + 1.0))
        )
        if len(members) > 1
        else float("nan"),
        "factorization_residual": abs(total - product) / scale,
        "closure_solve_residual": sample.solve_residual,
        "t0_condition": sample.condition,
    }
    if len(members) > 1:
        refined = refine_closure(space, omega, half_width=0.05 * omega, tolerance=1e-9)
        out["closure_min_refined"] = refined.margin
        out["freq_port_hz"] = refined.frequency_hz
    out["port_visible"] = bool(out["port_det_ratio_to_band_median"] < 1e-3)
    return out
