"""Journal gate 3 - IEEE 68-bus / NETS-NYPS replication, as preregistered.

Protocol: configs/ieee68/G3_preregistration.yaml (commit d9fa097e, frozen before
any 68-bus code existed). Data: configs/ieee68/ieee68_network.json
(experiments/G3_import_ieee68.py). Model: src/ibr_cycles/models/ieee68_devices.py.

Steps
    1. Eligibility: power flow against the report's Table 1; base modes against
       Table 4; base stability; distance to the documented controller limits.
    2. Map G68A (g x k, 31 x 31): H_RHP (primary) and H_IA (secondary) over the
       16 portfolios of the candidates {9, 6, 3, 4}; discrepancy classes as in G1.
    3. RHP boundaries on every pure-policy edge whose H_RHP differs, classified
       as an oscillatory imaginary-axis crossing or an aperiodic crossing through
       the origin.
    4. Resolution check: 60 random off-grid points on the direct path.
    5. Secondary family: the 12 physical plants G1-G12 (4 096 portfolios) at
       k = 1 and g = 0, 0.25, 1.

Exact assembly. The policy coordinate g enters only the converter rows of A_red,
and the regulator gain k enters only the DC4B regulator rows and the ST1A field
equation. Both enter linearly, and the operating point is common to every
subset and every (g, k) (matched dispatch). So

    A(g, k) = A(0, 1) + g [A(1, 1) - A(0, 1)] + (k - 1) [A(0, 2) - A(0, 1)]

holds exactly, from three direct solves per subset. The resolution check tests
this against the direct path at random points.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS, ROOT
from _gates import GateExperiment
from _overnight import choose_workers, pin_blas_threads
from G1_whole_rhp import classify_mode, rhp_modes
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_order,
    incompatibility_hypergraph,
    parse_hypergraph_label,
)
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

OUT = RESULTS / "G3"
NETWORK = ROOT / "configs" / "ieee68" / "ieee68_network.json"
CANDIDATES = (9, 6, 3, 4)
PHYSICAL = tuple(range(1, 13))
SUBSETS = [
    tuple(sorted(c)) for n in range(5) for c in itertools.combinations(CANDIDATES, n)
]
LEAK = 0.05
N = 31
G_GRID = (np.arange(N) / (N - 1)) ** 2
K_GRID = np.linspace(0.5, 2.0, N)
BAND = (0.3, 1.5)
ZERO = 1e-3
SEED = 20260910
SPOTS = 60
ITER = 40

# Table 1 of the report: bus voltage magnitude (pu) and angle (deg)
TABLE1 = {
    1: (1.045, -8.9563),
    2: (0.98, -0.9835),
    3: (0.983, 1.6129),
    4: (0.997, 1.6683),
    5: (1.011, -0.6276),
    6: (1.05, 3.8425),
    7: (1.063, 6.0307),
    8: (1.03, -2.841),
    9: (1.025, 2.6524),
    10: (1.01, -9.6439),
    11: (1.0, -7.2245),
    12: (1.0156, -22.6313),
    13: (1.011, -28.6539),
    14: (1.0, 10.962),
    15: (1.0, 0.0168),
    16: (1.0, 0.0),
    17: (0.9499, -36.0269),
    18: (1.0023, -5.8054),
    19: (0.932, -4.2634),
    20: (0.9806, -5.8744),
    21: (0.9602, -7.0539),
    22: (0.9937, -1.801),
    23: (0.9961, -2.1606),
    24: (0.9587, -9.8767),
    25: (0.9981, -9.9995),
    26: (0.9869, -11.0194),
    27: (0.9679, -12.8571),
    28: (0.9897, -7.4968),
    29: (0.9921, -4.5464),
    30: (0.9762, -19.71),
    31: (0.9838, -17.464),
    32: (0.9699, -15.2375),
    33: (0.9738, -19.758),
    34: (0.98, -26.1159),
    35: (1.043, -27.0886),
    36: (0.9606, -28.8273),
    37: (0.9555, -11.788),
    38: (0.989, -18.7593),
    39: (0.9915, -39.2902),
    40: (1.0442, -13.64),
    41: (0.9996, 9.4272),
    42: (0.999, -0.8435),
    43: (0.9765, -37.9088),
    44: (0.9775, -37.9863),
    45: (1.0471, -29.3626),
    46: (0.9903, -20.5314),
    47: (1.0184, -19.4976),
    48: (1.0337, -18.3761),
    49: (0.9936, -19.8196),
    50: (1.0602, -19.0521),
    51: (1.0634, -27.2882),
    52: (0.9545, -12.8334),
    53: (0.9863, -18.9373),
    54: (0.9857, -11.537),
    55: (0.9571, -13.2179),
    56: (0.9208, -11.9593),
    57: (0.9102, -11.2129),
    58: (0.909, -10.4023),
    59: (0.9037, -13.311),
    60: (0.9062, -14.0368),
    61: (0.9556, -23.2222),
    62: (0.9121, -7.3117),
    63: (0.9096, -8.37),
    64: (0.8367, -8.377),
    65: (0.9128, -8.1847),
    66: (0.9194, -10.1909),
    67: (0.928, -11.4299),
    68: (0.9483, -10.0712),
}
# Table 4 of the report (PSS on G1-G12): damping ratio (%), frequency (Hz)
TABLE4 = [
    (33.537, 0.314),
    (3.621, 0.52),
    (9.625, 0.591),
    (3.381, 0.779),
    (27.136, 0.972),
    (18.566, 1.08),
    (23.607, 0.939),
    (30.111, 1.078),
    (28.306, 1.136),
    (13.426, 1.278),
    (18.83, 1.188),
    (32.062, 1.292),
    (39.542, 1.288),
    (33.119, 1.367),
    (23.87, 5.075),
]


def net():
    return load_network(NETWORK)


def solve(members, g, k):
    return solve_case(
        ReplacementPlan.of({b: 1.0 for b in members}),
        network=net(),
        converter=ConverterParameters(
            voltage_control=True, voltage_gain=float(g), voltage_leak=LEAK
        ),
        machine_scaling={"ka": float(k)},
    )


# ------------------------------------------------------------- eligibility ----


def eligibility() -> dict:
    network = net()
    pf = solve_power_flow(network)
    dv = max(abs(abs(pf.at(b)) - v) for b, (v, _) in TABLE1.items())
    da = max(abs(math.degrees(np.angle(pf.at(b))) - a) for b, (_, a) in TABLE1.items())
    base = solve((), 0.0, 1.0)
    values = np.linalg.eigvals(base.system.A)
    keep = values[(np.abs(values) > ZERO) & (values.imag > 0)]
    rows = []
    for zeta, f in TABLE4:
        target_z, target_f = zeta / 100.0, f
        fz = keep.imag / (2 * math.pi)
        zz = -keep.real / np.abs(keep)
        i = int(np.argmin(np.hypot((fz - target_f) / 0.05, (zz - target_z) / 0.01)))
        rows.append(
            {
                "table_zeta_pct": zeta,
                "table_f_hz": f,
                "ours_zeta_pct": 100 * zz[i],
                "ours_f_hz": fz[i],
                "d_zeta_pct": 100 * zz[i] - zeta,
                "d_f_hz": fz[i] - f,
            }
        )
    modes = pd.DataFrame(rows)
    inter = modes.iloc[:4]
    return {
        "power_flow": {
            "max_dV_pu": dv,
            "max_dangle_deg": da,
            "converged": pf.converged,
            "verdict": "REPRODUCED" if dv < 1e-3 and da < 0.05 else "NOT REPRODUCED",
        },
        "base_modes": modes.to_dict("records"),
        "interarea_verdict": "REPRODUCED"
        if (inter.d_f_hz.abs() <= 0.02).all() and (inter.d_zeta_pct.abs() <= 1.0).all()
        else "NOT REPRODUCED",
        "max_abs_d_f_hz_all15": float(modes.d_f_hz.abs().max()),
        "max_abs_d_zeta_pct_all15": float(modes.d_zeta_pct.abs().max()),
        "base_abscissa": float(values[np.abs(values) > ZERO].real.max()),
        "reference_zeros": int(np.count_nonzero(np.abs(values) <= ZERO)),
        "n_states_base": int(base.n_states),
        "min_limit_margin": float(
            min(
                s.device.limit_margin
                for s in base.dae.slots
                if hasattr(s.device, "limit_margin")
            )
        ),
    }


# ------------------------------------------------------------ exact assembly --


def parts(members):
    a00 = solve(members, 0.0, 1.0).system.A
    a10 = solve(members, 1.0, 1.0).system.A
    a02 = solve(members, 0.0, 2.0).system.A
    return members, a00, a10 - a00, a02 - a00


_PARTS: dict = {}


def matrix(members, g, k):
    a00, dg, dk = _PARTS[members]
    return a00 + g * dg + (k - 1.0) * dk


def counts_at(g, k, subsets):
    rhp, ia, vals = {}, {}, {}
    for m in subsets:
        v = np.linalg.eigvals(matrix(m, g, k))
        vals[m] = v
        r = rhp_modes(v)
        f = np.abs(r.imag) / (2 * math.pi)
        rhp[frozenset(m)] = int(r.size)
        ia[frozenset(m)] = int(np.count_nonzero((f >= BAND[0]) & (f <= BAND[1])))
    return rhp, ia, vals


def node(task):
    i, j = task
    g, k = float(G_GRID[i]), float(K_GRID[j])
    rhp, ia, vals = counts_at(g, k, SUBSETS)
    if rhp[frozenset()] > 0:
        return {
            "i": i,
            "j": j,
            "g": g,
            "k": k,
            "H_RHP": "BASE_UNSTABLE",
            "H_IA": "BASE_UNSTABLE",
        }
    h_rhp = hypergraph_label(incompatibility_hypergraph(rhp))
    h_ia = hypergraph_label(incompatibility_hypergraph(ia))
    classes = set()
    extra = [
        e
        for e in parse_hypergraph_label(h_rhp)
        if e not in set(parse_hypergraph_label(h_ia))
    ]
    for e in extra:
        for v in rhp_modes(vals[tuple(sorted(e))]):
            c = classify_mode(v)
            if c != "IN_BAND":
                classes.add(c)
    return {
        "i": i,
        "j": j,
        "g": g,
        "k": k,
        "H_RHP": h_rhp,
        "H_IA": h_ia,
        "kappa_RHP": hypergraph_order(parse_hypergraph_label(h_rhp)),
        "kappa_IA": hypergraph_order(parse_hypergraph_label(h_ia)),
        "indicator": "".join(str(int(rhp[frozenset(m)] > 0)) for m in SUBSETS),
        "discrepancy": "NO_DISCREPANCY"
        if h_rhp == h_ia
        else ("+".join(sorted(classes)) or "IN_BAND_ONLY"),
    }


def edge(task):
    """Bisect every subset whose RHP indicator flips between two policy nodes."""

    k, g0, g1, ind0, ind1 = task
    out = []
    for pos, m in enumerate(SUBSETS):
        if ind0[pos] == ind1[pos] or not m:
            continue
        lo, hi = g0, g1
        unsafe_lo = ind0[pos] == "1"
        for _ in range(ITER):
            mid = 0.5 * (lo + hi)
            n = rhp_modes(np.linalg.eigvals(matrix(m, mid, k))).size
            if (n > 0) == unsafe_lo:
                lo = mid
            else:
                hi = mid
        v_u = np.linalg.eigvals(matrix(m, lo if unsafe_lo else hi, k))
        crit = v_u[(np.abs(v_u) > ZERO)]
        lam = crit[np.argmax(crit.real)]
        f = abs(lam.imag) / (2 * math.pi)
        out.append(
            {
                "k": k,
                "witness": "+".join(map(str, m)),
                "g_star": 0.5 * (lo + hi),
                "freq_hz": f,
                "lambda_re": float(lam.real),
                "crossing": "APERIODIC_THROUGH_ORIGIN"
                if f < 1e-3
                else "OSCILLATORY_AXIS_CROSSING",
                "direction": "LEAVES_UNSAFE" if unsafe_lo else "ENTERS_UNSAFE",
            }
        )
    return out


def spot(task):
    g, k, nearest = task
    counts = {}
    for m in SUBSETS:
        v = np.linalg.eigvals(solve(m, g, k).system.A)
        counts[frozenset(m)] = int(rhp_modes(v).size)
    direct = (
        hypergraph_label(incompatibility_hypergraph(counts))
        if counts[frozenset()] == 0
        else "BASE_UNSTABLE"
    )
    rhp, _, _ = counts_at(g, k, SUBSETS)
    assembled = (
        hypergraph_label(incompatibility_hypergraph(rhp))
        if rhp[frozenset()] == 0
        else "BASE_UNSTABLE"
    )
    return {
        "g": g,
        "k": k,
        "direct": direct,
        "assembled": assembled,
        "nearest_node": nearest,
    }


# -------------------------------------------------------- secondary family ----


def family_task(members):
    a0 = solve(members, 0.0, 1.0).system.A
    a1 = solve(members, 1.0, 1.0).system.A
    return members, {
        g: int(rhp_modes(np.linalg.eigvals(a0 + g * (a1 - a0))).size)
        for g in (0.0, 0.25, 1.0)
    }


def init():
    pin_blas_threads()


def init_parts(table):
    pin_blas_threads()
    _PARTS.update(table)


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    workers = choose_workers()
    exp = GateExperiment(
        name="G3_ieee68",
        question="Does policy-dependent composability replicate on the 68-bus system?",
        config={
            "prereg": "configs/ieee68/G3_preregistration.yaml",
            "candidates": CANDIDATES,
            "grid": [N, N],
            "g_grid": G_GRID.tolist(),
            "k_grid": K_GRID.tolist(),
        },
        workers=workers,
    )
    started = time.time()
    summary: dict = {"eligibility": eligibility()}
    print(json.dumps(summary["eligibility"], indent=1, default=str)[:2500], flush=True)

    with Pool(workers, initializer=init) as pool:
        table = {m: (a, dg, dk) for m, a, dg, dk in pool.map(parts, SUBSETS)}
    _PARTS.update(table)
    with Pool(workers, initializer=init_parts, initargs=(table,)) as pool:
        nodes = pd.DataFrame(
            pool.map(node, [(i, j) for i in range(N) for j in range(N)], chunksize=8)
        )
        nodes.to_csv(OUT / "G3_nodes.csv", index=False)
        ok = nodes[nodes.H_RHP != "BASE_UNSTABLE"]
        edges_tasks = []
        for j in range(N):
            line = nodes[nodes.j == j].sort_values("i")
            for a, b in zip(
                line.itertuples(), line.iloc[1:].itertuples(), strict=False
            ):
                if (
                    a.H_RHP == "BASE_UNSTABLE"
                    or b.H_RHP == "BASE_UNSTABLE"
                    or a.indicator == b.indicator
                ):
                    continue
                edges_tasks.append(
                    (float(K_GRID[j]), a.g, b.g, a.indicator, b.indicator)
                )
        edges = pd.DataFrame([r for rows in pool.map(edge, edges_tasks) for r in rows])
        edges.to_csv(OUT / "G3_rhp_boundaries.csv", index=False)
        rng = np.random.default_rng(SEED)
        spots = []
        for _ in range(SPOTS):
            g, k = float(rng.uniform(0, 1)), float(rng.uniform(0.5, 2.0))
            i = int(np.argmin(np.abs(G_GRID - g)))
            j = int(np.argmin(np.abs(K_GRID - k)))
            spots.append((g, k, nodes[(nodes.i == i) & (nodes.j == j)].H_RHP.iloc[0]))
        spot_rows = pd.DataFrame(pool.map(spot, spots, chunksize=1))
        spot_rows.to_csv(OUT / "G3_spot_checks.csv", index=False)
    phys = [
        tuple(sorted(c))
        for n in range(len(PHYSICAL) + 1)
        for c in itertools.combinations(PHYSICAL, n)
    ]
    with Pool(workers, initializer=init) as pool:
        fam = dict(pool.map(family_task, phys, chunksize=16))
    family = {}
    for g in (0.0, 0.25, 1.0):
        counts = {frozenset(m): c[g] for m, c in fam.items()}
        if counts[frozenset()] > 0:
            family[str(g)] = {"H_RHP": "BASE_UNSTABLE"}
            continue
        edges12 = incompatibility_hypergraph(counts)
        outside = [e for e in edges12 if not e <= set(CANDIDATES)]
        restricted = hypergraph_label([e for e in edges12 if e <= set(CANDIDATES)])
        family[str(g)] = {
            "n_hyperedges": len(edges12),
            "kappa": hypergraph_order(edges12),
            "label": hypergraph_label(edges12),
            "restricted_to_candidates": restricted,
            "map_label_at_k1": nodes[
                (nodes.j == int(np.argmin(np.abs(K_GRID - 1.0))))
                & (np.isclose(nodes.g, g))
            ].H_RHP.tolist(),
            "hyperedges_outside_candidates": len(outside),
            "smallest_outside": sorted(len(e) for e in outside)[:5],
        }
    pd.DataFrame(
        [
            {"subset": "+".join(map(str, m)), **{f"rhp_g{g}": c[g] for g in c}}
            for m, c in fam.items()
        ]
    ).to_csv(OUT / "G3_family12.csv.gz", index=False)

    lines_change = int(sum(ok[ok.j == j].H_RHP.nunique() > 1 for j in range(N)))
    summary.update(
        {
            "nodes": int(len(nodes)),
            "base_unstable_nodes": int((nodes.H_RHP == "BASE_UNSTABLE").sum()),
            "H_RHP_counts": ok.H_RHP.value_counts().to_dict(),
            "H_IA_counts": ok.H_IA.value_counts().to_dict(),
            "H_equal": float((ok.H_RHP == ok.H_IA).mean()) if len(ok) else None,
            "kappa_RHP_values": sorted(int(x) for x in ok.kappa_RHP.dropna().unique()),
            "H_RHP_EMPTY_nodes": int((ok.H_RHP == "EMPTY").sum()),
            "discrepancy_classes": ok.discrepancy.value_counts().to_dict(),
            "policy_lines_H_RHP_changes": lines_change,
            "policy_lines": N,
            "rhp_boundaries": {
                "located": int(len(edges)),
                "by_type": edges.crossing.value_counts().to_dict()
                if len(edges)
                else {},
                "freq_range_hz": [
                    float(edges.freq_hz.min()),
                    float(edges.freq_hz.max()),
                ]
                if len(edges)
                else None,
                "witnesses": edges.witness.value_counts().to_dict()
                if len(edges)
                else {},
            },
            "spot_check": {
                "n": int(len(spot_rows)),
                "assembled_equals_direct": int(
                    (spot_rows.direct == spot_rows.assembled).sum()
                ),
                "nearest_node_equals_direct": int(
                    (spot_rows.direct == spot_rows.nearest_node).sum()
                ),
            },
            "family12": family,
        }
    )
    oscill = len(edges) and (edges.crossing == "OSCILLATORY_AXIS_CROSSING").any()
    summary["decision"] = (
        "REPRODUCED" if lines_change > 0 and len(edges) else "NOT REPRODUCED"
    )
    summary["decision_detail"] = {
        "lines_with_H_RHP_change": lines_change,
        "boundaries_located": int(len(edges)),
        "any_oscillatory": bool(oscill),
    }
    (OUT / "G3_summary.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8"
    )
    exp.finish("COMPUTED", elapsed_s=time.time() - started)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "eligibility"},
            indent=1,
            default=str,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
