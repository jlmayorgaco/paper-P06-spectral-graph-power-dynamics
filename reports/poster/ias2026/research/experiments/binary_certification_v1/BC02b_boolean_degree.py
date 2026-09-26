"""BC02b: exact Boolean dependence and polynomial degree of the binary family.

On the common realization C2 (filler variant, BC02), for every subset T of the
candidates the Moebius coefficient

    mu_T[F] = sum_{R subset T} (-1)^{|T|-|R|} F(R)

is computed for four objects:

  1. the descriptor Jacobian J(delta)        -> degree 1 expected (affine)
  2. the reduced matrix A(delta)             -> rational in delta; Moebius
                                                spectrum measured by order
  3. the device factor h_S(s) = det(sI - f_x(delta)): a product of per-device
     factors, so log h_S is ADDITIVE (order >= 2 coefficients of log h vanish)
  4. the port ratio r_S(s) = det T_S(s) / det T_0(s) = det(I + M_SS(s)),
     M = blkdiag(dY_i) E^T T_0^-1 E. Its Moebius coefficients are exactly the
     principal-minor sums  mu_T = sum_{U: blocks(U) = T} det M_UU  (checked
     against the vertex values), so its degree is the largest |T| with a
     non-vanishing fully-touching minor sum.

Hence det P_S(s) = h_0 * prod_{i in S} kappa_i(s) * det T_0 * det(I + M_SS): all
Boolean interaction of the characteristic function sits in det(I + M_SS).
Families: IEEE-39 core {30,33,35,37} at P4, Kundur {2,3,4} at (g,k) = (0.08,
1.25), IEEE-68 {3,4,6,9} at g = 0 (the BC02 families; nothing new is tuned).
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
from _bc import WORKERS, BCExperiment, out_dir, write_json
from BC02_binary_representation import RES, _spec, port_t  # noqa: E402

from _overnight import pin_blas_threads  # noqa: E402
from ibr_cycles.certification.binary import (  # noqa: E402
    all_vertices,
    build_common_realization,
    reduced,
    stacked,
)

OUT = out_dir("BC02b")
S_VALUES = (0.5 + 2.0j, 0.1 + 4.0j, 2.0 + 0.5j, 0.05 + 3.46j)


def subsets(items):
    for r in range(len(items) + 1):
        yield from combinations(items, r)


def mobius(values: dict, items) -> dict:
    out = {}
    for t in subsets(items):
        acc = 0
        for r in subsets(t):
            acc = acc + (-1) ** (len(t) - len(r)) * values[r]
        out[t] = acc
    return out


def by_order(coefficients: dict, scale: float) -> dict:
    orders = {}
    for t, c in coefficients.items():
        v = float(np.linalg.norm(c)) if np.ndim(c) else float(abs(c))
        orders[len(t)] = max(orders.get(len(t), 0.0), v / scale)
    return orders


def family(task):
    name, candidates, spec = task
    builder, kwargs = _spec(spec)
    started = time.time()
    cr = build_common_realization(candidates, **kwargs)
    items = tuple(candidates)
    jac = {}
    for delta, s in all_vertices(candidates):
        j, _, _ = cr.jacobian(delta, "C2")
        jac[tuple(sorted(s))] = j
    net = cr.case.dae.network
    channels = {b: [2 * net.position(b), 2 * net.position(b) + 1] for b in items}
    out = {"family": name, "candidates": list(items)}
    j_vals = {t: stacked(jac[t]) for t in jac}
    a_vals = {t: reduced(jac[t]) for t in jac}
    out["J_moebius_rel_by_order"] = by_order(
        mobius(j_vals, items), np.linalg.norm(j_vals[()])
    )
    out["A_moebius_rel_by_order"] = by_order(
        mobius(a_vals, items), np.linalg.norm(a_vals[()])
    )
    first = out["A_moebius_rel_by_order"][1]
    out["A_moebius_by_order_over_order1"] = {
        k: v / first for k, v in out["A_moebius_rel_by_order"].items() if k >= 1
    }
    per_s = []
    for s in S_VALUES:
        n = {t: jac[t].fx.shape[0] for t in jac}
        log_h = {t: np.linalg.slogdet(s * np.eye(n[t]) - jac[t].fx) for t in jac}
        log_h = {t: np.log(sign) + ld for t, (sign, ld) in log_h.items()}
        ports = {t: port_t(jac[t], s) for t in jac}
        t0 = ports[()]
        sign0, ld0 = np.linalg.slogdet(t0)
        ratio = {}
        for t, p in ports.items():
            sg, ld = np.linalg.slogdet(p)
            ratio[t] = (sg / sign0) * np.exp(ld - ld0)
        # port increments and the principal-minor prediction
        idx = [c for b in items for c in channels[b]]
        k_mat = np.linalg.solve(t0, np.eye(t0.shape[0]))[np.ix_(idx, idx)]
        d_blocks = []
        for b in items:
            ch = channels[b]
            d_blocks.append((ports[(b,)] - t0)[np.ix_(ch, ch)])
        dmat = np.zeros((2 * len(items), 2 * len(items)), complex)
        for i, blk in enumerate(d_blocks):
            dmat[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = blk
        m_mat = dmat @ k_mat
        block_of = {2 * i + c: items[i] for i in range(len(items)) for c in (0, 1)}
        minor_sum = {}
        for t in subsets(items):
            acc = 0j if t else 1.0 + 0j
            if t:
                cols = [2 * items.index(b) + c for b in t for c in (0, 1)]
                for r in range(1, len(cols) + 1):
                    for u in combinations(cols, r):
                        if {block_of[x] for x in u} == set(t):
                            acc += np.linalg.det(m_mat[np.ix_(u, u)])
            minor_sum[t] = acc
        mu_ratio = mobius(ratio, items)
        predicted_vertex = {}
        for t in subsets(items):
            cols = [2 * items.index(b) + c for b in t for c in (0, 1)]
            predicted_vertex[t] = (
                np.linalg.det(np.eye(len(cols)) + m_mat[np.ix_(cols, cols)])
                if cols
                else 1.0
            )
        # the complex log is defined modulo 2 pi i: reduce the coefficients there
        mu_logh = {
            t: c.real + 1j * (((c.imag + np.pi) % (2 * np.pi)) - np.pi)
            for t, c in mobius(log_h, items).items()
        }
        per_s.append(
            {
                "s": str(s),
                "ratio_vertex_vs_det(I+M_SS)_max_abs": float(
                    max(abs(ratio[t] - predicted_vertex[t]) for t in ratio)
                ),
                "moebius_vs_principal_minor_sum_max_abs": float(
                    max(abs(mu_ratio[t] - minor_sum[t]) for t in mu_ratio)
                ),
                "ratio_moebius_abs_by_order": by_order(mu_ratio, 1.0),
                "ratio_moebius_order_m": {
                    "+".join(map(str, t)): [
                        float(mu_ratio[t].real),
                        float(mu_ratio[t].imag),
                    ]
                    for t in mu_ratio
                    if len(t) == len(items)
                },
                "log_h_moebius_abs_by_order": by_order(mu_logh, 1.0),
            }
        )
    out["per_s"] = per_s
    out["elapsed_s"] = round(time.time() - started, 1)
    return out


def tasks():
    f8 = json.loads((RES / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    return [
        ("IEEE-39 core at P4", (30, 33, 35, 37), ("39", f8["P4"])),
        ("Kundur K12A g=0.08 k=1.25", (2, 3, 4), ("K", 0.08, 1.25)),
        ("IEEE-68 cand g=0", (3, 4, 6, 9), ("68", 0.0)),
    ]


def main(argv) -> int:
    exp = BCExperiment(
        name="BC02b_boolean_degree",
        question="What is the exact Boolean degree of each object in delta?",
        config={"s_values": [str(s) for s in S_VALUES], "variant": "C2"},
    )
    started = time.time()
    with Pool(min(WORKERS, 3), initializer=pin_blas_threads) as pool:
        results = pool.map(family, tasks(), chunksize=1)
    summary = {r["family"]: r for r in results}
    write_json(OUT / "BC02b_summary.json", summary)
    exp.finish("COMPUTED", elapsed_s=time.time() - started, families=list(summary))
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
