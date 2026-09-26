# ruff: noqa: E501  -- formulas in docstrings and test labels kept on one line
"""CC02 - Phase 2: toy falsification tests for the finite-amplitude connected cumulants.

Synthetic only (no power-system model). Residuals are reported as maxima over fixed-
seed random draws; nothing is tuned except the one declared C4 (d) example.

T1  two independent pairs: full Boolean Moebius term != 0, chi_1234 = 0
T2  directed 4-ring: chi_1234 = spanning-cycle expression
T3  same kappa = 4: composite (symmetric path, chi_V = 0) vs connected (ring)
T4  random scalar Q, n <= 6: partition formula vs explicit cycle enumeration
T5  random 2x2-block Q: inversion, factorization annihilation, block expansion,
    two-block closed form, naive holonomy formula refuted
C4  (kappa, d_alg, d_conn) of five constructions
Thm connectivity of minimal coalitions on random toys
"""

from __future__ import annotations

import time
from itertools import combinations

import numpy as np
import pandas as pd
from _cc import out_dir, write_json

from ibr_cycles.cycles.connected import (
    algebraic_order,
    block_connected_expansion,
    block_cycle_trace_sum,
    boolean_mobius,
    characteristic_values,
    composite_part,
    connected_order,
    connected_share,
    cumulant_support_connected,
    cumulants,
    cumulants_direct,
    moments,
    spanning_cycle_sum,
    subsets,
)
from ibr_cycles.models.lag_network_toy import (
    LagNetworkToy,
    path,
    ring,
    ring_plus_pair,
    scaled_to_boundary,
    star,
    tuned_mobius_zero,
    two_pairs,
)

OUT = out_dir("CC02")
SEED = 20260912
PROBES = [0.1 + 0.4j, 0.0 + 1.3j, 0.5 - 0.2j, 2.0 + 0.0j]


def rand_c(rng, *shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def zdiag(q, p):
    q = q.copy()
    for b in range(q.shape[0] // p):
        q[p * b : p * b + p, p * b : p * b + p] = 0.0
    return q


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


def main() -> int:
    started = time.time()
    rng = np.random.default_rng(SEED)
    rows = []

    def record(test, quantity, value, gate, n=None):
        rows.append(
            {
                "test": test,
                "quantity": quantity,
                "value": float(value),
                "gate": gate,
                "pass": bool(value <= gate) if isinstance(gate, float) else None,
                "n": n,
            }
        )

    # ---- C1 / recursion on random set functions -------------------------------
    inv, rec = 0.0, 0.0
    for n in range(1, 7):
        for _ in range(20):
            f = {
                frozenset(s): complex(rand_c(rng))
                for s in subsets(range(n), nonempty=True)
            }
            f[frozenset()] = 1 + 0j
            chi = cumulants(f, range(n))
            chi_d = cumulants_direct(f, range(n))
            back = moments(chi, range(n))
            scale = max(abs(v) for v in f.values())
            inv = max(inv, max(abs(back[k] - f[k]) for k in f) / scale)
            rec = max(
                rec,
                max(abs(chi[k] - chi_d[k]) for k in chi)
                / max(1, max(abs(v) for v in chi_d.values())),
            )
    record(
        "C1", "max rel |F - moments(cumulants(F))|, n<=6, 120 draws", inv, 1e-10, 120
    )
    record("C1", "max rel |recursion - partition definition|", rec, 1e-10, 120)

    # ---- C2 factorization annihilation -------------------------------------------
    worst = 0.0
    for n1, n2 in [(1, 1), (1, 3), (2, 2), (2, 3), (3, 3), (2, 4)]:
        for _ in range(10):
            g1 = {
                frozenset(s): complex(rand_c(rng))
                for s in subsets(range(n1), nonempty=True)
            }
            g2 = {
                frozenset(s): complex(rand_c(rng))
                for s in subsets(range(n2), nonempty=True)
            }
            g1[frozenset()] = g2[frozenset()] = 1 + 0j
            f = {}
            for a in subsets(range(n1 + n2)):
                a1 = frozenset(x for x in a if x < n1)
                a2 = frozenset(x - n1 for x in a if x >= n1)
                f[frozenset(a)] = g1[a1] * g2[a2]
            chi = cumulants(f, range(n1 + n2))
            scale = max(abs(v) for v in chi.values())
            worst = max(
                worst,
                max(abs(v) for k, v in chi.items() if min(k) < n1 <= max(k)) / scale,
            )
    record(
        "C2", "max |mixed chi| / max |chi| under exact factorization", worst, 1e-12, 60
    )

    # ---- C3 scalar cycle formula (T4) ---------------------------------------------
    worst_zero, worst_diag = 0.0, 0.0
    for n in range(2, 7):
        for _ in range(10):
            for zero in (True, False):
                q = rand_c(rng, n, n) * 0.7
                if zero:
                    np.fill_diagonal(q, 0.0)
                chi = cumulants_direct(characteristic_values(q, [1] * n), range(n))
                for s in subsets(range(n)):
                    if len(s) >= 2:
                        r = rel(
                            chi[frozenset(s)],
                            (-1) ** (len(s) - 1) * spanning_cycle_sum(q, s),
                        )
                        if zero:
                            worst_zero = max(worst_zero, r)
                        else:
                            worst_diag = max(worst_diag, r)
    record(
        "T4/C3",
        "scalar Q_ii = 0: max rel |chi_S - (-1)^(|S|-1) sum spanning cycles|, n<=6",
        worst_zero,
        1e-10,
        50,
    )
    record("T4/C3", "scalar Q_ii != 0 (|S|>=2): same", worst_diag, 1e-10, 50)

    # ---- T5 blocks --------------------------------------------------------------
    inv_b, fac_b, exp_b, two_b, naive_gap_min, gauge = 0.0, 0.0, 0.0, 0.0, np.inf, 0.0
    for m in (2, 3, 4):
        for trial in range(10):
            q = zdiag(rand_c(rng, 2 * m, 2 * m) * 0.6, 2)
            f = characteristic_values(q, [2] * m)
            chi = cumulants(f, range(m))
            back = moments(chi, range(m))
            inv_b = max(inv_b, max(abs(back[k] - f[k]) for k in f))
            qd = q.copy()
            qd[0:2, 2:] = 0.0
            qd[2:, 0:2] = 0.0
            chi_d = cumulants(characteristic_values(qd, [2] * m), range(m))
            fac_b = max(
                fac_b, max(abs(v) for k, v in chi_d.items() if 0 in k and len(k) > 1)
            )
            if m <= 3 or trial < 2:
                for s in subsets(range(m)):
                    if len(s) >= 2:
                        exp_b = max(
                            exp_b,
                            rel(
                                chi[frozenset(s)],
                                block_connected_expansion(q, [2] * m, s),
                            ),
                        )
            for s in subsets(range(m)):
                if len(s) >= 2:
                    naive = (-1) ** (len(s) - 1) * block_cycle_trace_sum(q, [2] * m, s)
                    naive_gap_min = min(
                        naive_gap_min,
                        abs(chi[frozenset(s)] - naive)
                        / max(abs(chi[frozenset(s)]), 1e-300),
                    )
            x, y = q[0:2, 2:4], q[2:4, 0:2]
            two_b = max(
                two_b,
                rel(
                    chi[frozenset({0, 1})],
                    -np.trace(x @ y) + np.linalg.det(x) * np.linalg.det(y),
                ),
            )
            sb = np.zeros((2 * m, 2 * m), complex)
            for b in range(m):
                sb[2 * b : 2 * b + 2, 2 * b : 2 * b + 2] = rand_c(rng, 2, 2)
            chi_g = cumulants(
                characteristic_values(np.linalg.solve(sb, q @ sb), [2] * m), range(m)
            )
            gauge = max(gauge, max(rel(chi_g[k], chi[k]) for k in chi))
    record("T5", "2x2 blocks: max |F - moments(cumulants(F))|", inv_b, 1e-10, 30)
    record(
        "T5", "2x2 blocks: max |mixed chi| after decoupling block 0", fac_b, 1e-12, 30
    )
    record(
        "T5/C3",
        "2x2 blocks: max rel |chi_S - block-connected permutation sum|",
        exp_b,
        1e-10,
        30,
    )
    record(
        "T5/C3",
        "2x2 blocks: max rel |chi_ab - (-tr(Q_ab Q_ba) + det Q_ab det Q_ba)|",
        two_b,
        1e-12,
        30,
    )
    record(
        "T5/C3",
        "2x2 blocks: MIN relative gap |chi_S - naive holonomy formula| / |chi_S| (refutation if >> 0)",
        naive_gap_min,
        "report",
        30,
    )
    record(
        "T5",
        "2x2 blocks: max rel gauge drift of chi under block similarity",
        gauge,
        1e-9,
        30,
    )

    # ---- T1 / T2 --------------------------------------------------------------------
    t1_mu, t1_chi, t2 = np.inf, 0.0, 0.0
    for s in PROBES:
        f = characteristic_values(two_pairs().q(s), [1] * 4)
        mu = boolean_mobius(f, range(4))[frozenset(range(4))]
        chi = cumulants(f, range(4))[frozenset(range(4))]
        t1_mu, t1_chi = min(t1_mu, abs(mu)), max(t1_chi, abs(chi))
        q = ring().q(s)
        chi_r = cumulants(characteristic_values(q, [1] * 4), range(4))[
            frozenset(range(4))
        ]
        t2 = max(t2, rel(chi_r, -q[0, 1] * q[1, 2] * q[2, 3] * q[3, 0]))
    record(
        "T1",
        "two independent pairs: MIN |mu_1234| over probes (must be > 0)",
        t1_mu,
        "report",
        4,
    )
    record("T1", "two independent pairs: MAX |chi_1234| over probes", t1_chi, 1e-12, 4)
    record("T2", "4-ring: max rel |chi_1234 - (-q12 q23 q34 q41)|", t2, 1e-12, 4)

    # ---- T3: same kappa = 4 --------------------------------------------------------
    t3 = []
    for name, maker in (
        ("path (composite)", path),
        ("ring (connected)", ring),
        ("star (composite)", star),
    ):
        at = scaled_to_boundary(maker(), (0, 1, 2, 3))
        lam = at.critical((0, 1, 2, 3))
        f = characteristic_values(at.q(lam), [1] * 4)
        chi = cumulants(f, range(4))
        pushed = LagNetworkToy(1.02 * at.g, at.a)
        proper = max(at.alpha(s) for r in (1, 2, 3) for s in combinations(range(4), r))
        t3.append(
            {
                "system": name,
                "kappa_pushed_1.02": pushed.kappa(),
                "H_pushed": str(pushed.hypergraph()),
                "boundary_eigenvalue": str(lam),
                "max_alpha_proper_at_boundary": proper,
                "abs_F_V_at_zero": abs(f[frozenset(range(4))]),
                "abs_chi_V": abs(chi[frozenset(range(4))]),
                "abs_composite": abs(composite_part(chi, range(4))),
                "nu_V": connected_share(chi, range(4)),
                "support_connected": cumulant_support_connected(chi, range(4), 1e-10),
                "pair_cumulants": {
                    "".join(map(str, k)): round(abs(v), 6)
                    for k, v in chi.items()
                    if len(k) == 2 and abs(v) > 1e-12
                },
            }
        )
    pd.DataFrame(t3).to_csv(OUT / "CC02_T3_same_kappa.csv", index=False)

    # ---- C4 ---------------------------------------------------------------------
    def orders(toy):
        da = dc = 0
        for s in PROBES:
            f = characteristic_values(toy.q(s), [1] * toy.m)
            da = max(da, algebraic_order(boolean_mobius(f, range(toy.m)), 1e-9))
            dc = max(dc, connected_order(cumulants(f, range(toy.m)), 1e-9))
        return da, dc

    c4 = []
    cases = [
        (
            "(a) directed ring",
            LagNetworkToy(1.02 * scaled_to_boundary(ring(), (0, 1, 2, 3)).g, 1.0),
        ),
        (
            "(b) symmetric path",
            LagNetworkToy(1.02 * scaled_to_boundary(path(), (0, 1, 2, 3)).g, 1.0),
        ),
        (
            "(c) symmetric star",
            LagNetworkToy(1.02 * scaled_to_boundary(star(), (0, 1, 2, 3)).g, 1.0),
        ),
        ("(d) tuned det G = 0", tuned_mobius_zero()),
        ("(e) ring + unstable pair", ring_plus_pair()),
    ]
    for name, toy in cases:
        da, dc = orders(toy)
        c4.append(
            {
                "example": name,
                "kappa": toy.kappa(),
                "d_alg": da,
                "d_conn": dc,
                "H": str(toy.hypergraph()),
            }
        )
    pd.DataFrame(c4).to_csv(OUT / "CC02_C4_orders.csv", index=False)

    # ---- connectivity theorem on random toys ----------------------------------------
    checked, violations, top_zero = 0, 0, 0
    trng = np.random.default_rng(SEED + 1)
    for _ in range(500):
        g = trng.normal(size=(5, 5)) * (trng.random((5, 5)) < 0.45)
        np.fill_diagonal(g, 0.0)
        toy = LagNetworkToy(g.astype(complex), 1.0)
        for h in toy.hypergraph():
            if len(h) < 2:
                continue
            sub = LagNetworkToy(g[np.ix_(h, h)].astype(complex), 1.0)
            lam = toy.critical(h)
            f = characteristic_values(sub.q(lam), [1] * len(h))
            chi = cumulants(f, range(len(h)))
            checked += 1
            violations += not cumulant_support_connected(chi, range(len(h)), 1e-9)
            top_zero += abs(chi[frozenset(range(len(h)))]) < 1e-9
    record(
        "Thm 2.3",
        "minimal coalitions whose cumulant support is NOT connected (random toys)",
        violations,
        0.0,
        checked,
    )
    rows.append(
        {
            "test": "Thm 2.3",
            "quantity": "minimal coalitions (|H|>=2) with chi_H = 0 at the zero (composite)",
            "value": float(top_zero),
            "gate": "report",
            "pass": None,
            "n": checked,
        }
    )

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "CC02_residuals.csv", index=False)
    summary = {
        "seed": SEED,
        "all_gated_pass": bool(table[table["pass"].notna()]["pass"].all()),
        "T3": t3,
        "C4": c4,
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "CC02_summary.json", summary)
    print(table.to_string(index=False))
    print(pd.DataFrame(t3).drop(columns=["pair_cumulants"]).to_string(index=False))
    print(pd.DataFrame(c4).to_string(index=False))
    print("all gated pass:", summary["all_gated_pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
