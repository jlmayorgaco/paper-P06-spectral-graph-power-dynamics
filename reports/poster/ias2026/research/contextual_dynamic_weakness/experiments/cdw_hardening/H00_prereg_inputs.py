# ruff: noqa: E501
"""H00: generate the frozen hardening design inputs. NO model evaluation happens here.

Outputs (results/hardening/prereg_inputs/):
- hardening_policies.json: HARDENING_H01..H24, deterministic maximin Latin hypercube over the
  validated TX4 policy domain, u in [0,1] (g = u^2), k in [0.5, 2.3], t in [0.5, 3.0],
  h in [0, 2] (the domain of the old holdout H01-H24, CDW00_prereg_inputs.py).
  Maximin: 2000 candidate LHS designs (seed SEED_POLICY + c), keep the one with the largest
  minimum pairwise distance in the unit cube; ties -> lowest candidate index. No model is
  evaluated; points are never selected for expected stability.
- hardening_envelope_draws.json: 10 fresh draws per TX4 envelope (EM-f, EM-u, EC, EMC) with
  the CDW00 generator (identical group bounds) and a NEW seed SEED_ENV (+ envelope index).
- corridor_null_groups.json: size-matched random branch groups, 500 per size and family
  (A = arbitrary, B = connected), SEED_NULL; complete enumeration when the population is
  smaller than 500.
- manifest.json: seeds, bounds and sha256 of every file.
"""

from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
RESEARCH = PROJECT.parent
OUT = PROJECT / "results" / "hardening" / "prereg_inputs"

SEED_POLICY = 20260929
N_CAND = 2000
N_POLICY = 24
SEED_ENV = 20260930
N_ENV = 10
SEED_NULL = 20260931
N_NULL = 500
BOUNDS = {"u": (0.0, 1.0), "k": (0.5, 2.3), "t": (0.5, 3.0), "h": (0.0, 2.0)}
NULL_SIZES = (1, 2, 3, 4, 8, 12)  # every corridor size in results/CDW_E6_corridor_definitions.json


def sc(u, lo, hi):
    return round(float(lo + u * (hi - lo)), 10)


def maximin_design():
    best, best_d, best_c = None, -1.0, -1
    for c in range(N_CAND):
        x = qmc.LatinHypercube(d=4, rng=np.random.default_rng(SEED_POLICY + c)).random(N_POLICY)
        d = float(pdist(x).min())
        if d > best_d + 1e-15:
            best, best_d, best_c = x, d, c
    return best, best_d, best_c


def policies():
    x, dmin, cidx = maximin_design()
    rows = []
    for i, (u, k, t, h) in enumerate(x):
        rows.append({"id": f"HARDENING_H{i + 1:02d}", "u": round(float(u), 10), "g": round(float(u) ** 2, 10),
                     "k": sc(k, *BOUNDS["k"]), "t": sc(t, *BOUNDS["t"]), "h": sc(h, *BOUNDS["h"])})
    return rows, {"min_pairwise_distance_unit_cube": dmin, "selected_candidate": cidx}


def draws():
    import sys

    sys.path.insert(0, str(PROJECT / "experiments" / "cdw"))
    import CDW00_prereg_inputs as P  # the frozen CDW00 generator (bounds, groups, sc)

    mg, cg = list(P.MACHINE_GROUPS), list(P.CONV_GROUPS)
    nu, nc = len(P.SG_BUSES) * len(mg), len(P.CONV_BUSES) * len(cg)

    def unit(r):
        return {str(b): {g: P.sc(r[i * len(mg) + j], *P.MACHINE_GROUPS[g][1:]) for j, g in enumerate(mg)}
                for i, b in enumerate(P.SG_BUSES)}

    def conv(r):
        return {str(b): {g: P.sc(r[i * len(cg) + j], *P.CONV_GROUPS[g][1:]) for j, g in enumerate(cg)}
                for i, b in enumerate(P.CONV_BUSES)}

    out = []
    for i, r in enumerate(P.lhs(len(mg), SEED_ENV + 0, N_ENV)):
        out.append({"envelope": "EM-f", "draw": i, "fleet": {g: P.sc(r[j], *P.MACHINE_GROUPS[g][1:]) for j, g in enumerate(mg)}, "unit": {}, "conv": {}})
    for i, r in enumerate(P.lhs(nu, SEED_ENV + 1, N_ENV)):
        out.append({"envelope": "EM-u", "draw": i, "fleet": {}, "unit": unit(r), "conv": {}})
    for i, r in enumerate(P.lhs(nc, SEED_ENV + 2, N_ENV)):
        out.append({"envelope": "EC", "draw": i, "fleet": {}, "unit": {}, "conv": conv(r)})
    for i, r in enumerate(P.lhs(nu + nc, SEED_ENV + 3, N_ENV)):
        out.append({"envelope": "EMC", "draw": i, "fleet": {}, "unit": unit(r[:nu]), "conv": conv(r[nu:])})
    return out


def branch_graph():
    payload = json.loads((RESEARCH / "configs/ias2026/ieee39_network.json").read_text(encoding="utf-8"))
    ends = [(int(ln["bus1"]), int(ln["bus2"])) for ln in payload["lines"]]
    adj = {e: set() for e in range(len(ends))}
    for a, b in combinations(range(len(ends)), 2):
        if set(ends[a]) & set(ends[b]):
            adj[a].add(b)
            adj[b].add(a)
    return ends, adj


def is_connected_edges(group, adj):
    group = set(group)
    start = next(iter(group))
    seen, stack = {start}, [start]
    while stack:
        e = stack.pop()
        for f in adj[e] & group:
            if f not in seen:
                seen.add(f)
                stack.append(f)
    return seen == group


def null_groups():
    ends, adj = branch_graph()
    n = len(ends)
    rng = np.random.default_rng(SEED_NULL)
    out = {"A": {}, "B": {}}
    for k in NULL_SIZES:
        # family A: arbitrary size-k branch sets
        if k == 1:
            out["A"][k] = [[e] for e in range(n)]
        else:
            seen, lst = set(), []
            while len(lst) < N_NULL:
                g = tuple(sorted(int(x) for x in rng.choice(n, size=k, replace=False)))
                if g not in seen:
                    seen.add(g)
                    lst.append(list(g))
            out["A"][k] = lst
        # family B: connected size-k branch subgraphs (random growth from a uniform edge)
        if k == 1:
            out["B"][k] = [[e] for e in range(n)]
            continue
        if k <= 4:  # exact population of connected k-edge subgraphs; used whole if <= N_NULL
            population = [list(g) for g in combinations(range(n), k) if is_connected_edges(g, adj)]
            if len(population) <= N_NULL:
                out["B"][k] = population
                continue
        seen, lst, attempts = set(), [], 0
        while len(lst) < N_NULL and attempts < 400000:
            attempts += 1
            g = {int(rng.integers(n))}
            while len(g) < k:
                frontier = sorted(set().union(*(adj[e] for e in g)) - g)
                g.add(int(frontier[int(rng.integers(len(frontier)))]))
            t = tuple(sorted(g))
            if t not in seen:
                seen.add(t)
                lst.append(list(t))
        assert all(is_connected_edges(g, adj) for g in lst)
        out["B"][k] = lst
    return {fam: {str(k): v for k, v in d.items()} for fam, d in out.items()}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pol, pinfo = policies()
    files = {"hardening_policies.json": pol, "hardening_envelope_draws.json": draws(), "corridor_null_groups.json": null_groups()}
    man = {"seeds": {"policy": SEED_POLICY, "policy_candidates": N_CAND, "envelope": SEED_ENV, "null": SEED_NULL},
           "policy_bounds": {k: list(v) for k, v in BOUNDS.items()}, "policy_map": "g = u^2", "policy_design": pinfo,
           "n_envelope_draws_per_envelope": N_ENV, "null_sizes": list(NULL_SIZES), "n_null_per_size_family": N_NULL, "sha256": {}}
    for name, payload in files.items():
        p = OUT / name
        p.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        man["sha256"][name] = hashlib.sha256(p.read_bytes()).hexdigest()
    nulls = files["corridor_null_groups.json"]
    man["null_counts"] = {fam: {k: len(v) for k, v in d.items()} for fam, d in nulls.items()}
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1), encoding="utf-8")
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
