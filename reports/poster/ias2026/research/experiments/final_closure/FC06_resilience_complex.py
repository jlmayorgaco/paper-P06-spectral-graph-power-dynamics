"""FC06 (steps 8-9): K_NL(rho), H_NL(rho), kappa_NL(rho), R_k and the topology checks.

Input: FC05_thresholds.csv. A portfolio S tolerates rho iff rho < r_S (every
amplitude up to rho recovers). The complex is built only for
rho < rho_scope(point, family) = min over S of the scope-censored thresholds
(type OUTSIDE or NUMERICAL); beyond that some portfolio's status is not
observable in the model.

Checks, per (point, family, interval between consecutive critical radii):
  1. K_NL(rho) is downward closed (a simplicial complex).
  2. minimal non-faces of K_NL == minimal non-tolerating portfolios == H_NL.
  3. for each H in H_NL with |H| = k >= 2: K_NL[H] = boundary of Delta^{k-1},
     and its reduced homology is the field in degree k-2 and zero elsewhere
     (computed over the rationals).
  4. filtration: rho2 > rho1 implies K(rho2) subset K(rho1); R_{k+1} <= R_k.
  5. rho -> 0+: H_NL(0+) = H_RHP_perp.
The global reduced Betti numbers of K_NL are reported but are NOT claimed to
equal the local ones.
"""

from __future__ import annotations

import json
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from _fc import CORE, FCExperiment, out_dir, write_json

OUT = out_dir("FC06_resilience_complex")
SUBSETS = [tuple(sorted(s)) for k in range(5) for s in combinations(CORE, k)]


def lab(s):
    return "+".join(map(str, s)) or "BASE"


def down_closed(K):
    return all(
        tuple(sorted(r)) in K
        for s in K
        for k in range(len(s))
        for r in combinations(s, k)
    )


def minimal_nonfaces(K):
    if () not in K:  # the base itself does not tolerate rho: the void complex
        return [()]
    non = [s for s in SUBSETS if s not in K]
    return sorted(
        [
            s
            for s in non
            if all(tuple(sorted(r)) in K for r in combinations(s, len(s) - 1))
            and all(
                tuple(sorted(r)) in K for k in range(len(s)) for r in combinations(s, k)
            )
        ],
        key=lambda s: (len(s), s),
    )


def reduced_betti(faces):
    """Reduced Betti numbers over Q of a simplicial complex given by its faces."""

    faces = set(faces)
    if not faces:
        return {}
    by_dim = {}
    for f in faces:
        by_dim.setdefault(len(f) - 1, []).append(f)
    top = max(by_dim)
    index = {
        d: {f: i for i, f in enumerate(sorted(by_dim.get(d, [])))}
        for d in range(-1, top + 1)
    }
    index[-1] = {(): 0}
    ranks = {}
    for d in range(0, top + 1):
        rows, cols = index[d - 1], index[d]
        if not rows or not cols:
            ranks[d] = 0
            continue
        mat = np.zeros((len(rows), len(cols)))
        for f, j in cols.items():
            for i_drop in range(len(f)):
                face = f[:i_drop] + f[i_drop + 1 :]
                mat[rows[face], j] = (-1) ** i_drop
        ranks[d] = int(np.linalg.matrix_rank(mat))
    betti = {}
    for d in range(-1, top + 1):
        n_d = len(index.get(d, {}))
        rank_out = ranks.get(d, 0)  # boundary d -> d-1
        rank_in = ranks.get(d + 1, 0)  # boundary d+1 -> d
        betti[d] = n_d - rank_out - rank_in
    return {d: b for d, b in betti.items() if b}


def analyse(block, h_rhp_perp):
    r = {
        tuple(sorted(int(b) for b in s.split("+"))) if s != "BASE" else (): v
        for s, v in zip(block.subset, block.r, strict=True)
    }
    typ = dict(zip(block.subset, block.type, strict=True))
    # a transversely unstable portfolio has r = 0 by the frozen definition; its TDS
    # label at the smallest amplitude is corroboration only and never censors
    unstable = dict(zip(block.subset, block.alpha_perp > 0, strict=True))
    censored = [
        v
        for s, v in zip(block.subset, block.r, strict=True)
        if typ[s] in ("OUTSIDE_MODEL_SCOPE", "NUMERICAL_FAILURE") and not unstable[s]
    ]
    rho_scope = min(censored) if censored else float("inf")
    crit = sorted({v for v in r.values() if np.isfinite(v) and v < rho_scope} | {0.0})
    edges = crit + [
        rho_scope if np.isfinite(rho_scope) else (crit[-1] * 2 + 1 if crit else 1.0)
    ]
    intervals = []
    prev_k = None
    checks = {
        "down_closed": True,
        "nonfaces_eq_minimal_unsafe": True,
        "local_sphere": True,
        "filtration": True,
    }
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        if hi <= lo:
            continue
        rho = 0.5 * (lo + hi) if np.isfinite(hi) else lo + 1.0
        tol = {s: rho < v for s, v in r.items()}
        K = {
            s
            for s in SUBSETS
            if all(
                tol[tuple(sorted(q))]
                for k in range(len(s) + 1)
                for q in combinations(s, k)
            )
        }
        checks["down_closed"] &= down_closed(K)
        H = minimal_nonfaces(K)
        unsafe = [s for s in SUBSETS if not tol[s]]
        mins = sorted(
            [s for s in unsafe if not any(set(q) < set(s) for q in unsafe)],
            key=lambda s: (len(s), s),
        )
        checks["nonfaces_eq_minimal_unsafe"] &= H == mins
        local = []
        for h in H:
            if len(h) >= 2:
                sub = {s for s in K if set(s) <= set(h)}
                boundary = {
                    tuple(sorted(q)) for k in range(len(h)) for q in combinations(h, k)
                }
                ok = sub == boundary
                betti = reduced_betti(sub)
                sphere = betti == {len(h) - 2: 1}
                checks["local_sphere"] &= ok and sphere
                local.append(
                    {"H": lab(h), "boundary_of_simplex": ok, "reduced_betti": betti}
                )
        if prev_k is not None:
            checks["filtration"] &= K <= prev_k
        prev_k = K
        glob = reduced_betti(K) if K else {}
        kappa = min((len(h) for h in H), default=-1)
        intervals.append(
            {
                "rho_lo": lo,
                "rho_hi": hi,
                "rho_rep": rho,
                "faces": len(K),
                "H_NL": "|".join(lab(h) for h in H) or "EMPTY",
                "kappa_NL": kappa,
                "local": local,
                "global_reduced_betti": glob,
            }
        )
    # Limiter-activation complex (NOT composability): the first event of ANY kind
    # (in-scope failure or activation of an omitted limiter / ride-through band).
    # It says which portfolios drive the plant out of the model's validity first;
    # it says nothing about stability beyond that point.
    any_int = []
    crit_any = sorted({v for v in r.values() if np.isfinite(v)} | {0.0})
    edges_any = crit_any + [crit_any[-1] * 1.5 + 1.0]
    for lo, hi in zip(edges_any[:-1], edges_any[1:], strict=True):
        if hi <= lo:
            continue
        rho = 0.5 * (lo + hi)
        tol = {s: rho < v for s, v in r.items()}
        K = {
            s
            for s in SUBSETS
            if all(
                tol[tuple(sorted(q))]
                for k in range(len(s) + 1)
                for q in combinations(s, k)
            )
        }
        H = minimal_nonfaces(K)
        any_int.append(
            {
                "rho_lo": lo,
                "rho_hi": hi,
                "H_first_event": "|".join(lab(h) for h in H) or "EMPTY",
                "kappa_first_event": min((len(h) for h in H), default=-1),
            }
        )
    finite = {s: v for s, v in r.items() if np.isfinite(v)}
    rk = {
        k: min((v for s, v in r.items() if 1 <= len(s) <= k), default=float("inf"))
        for k in (1, 2, 3, 4)
    }
    checks["R_k_monotone"] = all(rk[k + 1] <= rk[k] for k in (1, 2, 3))
    zero = intervals[0]["H_NL"] if intervals else ""
    return {
        "rho_scope": rho_scope,
        "critical_radii": crit,
        "intervals": intervals,
        "R_k": rk,
        "first_event_intervals": any_int,
        "checks": checks,
        "H_NL_0plus": zero,
        "H_RHP_perp": h_rhp_perp,
        "H_NL_0plus_equals_H_RHP_perp": zero == h_rhp_perp,
        "n_finite_thresholds": len(finite),
    }


def main(argv) -> int:
    exp = FCExperiment(
        name="FC06_resilience_complex",
        question=(
            "What are K_NL, H_NL, kappa_NL, R_k, and do the topology statements hold?"
        ),
    )
    src = out_dir("FC05_nonlinear_thresholds") / "FC05_thresholds.csv"
    frame = pd.read_csv(src)
    import yaml
    from _fc import CONFIGS

    cfg = yaml.safe_load(
        (CONFIGS / "final_nonlinear_composability_v1.yaml").read_text(encoding="utf-8")
    )
    res, rows = {}, []
    for (pt, fam), blk in frame.groupby(["point", "family"]):
        h_rhp = cfg["policy_points"][pt]["H_RHP_perp"]
        a = analyse(blk, h_rhp)
        res[f"{pt}|{fam}"] = a
        for iv in a["intervals"]:
            rows.append(
                {
                    "point": pt,
                    "family": fam,
                    "rho_lo": iv["rho_lo"],
                    "rho_hi": iv["rho_hi"],
                    "H_NL": iv["H_NL"],
                    "kappa_NL": iv["kappa_NL"],
                    "kappa_RHP_perp": cfg["policy_points"][pt]["kappa_RHP_perp"],
                    "rho_scope": a["rho_scope"],
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "FC06_kappa_NL.csv", index=False)
    fe = pd.DataFrame(
        [
            {"point": k.split("|")[0], "family": k.split("|")[1], **iv}
            for k, a in res.items()
            for iv in a["first_event_intervals"]
        ]
    )
    fe.to_csv(OUT / "FC06_first_event_complex.csv", index=False)
    rk = pd.DataFrame(
        [
            {"point": k.split("|")[0], "family": k.split("|")[1], "k": kk, "R_k": v}
            for k, a in res.items()
            for kk, v in a["R_k"].items()
        ]
    )
    rk.to_csv(OUT / "FC06_R_k.csv", index=False)
    write_json(OUT / "FC06_summary.json", res)
    exp.finish(
        "COMPUTED",
        **{
            k: {
                kk: v[kk]
                for kk in ("rho_scope", "checks", "H_NL_0plus_equals_H_RHP_perp", "R_k")
            }
            for k, v in res.items()
        },
    )
    print(table.to_string(index=False))
    print(json.dumps({k: v["checks"] for k, v in res.items()}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
