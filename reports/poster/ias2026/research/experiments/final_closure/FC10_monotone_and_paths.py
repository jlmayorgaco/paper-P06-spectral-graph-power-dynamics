"""FC10 (steps 13-14): monotone-class diagnostic and safe-path lattice on the census.

1. The frozen E12 census (512 portfolios of the 9 candidates, the E12 plan: rho = 1,
   matched reactive policy, default converter) is re-labelled with the TRANSVERSE
   criterion on the direct path (quotient by C = span{R_x, w}; near-axis cases
   through the four-state classifier).
2. Monotone-class test (theory: FINAL_TRANSACTION_THEORY_AND_EVIDENCE A.6). Robust
   counterexamples S subset T = S + {i} with alpha(S) > +0.01 and
   alpha(T) < -0.01 at the same policy rule out EVERY common order-preserving
   Metzler realization. The same test is run on the 16-subset lattices of the six
   F8 points (FC01 structure table).
3. Safe paths. Lattice graph: a node is a portfolio, an edge is S -> S + {i};
   both endpoints must be transversely stable. For each target T:
       A  final-state stable
       B  a one-at-a-time stable monotone path from the empty set exists
       C  every subset stable (any-order safe; no hyperedge contained)
   Counts and examples of A and not C, and of A and not B.
4. Spectral planning on the census (P1): maximize replaced_pg_mw subject to A,
   B or C; replaced_sn_mva is reported on its own axis.
"""

from __future__ import annotations

import sys
import time
from collections import deque
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _fc import RESULTS, WORKERS, FCExperiment, out_dir, write_json

from _overnight import pin_blas_threads  # noqa: E402
from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.physical import physical_matrices  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.units import measure  # noqa: E402

OUT = out_dir("FC10_monotone_and_paths")
ROBUST = 0.01


def lab(s):
    return "+".join(map(str, s)) or "BASE"


def task(members):
    case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}, q_policy="matched"))
    a, d, _ = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    tr = transverse_operator(a, r_x, frequency_partner(case.dae).w)
    ev = np.linalg.eigvals(tr.a_perp)
    if np.abs(ev.real).min() < 2e-3:
        st = classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
    else:
        st = "UNSTABLE" if (ev.real > 0).any() else "STABLE"
    q = measure(case)
    return {
        "members": lab(members),
        "size": len(members),
        "alpha_perp": float(ev.real.max()),
        "status": st,
        "replaced_pg_mw": q.replaced_pg_mw,
        "replaced_sn_mva": q.replaced_sn_mva,
    }


def lattice(stat, candidates):
    nodes = [
        tuple(s)
        for k in range(len(candidates) + 1)
        for s in combinations(candidates, k)
    ]
    stable = {s: stat[s] == "STABLE" for s in nodes}
    reach = {(): stable[()]}
    queue = deque([()] if stable[()] else [])
    while queue:
        s = queue.popleft()
        for c in candidates:
            if c in s:
                continue
            t = tuple(sorted(s + (c,)))
            if stable[t] and t not in reach:
                reach[t] = True
                queue.append(t)
    rows = []
    for t in nodes:
        a_ok = stable[t]
        b_ok = bool(reach.get(t, False))
        c_ok = all(
            stable[tuple(sorted(q))]
            for k in range(len(t) + 1)
            for q in combinations(t, k)
        )
        rows.append(
            {
                "members": lab(t),
                "size": len(t),
                "A_final_stable": a_ok,
                "B_safe_path": b_ok,
                "C_any_order_safe": c_ok,
            }
        )
    return pd.DataFrame(rows)


def monotone_examples(frame, candidates):
    alpha = {
        tuple(int(b) for b in m.split("+")) if m != "BASE" else (): a
        for m, a in zip(frame.members, frame.alpha_perp, strict=True)
    }
    ex, weak = [], 0
    for s, a_s in alpha.items():
        for c in candidates:
            if c in s:
                continue
            t = tuple(sorted(s + (c,)))
            if t not in alpha:
                continue
            if a_s > 0 and alpha[t] < 0:
                weak += 1
                if a_s > ROBUST and alpha[t] < -ROBUST:
                    ex.append(
                        {
                            "S": lab(s),
                            "alpha_S": a_s,
                            "added": c,
                            "T": lab(t),
                            "alpha_T": alpha[t],
                        }
                    )
    return ex, weak


def main(argv) -> int:
    pin_blas_threads()
    exp = FCExperiment(
        name="FC10_monotone_and_paths",
        question=(
            "Is IEEE-39 in an order-preserving class; which targets are A / B / C?"
        ),
        config={"robust": ROBUST},
        workers=WORKERS,
    )
    started = time.time()
    census = pd.read_csv(RESULTS / "tables" / "E12_compatibility_census_portfolios.csv")
    cands = sorted(
        {int(b) for m in census.members if m != "BASE" for b in m.split("+")}
    )
    members = [
        tuple(sorted(int(b) for b in m.split("+"))) if m != "BASE" else ()
        for m in census.members
    ]
    with Pool(WORKERS, initializer=pin_blas_threads) as pool:
        rows = pool.map(task, members, chunksize=4)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "FC10_census_transverse.csv", index=False)
    stat = {
        tuple(int(b) for b in m.split("+")) if m != "BASE" else (): s
        for m, s in zip(frame.members, frame.status, strict=True)
    }
    lat = lattice(stat, cands).merge(frame, on=["members", "size"])
    lat.to_csv(OUT / "FC10_census_lattice.csv", index=False)
    ex, weak = monotone_examples(frame, cands)
    # F8 16-subset lattices (FC01 structure table)
    s1 = pd.read_csv(out_dir("FC01_transverse_quotient") / "FC01_structure.csv")
    f8_ex = []
    for pt, blk in s1[s1.benchmark == "IEEE-39"].groupby("point"):
        e, _ = monotone_examples(
            blk.rename(columns={"subset": "members"}), [30, 33, 35, 37]
        )
        f8_ex += [{**x, "point": pt} for x in e]

    def best(mask):
        sub = lat[mask]
        if sub.empty:
            return None
        r = sub.loc[sub.replaced_pg_mw.idxmax()]
        return {
            "members": r.members,
            "replaced_pg_mw": float(r.replaced_pg_mw),
            "replaced_sn_mva": float(r.replaced_sn_mva),
        }

    summary = {
        "census_status": frame.status.value_counts().to_dict(),
        "unresolved": int((~frame.status.isin(["STABLE", "UNSTABLE"])).sum()),
        "monotone_counterexamples_robust": len(ex),
        "monotone_counterexamples_any": weak,
        "examples": ex[:15],
        "f8_lattice_counterexamples_robust": f8_ex,
        "A": int(lat.A_final_stable.sum()),
        "B": int(lat.B_safe_path.sum()),
        "C": int(lat.C_any_order_safe.sum()),
        "A_not_C": int((lat.A_final_stable & ~lat.C_any_order_safe).sum()),
        "A_not_B": int((lat.A_final_stable & ~lat.B_safe_path).sum()),
        "B_not_C": int((lat.B_safe_path & ~lat.C_any_order_safe).sum()),
        "examples_A_not_B": lat[lat.A_final_stable & ~lat.B_safe_path]
        .nlargest(5, "replaced_pg_mw")[["members", "replaced_pg_mw", "replaced_sn_mva"]]
        .to_dict("records"),
        "examples_A_not_C_high_mw": lat[lat.A_final_stable & ~lat.C_any_order_safe]
        .nlargest(5, "replaced_pg_mw")[
            ["members", "replaced_pg_mw", "replaced_sn_mva", "B_safe_path"]
        ]
        .to_dict("records"),
        "P1_best_C": best(lat.C_any_order_safe),
        "P1_best_B": best(lat.B_safe_path),
        "P1_best_A": best(lat.A_final_stable),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC10_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    import json

    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
