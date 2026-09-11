"""FC01 (amendment A): transverse stability quotient and the re-audit of every label.

Part 1, structure (direct path): at the six F8 points x 16 flagship subsets, one
Kundur node x 8 subsets and two IEEE-68 cases, report dim C, the Jordan block of
A on C, the invariance, rotation and left-basis residuals, the coupling
residual Z^T A U, the spectral identity sigma(A) = {0, 0} U sigma(A_perp), the
old-cutoff RHP count (|lambda| > 1e-3) against N_RHP^perp, and every mode with
|lambda| < 1e-2 in A and in A_perp.

Part 2, F7 re-audit (fast path, exact parametric assembly _f7_common.FastModel,
validated 0/1 153 against the direct path in F7): all 345 229 frozen F7 points x
16 subsets. Transverse spectrum -> N_RHP^perp. A subset whose transverse spectrum
has min |Re| < MARGIN is classified with the frozen four-state classifier
(d_axis against 10 eps; anchor error matrices transformed to the quotient);
otherwise the sign count is exact. H_perp, kappa_perp per point, compared with
the frozen G1 H_RHP.

Part 3: every point whose H_perp is unresolved or differs from G1 and has
g <= 0.005 (the strip BC01 found non-certifiable on the fast path) is recomputed
on the direct path (direct Jacobians, same classifier).
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _fc import RESULTS, WORKERS, FCExperiment, out_dir, publish, write_json

from _f7_common import SUBSETS, Theta, fast_model, solve_subset  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from G3_ieee68 import solve as s68  # noqa: E402
from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.physical import physical_matrices  # noqa: E402
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import center_subspace, transverse_operator  # noqa: E402
from ibr_cycles.diagnosis.composability import hypergraph_label  # noqa: E402

OUT = out_dir("FC01_transverse_quotient")
MARGIN = 2e-3  # 5x the largest classifier threshold observed in BC01 (4.3e-4)
CHUNK = 400


def lab(m):
    return "+".join(map(str, m)) or "BASE"


def h_of(statuses: dict) -> dict:
    base = statuses[()]
    if base != "STABLE":
        return {"H": f"BASE_{base}", "kappa": np.nan, "exact": False, "unresolved": ""}
    unsafe = [s for s, v in statuses.items() if v == "UNSTABLE"]
    unres = [s for s, v in statuses.items() if v == "BOUNDARY_OR_UNRESOLVED"]
    minimal = [s for s in unsafe if not any(set(r) < set(s) for r in unsafe)]
    sizes = [len(e) for e in minimal]
    return {
        "H": hypergraph_label([frozenset(e) for e in minimal]),
        "kappa": min(sizes) if sizes else -1,
        "exact": not unres,
        "unresolved": "|".join(lab(u) for u in unres),
    }


# ---------------------------------------------------------------- part 1 --


def structure_row(bench, point, members, case):
    a, d, _ = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    partner = frequency_partner(case.dae)
    cs = center_subspace(a, r_x, partner.w)
    tr = transverse_operator(a, r_x, partner.w)
    full = np.linalg.eigvals(a)
    perp = np.linalg.eigvals(tr.a_perp)
    # spectral identity: remove the two eigenvalues of A closest to zero, match the rest
    order = np.argsort(np.abs(full))
    rest_full = np.sort_complex(full[order[cs.dimension :]])
    from scipy.optimize import linear_sum_assignment

    cmat = np.abs(rest_full[:, None] - perp[None, :])
    r, c = linear_sum_assignment(cmat)
    old = int(((full.real > 0) & (np.abs(full) > 1e-3)).sum())
    return {
        "benchmark": bench,
        "point": point,
        "subset": lab(members),
        "n_x": a.shape[0],
        "dim_C": cs.dimension,
        "jordan_11": cs.jordan[0, 0],
        "jordan_12": cs.jordan[0, 1] if cs.dimension == 2 else np.nan,
        "jordan_21": cs.jordan[1, 0] if cs.dimension == 2 else np.nan,
        "jordan_22": cs.jordan[1, 1] if cs.dimension == 2 else np.nan,
        "jordan_chain": cs.is_jordan_chain,
        "invariance_residual": cs.invariance_residual,
        "rotation_residual": cs.rotation_residual,
        "left_residual": cs.left_residual,
        "left_condition": cs.left_condition,
        "coupling_residual": tr.coupling_residual,
        "spectral_identity_maxdist": float(cmat[r, c].max()),
        "two_smallest_abs_eig_A": ";".join(f"{v:.2e}" for v in full[order[:2]]),
        "rhp_old_cutoff": old,
        "rhp_perp": int((perp.real > 0).sum()),
        "alpha_perp": float(perp.real.max()),
        "modes_A_abs_lt_1e-2": ";".join(f"{v:.3e}" for v in full[np.abs(full) < 1e-2]),
        "modes_perp_abs_lt_1e-2": ";".join(f"{v:.3e}" for v in perp[np.abs(perp) < 1e-2]),
        "partner": partner.reason,
    }


def part1():
    rows = []
    f8 = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    rnone = {c["name"]: c for c in r_configs()}["R_none"]
    for name, p in f8.items():
        th = Theta(p["g"], p["k"], p["t"], p["h"])
        for m in SUBSETS:
            rows.append(structure_row("IEEE-39", name, m, solve_config(m, th, rnone)))
    for m in ((), (2,), (3,), (4,), (2, 3), (2, 4), (3, 4), (2, 3, 4)):
        rows.append(structure_row("Kundur", "K12A g=0.08 k=1.25", m,
                                  k_solve(m, {"g": 0.08, "k": 1.25, "t": 1.0})))
    for m in ((), (3, 4, 6, 9)):
        rows.append(structure_row("IEEE-68", "g=0 k=1", m, s68(m, 0.0, 1.0)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- part 2 --

_W: dict = {}


def _init():
    pin_blas_threads()
    fm = fast_model()
    for m in SUBSETS:
        errs, first = [], None
        for th in (Theta(0.0), Theta(1.0), Theta(0.0, k=2.0)):
            case = solve_subset(m, th)
            _, d, _ = physical_matrices(case)
            errs.append(np.abs(d))
            first = first or case
        r_x, _ = rotation_generator(first.dae, first.equilibrium.z)
        w = frequency_partner(first.dae).w
        tr = transverse_operator(fm.matrix(m, Theta(0.3, 1.2)), r_x, w)
        err = errs[0] + errs[1] + errs[2]
        _W[m] = (tr.z, tr.z.T @ err @ tr.z, r_x, w)


def chunk_task(chunk):
    fm = fast_model()
    out = []
    for r in chunk:
        th = Theta(r["g"], r["k"], r["t"], r["h"])
        statuses, near, alphas, resid = {}, 0, {}, 0.0
        for m in SUBSETS:
            z, err_q, r_x, w = _W[m]
            a = fm.matrix(m, th)
            resid = max(resid, float(np.linalg.norm(a @ r_x) / np.linalg.norm(a, 1)))
            ap = z.T @ a @ z
            ev = np.linalg.eigvals(ap)
            alphas[m] = float(ev.real.max())
            if np.abs(ev.real).min() < MARGIN:
                near += 1
                statuses[m] = classify_spectrum(ap, err_q, SAFETY).status
            else:
                statuses[m] = "UNSTABLE" if (ev.real > 0).any() else "STABLE"
        info = h_of(statuses)
        out.append(
            {
                "idx": r["idx"],
                "map": r["map"],
                "g": r["g"],
                "k": r["k"],
                "t": r["t"],
                "h": r["h"],
                "H_perp": info["H"],
                "kappa_perp": info["kappa"],
                "exact": info["exact"],
                "unresolved": info["unresolved"],
                "near_axis_subsets": near,
                "alpha_perp_flagship": alphas[tuple(sorted((30, 33, 35, 37)))],
                "max_rotation_residual": resid,
                "H_RHP_G1": r["H_RHP"],
                "kappa_RHP_G1": r["kappa_RHP"],
            }
        )
    return out


def part2(workers):
    g1 = pd.read_csv(RESULTS / "G1" / "G1_ieee39_points.csv.gz")
    g1["idx"] = np.arange(len(g1))
    recs = g1[["idx", "map", "g", "k", "t", "h", "H_RHP", "kappa_RHP"]].to_dict("records")
    chunks = [recs[i : i + CHUNK] for i in range(0, len(recs), CHUNK)]
    rows = []
    with Pool(workers, initializer=_init) as pool:
        for k, out in enumerate(pool.imap_unordered(chunk_task, chunks, chunksize=1)):
            rows += out
            if k % 50 == 0:
                print(f"  part2 {k}/{len(chunks)}", flush=True)
    return pd.DataFrame(rows).sort_values("idx")


# ---------------------------------------------------------------- part 3 --

_R: dict = {}


def _init3():
    pin_blas_threads()
    _R["cfg"] = {c["name"]: c for c in r_configs()}["R_none"]


def direct_task(r):
    th = Theta(r["g"], r["k"], r["t"], r["h"])
    statuses = {}
    for m in SUBSETS:
        case = solve_config(m, th, _R["cfg"])
        a, d, _ = physical_matrices(case)
        r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
        w = frequency_partner(case.dae).w
        tr = transverse_operator(a, r_x, w)
        statuses[m] = classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
    info = h_of(statuses)
    return {"idx": r["idx"], "H_perp_direct": info["H"], "kappa_perp_direct": info["kappa"],
            "exact_direct": info["exact"], "unresolved_direct": info["unresolved"]}


def main(argv) -> int:
    exp = FCExperiment(
        name="FC01_transverse_quotient",
        question="Do the frozen F7 / H / kappa labels survive the exact transverse quotient?",
        config={"margin": MARGIN, "safety": SAFETY, "chunk": CHUNK, "workers": WORKERS},
        workers=WORKERS,
    )
    started = time.time()
    p1 = part1()
    p1.to_csv(OUT / "FC01_structure.csv", index=False)
    print(p1.groupby("benchmark")[["dim_C", "invariance_residual", "coupling_residual",
                                   "spectral_identity_maxdist", "left_condition"]].max())
    p2 = part2(WORKERS)
    p2["same_as_G1"] = p2.H_perp == p2.H_RHP_G1
    redo = p2[(~p2.exact | ~p2.same_as_G1) & (p2.g <= 0.005)]
    print(f"part3: {len(redo)} low-g points recomputed on the direct path", flush=True)
    if len(redo):
        with Pool(WORKERS, initializer=_init3) as pool:
            p3 = pd.DataFrame(pool.map(direct_task, redo.to_dict("records"), chunksize=1))
        p2 = p2.merge(p3, on="idx", how="left")
    p2.to_csv(OUT / "FC01_points.csv.gz", index=False)
    final_h = p2.get("H_perp_direct", pd.Series(index=p2.index, dtype=object)).fillna(p2.H_perp)
    final_exact = p2.get("exact_direct", pd.Series(index=p2.index, dtype=object)).fillna(p2.exact)
    p2["H_perp_final"], p2["exact_final"] = final_h, final_exact.astype(bool)
    p2["same_final"] = p2.H_perp_final == p2.H_RHP_G1
    diff = p2[~p2.same_final | ~p2.exact_final]
    diff.to_csv(OUT / "TSQ_ieee39_reaudit.csv", index=False)
    publish(OUT / "TSQ_ieee39_reaudit.csv")
    by_map = {}
    for mp, grp in p2.groupby("map"):
        by_map[mp] = {
            "points": int(len(grp)),
            "exact_final": int(grp.exact_final.sum()),
            "H_perp_equals_G1": int(grp.same_final.sum()),
            "differs_and_exact": int((~grp.same_final & grp.exact_final).sum()),
            "unresolved_final": int((~grp.exact_final).sum()),
            "unresolved_final_g_le_0.005": int((~grp.exact_final & (grp.g <= 0.005)).sum()),
            "near_axis_subset_evaluations": int(grp.near_axis_subsets.sum()),
            "max_rotation_residual": float(grp.max_rotation_residual.max()),
        }
    summary = {
        "structure": {
            "dim_C": sorted(p1.dim_C.unique().tolist()),
            "jordan_chain_all": bool(p1.jordan_chain.all()),
            "jordan_12_range": [float(p1.jordan_12.min()), float(p1.jordan_12.max())],
            "max_invariance_residual": float(p1.invariance_residual.max()),
            "max_coupling_residual": float(p1.coupling_residual.max()),
            "max_spectral_identity_dist": float(p1.spectral_identity_maxdist.max()),
            "max_left_condition": float(p1.left_condition.max()),
            "rhp_old_equals_perp": int((p1.rhp_old_cutoff == p1.rhp_perp).sum()),
            "rows": int(len(p1)),
        },
        "f7": by_map,
        "total_points": int(len(p2)),
        "total_equal_G1": int(p2.same_final.sum()),
        "total_unresolved_final": int((~p2.exact_final).sum()),
        "total_differs_exact": int((~p2.same_final & p2.exact_final).sum()),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC01_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
