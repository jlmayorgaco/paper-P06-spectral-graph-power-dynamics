"""Journal gate 1 - whole right half plane as the primary stability region.

    Gamma_RHP = { Re s > 0 }                         primary (small-signal safety)
    Gamma_IA  = { Re s > 0, 0.3 <= f <= 1.5 Hz }     secondary (inter-area mechanism)

The reference double zero of the angle-reference DAE is excluded by |s| > 1e-3
in both. For every point of the frozen F7 maps (IEEE-39) and the frozen F12 grids
(Kundur): H_RHP, kappa_RHP against H_IA, kappa_IA, and every discrepancy
classified, with thresholds declared here before the run:

    APERIODIC           an RHP eigenvalue with |f| < 1e-3 Hz (real)
    BAND_EDGE_ARTIFACT  an oscillatory RHP eigenvalue within 0.05 Hz outside a band
                        edge (0.25 <= f < 0.30 or 1.50 < f <= 1.55)
    OSCILLATORY         any other oscillatory RHP eigenvalue outside the band

classified on the out-of-band RHP modes of the hyperedges in H_RHP \\ H_IA.
Also: the F8 service factorial recomputed under Gamma_RHP (planning claims),
and RHP boundaries on pure-policy edges localized, classified (oscillatory
crossing at i omega, or aperiodic crossing through the origin) and, where
oscillatory, tested against the port closure.
"""

from __future__ import annotations

import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _f7_common import LABELS as F7_LABELS
from _f7_common import SUBSETS as F7_SUBSETS
from _f7_common import Theta, eigenvalues, fast_model, port_diagnostics
from _gates import GateExperiment
from _overnight import choose_workers, pin_blas_threads
from F7_policy_hypergraph import MAPS as F7_MAPS
from F7_report import points_file
from F8_service_attribution import a_configs, r_configs, solve_config
from F12_kundur import G_GRID, Y_GRID, theta_of
from F12_kundur import SUBSETS as K_SUBSETS
from F12_kundur import solve as k_solve
from ibr_cycles.diagnosis.composability import (
    hypergraph_label,
    hypergraph_order,
    incompatibility_hypergraph,
    parse_hypergraph_label,
)
from ibr_cycles.models.port_admittance import build_action_space

OUT = RESULTS / "G1"
ZERO = 1e-3
F_REAL = 1e-3
EDGE = 0.05
BAND = (0.3, 1.5)
BAD = {"BASE_UNSTABLE", "INFEASIBLE"}


def rhp_modes(values):
    return values[(values.real > 0) & (np.abs(values) > ZERO) & (values.imag >= 0)]


def classify_mode(v) -> str:
    f = abs(v.imag) / (2 * math.pi)
    if f < F_REAL:
        return "APERIODIC"
    if BAND[0] <= f <= BAND[1]:
        return "IN_BAND"
    if BAND[0] - EDGE <= f < BAND[0] or BAND[1] < f <= BAND[1] + EDGE:
        return "BAND_EDGE_ARTIFACT"
    return "OSCILLATORY"


def label_of(counts):
    return hypergraph_label(incompatibility_hypergraph(counts))


# ------------------------------------------------------------- IEEE-39 maps ----


def f7_task(task):
    name, x, y, h_rhp, h_ia = task
    theta = F7_MAPS[name].theta_xy(x, y)
    extra = [
        e
        for e in parse_hypergraph_label(h_rhp)
        if e not in set(parse_hypergraph_label(h_ia))
    ]
    classes, freqs = set(), []
    for e in extra:
        for v in rhp_modes(eigenvalues(tuple(sorted(e)), theta)):
            c = classify_mode(v)
            if c != "IN_BAND":
                classes.add(c)
                freqs.append(abs(v.imag) / (2 * math.pi))
    return {"classes": "+".join(sorted(classes)) or "NONE", "freqs": freqs}


def ieee39_maps(pool) -> dict:
    report, frames = {}, []
    for name in F7_MAPS:
        p = pd.read_csv(points_file(name), low_memory=False)
        p = p[~p.label.isin(BAD)].copy()
        rhp = p[[f"rhp_{lab}" for lab in F7_LABELS]].to_numpy()
        labels = [
            label_of({frozenset(s): int(n) for s, n in zip(F7_SUBSETS, r, strict=True)})
            for r in rhp
        ]
        p["H_RHP"] = labels
        p["kappa_RHP"] = [hypergraph_order(parse_hypergraph_label(s)) for s in labels]
        p["H_IA"], p["kappa_IA"] = p.label, p.kappa
        diff = p[p.H_RHP != p.H_IA]
        res = (
            pool.map(
                f7_task,
                [(name, r.x, r.y, r.H_RHP, r.H_IA) for r in diff.itertuples()],
                chunksize=8,
            )
            if len(diff)
            else []
        )
        p["discrepancy"] = "NO_DISCREPANCY"
        if len(diff):
            p.loc[diff.index, "discrepancy"] = [r["classes"] for r in res]
        frames.append(
            p[
                [
                    "x",
                    "y",
                    "g",
                    "k",
                    "t",
                    "h",
                    "H_IA",
                    "kappa_IA",
                    "H_RHP",
                    "kappa_RHP",
                    "discrepancy",
                ]
            ].assign(map=name)
        )
        report[name] = {
            "points": int(len(p)),
            "H_equal": float((p.H_RHP == p.H_IA).mean()),
            "kappa_equal": float((p.kappa_RHP == p.kappa_IA).mean()),
            "discrepancy_classes": p.discrepancy.value_counts().to_dict(),
            "distinct_H_RHP": int(p.H_RHP.nunique()),
            "distinct_H_IA": int(p.H_IA.nunique()),
            "kappa_RHP_values": sorted(int(k) for k in p.kappa_RHP.unique()),
            "H_RHP_EMPTY_points": int((p.H_RHP == "EMPTY").sum()),
        }
    pd.concat(frames).to_csv(OUT / "G1_ieee39_points.csv.gz", index=False)
    return report


# -------------------------------------------------------------------- Kundur ----


def k_task(task):
    m, i, j = task
    th = theta_of(m, G_GRID[i], Y_GRID[j])
    counts_rhp, counts_ia, modes = {}, {}, {}
    for s in K_SUBSETS:
        v = np.linalg.eigvals(k_solve(s, th).system.A)
        r = rhp_modes(v)
        modes[s] = r
        counts_rhp[frozenset(s)] = int(r.size)
        counts_ia[frozenset(s)] = int(sum(classify_mode(x) == "IN_BAND" for x in r))
    h_rhp, h_ia = label_of(counts_rhp), label_of(counts_ia)
    classes = set()
    for e in parse_hypergraph_label(h_rhp):
        if e in set(parse_hypergraph_label(h_ia)):
            continue
        for x in modes[tuple(sorted(e))]:
            c = classify_mode(x)
            if c != "IN_BAND":
                classes.add(c)
    return {
        "map": m,
        "i": i,
        "j": j,
        "H_RHP": h_rhp,
        "H_IA": h_ia,
        "kappa_RHP": hypergraph_order(parse_hypergraph_label(h_rhp)),
        "kappa_IA": hypergraph_order(parse_hypergraph_label(h_ia)),
        "discrepancy": "NO_DISCREPANCY"
        if h_rhp == h_ia
        else ("+".join(sorted(classes)) or "NONE"),
        "indicator_rhp": "".join(
            str(int(counts_rhp[frozenset(s)] > 0)) for s in K_SUBSETS
        ),
    }


def kundur_maps(pool) -> tuple[dict, pd.DataFrame]:
    nodes = pd.read_csv(RESULTS / "F12" / "F12_nodes.csv")
    ok = nodes[~nodes.label.isin(BAD)]
    rows = pd.DataFrame(
        pool.map(
            k_task, [(r.map, int(r.i), int(r.j)) for r in ok.itertuples()], chunksize=8
        )
    )
    rows = rows.merge(
        ok[["map", "i", "j", "g", "k", "t", "label"]], on=["map", "i", "j"]
    )
    rows.to_csv(OUT / "G1_kundur_nodes.csv", index=False)
    report = {}
    for m in ("K12A", "K12B"):
        f = rows[rows["map"] == m]
        lines = f.groupby("j").H_RHP.nunique()
        report[m] = {
            "nodes": int(len(f)),
            "H_equal": float((f.H_RHP == f.H_IA).mean()),
            "band_label_reproduced": float((f.H_IA == f.label).mean()),
            "kappa_equal": float((f.kappa_RHP == f.kappa_IA).mean()),
            "discrepancy_classes": f.discrepancy.value_counts().to_dict(),
            "distinct_H_RHP": int(f.H_RHP.nunique()),
            "H_RHP_counts": f.H_RHP.value_counts().to_dict(),
            "kappa_RHP_values": sorted(int(k) for k in f.kappa_RHP.unique()),
            "policy_lines_H_RHP_changes": int((lines > 1).sum()),
            "policy_lines": int(len(lines)),
            "H_RHP_EMPTY_nodes": int((f.H_RHP == "EMPTY").sum()),
        }
    return report, rows


# ------------------------------------------------ RHP boundaries on g-edges ----


def _k_count(members, th):
    return rhp_modes(np.linalg.eigvals(k_solve(members, th).system.A))


def k_edge(task):
    m, i, j, members = task
    ta, tb = theta_of(m, G_GRID[i], Y_GRID[j]), theta_of(m, G_GRID[i + 1], Y_GRID[j])
    lerp = lambda s: {k: ta[k] + s * (tb[k] - ta[k]) for k in ta}  # noqa: E731
    ra, rb = _k_count(members, ta), _k_count(members, tb)
    ua = ra.size > 0
    if ua == (rb.size > 0):
        return None
    lo, hi, vlo, vhi = 0.0, 1.0, ra, rb
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        r = _k_count(members, lerp(mid))
        if (r.size > 0) == ua:
            lo, vlo = mid, r
        else:
            hi, vhi = mid, r
    unsafe = vlo if ua else vhi
    th = lerp(0.5 * (lo + hi))
    full = np.linalg.eigvals(k_solve(members, th).system.A)
    lam = unsafe[np.argmin(np.abs(unsafe.real))] if unsafe.size else complex("nan")
    kind = (
        "APERIODIC_THROUGH_ORIGIN"
        if abs(lam.imag) / (2 * math.pi) < F_REAL
        else "OSCILLATORY_AXIS_CROSSING"
    )
    row = {
        "benchmark": "Kundur",
        "map": m,
        "witness": "+".join(map(str, members)),
        "g_star": th["g"],
        "y_star": th["k"] if m == "K12A" else th["t"],
        "freq_hz": abs(lam.imag) / (2 * math.pi),
        "lambda_re": float(lam.real),
        "crossing": kind,
    }
    if kind == "OSCILLATORY_AXIS_CROSSING":
        near = full[np.argmin(np.abs(full - complex(0, abs(lam.imag))))]
        space = build_action_space(k_solve((), th), k_solve(members, th), members)
        omega = abs(near.imag)
        split = space.split(complex(0.0, omega))
        ref = np.median(
            [
                abs(complex(space.split(complex(0, w))["full"]))
                for w in np.linspace(0.3 * omega, 3 * omega, 25)
            ]
        )
        row.update(
            port_visible=bool(abs(complex(split["full"])) / ref < 1e-3),
            closure_distance=float(
                np.min(np.abs(np.asarray(split["q_eigenvalues"]) + 1.0))
            )
            if len(members) > 1
            else float("nan"),
            individual_min=float(min(abs(d) for d in split["diagonals"])),
        )
    return row


def kundur_edges(pool, rows) -> pd.DataFrame:
    tasks = []
    for m in ("K12A", "K12B"):
        f = rows[rows["map"] == m].set_index(["i", "j"])
        for (i, j), r in f.iterrows():
            if (i + 1, j) not in f.index:
                continue
            a, b = r.indicator_rhp, f.loc[(i + 1, j)].indicator_rhp
            for s, x, y in zip(K_SUBSETS, a, b, strict=True):
                if s and x != y:
                    tasks.append((m, int(i), int(j), s))
    return pd.DataFrame([r for r in pool.map(k_edge, tasks, chunksize=2) if r])


def f7b_edge(task):
    ka, kb, members = task
    spec = F7_MAPS["F7B"]
    a, b = spec.theta(ka), spec.theta(kb)
    count = lambda th: rhp_modes(eigenvalues(members, th))  # noqa: E731
    ra, rb = count(a), count(b)
    ua = ra.size > 0
    if ua == (rb.size > 0):
        return None
    lo, hi, vlo, vhi = 0.0, 1.0, ra, rb
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        r = count(a.lerp(b, mid))
        if (r.size > 0) == ua:
            lo, vlo = mid, r
        else:
            hi, vhi = mid, r
    unsafe = vlo if ua else vhi
    lam = unsafe[np.argmin(np.abs(unsafe.real))]
    th = a.lerp(b, 0.5 * (lo + hi))
    kind = (
        "APERIODIC_THROUGH_ORIGIN"
        if abs(lam.imag) / (2 * math.pi) < F_REAL
        else "OSCILLATORY_AXIS_CROSSING"
    )
    row = {
        "benchmark": "IEEE-39",
        "map": "F7B",
        "witness": "+".join(map(str, members)),
        "g_star": th.g,
        "y_star": th.t,
        "freq_hz": abs(lam.imag) / (2 * math.pi),
        "lambda_re": float(lam.real),
        "crossing": kind,
    }
    if kind == "OSCILLATORY_AXIS_CROSSING":
        d = port_diagnostics(members, th, abs(lam.imag))
        row.update(
            port_visible=d["port_visible"],
            closure_distance=d["closure_distance"],
            individual_min=d["port_min_individual"],
        )
    return row


def f7b_edges(pool) -> pd.DataFrame:
    p = pd.read_csv(points_file("F7B"), low_memory=False)
    p = p[~p.label.isin(BAD)].set_index(["I", "J"])
    tasks = []
    for (i, j), r in p.iterrows():
        if (i + 1, j) not in p.index:
            continue
        q = p.loc[(i + 1, j)]
        for s, lab in zip(F7_SUBSETS, F7_LABELS, strict=True):
            if not s:
                continue
            rhp_change = (r[f"rhp_{lab}"] > 0) != (q[f"rhp_{lab}"] > 0)
            ia_change = (r[f"N_{lab}"] > 0) != (q[f"N_{lab}"] > 0)
            if rhp_change and not ia_change:
                tasks.append(((i, j), (i + 1, j), s))
    return pd.DataFrame([r for r in pool.map(f7b_edge, tasks, chunksize=2) if r])


# --------------------------------------------------------------- F8 under RHP ----


def f8_task(task):
    point, theta_d, config = task
    theta = Theta(**theta_d)
    counts = {}
    for members in F7_SUBSETS:
        try:
            v = np.linalg.eigvals(solve_config(members, theta, config).system.A)
        except Exception as error:  # noqa: BLE001
            return {
                "point": point,
                "config": config["name"],
                "H_RHP": "INFEASIBLE",
                "err": str(error)[:60],
            }
        if not members and v[np.abs(v) > ZERO].real.max() >= 0:
            return {"point": point, "config": config["name"], "H_RHP": "BASE_UNSTABLE"}
        counts[frozenset(members)] = int(rhp_modes(v).size)
    return {"point": point, "config": config["name"], "H_RHP": label_of(counts)}


def f8_rhp(pool) -> tuple[dict, pd.DataFrame]:
    pts = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    frozen = pd.read_csv(RESULTS / "F8" / "F8_interventions.csv", low_memory=False)
    configs = r_configs() + a_configs()
    tasks = [
        (n, {k: p[k] for k in ("g", "k", "t", "h")}, c)
        for n, p in pts.items()
        for c in configs
    ]
    rows = pd.DataFrame(pool.map(f8_task, tasks, chunksize=2))
    rows = rows.merge(frozen[["point", "config", "label"]], on=["point", "config"])
    rows.to_csv(OUT / "G1_f8_rhp.csv", index=False)
    same = rows.H_RHP == rows.label
    key = rows[
        rows.config.isin(
            [
                "R_none",
                "R_0000000",
                "R_1000000",
                "R_virtual_inertia",
                "A_00111",
                "A_20111",
                "A_11011",
            ]
        )
    ]
    return {
        "cases": int(len(rows)),
        "H_RHP_equals_frozen_H_IA": float(same.mean()),
        "differences": rows[~same][["point", "config", "label", "H_RHP"]].to_dict(
            "records"
        )[:40],
        "key": key[["point", "config", "label", "H_RHP"]].to_dict("records"),
    }, rows


def main(argv) -> int:
    pin_blas_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    exp = GateExperiment(
        name="G1_whole_rhp",
        question="Is whole-RHP composability what the claims say?",
        config={"F_REAL_hz": F_REAL, "BAND_EDGE_hz": EDGE, "zero_floor": ZERO},
        workers=choose_workers(),
    )
    started = time.time()
    with Pool(choose_workers(), initializer=_init) as pool:
        ieee = ieee39_maps(pool)
        kund, rows = kundur_maps(pool)
        k_edges = kundur_edges(pool, rows)
        f_edges = f7b_edges(pool)
        f8, _ = f8_rhp(pool)
    edges = pd.concat([k_edges, f_edges], ignore_index=True)
    edges.to_csv(OUT / "G1_rhp_boundaries.csv", index=False)
    osc = edges[edges.crossing == "OSCILLATORY_AXIS_CROSSING"]
    summary = {
        "IEEE39": ieee,
        "Kundur": kund,
        "F8_under_RHP": f8,
        "rhp_boundaries": {
            "located": int(len(edges)),
            "by_benchmark_and_type": edges.groupby(["benchmark", "crossing"])
            .size()
            .to_dict()
            if len(edges)
            else {},
            "oscillatory_port_visible": int(
                osc.port_visible.fillna(False).astype(bool).sum()
            )
            if len(osc)
            else 0,
            "oscillatory": int(len(osc)),
            "closure_max_multi": float(
                osc[osc.witness.str.contains(r"\+")].closure_distance.max()
            )
            if len(osc) and osc.witness.str.contains(r"\+").any()
            else None,
            "individual_max_single": float(
                osc[~osc.witness.str.contains(r"\+")].individual_min.max()
            )
            if len(osc) and (~osc.witness.str.contains(r"\+")).any()
            else None,
            "freq_range_hz": [float(edges.freq_hz.min()), float(edges.freq_hz.max())]
            if len(edges)
            else None,
        },
    }
    summary["rhp_boundaries"]["by_benchmark_and_type"] = {
        f"{a}|{b}": v
        for (a, b), v in summary["rhp_boundaries"]["by_benchmark_and_type"].items()
    }
    (OUT / "G1_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )
    exp.finish("COMPUTED", elapsed_s=time.time() - started)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "F8_under_RHP"},
            indent=1,
            default=str,
        )
    )
    print(
        json.dumps({k: v for k, v in f8.items() if k != "key"}, indent=1, default=str)[
            :3000
        ]
    )
    return 0


def _init():
    pin_blas_threads()
    fast_model()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
