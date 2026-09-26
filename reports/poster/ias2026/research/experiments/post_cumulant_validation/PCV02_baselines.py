# ruff: noqa: E501  -- table labels kept on one line
"""PCV02 - Phase 2/3: baseline recomputation with the current transverse truth labels.

Preregistration 5d0b1986, sections 2-5 and stopping rules S1/S3.

Core lattices (16 subsets of H4) at P4, G_S, G_S2, P_inf and at the 52 + 36 frozen F10
points: B0 truth (FC01 direct path), B3 modal sensitivity, B4/B5/B5b Moebius
truncations of the exact alpha_perp, B8/B9 closure (non-oracle / oracle), static B1/B2/B7.
Census (E12 policy): frozen FC10 transverse truth; B1, B2, B3, B4, B5, B5b, B6, B7, B8, B9.
Tasks A-D metrics. Task E is PCV04.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _pcv import (
    CORE,
    F10_BASELINES,
    FC03_POINTS,
    FC03_SUBSETS,
    FC10_CENSUS,
    POINTS,
    SUBSETS4,
    Realization,
    Static,
    b3_prediction,
    band_critical,
    closure_scan,
    hypergraph,
    label,
    mobius_truncation,
    modal_sensitivity,
    out_dir,
    participation,
    predicted_hypergraph,
    truth_task,
    write_json,
)
from scipy.stats import spearmanr

from _overnight import pin_blas_threads

OUT = out_dir("PCV02")
SEED = 20260916
PERM = 10000
CENSUS_CANDIDATES = (30, 31, 32, 33, 34, 35, 36, 37, 38)


# ------------------------------------------------------------------ workers --
def _truth(task):
    pin_blas_threads()
    return truth_task(task)


def _point_extras(task):
    """Per point: B3 modal sensitivity and the closure scan (needs the B0 critical freq)."""

    pin_blas_threads()
    pid, theta, crit_hz = task
    ms = modal_sensitivity(theta, CORE)
    r = Realization(CORE, theta)
    closure = closure_scan(r, CORE, SUBSETS4, crit_hz)
    return pid, {"lam0": ms["lam0"], "dlam": ms["dlam"]}, closure


def _census_task(task):
    pin_blas_threads()
    name, crit_hz = task
    members = tuple(int(b) for b in name.split("+"))
    r = Realization(members, None)
    return name, closure_scan(r, members, [members], {name: crit_hz})[name]


def _census_crit(name):
    """Band-critical frequency of a census portfolio from its C2 vertex spectrum."""

    pin_blas_threads()
    members = tuple(int(b) for b in name.split("+"))
    r = Realization(members, None)
    lam = band_critical(r.vertex_spectrum(members))
    return name, lam.imag / (2 * np.pi), lam.real


# ------------------------------------------------------------------ metrics --
def auc(score, y):
    pos, neg = score[y], score[~y]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float(
        (
            (pos[:, None] > neg[None, :]).sum()
            + 0.5 * (pos[:, None] == neg[None, :]).sum()
        )
        / (len(pos) * len(neg))
    )


def average_precision(score, y):
    order = np.argsort(-score, kind="mergesort")
    ys = y[order]
    if ys.sum() == 0:
        return float("nan")
    hits = np.cumsum(ys)
    precision = hits / np.arange(1, len(ys) + 1)
    return float((precision * ys).sum() / ys.sum())


def task_a(df, score_cols, truth_col="unstable", alpha_col="alpha", rng_seed=SEED):
    y = df[truth_col].to_numpy(bool)
    rng = np.random.default_rng(rng_seed)
    perms = [rng.permutation(y) for _ in range(PERM)]
    rows = []
    for name, sign in score_cols.items():
        s = sign * df[name].to_numpy(float)
        ok = np.isfinite(s)
        a = auc(s[ok], y[ok])
        null = np.array([auc(s[ok], p[ok]) for p in perms])
        rows.append(
            {
                "method": name,
                "n": int(ok.sum()),
                "n_unstable": int(y[ok].sum()),
                "roc_auc": a,
                "average_precision": average_precision(s[ok], y[ok]),
                "spearman_with_alpha": float(
                    spearmanr(s[ok], df[alpha_col].to_numpy(float)[ok]).statistic
                ),
                "permutation_p": float((np.sum(null >= a) + 1) / (PERM + 1)),
            }
        )
    return pd.DataFrame(rows)


def task_b(true_h, pred_h):
    t_edges = (
        []
        if true_h in ("EMPTY",) or true_h.startswith("BASE")
        else [set(e.split("+")) for e in true_h.split("|")]
    )
    p_edges = (
        []
        if pred_h in ("EMPTY",) or pred_h.startswith("BASE")
        else [set(e.split("+")) for e in pred_h.split("|")]
    )
    tu = set().union(*t_edges) if t_edges else set()
    pu = set().union(*p_edges) if p_edges else set()
    kt = min((len(e) for e in t_edges), default=np.inf)
    kp = min((len(e) for e in p_edges), default=np.inf)
    return {
        "H_exact": true_h == pred_h,
        "kappa_exact": kt == kp,
        "kappa_error": abs(kt - kp)
        if np.isfinite(kt) and np.isfinite(kp)
        else (0.0 if kt == kp else np.inf),
        "false_units": len(pu - tu),
        "missed_units": len(tu - pu),
    }


# --------------------------------------------------------------------- main --
def main(argv) -> int:
    started = time.time()
    static = Static()

    # ---------------- point list -------------------------------------------------
    f10 = pd.read_csv(F10_BASELINES)
    points = dict(POINTS)
    meta = {pid: {"set": "prereg", "H_F10_old": None} for pid in POINTS}
    for i, r in f10.iterrows():
        pid = f"F10_{i:02d}_{r['set']}"
        points[pid] = (float(r.g), float(r.k), float(r.t), float(r.h))
        meta[pid] = {
            "set": r["set"],
            "H_F10_old": r["H_true"],
            "B5_old": r.get("B5_H_pred"),
            "B5b_add_old": r.get("B5b_add_H_pred"),
            "B5b_pair_old": r.get("B5b_pair_H_pred"),
        }

    # ---------------- B0 truth ---------------------------------------------------
    tasks = [(pid, th, s) for pid, th in points.items() for s in SUBSETS4]
    with Pool(16) as pool:
        truth_rows = pool.map(_truth, tasks, chunksize=4)
    truth = pd.DataFrame(truth_rows)
    truth.to_csv(OUT / "PCV02_truth_core.csv", index=False)

    # ---------------- S1: frozen FC03 check --------------------------------------
    fc03 = pd.read_csv(FC03_SUBSETS)
    s1 = []
    for pid, tag in (("P4", "F8 P4"), ("P_inf", "F8 P_inf")):
        fz = fc03[fc03.tag == tag].set_index("subset")["alpha_frozen"]
        tr = truth[truth.point == pid].set_index("subset")["alpha"]
        s1.append(float((tr - fz.reindex(tr.index)).abs().max()))
    fc03p = pd.read_csv(FC03_POINTS)
    for pid, g in (("G_S", 0.25), ("G_S2", 1.0)):
        fz = fc03p[(fc03p.tag == "PATH k=1.425") & (np.isclose(fc03p.g, g))][
            "alpha_flag_frozen"
        ].iloc[0]
        tr = truth[(truth.point == pid) & (truth.subset == label(CORE))]["alpha"].iloc[
            0
        ]
        s1.append(abs(tr - fz))
    s1_max = float(max(s1))
    if s1_max > 1e-6:
        write_json(
            OUT / "PCV02_STOP.json", {"rule": "S1", "max_abs_alpha_difference": s1_max}
        )
        print("STOP S1: recomputed alpha differs from frozen FC03 by", s1_max)
        return 2

    # ---------------- per-point extras (B3, B8, B9) -------------------------------
    crit = {
        pid: dict(zip(d.subset, d.crit_hz, strict=True))
        for pid, d in truth.groupby("point")
    }
    with Pool(16) as pool:
        extras = pool.map(
            _point_extras,
            [(pid, th, crit[pid]) for pid, th in points.items()],
            chunksize=1,
        )
    ext = {pid: (ms, cl) for pid, ms, cl in extras}

    # ---------------- assemble core table ----------------------------------------
    rows, hrows = [], []
    for pid, th in points.items():
        d = truth[truth.point == pid].set_index("subset")
        statuses = {s: d.loc[label(s), "status"] for s in SUBSETS4}
        alpha = {s: float(d.loc[label(s), "alpha"]) for s in SUBSETS4}
        h_true = hypergraph(statuses)
        ms, cl = ext[pid]
        preds = {
            "B3": {s: b3_prediction(ms, s)[0] for s in SUBSETS4},
            "B4": mobius_truncation(alpha, 1),
            "B5": mobius_truncation(alpha, 2),
            "B5b": mobius_truncation(alpha, 3),
        }
        b3_verdict = {s: b3_prediction(ms, s)[1] for s in SUBSETS4}
        hr = {
            "point": pid,
            "set": meta[pid]["set"],
            "g": th[0],
            "k": th[1],
            "t": th[2],
            "h": th[3],
            "H_true": h_true["H"],
            "kappa_true": h_true["kappa"],
            "exact": h_true["exact"],
            "H_F10_old": meta[pid].get("H_F10_old"),
        }
        for mname, pa in preds.items():
            if mname == "B3":
                ph = hypergraph(
                    {s: ("UNSTABLE" if b3_verdict[s] else "STABLE") for s in SUBSETS4}
                )
            else:
                ph = predicted_hypergraph(pa)
            hr[f"{mname}_H"] = ph["H"]
            for k, v in task_b(h_true["H"], ph["H"]).items():
                hr[f"{mname}_{k}"] = v
        hrows.append(hr)
        for s in SUBSETS4:
            if not s:
                continue
            st = static.row(s)
            rec = {
                "point": pid,
                "set": meta[pid]["set"],
                "subset": label(s),
                "size": len(s),
                "status": statuses[s],
                "unstable": statuses[s] == "UNSTABLE",
                "alpha": alpha[s],
                "crit_hz": d.loc[label(s), "crit_hz"],
                **st,
            }
            for mname, pa in preds.items():
                rec[f"{mname}_alpha"] = pa[s]
            rec["B3_verdict_unstable"] = b3_verdict[s]
            rec.update(cl.get(label(s), {}))
            rows.append(rec)
    core = pd.DataFrame(rows)
    core.to_csv(OUT / "PCV02_core_portfolios.csv", index=False)
    hyper = pd.DataFrame(hrows)
    hyper.to_csv(OUT / "PCV02_core_hypergraphs.csv", index=False)

    # ---------------- census ------------------------------------------------------
    cen = pd.read_csv(FC10_CENSUS)
    cen["tuple"] = [
        tuple(sorted(int(b) for b in m.split("+"))) if m != "BASE" else ()
        for m in cen.members
    ]
    alpha_c = dict(zip(cen.tuple, cen.alpha_perp, strict=True))
    status_c = dict(zip(cen.tuple, cen.status, strict=True))
    h_c = hypergraph(
        {
            s: (
                "STABLE"
                if v == "STABLE"
                else "UNSTABLE"
                if v == "UNSTABLE"
                else "BOUNDARY_OR_UNRESOLVED"
            )
            for s, v in status_c.items()
        }
    )
    ms_c = modal_sensitivity(None, CENSUS_CANDIDATES)
    pf_c = participation(None, CENSUS_CANDIDATES)
    trunc = {r: mobius_truncation(alpha_c, r) for r in (1, 2, 3)}
    sizes = (4, 5, 6)
    names = [label(s) for s in cen.tuple if len(s) in sizes]
    with Pool(16) as pool:
        crit_c = {
            n: (hz, re) for n, hz, re in pool.map(_census_crit, names, chunksize=2)
        }
        clos_c = dict(
            pool.map(_census_task, [(n, crit_c[n][0]) for n in names], chunksize=2)
        )
    crows = []
    for s in cen.tuple:
        if len(s) not in sizes:
            continue
        n = label(s)
        rec = {
            "portfolio": n,
            "size": len(s),
            "status": status_c[s],
            "unstable": status_c[s] == "UNSTABLE",
            "alpha": alpha_c[s],
            **static.row(s),
            "pg_mw_frozen": float(cen.loc[cen.members == n, "replaced_pg_mw"].iloc[0]),
            "sn_mva_frozen": float(
                cen.loc[cen.members == n, "replaced_sn_mva"].iloc[0]
            ),
            "B3_alpha": b3_prediction(ms_c, s)[0],
            "B4_alpha": trunc[1][s],
            "B5_alpha": trunc[2][s],
            "B5b_alpha": trunc[3][s],
            "crit_hz_c2": crit_c[n][0],
            **clos_c[n],
        }
        crows.append(rec)
    census = pd.DataFrame(crows)
    census.to_csv(OUT / "PCV02_census.csv", index=False)
    # census minimum-set recovery (all 512)
    cen_b = {}
    for mname, pa in (
        ("B4", trunc[1]),
        ("B5", trunc[2]),
        ("B5b", trunc[3]),
        ("B3", {s: b3_prediction(ms_c, s)[0] for s in cen.tuple}),
    ):
        if mname == "B3":
            ph = hypergraph(
                {
                    s: ("UNSTABLE" if b3_prediction(ms_c, s)[1] else "STABLE")
                    for s in cen.tuple
                }
            )
        else:
            ph = predicted_hypergraph(pa)
        cen_b[mname] = {"H_pred": ph["H"], **task_b(h_c["H"], ph["H"])}
    true_edges = (
        [set(map(int, e.split("+"))) for e in h_c["H"].split("|")]
        if h_c["H"] not in ("EMPTY",)
        else []
    )
    top = sorted(CENSUS_CANDIDATES, key=lambda b: -pf_c[b])
    localization = {
        "participation": pf_c,
        "ranking": top,
        "hyperedges_equal_to_top_k_participation": sum(
            set(top[: len(e)]) == e for e in true_edges
        ),
        "n_hyperedges": len(true_edges),
    }

    # ---------------- Task A --------------------------------------------------------
    orient = {
        "pg_mw": 1,
        "sn_mva": 1,
        "penetration": 1,
        "gscr": -1,
        "min_scr": -1,
        "max_miif": 1,
        "compactness_score": 1,
        "B3_alpha": 1,
        "B4_alpha": 1,
        "B5_alpha": 1,
        "B5b_alpha": 1,
        "closure_nonoracle": -1,
        "closure_oracle": -1,
        "abs_chi_nonoracle": 1,
        "abs_chi_oracle": 1,
    }
    ta = []
    for size in sizes:
        t = task_a(census[census["size"] == size].reset_index(drop=True), orient)
        t.insert(0, "dataset", f"census_size{size}")
        ta.append(t)
    f10pool = core[core.set != "prereg"].reset_index(drop=True)
    t = task_a(f10pool, orient)
    t.insert(0, "dataset", "F10_points_pooled")
    ta.append(t)
    task_a_df = pd.concat(ta, ignore_index=True)
    task_a_df.to_csv(OUT / "PCV02_taskA.csv", index=False)

    # ---------------- Task B / C ---------------------------------------------------
    tb = []
    for dataset, sel in (
        ("prereg_points", hyper.set == "prereg"),
        ("F10_52", ~hyper.set.isin(["prereg", "policy_line"])),
        ("F10_line36", hyper.set == "policy_line"),
    ):
        h = hyper[sel & hyper.exact]
        for mname in ("B3", "B4", "B5", "B5b"):
            tb.append(
                {
                    "dataset": dataset,
                    "method": mname,
                    "n_points": len(h),
                    "H_exact": int(h[f"{mname}_H_exact"].sum()),
                    "kappa_exact": int(h[f"{mname}_kappa_exact"].sum()),
                    "mean_false_units": float(h[f"{mname}_false_units"].mean()),
                    "mean_missed_units": float(h[f"{mname}_missed_units"].mean()),
                    "n_unresolved_excluded": int(
                        (~hyper[sel].exact.astype(bool)).sum()
                    ),
                }
            )
    for mname, v in cen_b.items():
        tb.append(
            {
                "dataset": "census_512",
                "method": mname,
                "n_points": 1,
                "H_exact": int(v["H_exact"]),
                "kappa_exact": int(v["kappa_exact"]),
                "mean_false_units": v["false_units"],
                "mean_missed_units": v["missed_units"],
                "H_pred": v["H_pred"],
                "H_true": h_c["H"],
            }
        )
    task_b_df = pd.DataFrame(tb)
    task_b_df.to_csv(OUT / "PCV02_taskB.csv", index=False)
    tc = []
    for dataset, d in (
        ("prereg_points", core[core.set == "prereg"]),
        ("F10_52", core[~core.set.isin(["prereg", "policy_line"])]),
        ("census_4_5_6", census),
    ):
        for mname in ("B3", "B4", "B5", "B5b"):
            pred = d[f"{mname}_alpha"].to_numpy(float)
            tru = d["alpha"].to_numpy(float)
            err = pred - tru
            tc.append(
                {
                    "dataset": dataset,
                    "method": mname,
                    "n": len(d),
                    "sign_accuracy": float(
                        np.mean((pred > 0) == d["unstable"].to_numpy(bool))
                    ),
                    "mae": float(np.mean(np.abs(err))),
                    "rmse": float(np.sqrt(np.mean(err**2))),
                }
            )
    pd.DataFrame(tc).to_csv(OUT / "PCV02_taskC.csv", index=False)

    # ---------------- Task D (within-portfolio policy discrimination, F10 52) ---------
    td = []
    f10_52 = core[~core.set.isin(["prereg", "policy_line"])]
    for name, sign in orient.items():
        aucs = []
        for _, d in f10_52.groupby("subset"):
            y = d["unstable"].to_numpy(bool)
            if y.all() or (~y).all():
                continue
            s = sign * d[name].to_numpy(float)
            if not np.isfinite(s).all():
                continue
            aucs.append(auc(s, y))
        td.append(
            {
                "method": name,
                "subsets_with_both_classes": len(aucs),
                "mean_within_portfolio_auc": float(np.mean(aucs))
                if aucs
                else float("nan"),
            }
        )
    pd.DataFrame(td).to_csv(OUT / "PCV02_taskD.csv", index=False)

    # ---------------- F10 recount summary -------------------------------------------
    f10h = hyper[hyper.set != "prereg"]
    changed = f10h[f10h.H_true != f10h.H_F10_old]
    summary = {
        "prereg_commit": "5d0b1986",
        "S1_max_abs_alpha_difference_vs_FC03": s1_max,
        "points": len(points),
        "unresolved_points": int((~hyper.exact).sum()),
        "F10_truth_changed": int(len(changed)),
        "F10_truth_changed_list": changed[["point", "H_F10_old", "H_true"]].to_dict(
            "records"
        ),
        "F10_52_H_exact": {
            m: int(f10h[(f10h.set != "policy_line") & f10h.exact][f"{m}_H_exact"].sum())
            for m in ("B3", "B4", "B5", "B5b")
        },
        "F10_52_n_resolved": int(((f10h.set != "policy_line") & f10h.exact).sum()),
        "F10_line_H_exact": {
            m: int(f10h[(f10h.set == "policy_line") & f10h.exact][f"{m}_H_exact"].sum())
            for m in ("B3", "B4", "B5", "B5b")
        },
        "F10_line_n_resolved": int(((f10h.set == "policy_line") & f10h.exact).sum()),
        "census_H_true": h_c["H"],
        "census_minimum_set": cen_b,
        "census_localization": localization,
        "B3_base_mode_E12": [ms_c["lam0"].real, ms_c["lam0"].imag],
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "PCV02_summary.json", summary)
    pd.set_option("display.width", 250)
    print(hyper[hyper.set == "prereg"].to_string(index=False))
    print(task_a_df.to_string(index=False))
    print(task_b_df.to_string(index=False))
    print(pd.DataFrame(tc).to_string(index=False))
    print(pd.DataFrame(td).to_string(index=False))
    print({k: v for k, v in summary.items() if k != "F10_truth_changed_list"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
