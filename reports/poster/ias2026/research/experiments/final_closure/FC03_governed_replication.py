"""FC03 (Decision 2): robustness replication on the governed documented model.

Model version: ieee39_governed_documented_v1.

Question: do policy-dependent incompatibility and witness changes survive when the
source's own documented primary frequency control (TGOV1N, sheet TGOV1N of the
same workbook) is present? Nothing is tuned; D stays 0 as documented.

Cases (both the frozen model, transverse quotient by C = span{R_x, w}, and the
governed model, quotient by C = span{R_x}):
  A  the six F8 policy points x 16 flagship subsets (base, 15 proper, flagship)
  B  one F7 policy path: the pure-g line through P4 (k = 1.425, t = 1.5, h = 1)
  C  the tongue line k = 1.30 (the G2 line), same g grid
  D  a coarse F7A plane (t = 1.5, h = 1): g in linspace(0, 1, 21)^2, k in
     linspace(0.5, 2.3, 19)
  E  condenser mitigation (G1 damped condenser, D = 2, R_0010000) on the P4
     flagship at ratings 0, 1, 2, 2.48, 3, 5 percent of the retired rating
Classification: transverse spectrum; subsets with min |Re| < 2e-3 go through the
frozen four-state classifier with the two-scale Jacobian error.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _fc import RESULTS, WORKERS, FCExperiment, out_dir, write_json

from _f7_common import SUBSETS, Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from G1_f8_rhp_followup import _em_config  # noqa: E402
from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.diagnosis.composability import hypergraph_label  # noqa: E402
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.governed import MODEL_VERSION, govern  # noqa: E402

OUT = out_dir("FC03_governed_replication")
MARGIN = 2e-3
_CFG: dict = {}


def _init():
    pin_blas_threads()
    _CFG.update({c["name"]: c for c in r_configs()})


def _reduced(j):
    return j.fx - j.fz @ np.linalg.solve(j.gz, j.gx)


def status(dae, x, z, governed):
    j1 = central_difference_jacobians(dae, x, z, {})
    j2 = central_difference_jacobians(dae, x, z, {}, scale_x=2.0, scale_z=2.0)
    a = _reduced(j1)
    d = a - _reduced(j2)
    r_x, _ = rotation_generator(dae, z)
    w = None if governed else frequency_partner(dae).w
    # dead states (identically zero rows: the frozen EMFs of the condenser cases) are
    # deleted before the quotient, exactly as certification.physical does
    dead = np.flatnonzero(~np.any(a != 0.0, axis=1))
    if dead.size:
        keep = np.setdiff1d(np.arange(a.shape[0]), dead)
        assert np.all(r_x[dead] == 0.0)
        a, d, r_x = a[np.ix_(keep, keep)], d[np.ix_(keep, keep)], r_x[keep]
        w = None if w is None else w[keep]
    tr = transverse_operator(a, r_x, w)
    ev = np.linalg.eigvals(tr.a_perp)
    alpha = float(ev.real.max())
    if np.abs(ev.real).min() < MARGIN:
        st = classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
    else:
        st = "UNSTABLE" if (ev.real > 0).any() else "STABLE"
    return st, alpha


def point_task(task):
    tag, g, k, cond = task
    th = Theta(g, k, 1.5, 1.0)
    cfg = _CFG["R_none"] if not cond else _em_config(_CFG["R_0010000"], cond)
    out = []
    for m in SUBSETS:
        case = solve_config(m, th, cfg)
        st_f, al_f = status(case.dae, case.equilibrium.x, case.equilibrium.z, False)
        gdae, gx, gz = govern(case)
        st_g, al_g = status(gdae, gx, gz, True)
        out.append(
            {
                "tag": tag,
                "g": g,
                "k": k,
                "condenser": cond,
                "subset": "+".join(map(str, m)) or "BASE",
                "status_frozen": st_f,
                "alpha_frozen": al_f,
                "status_governed": st_g,
                "alpha_governed": al_g,
            }
        )
    return out


def h_of(block, col):
    st = {
        tuple(int(b) for b in s.split("+")) if s != "BASE" else (): v
        for s, v in zip(block.subset, block[col], strict=True)
    }
    if st[()] != "STABLE":
        return f"BASE_{st[()]}", np.nan
    uns = [s for s, v in st.items() if v == "UNSTABLE"]
    unr = [s for s, v in st.items() if v not in ("STABLE", "UNSTABLE")]
    mins = [s for s in uns if not any(set(r) < set(s) for r in uns)]
    lab = hypergraph_label([frozenset(e) for e in mins]) + ("?" if unr else "")
    return lab, (min(len(e) for e in mins) if mins else -1)


def main(argv) -> int:
    exp = FCExperiment(
        name="FC03_governed_replication",
        question=(
            "Does policy-dependent incompatibility survive documented primary "
            "frequency control?"
        ),
        config={"model": MODEL_VERSION, "margin": MARGIN},
        workers=WORKERS,
    )
    started = time.time()
    f8 = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    gline = np.linspace(0.0, 1.0, 41) ** 2
    tasks = [(f"F8 {n}", p["g"], p["k"], None) for n, p in f8.items()]
    tasks += [("PATH k=1.425", float(g), 1.425, None) for g in gline]
    tasks += [("TONGUE k=1.30", float(g), 1.30, None) for g in gline]
    tasks += [
        ("PLANE", float(g), float(k), None)
        for g in np.linspace(0.0, 1.0, 21) ** 2
        for k in np.linspace(0.5, 2.3, 19)
    ]
    tasks += [
        ("CONDENSER P4", f8["P4"]["g"], f8["P4"]["k"], r)
        for r in (None, 0.01, 0.02, 0.0248, 0.03, 0.05)
    ]
    rows = []
    if "--condenser-only" in argv:
        # rerun only the condenser block after the dead-state fix (no other task has
        # dead states, so their rows are unchanged by it)
        old = pd.read_csv(OUT / "FC03_subsets.csv.gz")
        rows = old[old.tag != "CONDENSER P4"].to_dict("records")
        tasks = [t for t in tasks if t[0] == "CONDENSER P4"]
    with Pool(min(WORKERS, len(tasks)), initializer=_init) as pool:
        for out in pool.imap_unordered(point_task, tasks, chunksize=1):
            rows += out
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "FC03_subsets.csv.gz", index=False)
    pts = []
    for (tag, g, k, cond), blk in frame.groupby(
        ["tag", "g", "k", "condenser"], dropna=False
    ):
        hf, kf = h_of(blk, "status_frozen")
        hg, kg = h_of(blk, "status_governed")
        flag = blk[blk.subset == "30+33+35+37"].iloc[0]
        pts.append(
            {
                "tag": tag,
                "g": g,
                "k": k,
                "condenser": cond,
                "H_frozen": hf,
                "kappa_frozen": kf,
                "H_governed": hg,
                "kappa_governed": kg,
                "alpha_flag_frozen": flag.alpha_frozen,
                "alpha_flag_governed": flag.alpha_governed,
            }
        )
    pts = pd.DataFrame(pts)
    pts.to_csv(OUT / "FC03_points.csv", index=False)
    plane = pts[pts.tag == "PLANE"]
    summary = {
        "model": MODEL_VERSION,
        "f8": pts[pts.tag.str.startswith("F8")][
            [
                "tag",
                "H_frozen",
                "H_governed",
                "alpha_flag_frozen",
                "alpha_flag_governed",
            ]
        ].to_dict("records"),
        "path_distinct_H_frozen": sorted(
            pts[pts.tag == "PATH k=1.425"].H_frozen.unique().tolist()
        ),
        "path_distinct_H_governed": sorted(
            pts[pts.tag == "PATH k=1.425"].H_governed.unique().tolist()
        ),
        "tongue_distinct_H_frozen": sorted(
            pts[pts.tag == "TONGUE k=1.30"].H_frozen.unique().tolist()
        ),
        "tongue_distinct_H_governed": sorted(
            pts[pts.tag == "TONGUE k=1.30"].H_governed.unique().tolist()
        ),
        "plane_points": int(len(plane)),
        "plane_distinct_H_frozen": int(plane.H_frozen.nunique()),
        "plane_distinct_H_governed": int(plane.H_governed.nunique()),
        "plane_H_governed_counts": plane.H_governed.value_counts().to_dict(),
        "plane_H_frozen_counts": plane.H_frozen.value_counts().to_dict(),
        "plane_same_H": int((plane.H_frozen == plane.H_governed).sum()),
        "condenser": pts[pts.tag == "CONDENSER P4"][
            [
                "condenser",
                "H_frozen",
                "H_governed",
                "alpha_flag_frozen",
                "alpha_flag_governed",
            ]
        ].to_dict("records"),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC03_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
