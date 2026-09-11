"""FC05: critical disturbance radii r_S^d for every flagship subset (frozen config).

configs/ias2026/final_nonlinear_composability_v1.yaml, committed before this run.
For each policy point (P4, P_inf), disturbance family (D1, D2, D3) and subset S of
{30, 33, 35, 37}:

  - continuation on the declared grid until the first amplitude whose label is
    not RECOVERS;
  - bisection between the last RECOVERS and that amplitude, to the declared
    tolerance;
  - the threshold type is the label of the upper end (FAILS, OUTSIDE, NUMERICAL).

A small-signal-unstable portfolio (alpha_perp > 0) has r = 0 by definition,
checked by one run at the smallest grid amplitude. D3 on the empty portfolio is
vacuous. The just-below and just-above traces are saved.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import yaml
from _fc import CONFIGS, WORKERS, FCExperiment, out_dir, write_json

from _f7_common import SUBSETS, Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from ibr_cycles.nonlinear.tds import RECOVERS, pulse_dae, simulate, transverse_alpha  # noqa: E402
from ibr_cycles.units import measure  # noqa: E402

CFG = yaml.safe_load((CONFIGS / "final_nonlinear_composability_v1.yaml").read_text(encoding="utf-8"))
OUT = out_dir("FC05_nonlinear_thresholds")
TRACES = out_dir("FC05_nonlinear_thresholds/traces")
_C: dict = {}


def _init():
    pin_blas_threads()
    _C["cfg"] = {c["name"]: c for c in r_configs()}["R_none"]


def tol(family, a):
    return 0.01 if family == "D3" else max(2.0, 0.02 * a)


def run_one(case, family, a, spec):
    t = time.time()
    out = simulate(case, pulse_dae(case, family, a, spec))
    s = out["summary"]
    return {"amplitude": a, "label": s["label"], "reason": s["reason"],
            "t_label": s["t_label"], "cpu_s": round(time.time() - t, 1)}, out["trace"]


def task(t):
    point, family, members = t
    p = CFG["policy_points"][point]
    spec = CFG["disturbances"][family]
    case = solve_config(members, Theta(p["g"], p["k"], p["t"], p["h"]), _C["cfg"])
    alpha = transverse_alpha(case)
    q = measure(case)
    base = {"point": point, "family": family, "subset": "+".join(map(str, members)) or "BASE",
            "size": len(members), "alpha_perp": alpha, "replaced_pg_mw": q.replaced_pg_mw}
    if family == "D3" and not members:
        return {**base, "r": float("inf"), "type": "VACUOUS", "runs": []}
    grid = spec["continuation_grid"]
    runs, last_ok, first_bad, bad_label = [], 0.0, None, None
    traces = {}
    candidates = grid[:1] if alpha > 0 else grid
    for a in candidates:
        rec, tr = run_one(case, family, float(a), spec)
        runs.append(rec)
        traces[a] = tr
        if rec["label"] != RECOVERS:
            first_bad, bad_label = float(a), rec["label"]
            break
        last_ok = float(a)
    if first_bad is None:
        return {**base, "r": float("inf"), "type": "NONE_IN_RANGE", "r_lower": last_ok,
                "runs": runs}
    if alpha > 0:
        out = {**base, "r": 0.0, "type": bad_label, "bracket": [0.0, first_bad], "runs": runs}
    else:
        lo, hi, lo_tr, hi_tr = last_ok, first_bad, traces.get(last_ok), traces[first_bad]
        while hi - lo > tol(family, hi):
            mid = 0.5 * (lo + hi)
            rec, tr = run_one(case, family, mid, spec)
            runs.append(rec)
            if rec["label"] == RECOVERS:
                lo, lo_tr = mid, tr
            else:
                hi, hi_tr, bad_label = mid, tr, rec["label"]
        out = {**base, "r": 0.5 * (lo + hi), "type": bad_label, "bracket": [lo, hi], "runs": runs}
        stem = f"{point}_{family}_{base['subset']}".replace("+", "p")
        if lo_tr is not None:
            lo_tr.to_csv(TRACES / f"{stem}_below.csv.gz", index=False)
        hi_tr.to_csv(TRACES / f"{stem}_above.csv.gz", index=False)
        return out
    stem = f"{point}_{family}_{base['subset']}".replace("+", "p")
    traces[first_bad].to_csv(TRACES / f"{stem}_above.csv.gz", index=False)
    return out


def main(argv) -> int:
    exp = FCExperiment(
        name="FC05_nonlinear_thresholds",
        question="What finite disturbance does each flagship subset tolerate?",
        config={"config": CFG},
        workers=WORKERS,
    )
    started = time.time()
    tasks = [(pt, fam, m) for pt in ("P4", "P_inf") for fam in ("D1", "D2", "D3") for m in SUBSETS]
    tasks.sort(key=lambda t: -len(t[2]))  # the long runs first
    rows = []
    with Pool(WORKERS, initializer=_init) as pool:
        for k, out in enumerate(pool.imap_unordered(task, tasks, chunksize=1)):
            rows.append(out)
            print(f"  {k + 1}/{len(tasks)} {out['point']} {out['family']} {out['subset']}: "
                  f"r={out['r']} {out['type']}", flush=True)
    frame = pd.DataFrame(rows)
    frame["runs"] = frame.runs.map(json.dumps)
    frame["bracket"] = frame.get("bracket", pd.Series(dtype=object)).map(
        lambda b: json.dumps(b) if isinstance(b, list) else "")
    frame.to_csv(OUT / "FC05_thresholds.csv", index=False)
    summary = {
        "tasks": int(len(frame)),
        "types": frame.groupby(["point", "family"]).type.value_counts().to_dict(),
        "total_runs": int(frame.runs.map(lambda r: len(json.loads(r))).sum()),
        "elapsed_s": round(time.time() - started, 1),
    }
    summary["types"] = {f"{a}|{b}|{c}": v for (a, b, c), v in summary["types"].items()}
    write_json(OUT / "FC05_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
