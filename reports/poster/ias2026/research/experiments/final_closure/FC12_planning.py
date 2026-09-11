"""FC12 (steps 14-15): spectral vs nonlinear safe planning, minimum synchronous support.

The declared objectives (final_nonlinear_composability_v1.yaml, `planning`,
committed before any threshold was computed):

    P1  maximize replaced_pg_mw s.t. no H_RHP_perp hyperedge contained
    P2  the same with, in addition, no H_NL(rho*) hyperedge contained
        (per family, and for all three families jointly)
    P3  minimum damped-condenser rating (G1 R_0010000, D = 2, common fraction at
        the four core buses; MVA = fraction x 4270.7) making the P4 flagship
        (a) transversely stable and (b) tolerant of rho* for each family

A / B / C (final stable / safe one-at-a-time path / any order safe) are
classified on the 16-subset lattice, spectrally and nonlinearly.
If rho* >= rho_scope(family), the nonlinear status is OUT_OF_MODEL_SCOPE and P2
is not evaluated for that family: it is never replaced by a guess.
MW = replaced_pg_mw [MW], sync MVA [MVA]; the axes are never merged.
"""

from __future__ import annotations

import json
import sys
import time
from collections import deque
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
import yaml
from _fc import CONFIGS, CORE, WORKERS, FCExperiment, out_dir, write_json

from _f7_common import Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from G1_f8_rhp_followup import _em_config  # noqa: E402
from ibr_cycles.certification.symmetry import rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.nonlinear.tds import RECOVERS, pulse_dae, simulate  # noqa: E402

CFG = yaml.safe_load(
    (CONFIGS / "final_nonlinear_composability_v1.yaml").read_text(encoding="utf-8")
)
OUT = out_dir("FC12_planning")
SUBSETS = [tuple(sorted(s)) for k in range(5) for s in combinations(CORE, k)]
SN_CORE = 4270.7


def lab(s):
    return "+".join(map(str, s)) or "BASE"


def classify(unsafe):
    """A / B / C for every target, given the unsafe set."""

    stable = {s: s not in unsafe for s in SUBSETS}
    reach, queue = {(): stable[()]}, deque([()] if stable[()] else [])
    while queue:
        s = queue.popleft()
        for c in CORE:
            if c not in s:
                t = tuple(sorted(s + (c,)))
                if stable[t] and t not in reach:
                    reach[t] = True
                    queue.append(t)
    out = {}
    for t in SUBSETS:
        out[t] = {
            "A": stable[t],
            "B": bool(reach.get(t)),
            "C": all(
                stable[tuple(sorted(q))]
                for k in range(len(t) + 1)
                for q in combinations(t, k)
            ),
        }
    return out


def best(classes, pg, key):
    ok = [t for t in SUBSETS if classes[t][key]]
    if not ok:
        return None
    t = max(ok, key=lambda s: pg[s])
    return {"portfolio": lab(t), "replaced_pg_mw": pg[t]}


# ------------------------------------------------------------------ P3 --

_C: dict = {}


def _init():
    pin_blas_threads()
    _C.update({c["name"]: c for c in r_configs()})


def condenser_case(fraction):
    p = CFG["policy_points"]["P4"]
    return solve_config(
        CORE,
        Theta(p["g"], p["k"], p["t"], p["h"]),
        _em_config(_C["R_0010000"], fraction),
    )


def alpha_c(case):
    a = case.system.A
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    # dead states (the condenser's frozen EMFs: identically zero rows) are deleted
    # before the quotient, as certification.physical does; D = 2 condenser: no
    # frequency partner, C = span{R_x}
    keep = np.flatnonzero(np.any(a != 0.0, axis=1))
    a, r_x = a[np.ix_(keep, keep)], r_x[keep]
    return float(np.linalg.eigvals(transverse_operator(a, r_x, None).a_perp).real.max())


def tolerates(case, family, rho):
    spec = CFG["disturbances"][family]
    for a in [g for g in spec["continuation_grid"] if g < rho] + [rho]:
        s = simulate(case, pulse_dae(case, family, float(a), spec))["summary"]
        if s["label"] != RECOVERS:
            return False, s["label"], float(a)
    return True, RECOVERS, float(rho)


def p3_task(task):
    kind, family = task
    lo, hi = 0.0, 0.25
    if kind == "spectral":
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if alpha_c(condenser_case(mid)) < 0:
                hi = mid
            else:
                lo = mid
        return {
            "requirement": "transversely stable",
            "family": "",
            "fraction": hi,
            "condenser_mva": hi * SN_CORE,
            "status": "BISECTED",
        }
    rho = CFG["planning"]["rho_star"][family]
    ok, label, _ = tolerates(condenser_case(hi), family, rho)
    if not ok:
        return {
            "requirement": f"tolerates rho*={rho}",
            "family": family,
            "fraction": np.nan,
            "condenser_mva": np.nan,
            "status": f"NOT_REACHED at 25% ({label})",
        }
    lo = task[2] if len(task) > 2 else 0.0
    labels = []
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        ok, label, a = tolerates(condenser_case(mid), family, rho)
        labels.append((mid, label, a))
        if ok:
            hi = mid
        else:
            lo = mid
    return {
        "requirement": f"tolerates rho*={rho}",
        "family": family,
        "fraction": hi,
        "condenser_mva": hi * SN_CORE,
        "status": "BISECTED",
        "trail": json.dumps(labels),
    }


def main(argv) -> int:
    exp = FCExperiment(
        name="FC12_planning",
        question="Does nonlinear composability change the planning decision?",
        config={"planning": CFG["planning"]},
        workers=WORKERS,
    )
    started = time.time()
    if "--p3-spectral-only" in argv:
        # recompute only the spectral P3 row after the dead-state fix in alpha_c
        # (the nonlinear P3 rows use the TDS and are unaffected)
        _init()
        row = p3_task(("spectral", ""))
        p3_df = pd.read_csv(OUT / "FC12_P3_support.csv")
        p3_df = pd.concat(
            [pd.DataFrame([row]), p3_df[p3_df.requirement != "transversely stable"]],
            ignore_index=True,
        )
        p3_df.to_csv(OUT / "FC12_P3_support.csv", index=False)
        summary = json.loads((OUT / "FC12_summary.json").read_text(encoding="utf-8"))
        summary["P3"] = p3_df.to_dict("records")
        summary["p3_spectral_rerun"] = (
            "dead-state deletion before the quotient (alpha_c fix)"
        )
        write_json(OUT / "FC12_summary.json", summary)
        exp.finish("COMPUTED", **summary)
        print(json.dumps(row, indent=1, default=str))
        return 0
    thr = pd.read_csv(out_dir("FC05_nonlinear_thresholds") / "FC05_thresholds.csv")
    res6 = json.loads(
        (out_dir("FC06_resilience_complex") / "FC06_summary.json").read_text(
            encoding="utf-8"
        )
    )
    plans, abc_rows = [], []
    for point in ("P4", "P_inf"):
        blk = thr[(thr.point == point) & (thr.family == "D1")]
        pg = {
            tuple(sorted(int(b) for b in s.split("+"))) if s != "BASE" else (): v
            for s, v in zip(blk.subset, blk.replaced_pg_mw, strict=True)
        }
        spec_unsafe = {
            tuple(sorted(int(b) for b in s.split("+"))) if s != "BASE" else ()
            for s, a in zip(blk.subset, blk.alpha_perp, strict=True)
            if a >= 0
        }
        cls = classify(spec_unsafe)
        plans.append(
            {
                "point": point,
                "constraint": "P1 spectral (H_RHP_perp)",
                **(best(cls, pg, "C") or {}),
                "final_stable_best": best(cls, pg, "A"),
            }
        )
        for t in SUBSETS:
            abc_rows.append(
                {
                    "point": point,
                    "criterion": "spectral",
                    "target": lab(t),
                    **cls[t],
                    "replaced_pg_mw": pg[t],
                }
            )
        joint_unsafe, joint_valid = set(spec_unsafe), True
        for fam in ("D1", "D2", "D3"):
            rho = CFG["planning"]["rho_star"][fam]
            a6 = res6[f"{point}|{fam}"]
            fb = thr[(thr.point == point) & (thr.family == fam)]
            if rho >= a6["rho_scope"]:
                plans.append(
                    {
                        "point": point,
                        "constraint": f"P2 nonlinear {fam} rho*={rho}",
                        "status": f"OUT_OF_MODEL_SCOPE (rho_scope={a6['rho_scope']})",
                    }
                )
                joint_valid = False
                continue
            unsafe = {
                tuple(sorted(int(b) for b in s.split("+"))) if s != "BASE" else ()
                for s, r in zip(fb.subset, fb.r, strict=True)
                if not (r > rho)
            }
            c2 = classify(unsafe)
            joint_unsafe |= unsafe
            plans.append(
                {
                    "point": point,
                    "constraint": f"P2 nonlinear {fam} rho*={rho}",
                    **(best(c2, pg, "C") or {"portfolio": "NONE"}),
                    "status": "EVALUATED",
                }
            )
            for t in SUBSETS:
                abc_rows.append(
                    {
                        "point": point,
                        "criterion": f"nonlinear {fam}",
                        "target": lab(t),
                        **c2[t],
                        "replaced_pg_mw": pg[t],
                    }
                )
        if joint_valid:
            cj = classify(joint_unsafe)
            plans.append(
                {
                    "point": point,
                    "constraint": "P2 nonlinear all families",
                    **(best(cj, pg, "C") or {"portfolio": "NONE"}),
                    "status": "EVALUATED",
                }
            )
    plan_df = pd.DataFrame(plans)
    plan_df.to_csv(OUT / "FC12_plans.csv", index=False)
    pd.DataFrame(abc_rows).to_csv(OUT / "FC12_abc.csv", index=False)
    with Pool(min(WORKERS, 4), initializer=_init) as pool:
        p3 = pool.map(
            p3_task,
            [("spectral", "")] + [("nonlinear", f) for f in ("D1", "D2", "D3")],
            chunksize=1,
        )
    p3_df = pd.DataFrame(p3)
    p3_df.to_csv(OUT / "FC12_P3_support.csv", index=False)
    summary = {
        "plans": plan_df.to_dict("records"),
        "P3": p3,
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC12_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
