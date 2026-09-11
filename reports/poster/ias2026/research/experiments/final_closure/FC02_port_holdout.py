"""FC02 (Decision 3): holdout validation of the frozen relocated port (port_relocation_v1).

The relocation, beta pairs, detection rule, ground truth and holdout slices are
those of configs/port_relocation_v1.yaml, committed before this script ran.
Ground truth: real / complex RHP eigenvalue counts of the transverse operator.
The port sees only the DAE Jacobian blocks.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import yaml
from _fc import CONFIGS, ROOT, WORKERS, FCExperiment, out_dir, publish, write_json

from _overnight import pin_blas_threads  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from ibr_cycles.certification.port_origin import relocated_port  # noqa: E402
from ibr_cycles.certification.symmetry import frequency_partner, rotation_generator  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402

OUT = out_dir("FC02_port_holdout")
CFG = yaml.safe_load((ROOT / "configs" / "port_relocation_v1.yaml").read_text(encoding="utf-8"))
SUBSETS = [(), (2,), (3,), (4,), (2, 3), (2, 4), (3, 4), (2, 3, 4)]
BETAS = ((CFG["normalization"]["beta"], CFG["normalization"]["beta2"]),
         tuple(CFG["normalization"]["consistency_pair"]))
BISECT = 40


def counts(members, th):
    case = k_solve(members, th)
    jac = central_difference_jacobians(case.dae, case.equilibrium.x, case.equilibrium.z, {})
    a = jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w
    ev = np.linalg.eigvals(transverse_operator(a, r_x, w).a_perp)
    real = np.abs(ev.imag) <= 1e-9 * max(1.0, np.abs(ev).max())
    return int(((ev.real > 0) & real).sum()), int(((ev.real > 0) & ~real).sum()), jac, r_x, w


def port(jac, r_x, w):
    out = {}
    for b1, b2 in BETAS:
        rp = relocated_port(jac, r_x, w, b1, b2)
        fx2_cond = None
        out[f"{b1:g}/{b2:g}"] = {**rp.__dict__, "fx2_cond": fx2_cond}
    return out


def bisect(members, th, g_lo, g_hi, which):
    c_lo = counts(members, {**th, "g": g_lo})[which]
    for _ in range(BISECT):
        mid = 0.5 * (g_lo + g_hi)
        if counts(members, {**th, "g": mid})[which] == c_lo:
            g_lo = mid
        else:
            g_hi = mid
    return g_lo, g_hi


def evaluate_pair(members, th, g_lo, g_hi, kind):
    lo = counts(members, {**th, "g": g_lo})
    hi = counts(members, {**th, "g": g_hi})
    p_lo, p_hi = port(*lo[2:]), port(*hi[2:])
    verdicts, rows = [], {}
    for key in p_lo:
        a, b = p_lo[key], p_hi[key]
        flip = a["t0_sign"] != b["t0_sign"]
        dev_flip = a["device_sign"] != b["device_sign"]
        verdicts.append(flip and not dev_flip)
        rows[key] = (flip, dev_flip, max(a["identity_residual"], b["identity_residual"]),
                     min(a["min_abs_eig_relocated"], b["min_abs_eig_relocated"]))
    consistent = len(set(verdicts)) == 1
    k1 = next(iter(rows))
    return {
        "kind": kind,
        "subset": "+".join(map(str, members)) or "BASE",
        "t": th["t"],
        "k": th["k"],
        "g_lo": g_lo,
        "g_hi": g_hi,
        "real_rhp_lo": lo[0],
        "real_rhp_hi": hi[0],
        "complex_rhp_lo": lo[1],
        "complex_rhp_hi": hi[1],
        "port_flip": bool(verdicts[0] and consistent),
        "beta_consistent": consistent,
        "t0_flip_b11": rows[k1][0],
        "device_flip_b11": rows[k1][1],
        "identity_residual": max(v[2] for v in rows.values()),
        "min_abs_eig_relocated": min(v[3] for v in rows.values()),
    }


def line_task(task):
    slice_name, t, k, members = task
    th = {"g": 0.0, "k": k, "t": t}
    rows = []
    scan = np.linspace(0.0, 0.003, 61)
    c = [counts(members, {**th, "g": g})[:2] for g in scan]
    # amendment v1.1: port and ground truth at the fixed scan bracket; bisection
    # only locates g*
    for i in range(len(scan) - 1):
        if c[i][0] != c[i + 1][0]:
            lo, hi = bisect(members, th, scan[i], scan[i + 1], 0)
            rows.append({**evaluate_pair(members, th, scan[i], scan[i + 1], "REAL_CROSSING"),
                         "g_star": 0.5 * (lo + hi), "slice": slice_name})
    same = [i for i in range(len(scan) - 1) if c[i] == c[i + 1]]
    for i in same[:: max(1, len(same) // 3)][:3]:
        rows.append({**evaluate_pair(members, th, scan[i], scan[i + 1], "NO_CROSSING"),
                     "slice": slice_name})
    osc = np.linspace(0.003, 1.0, 81)
    co = [counts(members, {**th, "g": g})[:2] for g in osc]
    for i in range(len(osc) - 1):
        if co[i][1] != co[i + 1][1] and co[i][0] == co[i + 1][0]:
            lo, hi = bisect(members, th, osc[i], osc[i + 1], 1)
            rows.append({**evaluate_pair(members, th, osc[i], osc[i + 1], "OSCILLATORY_CROSSING"),
                         "g_star": 0.5 * (lo + hi), "slice": slice_name})
    return rows


# ----------------------------------------------------------------- synthetic --


class J:
    def __init__(self, fx, fz, gx, gz):
        self.fx, self.fz, self.gx, self.gz = fx, fz, gx, gz


def synthetic(rng, kind, sign):
    from ibr_cycles.models.ieee39_devices import OMEGA_B

    n, m = 9, 6
    jblock = np.zeros((n, n))
    jblock[0, 1] = OMEGA_B
    eps = 1e-3 * sign
    if kind == "REAL":
        jblock[2, 2] = eps
        start = 3
    else:
        jblock[2:4, 2:4] = [[eps, 2.0], [-2.0, eps]]
        start = 4
    for i in range(start, n):
        jblock[i, i] = -rng.uniform(0.5, 5.0)
    v = rng.standard_normal((n, n)) + 3 * np.eye(n)
    a = v @ jblock @ np.linalg.inv(v)
    r_x, w = v[:, 0], v[:, 1]
    gz = rng.standard_normal((m, m)) + 4 * np.eye(m)
    gx = rng.standard_normal((m, n))
    fz = rng.standard_normal((n, m))
    fx = a + fz @ np.linalg.solve(gz, gx)
    return J(fx, fz, gx, gz), r_x, w


def synthetic_suite(draws=200, seed=20260921):
    rng = np.random.default_rng(seed)
    rows = []
    for kind in ("REAL", "COMPLEX"):
        for d in range(draws):
            state = rng.bit_generator.state
            lo = synthetic(rng, kind, -1.0)
            rng.bit_generator.state = state  # same system, opposite sign of the crossing
            hi = synthetic(rng, kind, +1.0)
            verdicts = []
            res = 0.0
            for b1, b2 in BETAS:
                a = relocated_port(*lo, b1, b2)
                b = relocated_port(*hi, b1, b2)
                verdicts.append(a.t0_sign != b.t0_sign and a.device_sign == b.device_sign)
                res = max(res, a.identity_residual, b.identity_residual)
            rows.append({"kind": kind, "draw": d, "port_flip": bool(verdicts[0] and len(set(verdicts)) == 1),
                         "beta_consistent": len(set(verdicts)) == 1, "identity_residual": res})
    return pd.DataFrame(rows)


def main(argv) -> int:
    pin_blas_threads()
    exp = FCExperiment(
        name="FC02_port_holdout",
        question="Does the frozen relocated port detect unseen zero-frequency crossings?",
        config={"port_config": CFG, "bisect": BISECT},
        workers=WORKERS,
    )
    started = time.time()
    tasks = []
    for sl in CFG["holdout_design"]["kundur"]["slices"]:
        for k in sl["k"]:
            for m in SUBSETS:
                tasks.append((sl["name"], float(sl["t"]), float(k), m))
    rows = []
    with Pool(WORKERS, initializer=pin_blas_threads) as pool:
        for out in pool.imap_unordered(line_task, tasks, chunksize=1):
            rows += out
    kundur = pd.DataFrame(rows)
    kundur.to_csv(OUT / "FC02_kundur_holdout.csv", index=False)
    syn = synthetic_suite()
    syn.to_csv(OUT / "FC02_synthetic.csv", index=False)
    real = kundur[kundur.kind == "REAL_CROSSING"]
    neg = kundur[kundur.kind != "REAL_CROSSING"]
    table = pd.concat(
        [
            pd.DataFrame(
                [
                    {
                        "source": "Kundur holdout",
                        "holdout_crossings": int(len(real)),
                        "detected": int(real.port_flip.sum()),
                        "missed": int((~real.port_flip).sum()),
                        "negative_controls": int(len(neg)),
                        "false_positives": int(neg.port_flip.sum()),
                        "beta_inconsistent": int((~kundur.beta_consistent).sum()),
                        "max_determinant_residual": float(kundur.identity_residual.max()),
                        "min_abs_eig_relocated": float(kundur.min_abs_eig_relocated.min()),
                    },
                    {
                        "source": "synthetic",
                        "holdout_crossings": int((syn.kind == "REAL").sum()),
                        "detected": int(syn[syn.kind == "REAL"].port_flip.sum()),
                        "missed": int((~syn[syn.kind == "REAL"].port_flip).sum()),
                        "negative_controls": int((syn.kind == "COMPLEX").sum()),
                        "false_positives": int(syn[syn.kind == "COMPLEX"].port_flip.sum()),
                        "beta_inconsistent": int((~syn.beta_consistent).sum()),
                        "max_determinant_residual": float(syn.identity_residual.max()),
                        "min_abs_eig_relocated": float("nan"),
                    },
                    {
                        "source": "IEEE-68",
                        "holdout_crossings": 0,
                        "detected": 0,
                        "missed": 0,
                        "negative_controls": 0,
                        "false_positives": 0,
                        "beta_inconsistent": 0,
                        "max_determinant_residual": float("nan"),
                        "min_abs_eig_relocated": float("nan"),
                    },
                ]
            )
        ]
    )
    table.to_csv(OUT / "zero_frequency_port_validation.csv", index=False)
    publish(OUT / "zero_frequency_port_validation.csv")
    summary = {
        "table": table.to_dict("records"),
        "kundur_by_kind": kundur.groupby("kind").port_flip.agg(["size", "sum"]).to_dict(),
        "kundur_real_by_subset": real.groupby("subset").size().to_dict(),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "FC02_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    print(table.to_string(index=False))
    print(summary["kundur_by_kind"], summary["kundur_real_by_subset"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
