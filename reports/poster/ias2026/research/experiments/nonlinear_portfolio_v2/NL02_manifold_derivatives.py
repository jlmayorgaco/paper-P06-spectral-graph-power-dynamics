"""NL02: derivatives of the algebraic manifold, verified against direct derivation.

For 11 declared points (equilibria and off-equilibrium states on the manifold) of
IEEE-39, Kundur and IEEE-68, the manifold formulas of
src/ibr_cycles/nonlinear/manifold.py are compared with direct differentiation of
the reduced field r(x) = f(x, psi(x)), where psi is re-solved by Newton at every
evaluation:

  - J_psi a against [psi(x + h a) - psi(x - h a)] / 2h;
  - H_r[a, b] against the second difference of r, for 9 step sizes
    (1e-1 ... 1e-5). The whole convergence curve is stored, including where
    rounding takes over;
  - B, C, D against direct derivatives of r and of the outputs with respect to
    the declared inputs (first load P and Q, first machine Pm) at fixed x.

Directions: seeded unit vectors (numpy default_rng 20260912), plus a = b.
Nothing is tuned: the formula uses h = 1e-4 throughout, and its sensitivity
to 1e-3 and 1e-5 is reported.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import _bootstrap  # noqa: E402,F401
from _bootstrap import ROOT  # noqa: E402
from _f7_common import Theta  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from G1_f8_rhp_followup import _em_config  # noqa: E402
from G3_ieee68 import solve as s68  # noqa: E402
from ibr_cycles.nonlinear.manifold import (  # noqa: E402
    chart,
    hessian_direct,
    hessian_reduced,
    input_output,
    psi,
)
from ibr_cycles.nonlinear.model import PhasorModel  # noqa: E402

OUT = ROOT / "results" / "NL" / "NL02"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260912
STEPS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
H_FORMULA = 1e-4


def points():
    cfg = {c["name"]: c for c in r_configs()}
    p4 = Theta(0.03625, 1.425, 1.5)
    flag = (30, 33, 35, 37)
    return [
        ("IEEE-39 base P4", lambda: solve_config((), p4, cfg["R_none"]), 0.0),
        ("IEEE-39 base P4 off-eq", lambda: solve_config((), p4, cfg["R_none"]), 1e-3),
        ("IEEE-39 flagship P4", lambda: solve_config(flag, p4, cfg["R_none"]), 0.0),
        (
            "IEEE-39 flagship P4 off-eq",
            lambda: solve_config(flag, p4, cfg["R_none"]),
            1e-3,
        ),
        (
            "IEEE-39 tongue gap g=0.042",
            lambda: solve_config(flag, Theta(0.042, 1.30, 1.5), cfg["R_none"]),
            0.0,
        ),
        (
            "IEEE-39 condenser D=2 3%",
            lambda: solve_config(flag, p4, _em_config(cfg["R_0010000"], 0.03)),
            0.0,
        ),
        (
            "Kundur 2+3 g=0 k=0.75",
            lambda: k_solve((2, 3), {"g": 0.0, "k": 0.75, "t": 1.0}),
            0.0,
        ),
        (
            "Kundur 2 g=0.08 k=1.25",
            lambda: k_solve((2,), {"g": 0.08, "k": 1.25, "t": 1.0}),
            0.0,
        ),
        (
            "Kundur 2 g=0.08 off-eq",
            lambda: k_solve((2,), {"g": 0.08, "k": 1.25, "t": 1.0}),
            1e-3,
        ),
        ("IEEE-68 base", lambda: s68((), 0.0, 1.0), 0.0),
        ("IEEE-68 3+4+6+9 off-eq", lambda: s68((3, 4, 6, 9), 0.0, 1.0), 1e-3),
    ]


def on_manifold(case, eps, rng):
    dae = case.dae
    x0, z0 = case.equilibrium.x, case.equilibrium.z
    if eps == 0.0:
        return x0.copy(), z0.copy()
    x = x0 + eps * rng.standard_normal(x0.size)
    ch0 = chart(dae, x0, z0)
    z = psi(dae, x, z0, ch0.lu)
    return x, z


def main() -> int:
    rng = np.random.default_rng(SEED)
    conv_rows, summary = [], []
    started = time.time()
    for name, build, eps in points():
        t0 = time.time()
        case = build()
        dae = case.dae
        x, z = on_manifold(case, eps, rng)
        ch = chart(dae, x, z)
        n = x.size
        a = rng.standard_normal(n)
        a /= np.linalg.norm(a)
        b = rng.standard_normal(n)
        b /= np.linalg.norm(b)
        # J_psi
        hj = 1e-6
        jpsi_fd = (psi(dae, x + hj * a, z, ch.lu) - psi(dae, x - hj * a, z, ch.lu)) / (
            2 * hj
        )
        jpsi_err = float(
            np.linalg.norm(ch.j_psi @ a - jpsi_fd) / np.linalg.norm(jpsi_fd)
        )
        # A against a direct first difference of r
        r_fd = (
            dae.f(x + hj * a, psi(dae, x + hj * a, z, ch.lu), {})
            - dae.f(x - hj * a, psi(dae, x - hj * a, z, ch.lu), {})
        ) / (2 * hj)
        a_err = float(np.linalg.norm(ch.a @ a - r_fd) / np.linalg.norm(r_fd))
        for tag, (p, q) in (("a,b", (a, b)), ("a,a", (a, a))):
            _, h_formula = hessian_reduced(dae, ch, p, q, H_FORMULA)
            sens = {
                hh: float(
                    np.linalg.norm(hessian_reduced(dae, ch, p, q, hh)[1] - h_formula)
                    / np.linalg.norm(h_formula)
                )
                for hh in (1e-3, 1e-5)
            }
            best = np.inf
            for h in STEPS:
                try:
                    direct = hessian_direct(dae, ch, p, q, h)
                    err = float(
                        np.linalg.norm(direct - h_formula) / np.linalg.norm(h_formula)
                    )
                    status = "OK"
                except RuntimeError as error:
                    err, status = float("nan"), f"CHART_FAIL: {str(error)[:40]}"
                conv_rows.append(
                    {
                        "point": name,
                        "pair": tag,
                        "h": h,
                        "rel_error": err,
                        "status": status,
                        "norm_H": float(np.linalg.norm(h_formula)),
                    }
                )
                if err == err:
                    best = min(best, err)
            summary.append(
                {
                    "point": name,
                    "pair": tag,
                    "n_x": n,
                    "n_z": z.size,
                    "gz_cond": float(np.linalg.cond(ch.gz)),
                    "jpsi_rel_error": jpsi_err,
                    "A_rel_error": a_err,
                    "best_H_rel_error": best,
                    "formula_step_sensitivity_1e-3": sens[1e-3],
                    "formula_step_sensitivity_1e-5": sens[1e-5],
                }
            )
        # B, C, D on three declared inputs at fixed x
        if eps == 0.0:
            model = PhasorModel(case)
            io = input_output(model, ch)
            nl = len(model.load_buses)
            pick = [0, nl, 2 * nl] if model.n_u > 2 * nl else [0, nl]
            du = 1e-6
            for j in pick:
                e = np.zeros(model.n_u)
                e[j] = du
                dp, dm = model._with_inputs(e), model._with_inputs(-e)
                zp, zm = psi(dp, x, z, ch.lu), psi(dm, x, z, ch.lu)
                b_fd = (dp.f(x, zp, {}) - dm.f(x, zm, {})) / (2 * du)
                d_fd = (
                    model.outputs(x, zp, e).values - model.outputs(x, zm, -e).values
                ) / (2 * du)
                summary.append(
                    {
                        "point": name,
                        "pair": f"input {model.input_names[j]}",
                        "B_rel_error": float(
                            np.linalg.norm(io["B"][:, j] - b_fd)
                            / max(np.linalg.norm(b_fd), 1e-300)
                        ),
                        "D_rel_error": float(
                            np.linalg.norm(io["D"][:, j] - d_fd)
                            / max(np.linalg.norm(d_fd), 1e-300)
                        ),
                        "D_norm": float(np.linalg.norm(io["D"][:, j])),
                    }
                )
            ca = io["C"] @ a
            ho = model.outputs(x + hj * a, psi(dae, x + hj * a, z, ch.lu)).values
            hm = model.outputs(x - hj * a, psi(dae, x - hj * a, z, ch.lu)).values
            c_fd = (ho - hm) / (2 * hj)
            summary.append(
                {
                    "point": name,
                    "pair": "C a",
                    "C_rel_error": float(
                        np.linalg.norm(ca - c_fd) / np.linalg.norm(c_fd)
                    ),
                }
            )
        print(f"{name}: {time.time() - t0:.1f}s", flush=True)
    conv = pd.DataFrame(conv_rows)
    conv.to_csv(OUT / "NL02_convergence.csv", index=False)
    summ = pd.DataFrame(summary)
    summ.to_csv(OUT / "NL02_summary.csv", index=False)
    out = {
        "points": int(conv.point.nunique()),
        "hessian_best_rel_error_max": float(summ.best_H_rel_error.max()),
        "hessian_best_rel_error_median": float(summ.best_H_rel_error.median()),
        "jpsi_rel_error_max": float(summ.jpsi_rel_error.max()),
        "A_rel_error_max": float(summ.A_rel_error.max()),
        "B_rel_error_max": float(summ.get("B_rel_error", pd.Series(dtype=float)).max()),
        "C_rel_error_max": float(summ.get("C_rel_error", pd.Series(dtype=float)).max()),
        "D_rel_error_max": float(summ.get("D_rel_error", pd.Series(dtype=float)).max()),
        "chart_failures": int((conv.status != "OK").sum()),
        # the direct difference at h = H_FORMULA shares the formula's own step,
        # so the independent agreement excludes it
        "hessian_best_independent_by_point": conv[conv.h != H_FORMULA]
        .groupby("point")
        .rel_error.min()
        .round(8)
        .to_dict(),
        "elapsed_s": round(time.time() - started, 1),
    }
    (OUT / "NL02_summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    pd.set_option("display.width", 220)
    print(summ.to_string())
    print(
        conv.pivot_table(
            index=["point", "pair"], columns="h", values="rel_error"
        ).to_string()
    )
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
