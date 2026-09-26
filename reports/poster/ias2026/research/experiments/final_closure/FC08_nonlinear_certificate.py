"""FC08 (step 11): nonlinear local certificate vs the TDS threshold.

Energy / passivity certificate: NOT APPLICABLE to the implemented model. Every
known transient-energy function needs a lossless reduced network with
constant EMFs and no controllers (classical model). The implemented L0 has:

- network losses and constant-power loads;
- two-axis flux dynamics;
- first-order AVRs with a power-input PSS;
- GFL PLL / PI / current loops.

No physically justified energy or storage function is available for that
combination, and none is invented here (FINAL_NONLINEAR_MODEL_TRACEABILITY §G).

Generic Lyapunov + curvature certificate (v1 T8), on the exact transverse
system `y' = Z^T r(x* + Z y)` (TRANSVERSE_STABILITY_QUOTIENT §2), in the model's
own per-unit coordinates, with no coordinate optimization:

    A_perp^T P + P A_perp = −I,   beta = 1 / (2 lambda_max(P))
    ||rem(y)|| <= L ||y||^2 on ||y|| <= r0,  rem = f_perp(y) − A_perp y
    r_cert = min(r0, beta lambda_min(P) / (||P|| L)),   c = lambda_min(P) r_cert^2
    amplitude a_cert: V(a v) = c, with v the transverse state at the pulse end
    of the LINEARIZED response to a unit pulse

`L` is estimated from samples of `rem` (200 random directions, 4 radii). A
sampled `L` is not a bound, so the result is labelled ESTIMATE, not CERTIFIED.
The forced phase is also linearized. The comparison with the FC05 TDS threshold
gives the conservatism ratio `a_cert / r_TDS`.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import yaml
from _fc import CONFIGS, FCExperiment, out_dir, write_json
from scipy.linalg import expm, solve_continuous_lyapunov

from _f7_common import Theta  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.nonlinear.manifold import chart, psi  # noqa: E402
from ibr_cycles.nonlinear.model import PhasorModel  # noqa: E402

OUT = out_dir("FC08_nonlinear_certificate")
CFG = yaml.safe_load(
    (CONFIGS / "final_nonlinear_composability_v1.yaml").read_text(encoding="utf-8")
)
CASES = [
    ("P4", ()),
    ("P4", (30,)),
    ("P4", (30, 33)),
    ("P4", (30, 33, 35)),
    ("P_inf", ()),
    ("P_inf", (30, 33, 35, 37)),
]
SEED = 20260922


def main(argv) -> int:
    exp = FCExperiment(
        name="FC08_nonlinear_certificate",
        question=(
            "What fraction of the TDS threshold does a local Lyapunov "
            "certificate reach?"
        ),
        config={"cases": [f"{p}:{m}" for p, m in CASES], "seed": SEED},
    )
    rng = np.random.default_rng(SEED)
    cfg = {c["name"]: c for c in r_configs()}["R_none"]
    thr_path = out_dir("FC05_nonlinear_thresholds") / "FC05_thresholds.csv"
    thr = pd.read_csv(thr_path) if thr_path.exists() else None
    rows = []
    for point, members in CASES:
        t0 = time.time()
        p = CFG["policy_points"][point]
        case = solve_config(members, Theta(p["g"], p["k"], p["t"], p["h"]), cfg)
        dae = case.dae
        x0, z0 = case.equilibrium.x, case.equilibrium.z
        r_x, _ = rotation_generator(dae, z0)
        w = frequency_partner(dae).w
        tr = transverse_operator(case.system.A, r_x, w)
        z_basis, ap = tr.z, tr.a_perp
        if np.linalg.eigvals(ap).real.max() >= 0:
            rows.append(
                {
                    "point": point,
                    "subset": "+".join(map(str, members)) or "BASE",
                    "status": "NOT_APPLICABLE (transversely unstable)",
                }
            )
            continue
        pmat = solve_continuous_lyapunov(ap.T, -np.eye(ap.shape[0]))
        ev = np.linalg.eigvalsh(pmat)
        lmin, lmax = float(ev.min()), float(ev.max())
        beta = 1.0 / (2.0 * lmax)
        ch = chart(dae, x0, z0)

        def fperp(y, dae=dae, x0=x0, z0=z0, z_basis=z_basis, ch=ch):
            x = x0 + z_basis @ y
            z = psi(dae, x, z0, ch.lu)
            return z_basis.T @ dae.f(x, z, {})

        f0 = fperp(np.zeros(ap.shape[0]))
        ratios = []
        for rad in (1e-3, 3e-3, 1e-2, 3e-2):
            for _ in range(50):
                e = rng.standard_normal(ap.shape[0])
                e *= rad / np.linalg.norm(e)
                try:
                    rem = fperp(e) - f0 - ap @ e
                except RuntimeError:
                    continue
                ratios.append(np.linalg.norm(rem) / rad**2)
        l_hat = float(np.max(ratios))
        r0 = 3e-2
        r_cert = min(r0, beta * lmin / (lmax * l_hat))
        c = lmin * r_cert**2
        # linearized transverse response to a unit D1 pulse (1 MW) at the pulse end
        model = PhasorModel(case)
        k = model.input_names.index(f"dP_load_{CFG['disturbances']['D1']['bus']}")
        jac = central_difference_jacobians(dae, x0, z0, {})
        u = np.zeros(model.n_u)
        u[k] = 1e-6
        gu = (
            model.residual_g(x0, z0, u).values - model.residual_g(x0, z0, -u).values
        ) / 2e-6
        b_full = -jac.fz @ np.linalg.solve(jac.gz, gu)
        b_perp = z_basis.T @ b_full * 0.01  # per MW
        tp = CFG["disturbances"]["D1"]["duration_s"]
        n = ap.shape[0]
        aug = np.zeros((n + 1, n + 1))
        aug[:n, :n] = ap * tp
        aug[:n, n] = b_perp * tp
        v = expm(aug)[:n, n]
        a_cert = float(np.sqrt(c / (v @ pmat @ v)))
        r_tds, r_type = np.nan, ""
        if thr is not None:
            sel = thr[
                (thr.point == point)
                & (thr.family == "D1")
                & (thr.subset == ("+".join(map(str, members)) or "BASE"))
            ]
            if len(sel):
                r_tds, r_type = float(sel.r.iloc[0]), str(sel.type.iloc[0])
        rows.append(
            {
                "point": point,
                "subset": "+".join(map(str, members)) or "BASE",
                "status": "ESTIMATE (sampled curvature, linearized forced phase)",
                "beta": beta,
                "lambda_min_P": lmin,
                "lambda_max_P": lmax,
                "L_hat": l_hat,
                "r_cert_state": r_cert,
                "a_cert_MW": a_cert,
                "r_TDS_MW": r_tds,
                "r_TDS_type": r_type,
                "conservatism_ratio": a_cert / r_tds
                if np.isfinite(r_tds) and r_tds > 0
                else np.nan,
                "cpu_s": round(time.time() - t0, 1),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "FC08_certificates.csv", index=False)
    summary = {
        "rows": rows,
        "energy_certificate": "NOT_APPLICABLE (see module docstring)",
    }
    write_json(OUT / "FC08_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    pd.set_option("display.width", 250)
    print(frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
