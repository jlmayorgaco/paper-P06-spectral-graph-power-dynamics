"""FC13 (step 16, optional diagnostic): first Lyapunov coefficient at three boundaries.

Standard Hopf theory applies to the exact transverse reduced system
`y' = f_perp(y) = Z^T r(x* + Z y)` (TRANSVERSE_STABILITY_QUOTIENT §2), where the
structural zero modes are quotiented out exactly. At each boundary:

- a simple pair `±i omega` of `A_perp` is located by bisection on the policy
  parameter;
- no other eigenvalue lies on the axis (checked);
- transversality is checked (`d alpha / d parameter`).

Then (Kuznetsov, *Elements of Applied Bifurcation Theory*, eq. 3.20):

    l1 = (1 / (2 omega)) Re < p, C(q, q, qbar) − 2 B(q, A^-1 B(q, qbar))
                                 + B(qbar, (2 i omega − A)^-1 B(q, q)) >,
    A q = i omega q,   A^T p = −i omega p,   <p, q> = pbar^T q = 1.

`B` and `C` are the second and third derivatives of `f_perp`. They are
obtained by polarization of central differences of `f_perp`, with `psi` solved by
Newton at every evaluation, and computed at two step sizes to show convergence.
Only established boundaries are used. There is no search for an interesting
sign.

    (1) flagship boundary on the P4 policy line (k = 1.425, t = 1.5), parameter g
    (2) the F7 tongue boundary on k = 1.30 (G2 line), the crossing with
        g in [0.05, 0.15]
    (3) the policy-repaired boundary: P4 flagship plus the damped condenser (D = 2),
        parameter = rating (G1 threshold about 2.48 %)
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
from _fc import FCExperiment, out_dir, write_json

from _f7_common import Theta  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from G1_f8_rhp_followup import _em_config  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402
from ibr_cycles.nonlinear.manifold import chart, psi  # noqa: E402

OUT = out_dir("FC13_hopf_l1")
CORE = (30, 33, 35, 37)
_CFG = {c["name"]: c for c in r_configs()}


def build(kind, p):
    if kind == "cond":
        return solve_config(
            CORE, Theta(0.03625, 1.425, 1.5, 1.0), _em_config(_CFG["R_0010000"], p)
        )
    k = 1.425 if kind == "p4line" else 1.30
    return solve_config(CORE, Theta(p, k, 1.5, 1.0), _CFG["R_none"])


class _Perp:
    def __init__(self, a_perp, z):
        self.a_perp, self.z = a_perp, z


def perp(case):
    """Transverse operator; dead states (identically zero rows: service switches at 0)
    are exact zeros, held fixed and removed exactly, as in the frozen classifier."""

    from scipy.linalg import null_space

    from ibr_cycles.certification.symmetry import orthonormal

    a = case.system.A
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    w = frequency_partner(case.dae).w  # None when a machine has D != 0 (condenser)
    dead = np.flatnonzero(~np.any(a != 0.0, axis=1))
    if dead.size == 0:
        tr = transverse_operator(a, r_x, w)
        return _Perp(tr.a_perp, tr.z)
    cols = [r_x] + ([] if w is None else [w]) + [np.eye(a.shape[0])[:, k] for k in dead]
    u = orthonormal(np.column_stack(cols))
    z = null_space(u.T)
    return _Perp(z.T @ a @ z, z)


def osc_alpha(case):
    ev = np.linalg.eigvals(perp(case).a_perp)
    osc = ev[np.abs(ev.imag) > 0.5]
    return float(osc.real.max()), ev


def locate(kind, lo, hi):
    a_lo = osc_alpha(build(kind, lo))[0]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        a = osc_alpha(build(kind, mid))[0]
        if np.sign(a) == np.sign(a_lo):
            lo, a_lo = mid, a
        else:
            hi = mid
    return 0.5 * (lo + hi)


def l1_at(case, hsteps=(2e-3, 1e-3)):
    dae = case.dae
    x0, z0 = case.equilibrium.x, case.equilibrium.z
    tr = perp(case)
    zb, a = tr.z, tr.a_perp
    ch = chart(dae, x0, z0)

    def f(y):
        x = x0 + zb @ y
        return zb.T @ dae.f(x, psi(dae, x, z0, ch.lu), {})

    ev, vr = np.linalg.eig(a)
    k = int(np.argmin(np.abs(ev.real) + 1e3 * (ev.imag <= 0)))
    om = float(ev[k].imag)
    q = vr[:, k] / np.linalg.norm(vr[:, k])
    evl, vl = np.linalg.eig(a.T)
    kl = int(np.argmin(np.abs(evl - np.conj(ev[k]))))
    p = vl[:, kl]
    p = p / np.conj(np.vdot(p, q))  # pbar^T q = 1
    others = np.delete(ev, [k, int(np.argmin(np.abs(ev - np.conj(ev[k]))))])
    gap = float(np.abs(others.real).min())
    out = {
        "omega": om,
        "freq_hz": om / (2 * np.pi),
        "alpha_at_point": float(ev[k].real),
        "min_abs_re_other": gap,
    }
    n = a.shape[0]
    for h in hsteps:

        def d2(u, h=h):
            return (f(h * u) - 2 * f(np.zeros(n)) + f(-h * u)) / h**2

        def d3(u, h=h):
            return (f(2 * h * u) - 2 * f(h * u) + 2 * f(-h * u) - f(-2 * h * u)) / (
                2 * h**3
            )

        def B(u, v):  # real bilinear by polarization
            return (d2(u + v) - d2(u - v)) / 4.0

        def Bc(u, v):  # complex arguments
            ur, ui, vr_, vi = u.real, u.imag, v.real, v.imag
            return (B(ur, vr_) - B(ui, vi)) + 1j * (B(ur, vi) + B(ui, vr_))

        def C3(u, v, w):  # real symmetric trilinear by polarization of the cubic form
            return (
                d3(u + v + w) - d3(u + v - w) - d3(u - v + w) + d3(u - v - w)
            ) / 24.0

        ar, ai = q.real, q.imag
        cqqq = (C3(ar, ar, ar) + C3(ar, ai, ai)) + 1j * (
            C3(ar, ar, ai) + C3(ai, ai, ai)
        )
        bqqb = Bc(q, np.conj(q))
        bqq = Bc(q, q)
        t1 = np.linalg.solve(a, bqqb)
        t2 = np.linalg.solve(2j * om * np.eye(n) - a, bqq)
        val = np.vdot(p, cqqq - 2 * Bc(q, t1) + Bc(np.conj(q), t2))
        out[f"l1_h{h:g}"] = float(np.real(val) / (2 * om))
    return out


def main(argv) -> int:
    exp = FCExperiment(
        name="FC13_hopf_l1",
        question="Are the established oscillatory boundaries super- or subcritical?",
    )
    started = time.time()
    rows = []
    for kind, lo, hi, label in (
        ("p4line", 0.03625, 0.30, "flagship boundary, P4 line k = 1.425"),
        ("tongue", 0.05, 0.15, "tongue boundary, k = 1.30"),
        ("cond", 0.02, 0.03, "policy-repaired boundary: damped condenser rating"),
    ):
        a_lo, a_hi = osc_alpha(build(kind, lo))[0], osc_alpha(build(kind, hi))[0]
        if np.sign(a_lo) == np.sign(a_hi):
            rows.append(
                {
                    "boundary": label,
                    "status": f"NO CROSSING in [{lo}, {hi}]",
                    "alpha_lo": a_lo,
                    "alpha_hi": a_hi,
                }
            )
            continue
        p_star = locate(kind, lo, hi)
        dp = 1e-4 * max(abs(p_star), 1e-3)
        dalpha = (
            osc_alpha(build(kind, p_star + dp))[0]
            - osc_alpha(build(kind, p_star - dp))[0]
        ) / (2 * dp)
        res = l1_at(build(kind, p_star))
        l1s = [v for k, v in res.items() if k.startswith("l1_h")]
        rows.append(
            {
                "boundary": label,
                "parameter": p_star,
                "d_alpha_d_parameter": dalpha,
                **res,
                "l1_step_spread": float(
                    abs(l1s[0] - l1s[1]) / max(abs(l1s[1]), 1e-300)
                ),
                "type": "supercritical" if np.mean(l1s) < 0 else "subcritical",
                "status": "COMPUTED",
            }
        )
        print(rows[-1], flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "FC13_hopf.csv", index=False)
    write_json(
        OUT / "FC13_summary.json",
        {"rows": rows, "elapsed_s": round(time.time() - started, 1)},
    )
    exp.finish("COMPUTED", rows=rows)
    print(frame.T.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
