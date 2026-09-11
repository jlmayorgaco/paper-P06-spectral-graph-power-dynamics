"""FC09 (BC03 gate; campaign steps 12-13): binary / cardinality certificates, IEEE-39.

Object (theory/PRINCIPAL_MINOR_PORTFOLIO_STRUCTURE.md): at a common operating
point (C2 common realization, BC02),

    det T_S / det T_0 = det(I + M_SS),  M = blkdiag(dY_i) E^T T_0^{-1} E,
    Q = blkdiag(I + M_ii)^{-1} (I + M) - I      (zero block diagonal).

Assumptions checked or reported:

- (A1) exact locality at the common point;
- (A2) N_h = 0 (BC02 ledger);
- (A3) base and every singleton transversely stable;
- (A4) Q bounded on the contour (the structural pole of T_0^{-1} at s = 0);
- (A5) continuous-frequency bounds. Only a grid is available here, so every
  verdict is SCREENING, never CERTIFIED.

Candidates, evaluated frequency-wise on s = j omega (log grid 1e-3 .. 1e3 rad/s,
plus the band 0.3-1.5 Hz refined):

    SG    small gain / Perron (H-matrix): rho(R(omega)) < 1 for all omega,
          R_ij = ||Q_ij||, R_ii = 0   -> all subsets safe
    TS_r  top-sum comparison: max_i TopSum_{r-1}{R_ij} < 1 -> kappa >= r + 1
    TSW_r the same with Perron weights d(omega)
    GB    block Gershgorin row sums (the unweighted TS at r = m)

The exact kappa comes from the transverse spectra of the 16 subsets (direct
path). The CPU time of each certificate includes the construction (realization,
ports, grid). Also recorded: the exact Boolean degree of A_red(delta), J(delta)
and r(S) at each point (step 13).
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd
from _fc import CORE, RESULTS, FCExperiment, out_dir, write_json
from BC02_binary_representation import _spec, port_t  # noqa: E402
from BC02b_boolean_degree import by_order, mobius  # noqa: E402

from _f7_common import SUBSETS, Theta  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from ibr_cycles.certification.binary import (  # noqa: E402
    all_vertices,
    build_common_realization,
    reduced,
    stacked,
)
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402

OUT = out_dir("FC09_bc03_certificates")
POINTS = ("P4", "P2", "P_inf")


def grid():
    g = np.logspace(-3, 3, 700)
    band = 2 * np.pi * np.linspace(0.3, 1.5, 300)
    return np.unique(np.concatenate([g, band]))


def exact_kappa(point):
    t = time.time()
    cfg = {c["name"]: c for c in r_configs()}["R_none"]
    th = Theta(point["g"], point["k"], point["t"], point["h"])
    st = {}
    for m in SUBSETS:
        case = solve_config(m, th, cfg)
        r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
        ev = np.linalg.eigvals(
            transverse_operator(
                case.system.A, r_x, frequency_partner(case.dae).w
            ).a_perp
        )
        st[m] = bool((ev.real > 0).any())
    uns = [s for s, u in st.items() if u]
    mins = [s for s in uns if not any(set(q) < set(s) for q in uns)]
    return (min((len(s) for s in mins), default=-1), mins, st, time.time() - t)


def topsum(row, r):
    vals = np.sort(row)[::-1]
    return float(vals[: max(r - 1, 0)].sum())


def certificates(point):
    started = time.time()
    builder, kwargs = _spec(("39", point))
    cr = build_common_realization(CORE, **kwargs)
    jac = {}
    for delta, s in all_vertices(CORE):
        j, _, _ = cr.jacobian(delta, "C2")
        jac[tuple(sorted(s))] = j
    build_s = time.time() - started
    net = cr.case.dae.network
    ch = {b: [2 * net.position(b), 2 * net.position(b) + 1] for b in CORE}
    idx = [c for b in CORE for c in ch[b]]
    m = len(CORE)
    om = grid()
    rho, ts, tsw, rmax, locality, qnorm_low = (
        [],
        {r: [] for r in range(2, m + 1)},
        {r: [] for r in range(2, m + 1)},
        np.zeros((m, m)),
        0.0,
        [],
    )
    for w in om:
        s = 1j * w
        t0 = port_t(jac[()], s)
        kmat = np.linalg.solve(t0, np.eye(t0.shape[0]))[np.ix_(idx, idx)]
        dmat = np.zeros((2 * m, 2 * m), complex)
        for i, b in enumerate(CORE):
            dmat[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = (port_t(jac[(b,)], s) - t0)[
                np.ix_(ch[b], ch[b])
            ]
        mm = dmat @ kmat
        dblk = np.zeros_like(mm)
        for i in range(m):
            dblk[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = (
                np.eye(2) + mm[2 * i : 2 * i + 2, 2 * i : 2 * i + 2]
            )
        q = np.linalg.solve(dblk, np.eye(2 * m) + mm) - np.eye(2 * m)
        rmat = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                if i != j:
                    rmat[i, j] = np.linalg.norm(
                        q[2 * i : 2 * i + 2, 2 * j : 2 * j + 2], 2
                    )
        rmax = np.maximum(rmax, rmat)
        ev, vec = np.linalg.eig(rmat)
        k = int(np.argmax(ev.real))
        rho.append(float(ev.real[k]))
        d = np.abs(vec[:, k].real) + 1e-12
        for r in range(2, m + 1):
            ts[r].append(max(topsum(rmat[i], r) for i in range(m)))
            scaled = rmat * d[None, :] / d[:, None]
            tsw[r].append(max(topsum(scaled[i], r) for i in range(m)))
        if w < 1e-2:
            qnorm_low.append(float(np.abs(q).max()))
        # locality check at this s for the full core
        full = tuple(CORE)
        tf = port_t(jac[full], s)
        pred = t0.copy()
        for b in CORE:
            pred += port_t(jac[(b,)], s) - t0
        locality = max(locality, float(np.linalg.norm(tf - pred) / np.linalg.norm(tf)))
    cpu = time.time() - started
    rho = np.array(rho)
    out = {
        "build_realization_s": round(build_s, 1),
        "certificate_cpu_s_incl_construction": round(cpu, 1),
        "grid_points": int(om.size),
        "max_rho_R_omega": float(rho.max()),
        "omega_at_max_rho": float(om[int(np.argmax(rho))]),
        "SG_all_subsets_screened": bool(rho.max() < 1),
        "locality_rel_max": locality,
        "Q_abs_max_below_1e-2_rad_s": float(max(qnorm_low)) if qnorm_low else None,
        "Q_blows_up_at_zero": bool(qnorm_low and qnorm_low[0] > 10 * qnorm_low[-1]),
        "sup_rho_of_Rbar": float(np.max(np.abs(np.linalg.eigvals(rmax)))),
    }
    cert_r = 0
    cert_rw = 0
    for r in range(2, m + 1):
        out[f"TS_{r}_max"] = float(np.max(ts[r]))
        out[f"TSW_{r}_max"] = float(np.max(tsw[r]))
        if np.max(ts[r]) < 1:
            cert_r = r
        if np.max(tsw[r]) < 1:
            cert_rw = r
    out["screened_max_subset_size_TS"] = cert_r if cert_r else 1
    out["screened_max_subset_size_TSW"] = cert_rw if cert_rw else 1
    out["kappa_lower_bound_screened"] = (
        max(cert_r, cert_rw, 1) + 1 if max(cert_r, cert_rw) else 2
    )
    # Boolean degree at this point
    j_vals = {t: stacked(jac[t]) for t in jac}
    a_vals = {t: reduced(jac[t]) for t in jac}
    jd = by_order(mobius(j_vals, CORE), np.linalg.norm(j_vals[()]))
    ad = by_order(mobius(a_vals, CORE), np.linalg.norm(a_vals[()]))
    out["J_moebius_rel_by_order"] = jd
    out["A_moebius_rel_by_order"] = ad
    out["J_exact_degree"] = max(k for k, v in jd.items() if v > 1e-10)
    out["A_exact_degree"] = max(k for k, v in ad.items() if v > 1e-10)
    return out, om, rho, ts, tsw


def main(argv) -> int:
    exp = FCExperiment(
        name="FC09_bc03_certificates",
        question=(
            "Does any cheap structured certificate certify all-subset safety or a "
            "kappa bound?"
        ),
        config={"points": POINTS, "grid": "logspace(-3,3,700) + band 300"},
    )
    started = time.time()
    f8 = json.loads((RESULTS / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    rows, curves = [], []
    for name in POINTS:
        p = f8[name]
        kap, mins, st, exact_s = exact_kappa(p)
        cert, om, rho, ts, tsw = certificates(p)
        singles_stable = all(not st[(b,)] for b in CORE) and not st[()]
        row = {
            "point": name,
            "exact_kappa": kap,
            "exact_H": "|".join("+".join(map(str, s)) for s in mins) or "EMPTY",
            "exact_cpu_s": round(exact_s, 1),
            "A3_base_and_singletons_stable": singles_stable,
            **{k: v for k, v in cert.items() if not isinstance(v, dict)},
        }
        row["gap"] = (
            (kap - row["kappa_lower_bound_screened"])
            if kap > 0
            else "n/a (kappa = inf)"
        )
        row["status"] = "SCREENING (sampled frequency grid; not a certificate)"
        rows.append(row)
        for w, r in zip(om, rho, strict=True):
            curves.append(
                {
                    "point": name,
                    "omega": w,
                    "rho_R": r,
                    **{
                        f"TS_{k}": ts[k][i]
                        for k in ts
                        for i in [int(np.searchsorted(om, w))]
                    },
                }
            )
        write_json(
            OUT / f"FC09_{name}_degree.json",
            {"J": cert["J_moebius_rel_by_order"], "A": cert["A_moebius_rel_by_order"]},
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "FC09_certificates.csv", index=False)
    pd.DataFrame(curves).to_csv(OUT / "FC09_rho_curves.csv.gz", index=False)
    summary = {"rows": rows, "elapsed_s": round(time.time() - started, 1)}
    write_json(OUT / "FC09_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    pd.set_option("display.width", 250)
    print(frame.T.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
