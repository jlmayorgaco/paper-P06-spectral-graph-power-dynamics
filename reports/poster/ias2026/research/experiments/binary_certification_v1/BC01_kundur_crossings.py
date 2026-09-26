"""BC01: the 194 aperiodic Kundur crossings of G1, with no disc around zero.

G1 located these boundaries where a real eigenvalue reached |lambda| = 1e-3
(the edge of the excluded disc), not where it crosses zero. For each one, with
the protocol frozen in configs/binary_certification_v1/BC00_config.yaml:

1. Bracket: evaluate the witness subset at g* (1 -/+ 0.25). Add g = 0 when
   0.75 g* < 0.001, and widen by factors of 2 (at most 6 times) until the real
   eigenvalue of A_qq closest to zero changes sign.
2. Localize the physical zero crossing g0 by 40 bisection steps.
3. Representations on both sides (original port kept for comparison):
   - full descriptor: finite generalized eigenvalues of the DAE pencil, checked
     against A;
   - quotient: A_q and A_qq, with the neutral mode in the ledger;
   - port: r_S(s) = det T_S(s) / det T_0(s), continued analytically through
     s = 0 by the mean value on |s| = rho. This is valid when the windings of
     det T_S and det T_0 on that circle are equal, i.e. the structural zeros
     at the origin cancel.
   The device factor h_S(0) = prod det(-A_i) is kept in the ledger
   (det P_S = h_S det T_S).
   - relocated port (src/ibr_cycles/certification/port_origin.py): the rotation
     and the neutral mode are moved to -beta and -beta2 by exact rank-one
     relocations in the device block. Then det T##(0) is regular, and a
     physical real crossing flips its sign unless the device ledger flips.
     beta, beta2 = (1, 1) and (2, 0.5), and the signs must agree.
4. The original port is EVALUABLE_AT_ZERO when the windings match on both
   sides, the mean-value estimate is stable under rho -> rho/2 (relative
   change < 1e-3), and sign(r_S(0)) changes across the crossing (the
   preregistered rule).
   The result is PORT_EVALUABLE_AT_ZERO when sign det T##(0) flips across the
   crossing, the device ledger does not flip, det(-A##) flips, and the signs
   are the same for both beta pairs. It is REQUIRES_DESCRIPTOR otherwise.
   No outcome count is targeted.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _bc import WORKERS, BCExperiment, out_dir, write_json
from scipy.linalg import eig

from _overnight import pin_blas_threads  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from F12_kundur import theta_of  # noqa: E402
from ibr_cycles.certification.physical import (  # noqa: E402
    classify_physical,
    physical_matrices,
)
from ibr_cycles.certification.port_origin import relocated_port  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    deflate_eigenvector,
    frequency_partner,
    quotient,
    rotation_generator,
)
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.port_admittance import build_port_operator  # noqa: E402

OUT = out_dir("BC01")
RES = OUT.parent.parent
RHO = 1e-5
N_CIRCLE = 64
ITER = 40


def _members(label: str) -> tuple[int, ...]:
    return tuple(sorted(int(b) for b in label.split("+")))


def rest_real_near_zero(case, *, with_jac: bool = False):
    """The real eigenvalue of A_qq closest to zero, and the physical report."""

    a, d, jac = physical_matrices(case)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    partner = frequency_partner(case.dae)
    rep = classify_physical(a, d, case.dae.labels, r_x, partner)
    q = quotient(a, r_x)
    rest = q.a_q
    if rep.neutral_verified:
        v = q.z.T @ partner.w
        rest = deflate_eigenvector(q.a_q, v / np.linalg.norm(v)).a_q
    vals = np.linalg.eigvals(rest)
    real = vals[np.abs(vals.imag) <= 1e-9 * max(1.0, np.abs(vals).max())].real
    lam = float(real[np.argmin(np.abs(real))]) if real.size else float("nan")
    if with_jac:
        return lam, rep, jac, r_x, partner
    return lam, rep


def descriptor_check(case) -> dict:
    jac = central_difference_jacobians(
        case.dae, case.equilibrium.x, case.equilibrium.z, {}
    )
    n_x = jac.fx.shape[0]
    j = np.block([[jac.fx, jac.fz], [jac.gx, jac.gz]])
    e = np.zeros_like(j)
    e[:n_x, :n_x] = np.eye(n_x)
    vals = eig(j, e, right=False)
    finite = vals[np.isfinite(vals) & (np.abs(vals) < 1e8)]
    a = case.system.A
    ref = np.linalg.eigvals(a)
    dist = (
        max(float(np.min(np.abs(ref - v))) for v in finite) if finite.size else np.nan
    )
    return {
        "descriptor_finite": int(finite.size),
        "n_x": int(n_x),
        "descriptor_vs_A_maxdist": dist,
        "descriptor_rhp_raw": int(np.count_nonzero(finite.real > 0)),
    }


def port_ratio(case_s, case_0, rho: float) -> dict:
    t_s = build_port_operator(case_s.dae, case_s.equilibrium.x, case_s.equilibrium.z)
    t_0 = build_port_operator(case_0.dae, case_0.equilibrium.x, case_0.equilibrium.z)
    ang = 2 * np.pi * np.arange(N_CIRCLE) / N_CIRCLE
    pts = rho * np.exp(1j * ang)
    ls = np.array([t_s.log_determinant(s) for s in pts])
    l0 = np.array([t_0.log_determinant(s) for s in pts])

    def winding(logs):
        phase = np.unwrap(np.concatenate([logs.imag, logs.imag[:1]]))
        return int(round((phase[-1] - phase[0]) / (2 * np.pi)))

    ratio = np.exp(ls - l0)
    h_s = np.prod([np.linalg.det(-p.a) for p in t_s.ports if p.n_states]).real
    h_0 = np.prod([np.linalg.det(-p.a) for p in t_0.ports if p.n_states]).real
    return {
        "wind_T_S": winding(ls),
        "wind_T_0": winding(l0),
        "r0": complex(ratio.mean()),
        "r_spread": float(
            np.abs(ratio - ratio.mean()).max() / max(abs(ratio.mean()), 1e-300)
        ),
        "h_ratio_sign": float(np.sign(h_s / h_0)),
    }


def side_eval(members, map_name, y, g) -> dict:
    th = theta_of(map_name, g, y)
    case_s, case_0 = k_solve(members, th), k_solve((), th)
    lam, rep, jac, r_x, partner = rest_real_near_zero(case_s, with_jac=True)
    reloc = {}
    if partner.w is not None:
        for tag, (b1, b2) in (("b11", (1.0, 1.0)), ("b205", (2.0, 0.5))):
            rp = relocated_port(jac, r_x, partner.w, b1, b2)
            reloc.update({f"reloc_{tag}_{k}": v for k, v in rp.__dict__.items()})
    p1 = port_ratio(case_s, case_0, RHO)
    p2 = port_ratio(case_s, case_0, RHO / 2)
    stable_mean = abs(p1["r0"] - p2["r0"]) <= 1e-3 * max(abs(p1["r0"]), 1e-300)
    return {
        "g": g,
        "lam_real_near_zero": lam,
        "status_rest": rep.status_rest,
        "n_positive_rest": rep.n_positive_rest,
        "status_physical": rep.status_physical,
        "neutral_verified": rep.neutral_verified,
        **{f"port_{k}": v for k, v in p1.items()},
        "port_mean_stable": bool(stable_mean),
        **reloc,
        **descriptor_check(case_s),
    }


def crossing(task):
    idx, row = task
    members = _members(row["witness"])
    map_name, y, g_star = row["map"], float(row["y_star"]), float(row["g_star"])

    def lam(g):
        return rest_real_near_zero(k_solve(members, theta_of(map_name, g, y)))[0]

    lo, hi = 0.75 * g_star, 1.25 * g_star
    if lo < 1e-3:
        lo = 0.0
    f_lo, f_hi = lam(lo), lam(hi)
    widen = 0
    while np.sign(f_lo) == np.sign(f_hi) and widen < 6:
        widen += 1
        lo = 0.0 if lo < 1e-3 else max(0.0, lo / 2)
        hi = min(1.0, hi * 2)
        f_lo, f_hi = lam(lo), lam(hi)
    out = {
        **row,
        "crossing_id": idx,
        "bracket_lo": lo,
        "bracket_hi": hi,
        "lam_lo": f_lo,
        "lam_hi": f_hi,
        "widened": widen,
    }
    if np.sign(f_lo) == np.sign(f_hi):
        out["result"] = "NO_PHYSICAL_ZERO_CROSSING_IN_DOMAIN"
        out["side_lo"] = side_eval(members, map_name, y, lo)
        return out
    a, b, fa = lo, hi, f_lo
    for _ in range(ITER):
        mid = 0.5 * (a + b)
        fm = lam(mid)
        if np.sign(fm) == np.sign(fa):
            a, fa = mid, fm
        else:
            b = mid
    g0 = 0.5 * (a + b)
    out["g0_physical"] = g0
    out["g0_minus_gstar"] = g0 - g_star
    side_a = side_eval(members, map_name, y, lo)
    side_b = side_eval(members, map_name, y, hi)
    out["side_lo"], out["side_hi"] = side_a, side_b
    windings_ok = (
        side_a["port_wind_T_S"] == side_a["port_wind_T_0"]
        and side_b["port_wind_T_S"] == side_b["port_wind_T_0"]
    )
    flips = np.sign(side_a["port_r0"].real) != np.sign(side_b["port_r0"].real)
    stable = side_a["port_mean_stable"] and side_b["port_mean_stable"]
    out["port_windings_match"] = bool(windings_ok)
    out["port_sign_flips"] = bool(flips)
    out["port_means_stable"] = bool(stable)
    out["original_port"] = (
        "EVALUABLE_AT_ZERO"
        if windings_ok and flips and stable
        else "NOT_EVALUABLE_AT_ZERO"
    )
    key = "reloc_b11_"
    if f"{key}t0_sign" not in side_a:
        out["result"] = "REQUIRES_DESCRIPTOR"
        out["reason"] = "no verified neutral partner"
        return out
    t_flip = side_a[f"{key}t0_sign"] != side_b[f"{key}t0_sign"]
    h_flip = side_a[f"{key}device_sign"] != side_b[f"{key}device_sign"]
    a_flip = side_a[f"{key}a_sign"] != side_b[f"{key}a_sign"]
    beta_ok = all(
        side[f"reloc_b11_{k}"] == side[f"reloc_b205_{k}"]
        for side in (side_a, side_b)
        for k in ("t0_sign", "device_sign", "a_sign")
    )
    out["reloc_t0_flips"], out["reloc_device_flips"] = bool(t_flip), bool(h_flip)
    out["reloc_a_flips"], out["reloc_beta_independent"] = bool(a_flip), bool(beta_ok)
    if t_flip and not h_flip and a_flip and beta_ok:
        out["result"] = "PORT_EVALUABLE_AT_ZERO"
    elif h_flip:
        out["result"] = "REQUIRES_DESCRIPTOR"
        out["reason"] = "crossing mode is device-internal (h ledger flips)"
    else:
        out["result"] = "REQUIRES_DESCRIPTOR"
        out["reason"] = "relocated port does not resolve the crossing"
    return out


def main(argv) -> int:
    edges = pd.read_csv(RES / "G1" / "G1_rhp_boundaries.csv")
    aper = edges[edges.crossing == "APERIODIC_THROUGH_ORIGIN"].reset_index(drop=True)
    exp = BCExperiment(
        name="BC01_kundur_crossings",
        question="Can the aperiodic Kundur crossings be analysed at s = 0?",
        config={
            "rho": RHO,
            "n_circle": N_CIRCLE,
            "iterations": ITER,
            "crossings": int(len(aper)),
        },
    )
    started = time.time()
    tasks = [(i, r._asdict()) for i, r in enumerate(aper.itertuples(index=False))]
    with Pool(WORKERS, initializer=pin_blas_threads) as pool:
        rows = pool.map(crossing, tasks, chunksize=1)
    flat = []
    for r in rows:
        base = {k: v for k, v in r.items() if not k.startswith("side_")}
        for side in ("side_lo", "side_hi"):
            if side in r:
                base.update({f"{side}:{k}": v for k, v in r[side].items()})
        flat.append(base)
    frame = pd.DataFrame(flat)
    frame.to_csv(OUT / "BC01_kundur_crossings.csv", index=False)
    summary = {
        "crossings": int(len(frame)),
        "results": frame.result.value_counts().to_dict(),
        "g0_below_domain_or_at_zero": int(
            (frame.get("g0_physical", pd.Series(dtype=float)) <= 1e-9).sum()
        ),
        "g0_minus_gstar_range": [
            float(frame.g0_minus_gstar.min()),
            float(frame.g0_minus_gstar.max()),
        ]
        if "g0_minus_gstar" in frame
        else None,
        "descriptor_vs_A_maxdist": float(
            np.nanmax(frame.filter(like="descriptor_vs_A_maxdist").to_numpy())
        ),
        "original_port": frame.get("original_port", pd.Series(dtype=str))
        .value_counts()
        .to_dict(),
        "reasons": frame.get("reason", pd.Series(dtype=str)).value_counts().to_dict(),
    }
    write_json(OUT / "BC01_kundur_crossings_summary.json", summary)
    exp.finish("COMPUTED", elapsed_s=time.time() - started, **summary)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
