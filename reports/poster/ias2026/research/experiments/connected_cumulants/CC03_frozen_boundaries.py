# ruff: noqa: E501  -- table labels kept on one line
"""CC03 - Phase 3: connected cumulants at the ten frozen FC18 boundaries.

Preregistered in configs/ias2026/connected_cumulants_prereg_v1.yaml (committed de09db3b
before this script was run). theta and s* are read verbatim from FC18_summary.json;
nothing is re-located or tuned.

Per boundary (H = the coalition whose label changes):
    F(B), mu_B, chi_B for every B subseteq H; |F(H)(s*)|; reconstruction residual
    nu_H and the preregistered label (primary 0.10; sensitivity 0.05 / 0.20)
    delta_B = |chi_B F(H minus B)| / |dF(H)/ds| for every B subseteq H, |B| >= 2
    Theorem 2.3 check (cumulant support connected and spanning H)
    naive single-cycle part and remainder R_H = chi_H - naive_H
    deletion test: zero of C_H = F(H) - chi_H nearest s* (Newton)
At FLAG: cross-check of F(B) against the PortActionSpace Q of the frozen F2c script.
"""

from __future__ import annotations

import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _cc import CORE, Realization, frozen_events, kwargs_of, label, out_dir, write_json

from _overnight import pin_blas_threads
from ibr_cycles.cycles.connected import (
    block_cycle_trace_sum,
    boolean_mobius,
    characteristic_values,
    composite_part,
    connected_share,
    cumulant_support_connected,
    cumulants,
    moments,
    subsets,
)

OUT = out_dir("CC03")
TAU, TAU_LO, TAU_HI = 0.10, 0.05, 0.20
H_S = 1e-5
NEWTON_RADIUS, MATERIAL_DISTANCE = 0.5, 0.05


def f_of(r: Realization, s: complex, h: tuple) -> dict:
    """F(B) for B subseteq h at s, keyed by frozensets of bus numbers."""

    q = r.q_matrix(s)[0]
    idx = [CORE.index(b) for b in h]
    sub = q[
        np.ix_(
            [2 * i + c for i in idx for c in (0, 1)],
            [2 * i + c for i in idx for c in (0, 1)],
        )
    ]
    vals = characteristic_values(sub, [2] * len(h))
    return {frozenset(h[i] for i in k): v for k, v in vals.items()}, sub


def chi_top(r, s, h):
    f, _ = f_of(r, s, h)
    chi = cumulants(f, h)
    return f, chi


def newton_composite(r, s0, h):
    """Zero of C_H(s) = F(H)(s) - chi_H(s) nearest s0 (Newton, central-difference slope)."""

    def c_of(s):
        f, chi = chi_top(r, s, h)
        return f[frozenset(h)] - chi[frozenset(h)]

    s = s0
    for _ in range(40):
        c = c_of(s)
        dc = (c_of(s + H_S) - c_of(s - H_S)) / (2 * H_S)
        if abs(dc) < 1e-14 * max(1.0, abs(c)):
            # |H| = 2: C_H = F(i) F(j) = 1 identically -> deleting chi_H removes the zero
            return None, "C_H constant (no zero)"
        step = c / dc
        s = s - step
        if abs(s - s0) > NEWTON_RADIUS:
            return None, "no zero within radius"
        if abs(step) < 1e-11:
            return complex(s), "converged"
    return complex(s), "not converged"


def cycle_score(r: Realization, s: complex, h: tuple) -> float:
    """E18 score: max spectral radius of the M-holonomy over simple directed cycles of h."""

    d, k = r.m_matrix(s)
    m = d @ k
    blocks = [CORE.index(b) for b in h]
    best = 0.0
    from itertools import permutations

    for size in range(2, len(blocks) + 1):
        for chosen in subsets(blocks):
            if len(chosen) != size:
                continue
            head, rest = chosen[0], chosen[1:]
            for tail in permutations(rest):
                cyc = (head, *tail)
                prod = np.eye(2, dtype=complex)
                for i, a in enumerate(cyc):
                    b = cyc[(i + 1) % len(cyc)]
                    prod = prod @ m[2 * a : 2 * a + 2, 2 * b : 2 * b + 2]
                best = max(best, float(np.abs(np.linalg.eigvals(prod)).max()))
    return best


def event_task(ev):
    pin_blas_threads()
    h = tuple(sorted(ev["subset"]))
    s_star = ev["s_star"]
    r = Realization(CORE, ev["theta"])
    f, sub = f_of(r, s_star, h)
    chi = cumulants(f, h)
    mu = boolean_mobius(f, h)
    back = moments(chi, h)
    recon = max(abs(back[k] - f[k]) for k in f)
    full = frozenset(h)
    nu = connected_share(chi, h)
    fp, _ = f_of(r, s_star + H_S, h)
    fm, _ = f_of(r, s_star - H_S, h)
    dfs = (fp[full] - fm[full]) / (2 * H_S)
    delta = {}
    for b in subsets(h):
        if len(b) >= 2:
            key = frozenset(b)
            delta[label(b)] = float(abs(chi[key] * f[full - key]) / abs(dfs))
    tol = 1e-8 * max(abs(v) for v in chi.values())
    naive = (-1) ** (len(h) - 1) * block_cycle_trace_sum(
        sub, [2] * len(h), tuple(range(len(h)))
    )
    remainder = chi[full] - naive
    z_del, status = newton_composite(r, s_star, h)
    dist = None if z_del is None else abs(z_del - s_star)
    if len(h) == 2:
        verdict = "CONNECTED (FORCED, |H| = 2)"
    else:
        verdict = "GENUINELY CONNECTED" if nu >= TAU else "COMPOSITE"
    material = None
    if len(h) >= 3:
        material = bool(nu >= TAU and (z_del is None or dist >= MATERIAL_DISTANCE))
    return {
        "event": ev["event"],
        "H": label(h),
        "size": len(h),
        "theta": ev["theta"],
        "s_star": [s_star.real, s_star.imag],
        "freq_hz": ev["freq_hz"],
        "abs_F_H_at_s_star": abs(f[full]),
        "reconstruction_residual": recon,
        "F": {label(k): [v.real, v.imag] for k, v in f.items() if k},
        "chi": {label(k): [v.real, v.imag, abs(v)] for k, v in chi.items()},
        "mu": {label(k): [v.real, v.imag, abs(v)] for k, v in mu.items() if k},
        "abs_chi_H": abs(chi[full]),
        "abs_composite_H": abs(composite_part(chi, h)),
        "nu_H": nu,
        "label_primary_0.10": verdict,
        "label_0.05": verdict
        if len(h) == 2
        else ("GENUINELY CONNECTED" if nu >= TAU_LO else "COMPOSITE"),
        "label_0.20": verdict
        if len(h) == 2
        else ("GENUINELY CONNECTED" if nu >= TAU_HI else "COMPOSITE"),
        "delta_B_per_s": delta,
        "dominant_support_delta": max(delta, key=delta.get),
        "abs_dF_ds": abs(dfs),
        "support_connected_thm23": cumulant_support_connected(chi, h, tol),
        "naive_single_cycle": [naive.real, naive.imag, abs(naive)],
        "remainder_R_H": [remainder.real, remainder.imag, abs(remainder)],
        "remainder_share": abs(remainder) / (abs(remainder) + abs(naive)),
        "deletion_zero": None if z_del is None else [z_del.real, z_del.imag],
        "deletion_status": status,
        "deletion_distance": dist,
        "materially_necessary": material,
        "abs_mu_H": abs(mu[full]),
        "cycle_score_M": cycle_score(r, s_star, h),
        "closure_distance": float(np.min(np.abs(np.linalg.eigvals(sub) + 1.0))),
    }


def flag_crosscheck(ev):
    """F(B) at FLAG from the PortActionSpace machinery (frozen F2c script) vs C2."""

    pin_blas_threads()
    from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
    from ibr_cycles.models.port_admittance import build_action_space

    kw = kwargs_of(ev["theta"])
    base = solve_case(ReplacementPlan.of({}), **kw)
    full = solve_case(ReplacementPlan.of({b: 1 for b in CORE}), **kw)
    space = build_action_space(base, full, CORE)
    m = space.m(ev["s_star"])
    eye = np.eye(m.shape[0], dtype=complex)
    total = eye + m
    sb = np.zeros_like(total)
    for k in range(4):
        sb[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = total[
            2 * k : 2 * k + 2, 2 * k : 2 * k + 2
        ]
    q = np.linalg.solve(sb, total) - eye
    vals = characteristic_values(q, [2] * 4)
    port = {frozenset(CORE[i] for i in k): v for k, v in vals.items()}
    r = Realization(CORE, ev["theta"])
    c2, _ = f_of(r, ev["s_star"], CORE)
    diff = max(abs(port[k] - c2[k]) for k in c2)
    chi_p = cumulants(port, CORE)
    chi_c = cumulants(c2, CORE)
    return {
        "max_abs_F_difference": diff,
        "chi_V_c2": abs(chi_c[frozenset(CORE)]),
        "chi_V_portspace": abs(chi_p[frozenset(CORE)]),
        "nu_V_c2": connected_share(chi_c, CORE),
        "nu_V_portspace": connected_share(chi_p, CORE),
    }


def main(argv) -> int:
    started = time.time()
    events = frozen_events()
    with Pool(min(11, len(events) + 1)) as pool:
        async_rows = pool.map_async(event_task, events, chunksize=1)
        flag = next(e for e in events if e["event"].startswith("FLAG"))
        cross = pool.apply_async(flag_crosscheck, (flag,))
        rows = async_rows.get()
        cross = cross.get()
    write_json(
        OUT / "CC03_boundaries.json",
        {
            "events": rows,
            "flag_crosscheck": cross,
            "elapsed_s": round(time.time() - started, 1),
        },
    )
    flat = []
    for e in rows:
        flat.append(
            {
                "event": e["event"],
                "H": e["H"],
                "freq_hz": round(e["freq_hz"], 4),
                "abs_F_H": e["abs_F_H_at_s_star"],
                "recon_resid": e["reconstruction_residual"],
                "abs_chi_H": e["abs_chi_H"],
                "abs_mu_H": e["abs_mu_H"],
                "nu_H": e["nu_H"],
                "label_0.10": e["label_primary_0.10"],
                "label_0.05": e["label_0.05"],
                "label_0.20": e["label_0.20"],
                "dominant_delta": e["dominant_support_delta"],
                "delta_H": e["delta_B_per_s"][e["H"]],
                "thm23_connected": e["support_connected_thm23"],
                "remainder_share": e["remainder_share"],
                "deletion_distance": e["deletion_distance"],
                "deletion_Re": None
                if e["deletion_zero"] is None
                else e["deletion_zero"][0],
                "deletion_status": e["deletion_status"],
                "materially_necessary": e["materially_necessary"],
                "cycle_score_M": e["cycle_score_M"],
                "closure_distance": e["closure_distance"],
            }
        )
    table = pd.DataFrame(flat)
    table.to_csv(OUT / "CC03_boundaries.csv", index=False)
    pd.set_option("display.width", 250)
    print(table.to_string(index=False))
    print("FLAG cross-check:", cross)
    for e in rows:
        print(
            e["event"], e["H"], {k: round(v, 5) for k, v in e["delta_B_per_s"].items()}
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
