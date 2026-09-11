"""MC02: Monte Carlo validation of the theorems on the frozen IEEE-39 model (P1-P4).

Pass rules: configs/ias2026/paper_mc_validation_v1.yaml (frozen, commit 7a772808).
Random policy points theta = (g, k, t, h) from the declared box and random
portfolios of the core {30, 33, 35, 37}; every case is solved on the direct path
(solve_config, R_none: matched dispatch, no condenser).

    P1  Theorem 1 on the physical model: center subspace, Jordan entry,
        spectral identity, and the exact NONLINEAR symmetries at random
        off-equilibrium states
    P2  Theorem 2 / Proposition 3: random policy segments; every change of the
        hypergraph is bracketed by a transverse eigenvalue crossing Re = 0
    P3  Theorem 6: det T_S / det T_0 = det(I + M_SS) and Moebius = minor sums at
        random theta and random s
    P4  nonlinear transverse recovery: TDS label vs sign of alpha_perp

Usage: python MC02_ieee39.py P1 [P2 P3 P4]
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd
import yaml
from _mc import CFG, CONFIGS, OMEGA_B, SEED, WORKERS, MCExperiment, out_dir, write_json
from BC02_binary_representation import _b39, port_t  # noqa: E402
from scipy.optimize import linear_sum_assignment

from _f7_common import SUBSETS, Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from ibr_cycles.certification.binary import (  # noqa: E402
    all_vertices,
    build_common_realization,
)
from ibr_cycles.certification.classify import SAFETY, classify_spectrum  # noqa: E402
from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import (  # noqa: E402
    center_subspace,
    transverse_operator,
)
from ibr_cycles.diagnosis.composability import hypergraph_label  # noqa: E402
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.nonlinear.manifold import chart, psi  # noqa: E402

OUT = out_dir("MC02_ieee39")
IEEE = CFG["ieee39"]
MARGIN = 2e-3
_C: dict = {}
CORE = (30, 33, 35, 37)
LO = np.array([0.0, 0.5, 0.5, 0.0])
HI = np.array([1.0, 2.3, 3.0, 2.0])


def _init():
    pin_blas_threads()
    _C["cfg"] = {c["name"]: c for c in r_configs()}["R_none"]


def draw_theta(rng):
    u = rng.random(4)
    g = u[0] ** 2
    return np.array([g, *(LO[1:] + u[1:] * (HI[1:] - LO[1:]))])


def theta_of(v):
    return Theta(float(v[0]), float(v[1]), float(v[2]), float(v[3]))


def lab(s):
    return "+".join(map(str, s)) or "BASE"


def solve(members, v):
    return solve_config(tuple(members), theta_of(v), _C["cfg"])


def transverse_status(case):
    """(status, alpha_perp, spectrum) with the frozen classifier near the axis."""

    dae, x, z = case.dae, case.equilibrium.x, case.equilibrium.z
    j1 = central_difference_jacobians(dae, x, z, {})
    j2 = central_difference_jacobians(dae, x, z, {}, scale_x=2.0, scale_z=2.0)
    a = j1.fx - j1.fz @ np.linalg.solve(j1.gz, j1.gx)
    d = a - (j2.fx - j2.fz @ np.linalg.solve(j2.gz, j2.gx))
    r_x, _ = rotation_generator(dae, z)
    tr = transverse_operator(a, r_x, frequency_partner(dae).w)
    ev = np.linalg.eigvals(tr.a_perp)
    if np.abs(ev.real).min() < MARGIN:
        st = classify_spectrum(tr.a_perp, tr.z.T @ d @ tr.z, SAFETY).status
    else:
        st = "UNSTABLE" if (ev.real > 0).any() else "STABLE"
    return st, float(ev.real.max()), ev


# ----------------------------------------------------------------------- P1 --
def rotate_z(z, phi):
    c, s = np.cos(phi), np.sin(phi)
    out = np.empty_like(z)
    out[0::2] = c * z[0::2] - s * z[1::2]
    out[1::2] = s * z[0::2] + c * z[1::2]
    return out


def p1_task(task):
    d, seed = task
    rng = np.random.default_rng(seed)
    v = draw_theta(rng)
    members = SUBSETS[int(rng.integers(len(SUBSETS)))]
    case = solve(members, v)
    dae, xs, zs = case.dae, case.equilibrium.x, case.equilibrium.z
    a = case.system.A
    r_x, _ = rotation_generator(dae, zs)
    w = frequency_partner(dae).w
    cs = center_subspace(a, r_x, w, left=False)
    tr = transverse_operator(a, r_x, w)
    ev_full = np.linalg.eigvals(a)
    # amendment v1.1: the structural Jordan pair splits numerically by
    # O(sqrt(eps ||A|| omega_b)); it is removed (the two eigenvalues nearest 0)
    # and its split is reported separately
    order = np.argsort(np.abs(ev_full))
    jordan_split = float(np.abs(ev_full[order[:2]]).max())
    ev_rest = ev_full[order[2:]]
    ev_perp = np.linalg.eigvals(tr.a_perp)
    c = np.abs(ev_rest[:, None] - ev_perp[None, :]) / (1 + np.abs(ev_perp[None, :]))
    r, cc = linear_sum_assignment(c)
    spec_id = float(c[r, cc].max())
    # nonlinear symmetries at a random off-equilibrium state
    ch = chart(dae, xs, zs)
    x = xs + 1e-2 * rng.standard_normal(xs.size)
    zx = psi(dae, x, zs, ch.lu)
    fx = dae.f(x, zx, {})
    phi = rng.uniform(-np.pi, np.pi)
    xr = x + phi * r_x
    zr = psi(dae, xr, rotate_z(zx, phi), ch.lu)
    rot_res = float(np.linalg.norm(dae.f(xr, zr, {}) - fx) / (1 + np.linalg.norm(fx)))
    cshift = rng.uniform(-0.02, 0.02)
    xc = x + cshift * w
    zc = psi(dae, xc, zx, ch.lu)
    drift_res = float(
        np.linalg.norm(dae.f(xc, zc, {}) - fx - cshift * OMEGA_B * r_x)
        / (1 + np.linalg.norm(fx))
    )
    return {
        "draw": d,
        "g": v[0],
        "k": v[1],
        "t": v[2],
        "h": v[3],
        "subset": lab(members),
        "dim_C": cs.dimension,
        "jordan_rel": abs(cs.jordan[0, 1] / OMEGA_B - 1.0),
        "invariance": cs.invariance_residual,
        "spectral_identity": spec_id,
        "jordan_split": jordan_split,
        "nonlinear_rotation": rot_res,
        "nonlinear_drift": drift_res,
        "phi": phi,
        "c": cshift,
        "alpha_perp": float(np.linalg.eigvals(tr.a_perp).real.max()),
    }


def run_p1():
    spec = IEEE["P1_structure"]
    seeds = np.random.SeedSequence(SEED + 1).spawn(spec["draws"])
    tasks = [(d, int(s.generate_state(1)[0])) for d, s in enumerate(seeds)]
    with Pool(WORKERS, initializer=_init) as pool:
        rows = pool.map(p1_task, tasks, chunksize=1)
    df = pd.DataFrame(rows)
    res = {
        "draws": len(df),
        "dim_C_2": int((df.dim_C == 2).sum()),
        "jordan_rel_max": float(df.jordan_rel.max()),
        "invariance_max": float(df.invariance.max()),
        "spectral_identity_max": float(df.spectral_identity.max()),
        "jordan_split_max": float(df.jordan_split.max()),
        "nonlinear_rotation_max": float(df.nonlinear_rotation.max()),
        "nonlinear_drift_max": float(df.nonlinear_drift.max()),
        "unstable_draws": int((df.alpha_perp > 0).sum()),
    }
    res["PASS"] = bool(
        res["dim_C_2"] == len(df)
        and res["jordan_rel_max"] <= 1e-6
        and res["invariance_max"] <= 1e-9
        and res["spectral_identity_max"] <= 1e-6
        and res["nonlinear_rotation_max"] <= 1e-8
        and res["nonlinear_drift_max"] <= 1e-8
    )
    return df, res


# ----------------------------------------------------------------------- P2 --
def lattice_classes(stable):
    from collections import deque

    reach, queue = {(): stable[()]}, deque([()] if stable[()] else [])
    while queue:
        s = queue.popleft()
        for c in CORE:
            if c not in s:
                t = tuple(sorted(s + (c,)))
                if stable[t] and t not in reach:
                    reach[t] = True
                    queue.append(t)
    return {
        t: (
            stable[t],
            bool(reach.get(t)),
            all(
                stable[tuple(sorted(q))]
                for k in range(len(t) + 1)
                for q in combinations(t, k)
            ),
        )
        for t in SUBSETS
    }


def h_label(status):
    if status[()] != "STABLE":
        return f"BASE_{status[()]}"
    uns = [s for s, v in status.items() if v == "UNSTABLE"]
    mins = [s for s in uns if not any(set(r) < set(s) for r in uns)]
    return hypergraph_label([frozenset(e) for e in mins]) + (
        "?" if any(v not in ("STABLE", "UNSTABLE") for v in status.values()) else ""
    )


def p2_task(task):
    seg, seed = task
    rng = np.random.default_rng(seed)
    va = draw_theta(rng)
    span = HI - LO
    vb = np.clip(va + rng.uniform(-0.15, 0.15, 4) * span, LO, HI)
    npts = IEEE["P2_segments"]["points_per_segment"]
    grid = np.linspace(0.0, 1.0, npts)
    status, alpha, gzc, labels_h, abc = [], [], [], [], []
    for sfrac in grid:
        v = va + sfrac * (vb - va)
        st, al = {}, {}
        worst = 0.0
        for m in SUBSETS:
            case = solve(m, v)
            st[m], al[m], _ = transverse_status(case)
            worst = max(worst, float(getattr(case, "gz_condition", np.nan) or np.nan))
        status.append(st)
        alpha.append(al)
        gzc.append(worst)
        labels_h.append(h_label(st))
        if st[()] == "STABLE" and all(x in ("STABLE", "UNSTABLE") for x in st.values()):
            cls = lattice_classes({m: st[m] == "STABLE" for m in SUBSETS})
            abc.append(
                (
                    sum(c[0] for c in cls.values()),
                    sum(c[0] and not c[1] for c in cls.values()),
                    sum(c[1] and not c[2] for c in cls.values()),
                )
            )
    events = []
    for j in range(npts - 1):
        if labels_h[j] == labels_h[j + 1]:
            continue
        changed = [
            m
            for m in SUBSETS
            if status[j][m] != status[j + 1][m]
            and {status[j][m], status[j + 1][m]} == {"STABLE", "UNSTABLE"}
        ]
        unresolved = any(
            "?" in lbl or lbl.startswith("BASE_")
            for lbl in (labels_h[j], labels_h[j + 1])
        )
        if not changed:
            events.append(
                {
                    "segment": seg,
                    "j": j,
                    "subset": "",
                    "bracketed": False,
                    "unresolved_endpoint": unresolved,
                }
            )
            continue
        for m in changed:
            lo, hi = grid[j], grid[j + 1]
            s_lo = status[j][m]
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if transverse_status(solve(m, va + mid * (vb - va)))[0] == s_lo:
                    lo = mid
                else:
                    hi = mid
            _, _, ev = transverse_status(solve(m, va + 0.5 * (lo + hi) * (vb - va)))
            k = int(np.argmin(np.abs(ev.real)))
            events.append(
                {
                    "segment": seg,
                    "j": j,
                    "subset": lab(m),
                    "bracketed": bool(abs(ev[k].real) <= 1e-6 and np.isfinite(ev[k])),
                    "re": float(ev[k].real),
                    "im": float(abs(ev[k].imag)),
                    "freq_hz": float(abs(ev[k].imag) / (2 * np.pi)),
                    "direction": f"{s_lo}->{status[j + 1][m]}",
                    "H_before": labels_h[j],
                    "H_after": labels_h[j + 1],
                    "unresolved_endpoint": unresolved,
                }
            )
    return {
        "segment": seg,
        "theta_a": va.tolist(),
        "theta_b": vb.tolist(),
        "H_sequence": labels_h,
        "distinct_H": len(set(labels_h)),
        "gz_condition_max": float(np.nanmax(gzc)) if np.isfinite(gzc).any() else None,
        "abc": abc,
        "events": events,
    }


def run_p2():
    spec = IEEE["P2_segments"]
    seeds = np.random.SeedSequence(SEED + 2).spawn(spec["segments"])
    tasks = [(d, int(s.generate_state(1)[0])) for d, s in enumerate(seeds)]
    with Pool(WORKERS, initializer=_init) as pool:
        segs = pool.map(p2_task, tasks, chunksize=1)
    ev = pd.DataFrame([e for s in segs for e in s["events"]])
    write_json(OUT / "P2_segments.json", segs)
    abc = (
        np.array([a for s in segs for a in s["abc"]])
        if any(s["abc"] for s in segs)
        else np.zeros((0, 3))
    )
    decided = ev[~ev.unresolved_endpoint] if len(ev) else ev
    res = {
        "segments": len(segs),
        "points": int(sum(len(s["H_sequence"]) for s in segs)),
        "H_changes": int(len(ev.drop_duplicates(["segment", "j"]))) if len(ev) else 0,
        "located_crossings": int(ev.subset.astype(bool).sum()) if len(ev) else 0,
        "exceptions_decided": int((~decided.bracketed).sum()) if len(decided) else 0,
        "exceptions_at_unresolved_endpoints": int(
            (~ev[ev.unresolved_endpoint].bracketed).sum()
        )
        if len(ev)
        else 0,
        "crossing_freq_hz_quantiles": (
            ev.loc[ev.subset.astype(bool), "freq_hz"]
            .quantile([0, 0.25, 0.5, 0.75, 1])
            .round(4)
            .tolist()
            if len(ev)
            else []
        ),
        "real_crossings": int(((ev.im < 1e-6) & ev.subset.astype(bool)).sum())
        if len(ev)
        else 0,
        "segments_with_H_change": int(sum(s["distinct_H"] > 1 for s in segs)),
        "distinct_H_total": len({h for s in segs for h in s["H_sequence"]}),
        "points_classified": int(len(abc)),
        "targets_A": int(abc[:, 0].sum()) if len(abc) else 0,
        "targets_A_not_B": int(abc[:, 1].sum()) if len(abc) else 0,
        "targets_B_not_C": int(abc[:, 2].sum()) if len(abc) else 0,
        "gz_condition_max": float(
            np.nanmax([s["gz_condition_max"] or np.nan for s in segs])
        ),
    }
    res["PASS"] = res["exceptions_decided"] == 0
    return ev, res


# ----------------------------------------------------------------------- P3 --
def p3_task(task):
    d, seed = task
    rng = np.random.default_rng(seed)
    v = draw_theta(rng)
    point = dict(zip("gkth", map(float, v), strict=True))
    _, kwargs = _b39(point)
    cr = build_common_realization(CORE, **kwargs)
    jac = {}
    for delta, s in all_vertices(CORE):
        j, _, _ = cr.jacobian(delta, "C2")
        jac[tuple(sorted(s))] = j
    net = cr.case.dae.network
    items = CORE
    channels = {b: [2 * net.position(b), 2 * net.position(b) + 1] for b in items}
    out = []
    for _ in range(IEEE["P3_ports"]["s_per_draw"]):
        s = complex(rng.uniform(0.02, 1.0), rng.uniform(0.5, 10.0))
        ports = {t: port_t(jac[t], s) for t in jac}
        t0 = ports[()]
        sign0, ld0 = np.linalg.slogdet(t0)
        ratio = {}
        for t, p in ports.items():
            sg, ld = np.linalg.slogdet(p)
            ratio[t] = (sg / sign0) * np.exp(ld - ld0)
        idx = [c for b in items for c in channels[b]]
        kmat = np.linalg.solve(t0, np.eye(t0.shape[0]))[np.ix_(idx, idx)]
        dmat = np.zeros((2 * len(items), 2 * len(items)), complex)
        for i, b in enumerate(items):
            ch = channels[b]
            dmat[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = (ports[(b,)] - t0)[
                np.ix_(ch, ch)
            ]
        mm = dmat @ kmat
        block_of = {2 * i + c: items[i] for i in range(len(items)) for c in (0, 1)}
        allcols = tuple(range(2 * len(items)))
        minors = {
            u: np.linalg.det(mm[np.ix_(u, u)]) if u else 1.0 + 0j
            for r in range(len(allcols) + 1)
            for u in combinations(allcols, r)
        }
        vert_err = mob_err = 0.0
        for t in jac:
            cols = [2 * items.index(b) + c for b in t for c in (0, 1)]
            pred = (
                np.linalg.det(np.eye(len(cols)) + mm[np.ix_(cols, cols)])
                if cols
                else 1.0
            )
            vert_err = max(vert_err, abs(ratio[t] - pred) / max(1.0, abs(ratio[t])))
            mu = sum(
                (-1) ** (len(t) - len(q)) * ratio[tuple(sorted(q))]
                for k in range(len(t) + 1)
                for q in combinations(t, k)
            )
            touch = sum(
                val for u, val in minors.items() if {block_of[c] for c in u} == set(t)
            )
            mob_err = max(mob_err, abs(mu - touch) / max(1.0, abs(ratio[t])))
        out.append(
            {
                "draw": d,
                "g": v[0],
                "k": v[1],
                "t": v[2],
                "h": v[3],
                "s": str(s),
                "vertex_identity_rel": vert_err,
                "moebius_minor_rel": mob_err,
                "full_order_coeff_abs": float(
                    abs(
                        sum(
                            val
                            for u, val in minors.items()
                            if {block_of[c] for c in u} == set(items)
                        )
                    )
                ),
            }
        )
    return out


def run_p3():
    spec = IEEE["P3_ports"]
    seeds = np.random.SeedSequence(SEED + 3).spawn(spec["draws"])
    tasks = [(d, int(s.generate_state(1)[0])) for d, s in enumerate(seeds)]
    with Pool(WORKERS, initializer=_init) as pool:
        rows = [r for out in pool.map(p3_task, tasks, chunksize=1) for r in out]
    df = pd.DataFrame(rows)
    res = {
        "draws": spec["draws"],
        "evaluations": len(df),
        "vertex_identity_rel_max": float(df.vertex_identity_rel.max()),
        "moebius_minor_rel_max": float(df.moebius_minor_rel.max()),
        "full_order_coeff_abs_range": [
            float(df.full_order_coeff_abs.min()),
            float(df.full_order_coeff_abs.max()),
        ],
    }
    res["PASS"] = bool(
        res["vertex_identity_rel_max"] <= 1e-8 and res["moebius_minor_rel_max"] <= 1e-8
    )
    return df, res


# ----------------------------------------------------------------------- P4 --
def p4_sample(task):
    d, seed = task
    rng = np.random.default_rng(seed)
    v = draw_theta(rng)
    m = SUBSETS[int(rng.integers(len(SUBSETS)))]
    st, al, _ = transverse_status(solve(m, v))
    return {
        "cand": d,
        "theta": v.tolist(),
        "subset": list(m),
        "status": st,
        "alpha_perp": al,
    }


def p4_task(row):
    from ibr_cycles.nonlinear.tds import pulse_dae, simulate

    spec = yaml.safe_load(
        (CONFIGS / "final_nonlinear_composability_v1.yaml").read_text(encoding="utf-8")
    )["disturbances"]["D1"]
    case = solve(tuple(row["subset"]), np.array(row["theta"]))
    out = simulate(case, pulse_dae(case, "D1", 5.0, spec))
    s, tr = out["summary"], out["trace"]
    rate = np.nan
    if len(tr):
        post = tr[tr.t >= 5.5]
        if len(post) > 10:
            from ibr_cycles.nonlinear.tds import _envelope_growth

            rate = _envelope_growth(post.t.to_numpy(), post.D.to_numpy())
    return {
        **row,
        "label": s["label"],
        "reason": s["reason"],
        "t_label": s["t_label"],
        "envelope_rate": rate,
    }


def run_p4():
    spec = IEEE["P4_tds"]
    seeds = np.random.SeedSequence(SEED + 4).spawn(4000)
    tasks = [(d, int(s.generate_state(1)[0])) for d, s in enumerate(seeds)]
    with Pool(WORKERS, initializer=_init) as pool:
        cands = []
        for out in pool.imap(p4_sample, tasks, chunksize=4):
            cands.append(out)
            stab = [
                c for c in cands if c["status"] == "STABLE" and c["alpha_perp"] <= -0.10
            ]
            unst = [
                c
                for c in cands
                if c["status"] == "UNSTABLE" and c["alpha_perp"] >= 0.10
            ]
            if len(stab) >= spec["draws"] // 2 and len(unst) >= spec["draws"] // 2:
                break
        pool.terminate()
    chosen = stab[: spec["draws"] // 2] + unst[: spec["draws"] // 2]
    with Pool(WORKERS, initializer=_init) as pool:
        rows = pool.map(p4_task, chosen, chunksize=1)
    df = pd.DataFrame(rows)
    df["expected"] = np.where(df.alpha_perp < 0, "RECOVERS", "not RECOVERS")
    decided = df[df.label != "NUMERICAL_FAILURE"]
    agree = np.where(
        decided.alpha_perp < 0, decided.label == "RECOVERS", decided.label != "RECOVERS"
    )
    res = {
        "candidates_sampled": len(cands),
        "runs": len(df),
        "numerical_failures": int((df.label == "NUMERICAL_FAILURE").sum()),
        "agreement": f"{int(agree.sum())}/{len(decided)}",
        "labels_stable": df[df.alpha_perp < 0].label.value_counts().to_dict(),
        "labels_unstable": df[df.alpha_perp > 0].label.value_counts().to_dict(),
        "rate_vs_alpha_corr": float(
            np.corrcoef(
                df.alpha_perp[np.isfinite(df.envelope_rate)],
                df.envelope_rate[np.isfinite(df.envelope_rate)],
            )[0, 1]
        ),
    }
    res["PASS"] = bool(agree.all())
    return df, res


RUNS = {"P1": run_p1, "P2": run_p2, "P3": run_p3, "P4": run_p4}


def main(argv) -> int:
    names = [a for a in argv if a in RUNS] or list(RUNS)
    exp = MCExperiment(
        name="MC02_ieee39",
        question=(
            "Do the theorems hold on the frozen IEEE-39 model at random policies "
            "and portfolios?"
        ),
        config={"blocks": names, "seed": SEED, "margin": MARGIN},
        workers=WORKERS,
    )
    summary_path = OUT / "MC02_summary.json"
    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {}
    )
    for name in names:
        t0 = time.time()
        df, res = RUNS[name]()
        df.to_csv(OUT / f"{name}.csv", index=False)
        res["elapsed_s"] = round(time.time() - t0, 1)
        summary[name] = res
        write_json(summary_path, summary)
        print(name, json.dumps(res, default=str), flush=True)
    exp.finish("COMPUTED", **summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
