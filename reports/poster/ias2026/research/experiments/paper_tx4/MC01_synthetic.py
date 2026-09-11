"""MC01: synthetic Monte Carlo validation of the paper's theorems (S1-S7).

Pass rules: configs/ias2026/paper_mc_validation_v1.yaml (frozen, commit 7a772808).
Each block draws random systems that satisfy exactly the hypotheses of one
theorem and checks its conclusion; controls that are expected to fail are
reported as such. Usage: python MC01_synthetic.py [S1 S2 ...] (default: all).
"""

from __future__ import annotations

import sys
import time
from collections import deque
from itertools import combinations
from types import SimpleNamespace

import numpy as np
import pandas as pd
from _mc import CFG, OMEGA_B, SEED, MCExperiment, out_dir, write_json
from scipy.optimize import linear_sum_assignment

from ibr_cycles.certification.port_origin import relocated_port  # noqa: E402
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402

OUT = out_dir("MC01_synthetic")
SYN = CFG["synthetic"]


# ------------------------------------------------------------------ helpers --
def subsets(items):
    for r in range(len(items) + 1):
        yield from combinations(items, r)


def haar(n, rng):
    q, r = np.linalg.qr(rng.standard_normal((n, n)))
    return q * np.sign(np.diag(r))


def well_conditioned(n, rng, spread=0.3):
    return haar(n, rng) @ (
        np.eye(n) + spread * rng.standard_normal((n, n)) / np.sqrt(n)
    )


def real_form(eigs):
    """Block-diagonal real matrix with the given eigenvalues (complex ones in pairs,
    listed once with positive imaginary part)."""

    blocks = []
    for lam in eigs:
        lam = complex(lam)
        if abs(lam.imag) > 0:
            blocks.append(np.array([[lam.real, lam.imag], [-lam.imag, lam.real]]))
        else:
            blocks.append(np.array([[lam.real]]))
    n = sum(b.shape[0] for b in blocks)
    d = np.zeros((n, n))
    i = 0
    for b in blocks:
        k = b.shape[0]
        d[i : i + k, i : i + k] = b
        i += k
    return d


def random_spectrum(n, rng, lo=-2.0, hi=1.0, keep_away=1e-2):
    """Eigenvalue list (pairs listed once) filling dimension n."""

    eigs, dim = [], 0
    while dim < n:
        if n - dim >= 2 and rng.random() < 0.5:
            re = rng.uniform(lo, hi)
            while abs(re) < keep_away:
                re = rng.uniform(lo, hi)
            eigs.append(complex(re, rng.uniform(0.5, 10.0)))
            dim += 2
        else:
            re = rng.uniform(lo, hi)
            while abs(re) < keep_away:
                re = rng.uniform(lo, hi)
            eigs.append(re)
            dim += 1
    return eigs


def dim_of(eigs):
    return sum(2 if abs(complex(e).imag) > 0 else 1 for e in eigs)


def planted_system(transverse_eigs, rng, v=None, p=None, x=None):
    """A = P [[J, X], [0, B]] P^-1 with an exact rotation / Jordan partner pair."""

    m = dim_of(transverse_eigs)
    n = m + 2
    v = well_conditioned(m, rng) if v is None else v
    b = v @ real_form(transverse_eigs) @ np.linalg.inv(v)
    p = well_conditioned(n, rng) if p is None else p
    x = 0.5 * rng.standard_normal((2, m)) if x is None else x
    top = np.zeros((n, n))
    top[0, 1] = OMEGA_B
    top[:2, 2:] = x
    top[2:, 2:] = b
    a = p @ top @ np.linalg.inv(p)
    return a, p[:, 0].copy(), p[:, 1].copy(), (v, p, x)


def matched_distance(a, b):
    c = np.abs(a[:, None] - b[None, :])
    r, cc = linear_sum_assignment(c)
    return float(c[r, cc].max())


def alpha(m):
    return float(np.linalg.eigvals(m).real.max())


# ----------------------------------------------------------------------- S1 --
def s1(rng):
    spec = SYN["S1_transverse_quotient"]
    rows = []
    for d in range(spec["draws"]):
        n = int(rng.integers(spec["n_states"][0], spec["n_states"][1] + 1))
        planted = None
        if d % 2 == 1:  # plant one physical near-zero real mode
            u = float(np.exp(rng.uniform(np.log(1e-5), np.log(5e-4))))
            planted = u if rng.random() < 0.5 else -u
            eigs = random_spectrum(n - 3, rng) + [planted]
        else:
            eigs = random_spectrum(n - 2, rng)
        a, r_x, w, _ = planted_system(eigs, rng)
        true_rhp = sum(
            (2 if abs(complex(e).imag) > 0 else 1) for e in eigs if complex(e).real > 0
        )
        tr = transverse_operator(a, r_x, w)
        ev_perp = np.linalg.eigvals(tr.a_perp)
        rhp_perp = int((ev_perp.real > 0).sum())
        norm = float(np.linalg.norm(a, 2))
        ident = []
        for _ in range(3):
            s = np.exp(rng.uniform(np.log(0.1), np.log(10.0))) * np.exp(
                1j * rng.uniform(-np.pi, np.pi)
            )
            s1_, l1 = np.linalg.slogdet(s * np.eye(n) - a)
            s2_, l2 = np.linalg.slogdet(s * np.eye(n - 2) - tr.a_perp)
            ratio = (s1_ / (s2_ * s**2 / abs(s) ** 2)) * np.exp(
                l1 - l2 - 2 * np.log(abs(s))
            )
            ident.append(abs(ratio - 1.0))
        o = haar(n - 2, rng)
        basis = matched_distance(np.linalg.eigvals(o.T @ tr.a_perp @ o), ev_perp) / norm
        ev_full = np.linalg.eigvals(a)
        keep = np.abs(ev_full) > 1e-3
        cutoff_rhp = int((ev_full[keep].real > 0).sum())
        rows.append(
            {
                "draw": d,
                "n": n,
                "planted": planted,
                "true_rhp": true_rhp,
                "rhp_perp": rhp_perp,
                "count_exact": rhp_perp == true_rhp,
                "det_identity_rel": float(max(ident)),
                "basis_rel": basis,
                "cutoff_rhp": cutoff_rhp,
                "cutoff_wrong": cutoff_rhp != true_rhp,
                "coupling": tr.coupling_residual,
            }
        )
    df = pd.DataFrame(rows)
    pos = df.planted.notna() & (df.planted > 0)
    res = {
        "draws": len(df),
        "count_exact": int(df.count_exact.sum()),
        "det_identity_rel_max": float(df.det_identity_rel.max()),
        "basis_rel_max": float(df.basis_rel.max()),
        "PASS": bool(
            df.count_exact.all()
            and df.det_identity_rel.max() <= 1e-8
            and df.basis_rel.max() <= 1e-8
        ),
        "control_cutoff_wrong_in_positive_planted": (
            f"{int(df.cutoff_wrong[pos].sum())}/{int(pos.sum())}"
        ),
        "control_cutoff_wrong_elsewhere": (
            f"{int(df.cutoff_wrong[~pos].sum())}/{int((~pos).sum())}"
        ),
    }
    return df, res


# ------------------------------------------------------- lattice utilities --
def classes(stable: dict, items):
    """A (final stable), B (safe one-at-a-time path), C (every subset stable)."""

    reach, queue = {(): stable[()]}, deque([()] if stable[()] else [])
    while queue:
        s = queue.popleft()
        for c in items:
            if c not in s:
                t = tuple(sorted(s + (c,)))
                if stable[t] and t not in reach:
                    reach[t] = True
                    queue.append(t)
    return {
        t: (stable[t], bool(reach.get(t)), all(stable[q] for q in subsets(t)))
        for t in stable
    }


def minimal_unsafe(stable: dict):
    uns = [s for s, v in stable.items() if not v]
    return [s for s in uns if not any(set(r) < set(s) for r in uns)]


# ----------------------------------------------------------------------- S2 --
def s2(rng):
    spec = SYN["S2_hypergraph_and_order_safety"]
    n = spec["n_states"]
    rows = []
    viol = {"antichain": 0, "free_implies_stable": 0, "C_iff_free": 0, "hierarchy": 0}
    for d in range(spec["draws"]):
        m = int(rng.integers(spec["candidates"][0], spec["candidates"][1] + 1))
        g = rng.standard_normal((n, n)) / np.sqrt(n)
        a0 = g - (alpha(g) - rng.uniform(-1.0, -0.2)) * np.eye(n)
        acts = [
            rng.uniform(0.3, 1.5)
            * (rng.standard_normal((n, 2)) @ rng.standard_normal((2, n)))
            / np.sqrt(n)
            for _ in range(m)
        ]
        items = tuple(range(m))
        stable = {
            s: alpha(a0 + sum((acts[i] for i in s), np.zeros((n, n)))) < 0
            for s in subsets(items)
        }
        h = minimal_unsafe(stable)
        cls = classes(stable, items)
        if any(set(e) < set(f) for e in h for f in h if e != f):
            viol["antichain"] += 1
        for t, (ca, cb, cc) in cls.items():
            free = not any(set(e) <= set(t) for e in h)
            if free and not ca:
                viol["free_implies_stable"] += 1
            if cc != free:
                viol["C_iff_free"] += 1
            if (cc and not cb) or (cb and not ca):
                viol["hierarchy"] += 1
        a_not_b = sum(1 for v in cls.values() if v[0] and not v[1])
        b_not_c = sum(1 for v in cls.values() if v[1] and not v[2])
        rows.append(
            {
                "draw": d,
                "m": m,
                "kappa": min((len(e) for e in h), default=-1),
                "n_hyperedges": len(h),
                "targets_A": sum(v[0] for v in cls.values()),
                "A_not_B": a_not_b,
                "B_not_C": b_not_c,
                "non_hereditary": a_not_b + b_not_c > 0,
            }
        )
    df = pd.DataFrame(rows)
    res = {
        "draws": len(df),
        "violations": viol,
        "PASS": all(v == 0 for v in viol.values()),
        "families_non_hereditary": int(df.non_hereditary.sum()),
        "families_with_A_not_B": int((df.A_not_B > 0).sum()),
        "families_with_B_not_C": int((df.B_not_C > 0).sum()),
        "targets_A_total": int(df.targets_A.sum()),
        "targets_A_not_B_total": int(df.A_not_B.sum()),
        "targets_B_not_C_total": int(df.B_not_C.sum()),
        "kappa_counts": df.kappa.value_counts().sort_index().to_dict(),
    }
    return df, res


# ----------------------------------------------------------------------- S3 --
def s3(rng):
    spec = SYN["S3_boundary_localization"]
    m, n, grid = spec["candidates"], 10, spec["grid"]
    items = tuple(range(m))
    ps = np.linspace(0.0, 1.0, grid)
    rows, exceptions = [], 0
    for d in range(spec["draws"]):
        g = rng.standard_normal((n, n)) / np.sqrt(n)
        a0 = g - (alpha(g) - rng.uniform(-1.0, -0.2)) * np.eye(n)
        acts = [
            rng.uniform(0.3, 1.5) * rng.standard_normal((n, n)) / np.sqrt(n)
            for _ in range(m)
        ]
        slopes = [
            rng.uniform(0.3, 1.5) * rng.standard_normal((n, n)) / np.sqrt(n)
            for _ in range(m)
        ]

        def mat(s, p, a0=a0, acts=acts, slopes=slopes):
            return a0 + sum((acts[i] + p * slopes[i] for i in s), np.zeros((n, n)))

        al = {s: np.array([alpha(mat(s, p)) for p in ps]) for s in subsets(items)}
        hs = [
            tuple(sorted(minimal_unsafe({s: al[s][j] < 0 for s in al})))
            for j in range(grid)
        ]
        for j in range(grid - 1):
            if hs[j] == hs[j + 1]:
                continue
            changed = [s for s in al if (al[s][j] < 0) != (al[s][j + 1] < 0)]
            if not changed:
                exceptions += 1
                continue
            for s in changed:
                lo, hi = ps[j], ps[j + 1]
                sign_lo = al[s][j] < 0
                for _ in range(50):
                    mid = 0.5 * (lo + hi)
                    if (alpha(mat(s, mid)) < 0) == sign_lo:
                        lo = mid
                    else:
                        hi = mid
                ev = np.linalg.eigvals(mat(s, 0.5 * (lo + hi)))
                k = int(np.argmin(np.abs(ev.real)))
                ok = abs(ev[k].real) <= 1e-9 * (1 + abs(ev[k])) and np.isfinite(ev[k])
                exceptions += 0 if ok else 1
                rows.append(
                    {
                        "draw": d,
                        "p": 0.5 * (lo + hi),
                        "subset": "+".join(map(str, s)) or "BASE",
                        "re": float(ev[k].real),
                        "im": float(abs(ev[k].imag)),
                        "kind": "real" if abs(ev[k].imag) < 1e-9 else "complex",
                        "ok": ok,
                    }
                )
    df = pd.DataFrame(rows)
    res = {
        "draws": spec["draws"],
        "located_crossings": len(df),
        "exceptions": exceptions,
        "PASS": exceptions == 0,
        "crossing_kinds": df.kind.value_counts().to_dict() if len(df) else {},
    }
    return df, res


# ----------------------------------------------------------------------- S4 --
def metzler_family(rng, m, n, signed):
    b0 = rng.uniform(0, 1, (n, n)) * (rng.random((n, n)) < 0.5)
    np.fill_diagonal(b0, 0.0)
    b0 = b0 - (alpha(b0) - rng.uniform(-2.0, -0.5)) * np.eye(n)
    ns = []
    for _ in range(m):
        nn = rng.uniform(0, 0.6, (n, n)) * (rng.random((n, n)) < 0.3)
        if signed:
            nn = nn * rng.choice([-1.0, 1.0], size=(n, n))
        ns.append(nn)
    while True:
        t = well_conditioned(n, rng, spread=1.0)
        if np.linalg.cond(t) <= 1e3:
            break
    ti = np.linalg.inv(t)
    return b0, ns, t, ti


def s4(rng):
    spec = SYN["S4_monotone_class"]
    m, n = spec["candidates"], spec["n_states"]
    items = tuple(range(m))
    out = {}
    rows = []
    for signed in (False, True):
        bad_pairs = families_bad = abc_viol = 0
        for d in range(spec["draws"]):
            b0, ns, t, ti = metzler_family(rng, m, n, signed)
            al = {
                s: alpha(t @ (b0 + sum((ns[i] for i in s), np.zeros((n, n)))) @ ti)
                for s in subsets(items)
            }
            bad = 0
            for s in al:
                for c in items:
                    if c not in s:
                        u = tuple(sorted(s + (c,)))
                        if al[s] > al[u] + 1e-9 * (1 + abs(al[s])):
                            bad += 1
            bad_pairs += bad
            families_bad += bad > 0
            if not signed:
                cls = classes({s: v < 0 for s, v in al.items()}, items)
                abc_viol += sum(1 for v in cls.values() if len(set(v)) > 1)
            rows.append(
                {"signed_control": signed, "draw": d, "non_monotone_edges": bad}
            )
        out["control_signed" if signed else "metzler"] = {
            "draws": spec["draws"],
            "non_monotone_edges": bad_pairs,
            "families_with_violation": families_bad,
            **({} if signed else {"A_B_C_disagreements": abc_viol}),
        }
    out["PASS"] = (
        out["metzler"]["non_monotone_edges"] == 0
        and out["metzler"]["A_B_C_disagreements"] == 0
    )
    return pd.DataFrame(rows), out


# ----------------------------------------------------------------------- S5 --
def s5(rng):
    spec = SYN["S5_complexity_reduction"]
    rows, mismatches = [], 0
    for d in range(spec["draws"]):
        n = int(rng.integers(spec["vertices"][0], spec["vertices"][1] + 1))
        p = rng.uniform(*spec["edge_probability"])
        k = int(rng.integers(spec["k"][0], spec["k"][1] + 1))
        w = np.triu((rng.random((n, n)) < p).astype(float), 1)
        w = w + w.T
        found = cliques = 0
        for size in range(1, k + 1):
            for s in combinations(range(n), size):
                dlt = np.zeros(n)
                dlt[list(s)] = 1.0
                mat = np.outer(dlt, dlt) * w - (k - 1) * np.eye(n)
                dest = np.linalg.eigvalsh(mat)[-1] >= -1e-9
                is_clique = size == k and all(w[i, j] for i, j in combinations(s, 2))
                found += dest
                cliques += is_clique
                mismatches += dest != is_clique
        rows.append(
            {
                "draw": d,
                "n": n,
                "p": p,
                "k": k,
                "destabilizing": found,
                "k_cliques": cliques,
            }
        )
    df = pd.DataFrame(rows)
    return df, {
        "draws": len(df),
        "mismatches": mismatches,
        "PASS": mismatches == 0,
        "graphs_with_k_clique": int((df.k_cliques > 0).sum()),
    }


# ----------------------------------------------------------------------- S6 --
def minor_table(mat, channels):
    return {
        u: (np.linalg.det(mat[np.ix_(u, u)]) if u else 1.0 + 0j)
        for u in subsets(tuple(channels))
    }


def s6(rng):
    spec = SYN["S6_principal_minors"]
    rows = []
    for d in range(spec["draws"]):
        m = int(rng.integers(spec["blocks"][0], spec["blocks"][1] + 1))
        p = spec["block_size"]
        n = m * p
        mat = (
            rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        ) / np.sqrt(2 * n)
        block_of = {c: c // p for c in range(n)}
        items = tuple(range(m))

        def ch(s, p=p):
            return [c for b in s for c in range(b * p, (b + 1) * p)]

        def objects(mm, n=n, items=items, block_of=block_of):
            minors = minor_table(mm, range(n))
            r = {
                s: np.linalg.det(np.eye(len(ch(s))) + mm[np.ix_(ch(s), ch(s))])
                if s
                else 1.0 + 0j
                for s in subsets(items)
            }
            r_minor = {
                s: sum(v for u, v in minors.items() if set(u) <= set(ch(s)))
                for s in subsets(items)
            }
            mu = {}
            for t in subsets(items):
                mu[t] = sum((-1) ** (len(t) - len(q)) * r[q] for q in subsets(t))
            touch = {
                t: sum(
                    v for u, v in minors.items() if {block_of[c] for c in u} == set(t)
                )
                for t in subsets(items)
            }
            return minors, r, r_minor, mu, touch

        minors, r, r_minor, mu, touch = objects(mat)
        sb = np.zeros((n, n))
        for b in items:
            blk = rng.standard_normal((p, p)) + 2.0 * np.eye(p)
            sb[b * p : (b + 1) * p, b * p : (b + 1) * p] = blk
        mat2 = np.linalg.solve(sb, mat @ sb)
        minors2, r2, _, mu2, _ = objects(mat2)
        scale = max(1.0, max(abs(v) for v in r.values()))
        split = [
            abs(minors2[u] - minors[u]) / max(abs(minors[u]), 1e-300)
            for u in minors
            if u and len({block_of[c] for c in u}) * p != len(u)
        ]
        rows.append(
            {
                "draw": d,
                "m": m,
                "vertex_identity_rel": max(abs(r[s] - r_minor[s]) for s in r) / scale,
                "moebius_touch_rel": max(abs(mu[t] - touch[t]) for t in mu) / scale,
                "invariance_rel": max(
                    max(abs(r[s] - r2[s]) for s in r),
                    max(abs(mu[t] - mu2[t]) for t in mu),
                )
                / scale,
                "split_minor_change_median": float(np.median(split))
                if split
                else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    return df, {
        "draws": len(df),
        "vertex_identity_rel_max": float(df.vertex_identity_rel.max()),
        "moebius_touch_rel_max": float(df.moebius_touch_rel.max()),
        "invariance_rel_max": float(df.invariance_rel.max()),
        "PASS": bool(
            df.vertex_identity_rel.max() <= 1e-10
            and df.moebius_touch_rel.max() <= 1e-10
            and df.invariance_rel.max() <= 1e-10
        ),
        "control_split_minor_change_median": float(
            df.split_minor_change_median.median()
        ),
    }


# ----------------------------------------------------------------------- S7 --
EVENTS = (
    "single_real_crossing",
    "double_real_crossing",
    "real_pair_coalescence",
    "complex_pair_crossing",
    "no_crossing",
)


def event_spectra(kind, rng):
    """(lo, hi) eigenvalue lists for the event block."""

    u1, u2 = rng.uniform(0.01, 0.2, 2)
    a, b = rng.uniform(0.02, 0.3), rng.uniform(0.5, 6.0)
    if kind == "single_real_crossing":
        return [-u1], [u1 * rng.uniform(0.5, 2.0)]
    if kind == "double_real_crossing":
        return [-u1, -u2], [u1, u2 * rng.uniform(0.5, 2.0)]
    if kind == "real_pair_coalescence":
        dd = rng.uniform(0.005, 0.05)
        return [a - dd, a + dd], [complex(a, dd)]
    if kind == "complex_pair_crossing":
        return [complex(-a, b)], [complex(a, b)]
    return [-u1], [-u1 * rng.uniform(0.5, 2.0)]


def s7(rng):
    spec = SYN["S7_zero_frequency_port"]
    rows = []
    for kind in EVENTS:
        for d in range(spec["draws_per_class"]):
            nx = int(rng.integers(10, 31))
            nz = int(rng.integers(6, 21))
            lo_ev, hi_ev = event_spectra(kind, rng)
            base = random_spectrum(nx - 2 - dim_of(lo_ev), rng, keep_away=0.05)
            # the rest of the spectrum moves slightly along the path, without
            # changing sign or type (multiplicative 5 % perturbation)
            base_hi = [
                complex(
                    complex(e).real * (1 + 0.05 * rng.standard_normal()),
                    complex(e).imag * (1 + 0.05 * rng.standard_normal()),
                )
                if abs(complex(e).imag) > 0
                else float(e) * (1 + 0.05 * rng.standard_normal())
                for e in base
            ]
            a_lo, r_x, w, (v, p, x) = planted_system(base + lo_ev, rng)
            m = dim_of(base + lo_ev)
            if dim_of(base_hi + hi_ev) != m:
                raise AssertionError("event changes the dimension")
            a_hi, _, _, _ = planted_system(base_hi + hi_ev, rng, v=v, p=p, x=x)
            gz = rng.standard_normal((nz, nz)) + 3.0 * np.eye(nz)
            gx = rng.standard_normal((nz, nx)) / np.sqrt(nx)
            fz = rng.standard_normal((nx, nz)) / np.sqrt(nz)
            corr = fz @ np.linalg.solve(gz, gx)
            verdicts = {}
            for tag, a in (("lo", a_lo), ("hi", a_hi)):
                jac = SimpleNamespace(fx=a + corr, fz=fz, gx=gx, gz=gz)
                for pair in ((1.0, 1.0), (2.0, 0.5)):
                    verdicts[(tag, pair)] = relocated_port(jac, r_x, w, *pair)
            flip = {
                pair: verdicts[("lo", pair)].t0_sign != verdicts[("hi", pair)].t0_sign
                for pair in ((1.0, 1.0), (2.0, 0.5))
            }
            dev = {
                pair: verdicts[("lo", pair)].device_sign
                != verdicts[("hi", pair)].device_sign
                for pair in ((1.0, 1.0), (2.0, 0.5))
            }

            def pos_real(ev):
                return sum(
                    1 for e in ev if abs(complex(e).imag) == 0 and complex(e).real > 0
                )

            parity_change = (
                pos_real(base + lo_ev) - pos_real(base_hi + hi_ev)
            ) % 2 == 1
            rows.append(
                {
                    "kind": kind,
                    "draw": d,
                    "nx": nx,
                    "nz": nz,
                    "parity_change": parity_change,
                    "flip": flip[(1.0, 1.0)],
                    "flip_pair2": flip[(2.0, 0.5)],
                    "device_flip": dev[(1.0, 1.0)] or dev[(2.0, 0.5)],
                    "identity_residual": max(
                        v_.identity_residual for v_ in verdicts.values()
                    ),
                }
            )
    df = pd.DataFrame(rows)
    kept = df[~df.device_flip]
    per = {}
    for kind in EVENTS:
        k = kept[kept.kind == kind]
        per[kind] = {
            "kept": len(k),
            "device_flip_excluded": int((df[df.kind == kind].device_flip).sum()),
            "flips": int(k.flip.sum()),
            "beta_inconsistent": int((k.flip != k.flip_pair2).sum()),
            "flip_equals_parity": int((k.flip == k.parity_change).sum()),
        }
    ok = all(
        v["flip_equals_parity"] == v["kept"] and v["beta_inconsistent"] == 0
        for v in per.values()
    )
    ok &= per["single_real_crossing"]["flips"] == per["single_real_crossing"]["kept"]
    return df, {
        "per_class": per,
        "identity_residual_max": float(df.identity_residual.max()),
        "PASS": bool(ok),
    }


BLOCKS = {"S1": s1, "S2": s2, "S3": s3, "S4": s4, "S5": s5, "S6": s6, "S7": s7}


def main(argv) -> int:
    names = [a for a in argv if a in BLOCKS] or list(BLOCKS)
    exp = MCExperiment(
        name="MC01_synthetic",
        question=(
            "Do random systems satisfying each theorem's hypotheses obey its "
            "conclusion?"
        ),
        config={"blocks": names, "seed": SEED},
    )
    summary = {}
    for i, name in enumerate(BLOCKS):
        if name not in names:
            continue
        t0 = time.time()
        df, res = BLOCKS[name](np.random.default_rng(SEED + 101 * (i + 1)))
        df.to_csv(OUT / f"{name}.csv", index=False)
        res["elapsed_s"] = round(time.time() - t0, 1)
        summary[name] = res
        print(name, res, flush=True)
    prev = OUT / "MC01_summary.json"
    if prev.exists():
        import json

        old = json.loads(prev.read_text(encoding="utf-8"))
        old.update(summary)
        summary = old
    write_json(prev, summary)
    exp.finish("COMPUTED", **summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
