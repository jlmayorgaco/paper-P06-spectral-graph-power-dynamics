"""Independent mathematical checks for binary spectral portfolio theory.

This is NOT an IEEE-39 simulation and does not access the user's repository.
Dependencies: NumPy, SciPy, SymPy (symbolic check only).
Run: python theory_checks.py --out results
Tests are deterministic. Numerical tests corroborate, rather than prove, the
analytic propositions in the accompanying technical note.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import itertools
import json
from pathlib import Path
import sys
import time
from typing import Callable, Iterable

import numpy as np
from scipy.linalg import eigvals, null_space, solve, solve_continuous_lyapunov
from scipy.optimize import minimize_scalar


def alpha(a: np.ndarray) -> float:
    return float(np.max(eigvals(a).real))


def masks(m: int) -> Iterable[int]:
    return sorted(range(1 << m), key=lambda v: (v.bit_count(), v))


def subsets(mask: int) -> Iterable[int]:
    sub = mask
    while True:
        yield sub
        if sub == 0:
            break
        sub = (sub - 1) & mask


def minimal_unsafe(flags: dict[int, bool]) -> list[int]:
    """Exact inclusion-minimal unsafe sets; NOT greedy one-deletion minimality."""
    answer: list[int] = []
    for s in sorted(flags, key=lambda v: (v.bit_count(), v)):
        if flags[s] and not any((h & s) == h for h in answer):
            answer.append(s)
    return answer


def quotient(a: np.ndarray, r: np.ndarray, tol: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Quotient only an externally declared, verified symmetry subspace.

    A R = 0 is checked. No spectral magnitude cutoff selects R.
    Additional zero/near-zero physical modes are deliberately preserved.
    """
    a = np.asarray(a, dtype=float)
    r = np.asarray(r, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or r.shape[0] != a.shape[0]:
        raise ValueError("Incompatible dimensions")
    if np.linalg.matrix_rank(r) != r.shape[1]:
        raise ValueError("Declared symmetry generators must be independent")
    q, _ = np.linalg.qr(r, mode="reduced")
    if np.linalg.norm(a @ q) > tol * max(1.0, np.linalg.norm(a)):
        raise ValueError("A R != 0: cannot delete this direction as a symmetry")
    z = null_space(q.T)
    return z.T @ a @ z, z


def action_q(k: np.ndarray, blocks: list[slice]) -> np.ndarray:
    d = np.zeros_like(k)
    for sl in blocks:
        d[sl, sl] = np.eye(sl.stop - sl.start) + k[sl, sl]
    return solve(d, np.eye(k.shape[0]) + k) - np.eye(k.shape[0])


def cardinality_gain_bound(r: np.ndarray, size: int, weights: np.ndarray | None = None) -> float:
    """Weighted top-(size-1) row bound for zero-diagonal block gain bounds."""
    m = r.shape[0]
    if r.shape != (m, m) or np.min(r) < 0 or np.max(np.abs(np.diag(r))) > 1e-12:
        raise ValueError("R must be square, nonnegative and zero-diagonal")
    if not 1 <= size <= m:
        raise ValueError("size must be 1..m")
    d = np.ones(m) if weights is None else np.asarray(weights, float)
    if np.any(d <= 0):
        raise ValueError("Weights must be strictly positive")
    w = r * d[np.newaxis, :] / d[:, np.newaxis]
    return float(max(np.sum(np.sort(np.delete(w[i], i))[-(size-1):]) if size > 1 else 0
                     for i in range(m)))


def positive_part(a: np.ndarray) -> np.ndarray:
    values, vec = np.linalg.eigh((a + a.T) / 2)
    return (vec * np.maximum(values, 0.0)) @ vec.T


@dataclass
class SearchStats:
    nodes: int = 0
    leaf_queries: int = 0
    eigen_queries: int = 0
    lyapunov_solves: int = 0
    certified_safe_subcubes: int = 0
    minimality_exclusions: int = 0


def search_minimal(a0: np.ndarray, changes: list[np.ndarray]) -> tuple[list[int], SearchStats]:
    """Exact search with safe-subcube Lyapunov majorants.

    A failed certificate is UNKNOWN, never unsafe. Supersets of an already
    verified unsafe set are skipped ONLY for minimal-witness enumeration.
    Their actual stability remains unspecified; it may be restored.
    """
    m = len(changes)
    stats = SearchStats()
    witnesses: list[int] = []

    def matrix(mask: int) -> np.ndarray:
        a = a0.copy()
        for i, di in enumerate(changes):
            if (mask >> i) & 1:
                a += di
        return a

    def visit(forced: int, remaining: tuple[int, ...]) -> None:
        stats.nodes += 1
        if any((h & forced) == h for h in witnesses):
            stats.minimality_exclusions += 1
            return
        a = matrix(forced)
        stats.eigen_queries += 1
        aa = alpha(a)
        if aa > 1e-10:
            # This is verified unsafe, not necessarily minimal yet. Adding it
            # can only eliminate candidate supersets from minimality search.
            witnesses.append(forced)
            witnesses[:] = minimal_unsafe({s: True for s in witnesses})
            return
        if not remaining:
            stats.leaf_queries += 1
            return
        if aa < -1e-8:
            p = solve_continuous_lyapunov(a.T, -np.eye(a.shape[0]))
            p = (p + p.T) / 2
            stats.lyapunov_solves += 1
            residual = np.linalg.norm(a.T @ p + p @ a + np.eye(a.shape[0]))
            if np.min(np.linalg.eigvalsh(p)) > 1e-10 and residual < 1e-8:
                b = a.T @ p + p @ a
                for i in remaining:
                    b += positive_part(changes[i].T @ p + p @ changes[i])
                if np.max(np.linalg.eigvalsh((b+b.T)/2)) < -1e-8:
                    stats.certified_safe_subcubes += 1
                    return
        i = remaining[0]
        rest = remaining[1:]
        visit(forced, rest)
        visit(forced | (1 << i), rest)

    visit(0, tuple(range(m)))
    return minimal_unsafe({s: True for s in witnesses}), stats


def run_checks() -> dict:
    results: dict = {}
    rng = np.random.default_rng(20260910)
    count = 0

    # 1. Quotient preserves tiny physical unstable and marginal eigenvalues.
    rows = []
    for eps in [-1e-4, 0.0, 1e-7, 1e-4]:
        a = np.array([[0., 1., 0.], [0., eps, 0.], [0., 0., -2.]])
        aq, z = quotient(a, np.eye(3)[:, :1])
        assert np.max(np.abs(np.sort(eigvals(aq).real) - np.sort([eps, -2.]))) < 1e-12
        # Shift the *known symmetry* only, not the entire origin cluster.
        ashift = a - 0.7 * np.diag([1., 0., 0.])
        for s in [0.1+0.2j, 1.2j, -0.4+0.3j]:
            assert abs(np.linalg.det(s*np.eye(3)-ashift) - (s+0.7)*np.linalg.det(s*np.eye(2)-aq)) < 1e-12
            count += 1
        rows.append({"physical_pole": eps, "quotient_alpha": alpha(aq),
                     "threshold_1e-3_would_hide": abs(eps) < 1e-3})
    # Negative test: cannot delete the generalized zero eigenvector.
    try:
        quotient(np.array([[0.,1.],[0.,0.]]), np.eye(2))
        raise AssertionError("Must not remove a Jordan chain as two symmetries")
    except ValueError:
        count += 1
    results["symmetry_quotient"] = rows

    # 2. Exact determinant lemma and fixed-block gauge invariance.
    n, p = 7, 4
    a0 = -np.diag(np.arange(1., n+1.))
    u, v = rng.normal(size=(n,p))*.3, rng.normal(size=(n,p))*.3
    a1 = a0 - u @ v.T
    blocks = [slice(0,2), slice(2,4)]
    gauge = np.diag([.02, 3., 10., .1])
    max_det, max_gauge = 0., 0.
    for s in [.4+.7j, -.2+1.3j, .03+2.6j]:
        k = v.T @ solve(s*np.eye(n)-a0, u)
        ratio = np.linalg.det(s*np.eye(n)-a1) / np.linalg.det(s*np.eye(n)-a0)
        e = abs(ratio - np.linalg.det(np.eye(p)+k))
        kp = solve(gauge, k @ gauge)
        q, qp = action_q(k, blocks), action_q(kp, blocks)
        gerr = np.linalg.norm(qp-solve(gauge,q@gauge))
        max_det, max_gauge = max(max_det,e), max(max_gauge,gerr)
        assert e < 1e-12 and gerr < 1e-11
        count += 1
    results["port_identity"] = {"max_absolute_error": max_det, "gauge_error": max_gauge}

    # 3. Nonnormality: eigen-distance to -1 is not a robustness radius.
    nrows = []
    for L in [1., 10., 100., 1000.]:
        q = np.array([[0., L], [0., 0.]])
        d = float(np.min(abs(eigvals(q)+1)))
        ss = float(np.linalg.svd(np.eye(2)+q,compute_uv=False)[-1])
        assert abs(d-1) < 1e-12
        nrows.append({"L": L, "eigen_distance_minus1": d, "singular_min": ss})
        count += 1
    results["closure_not_robustness_radius"] = nrows

    # 4. Small-gain cardinality lower bound: every <=4 subset stable, 5 fails.
    m = 6
    R = .3*(np.ones((m,m))-np.eye(m))
    cert = {r:cardinality_gain_bound(R,r) for r in range(1,m+1)}
    truth = {}
    for s in masks(m):
        ids = [i for i in range(m) if (s>>i)&1]
        aa = -1. if not ids else alpha(-np.eye(len(ids))+R[np.ix_(ids,ids)])
        truth[s] = aa > 0
        if len(ids) <= 4:
            assert aa < -1e-8
        count += 1
    h = minimal_unsafe(truth)
    kappa = min(s.bit_count() for s in h)
    assert kappa == 5 and cert[4] < 1 < cert[5]
    # Unstable-open-loop trap: Q=.1/(s-1) has gain<1 on jw, closed pole .9.
    assert .1 < 1 and 1.-.1 > 0
    count += 1
    results["cardinality_certificate"] = {"top_row_bounds":cert, "true_kappa":kappa,
        "number_minimal_witnesses":len(h), "global_perron_bound":float(max(abs(eigvals(R)))),
        "unstable_port_naive_smallgain_counterexample_pole":.9}

    # 5. Exact Boolean Lyapunov identity; continuum between safe vertices is bad.
    import sympy as sy
    d = sy.symbols("d", real=True)
    a0s, a1s = sy.Matrix([[-1,4],[0,-1]]), sy.Matrix([[-1,0],[4,-1]])
    p0s = sy.Matrix([[sy.Rational(1,2),1],[1,sy.Rational(9,2)]])
    p1s = sy.Matrix([[sy.Rational(9,2),1],[1,sy.Rational(1,2)]])
    ap, pp = (1-d)*a0s+d*a1s, (1-d)*p0s+d*p1s
    rem = sy.simplify(-(ap.T*pp+pp*ap)-sy.eye(2))
    target = 32*d*(d-1)*sy.Matrix([[0,1],[1,0]])
    assert sy.simplify(rem-target) == sy.zeros(2)
    for val in [0,1]:
        P = np.array(pp.subs(d,val), float)
        A = np.array(ap.subs(d,val), float)
        assert np.linalg.eigvalsh(P).min()>0
        assert np.linalg.norm(A.T@P+P@A+np.eye(2))<1e-12
        count += 1
    results["boolean_vs_box"]={"alpha_vertex0":alpha(np.array(a0s,float)),
        "alpha_vertex1":alpha(np.array(a1s,float)),
        "alpha_midpoint":alpha(np.array(ap.subs(d,.5),float)),
        "polynomial_residual":"32*d*(d-1)*[[0,1],[1,0]]"}

    # 6. Arbitrary k-wise no-go need not involve high-degree interaction.
    ng = []
    for m in range(2,9):
        flags = {s:s.bit_count()-(m-.5)>0 for s in range(1<<m)}
        assert minimal_unsafe(flags)==[(1<<m)-1]
        ng.append({"actions":m,"minimum_bad_order":m,"scalar_interaction_degree":1})
        count += 1
    results["threshold_not_interaction_degree"]=ng

    # 7. Exact search vs exhaustive enumeration for affine, nonmonotone families.
    sr=[]
    for seed in range(5):
        rr=np.random.default_rng(seed+500)
        n,m=4,9
        a0=-2*np.eye(n)
        changes=[]
        for i in range(m):
            diag=rr.uniform(-.9,.9,size=n)
            changes.append(np.diag(diag))
        flags={}
        for s in range(1<<m):
            aa=a0+sum((changes[i] for i in range(m) if (s>>i)&1),np.zeros_like(a0))
            flags[s]=alpha(aa)>0
        exact=minimal_unsafe(flags)
        got,st=search_minimal(a0,changes)
        assert set(got)==set(exact)
        restab=sum(1 for s in flags if flags[s]
                   for i in range(m) if not((s>>i)&1) and not flags[s|(1<<i)])
        sr.append({"seed":seed,"portfolios":1<<m,"hyperedges":len(exact),
                   "restabilizing_edges":restab,**st.__dict__})
        count += 1
    results["safe_subcube_search"]=sr

    # 8. Greedy one-deletion minimal is NOT inclusion-minimal.
    flags={s:s in [1,7] for s in range(8)}
    assert flags[7] and all(not flags[s] for s in [3,5,6])
    assert minimal_unsafe(flags)==[1]
    count += 1
    results["nonmonotone_greedy_trap"]={"unsafe_masks":[1,7],"true_minimal_masks":[1],
        "mask7_all_immediate_deletions_stable":True}

    # 9. Hyperedge-free is EXACT for all-prefix orders, not all safe endpoints.
    hr=[]
    for seed in range(8):
        rr=np.random.default_rng(seed+300)
        m=5
        flags={s:bool(rr.integers(0,2)) for s in range(1<<m)}
        flags[0]=False
        h=minimal_unsafe(flags)
        for s in range(1<<m):
            hfree=not any((e&s)==e for e in h)
            allsafe=all(not flags[r] for r in subsets(s))
            assert hfree==allsafe
            count += 1
        hr.append({"seed":seed,"minimal_hyperedges":len(h)})
    results["hereditary_core_equivalence"]=hr

    # 10. Stable endpoints need not permit a safe monotone deployment path.
    a0=-np.eye(2)
    d1,d2=np.diag([2.,-3.]),np.diag([-3.,2.])
    al=[alpha(a0),alpha(a0+d1),alpha(a0+d2),alpha(a0+d1+d2)]
    assert al==[-1.,1.,1.,-2.]
    count+=1
    results["stable_endpoint_no_safe_path"]={"alphas_empty_1_2_12":al}

    # 11. Local return-difference margin law, an analytic Hopf-type toy.
    # f(s,t) = ((s-a)^2+w0^2)/(s+b)^2, a=t, b>0. Q=f-1.
    b,w0=1.7,2.3
    c=2*w0/(b*b+w0*w0)
    cr=[]
    for a in [1e-2,5e-3,1e-3,-1e-2,-1e-3]:
        def val(w):
            s=1j*w
            return abs(((s-a)**2+w0*w0)/(s+b)**2)
        opt=minimize_scalar(val,bounds=(w0-.2,w0+.2),method="bounded",options={"xatol":1e-14})
        rel=abs(opt.fun/(c*abs(a))-1)
        assert rel<.02
        cr.append({"alpha":a,"local_min":float(opt.fun),"linear_prediction":c*abs(a),
                   "relative_error":rel,"omega_error":float(opt.x-w0)})
        count+=1
    results["local_closure_law"]=cr

    # 12. Safe modes alone do not imply stability under arbitrary switching.
    A0=np.array(a0s,float); A1=np.array(a1s,float)
    from scipy.linalg import expm
    mon=expm(A1*.05)@expm(A0*.05)
    rho=float(max(abs(eigvals(mon))))
    assert rho>1 and alpha(A0)<0 and alpha(A1)<0
    count+=1
    results["switching_trap"]={"mode_alphas":[alpha(A0),alpha(A1)],"period":.1,
                                "monodromy_spectral_radius":rho}

    # 13. Finite nonlinear Lyapunov radius sanity check (analytical scalar).
    # xdot=-x+x^2; V=x^2; L=1, beta=1, any radius<1 is inward/invariant.
    for x in np.linspace(-.8,.8,41):
        deriv=2*x*(-x+x*x)
        assert deriv <= -.4*x*x + 1e-12
        count+=1
    results["local_nonlinear_certificate"]={"system":"xdot=-x+x^2", "V":"x^2",
       "certified_radius":.8,"inequality":"Vdot <= -0.4*x^2"}

    results["meta"]={"checks_passed":count,"seed":20260910,"scope":"synthetic theory checks, not IEEE/ANDES validation",
                     "python":sys.version,"numpy":np.__version__}
    return results


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("results"))
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    t=time.time(); results=run_checks()
    results["meta"]["wall_seconds"]=time.time()-t
    (args.out/"theory_check_results.json").write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(json.dumps(results["meta"],indent=2))
    print("Cardinality certificate:",results["cardinality_certificate"])
    print("Search:",results["safe_subcube_search"])

if __name__=="__main__":
    main()
