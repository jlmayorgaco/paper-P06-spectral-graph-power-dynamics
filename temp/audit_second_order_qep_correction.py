from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


def connected_laplacian(n: int, rng: np.random.Generator, extra_p: float = 0.10):
    W = np.zeros((n, n))
    for i in range(1, n):
        j = int(rng.integers(0, i))
        w = float(np.exp(rng.uniform(math.log(0.4), math.log(4.0))))
        W[i, j] = W[j, i] = w
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j] == 0 and rng.random() < extra_p:
                w = float(np.exp(rng.uniform(math.log(0.2), math.log(6.0))))
                W[i, j] = W[j, i] = w
    return np.diag(W.sum(axis=1)) - W


def complete_laplacian(n: int, eps: float, rng: np.random.Generator):
    W = np.ones((n, n)) - np.eye(n)
    if eps:
        noise = rng.normal(0, eps, size=(n, n))
        noise = (noise + noise.T) / 2
        W = np.maximum(0.05, W * (1 + noise))
        np.fill_diagonal(W, 0)
    return np.diag(W.sum(axis=1)) - W


def modal_objects(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    Mh = np.diag(1 / np.sqrt(M))
    Lt = Mh @ L @ Mh
    Dt = np.diag(D / M)
    nu_all, Q_all = np.linalg.eigh(Lt)
    keep = nu_all > 1e-8
    nu = nu_all[keep]
    Q = Q_all[:, keep]
    Gamma = Q.T @ Dt @ Q
    return Lt, Dt, nu, Gamma


def damping_margin_from_poles(poles: np.ndarray):
    zetas = np.array(
        [
            -lam.real / abs(lam)
            for lam in poles
            if abs(lam.imag) > 1e-7 and abs(lam) > 1e-9 and lam.real < 1e-7
        ]
    )
    return float(np.min(zetas)) if len(zetas) else np.nan


def exact_margin(L: np.ndarray, M: np.ndarray, D: np.ndarray):
    n = len(M)
    A = np.block(
        [
            [np.zeros((n, n)), np.eye(n)],
            [-np.diag(1 / M) @ L, -np.diag(D / M)],
        ]
    )
    return damping_margin_from_poles(np.linalg.eigvals(A))


def diagonal_poles(nu: np.ndarray, delta: np.ndarray):
    poles = []
    modes = []
    for k, (vk, dk) in enumerate(zip(nu, delta)):
        disc = vk - 0.25 * dk * dk
        if disc <= 1e-12:
            continue
        wd = math.sqrt(disc)
        poles.extend([complex(-0.5 * dk, wd), complex(-0.5 * dk, -wd)])
        modes.extend([k, k])
    return np.array(poles, dtype=complex), np.array(modes, dtype=int)


def diagonal_margin(nu: np.ndarray, Gamma: np.ndarray):
    delta = np.diag(Gamma)
    poles, _ = diagonal_poles(nu, delta)
    return damping_margin_from_poles(poles)


def qep_margin_from_mats(Lr: np.ndarray, Dr: np.ndarray):
    r = Lr.shape[0]
    A = np.block([[np.zeros((r, r)), np.eye(r)], [-Lr, -Dr]])
    return damping_margin_from_poles(np.linalg.eigvals(A))


def full_modal_qep_margin(nu: np.ndarray, Gamma: np.ndarray):
    return qep_margin_from_mats(np.diag(nu), Gamma)


def reduced_qep_margin(nu: np.ndarray, Gamma: np.ndarray, r: int):
    zdiag = np.diag(Gamma) / (2 * np.sqrt(nu))
    idx = np.argsort(zdiag)[: min(r, len(nu))]
    return qep_margin_from_mats(np.diag(nu[idx]), Gamma[np.ix_(idx, idx)])


def second_order_corrected_poles(
    nu: np.ndarray,
    Gamma: np.ndarray,
    *,
    use_conjugate: bool = False,
    sign: float = 1.0,
):
    delta = np.diag(Gamma)
    E = Gamma - np.diag(delta)
    poles0, modes = diagonal_poles(nu, delta)
    corrected = []
    valid = []
    for s0, k in zip(poles0, modes):
        pc_prime = 2 * s0 + delta[k]
        if abs(pc_prime) <= 1e-12:
            valid.append(False)
            corrected.append(s0)
            continue
        acc = 0.0 + 0.0j
        ok = True
        for ell in range(len(nu)):
            if ell == k:
                continue
            p_ell = s0 * s0 + delta[ell] * s0 + nu[ell]
            if abs(p_ell) <= 1e-10:
                ok = False
                break
            product = E[k, ell] * E[ell, k]
            if use_conjugate:
                product = E[k, ell] * np.conjugate(E[ell, k])
            acc += product / p_ell
        ds = sign * (s0 * s0 / pc_prime) * acc
        corrected.append(s0 + ds if ok else s0)
        valid.append(ok and np.isfinite(ds.real) and np.isfinite(ds.imag))
    return np.array(corrected, dtype=complex), np.array(valid, dtype=bool)


def second_order_margin(nu: np.ndarray, Gamma: np.ndarray, **kwargs):
    poles, valid = second_order_corrected_poles(nu, Gamma, **kwargs)
    if not np.all(valid):
        return np.nan
    return damping_margin_from_poles(poles)


def make_case(name: str, rng: np.random.Generator):
    if name == "proportional":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        D = 0.12 * M
        return L, M, D
    if name == "weak_ibr":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        ratio = np.exp(rng.uniform(math.log(0.006), math.log(0.08), n))
        D = M * ratio
        return L, M, D
    if name == "mild_ibr":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        ratio = np.exp(rng.uniform(math.log(0.08), math.log(0.8), n))
        D = M * ratio
        return L, M, D
    if name == "strong_ibr":
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(0.25), math.log(12.0), n))
        ratio = np.exp(rng.uniform(math.log(0.02), math.log(3.0), n))
        D = M * ratio
        return L, M, D
    if name == "near_degenerate":
        n = 24
        L = complete_laplacian(n, eps=0.015, rng=rng)
        M = np.ones(n)
        ratio = np.exp(rng.uniform(math.log(0.03), math.log(1.5), n))
        D = M * ratio
        return L, M, D
    if name.startswith("rho_"):
        rho = float(name.split("_", 1)[1])
        n = 39
        L = connected_laplacian(n, rng, extra_p=0.08)
        M = np.exp(rng.uniform(math.log(2.0), math.log(10.0), n))
        base = 0.12
        hetero = np.exp(rng.uniform(math.log(0.01), math.log(1.2), n))
        ratio = (1 - rho) * base + rho * hetero
        D = M * ratio
        return L, M, D
    raise ValueError(name)


def rel_error(est: float, truth: float):
    if not np.isfinite(est) or not np.isfinite(truth) or truth <= 0:
        return np.nan
    return abs(est - truth) / truth


@dataclass
class Row:
    family: str
    truth: float
    diag: float
    second_order: float
    second_order_neg: float
    second_order_conj: float
    qep6: float
    qep20: float
    diag_err: float
    second_order_err: float
    second_order_neg_err: float
    second_order_conj_err: float
    qep6_err: float
    qep20_err: float


def evaluate_family(name: str, cases: int, rng: np.random.Generator):
    rows: list[Row] = []
    for _ in range(cases):
        L, M, D = make_case(name, rng)
        _Lt, _Dt, nu, Gamma = modal_objects(L, M, D)
        truth = exact_margin(L, M, D)
        if not np.isfinite(truth) or truth <= 1e-10:
            continue
        diag = diagonal_margin(nu, Gamma)
        so = second_order_margin(nu, Gamma, sign=1.0, use_conjugate=False)
        so_neg = second_order_margin(nu, Gamma, sign=-1.0, use_conjugate=False)
        so_conj = second_order_margin(nu, Gamma, sign=1.0, use_conjugate=True)
        q6 = reduced_qep_margin(nu, Gamma, r=6)
        q20 = reduced_qep_margin(nu, Gamma, r=20)
        rows.append(
            Row(
                family=name,
                truth=truth,
                diag=diag,
                second_order=so,
                second_order_neg=so_neg,
                second_order_conj=so_conj,
                qep6=q6,
                qep20=q20,
                diag_err=rel_error(diag, truth),
                second_order_err=rel_error(so, truth),
                second_order_neg_err=rel_error(so_neg, truth),
                second_order_conj_err=rel_error(so_conj, truth),
                qep6_err=rel_error(q6, truth),
                qep20_err=rel_error(q20, truth),
            )
        )
    return rows


def summarize(rows: list[Row]):
    def stats(field: str):
        x = np.array([getattr(r, field) for r in rows], dtype=float)
        x = x[np.isfinite(x)]
        return {
            "median": float(np.median(x)) if len(x) else math.nan,
            "p95": float(np.quantile(x, 0.95)) if len(x) else math.nan,
            "finite_fraction": float(len(x) / len(rows)) if rows else math.nan,
        }

    return {
        "cases": len(rows),
        "diag": stats("diag_err"),
        "second_order": stats("second_order_err"),
        "second_order_neg_sign": stats("second_order_neg_err"),
        "second_order_conjugate": stats("second_order_conj_err"),
        "qep6": stats("qep6_err"),
        "qep20": stats("qep20_err"),
        "so_beats_diag_fraction": float(
            np.nanmean(
                np.array([r.second_order_err < r.diag_err for r in rows], dtype=float)
            )
        )
        if rows
        else math.nan,
        "so_beats_qep20_fraction": float(
            np.nanmean(
                np.array([r.second_order_err < r.qep20_err for r in rows], dtype=float)
            )
        )
        if rows
        else math.nan,
    }


def main():
    rng = np.random.default_rng(20260608)
    families = [
        "proportional",
        "weak_ibr",
        "mild_ibr",
        "strong_ibr",
        "near_degenerate",
        "rho_0.20",
        "rho_0.50",
        "rho_0.80",
    ]
    report = {}
    samples = {}
    for fam in families:
        rows = evaluate_family(fam, 300, rng)
        report[fam] = summarize(rows)
        samples[fam] = [asdict(r) for r in rows[:5]]
    out = {"families": report, "samples": samples}
    out_path = Path("temp") / "second_order_qep_correction_report.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
