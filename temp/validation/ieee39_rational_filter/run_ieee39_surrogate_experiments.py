"""IEEE 39-derived experiments for the rational graph-filter damping paper.

This script runs two validation layers:

1. Full ANDES baseline:
   - load `ieee39_full.xlsx`,
   - run power flow and eigen-analysis,
   - report the critical oscillatory poles and model inventory.

2. IEEE39-derived second-order surrogate:
   - build the active-power stiffness Laplacian from the solved IEEE39 network,
   - Kron-reduce load buses to the ten generator buses,
   - assign SG/IBR-like inertia and damping scenarios,
   - compare diagonal, second-order, reduced-QEP, and full-QEP margins,
   - compute weak-node and weak-link/Braess-like diagnostic tables.

Important limitation
--------------------
The surrogate experiments are NOT full ANDES IBR simulations with REGCP/REGF
models. The included ANDES IEEE39 case is a synchronous-generator dynamic case
with no active REGCA/REGCP/REGF devices. These experiments are therefore
"IEEE39-derived surrogate" evidence for the math, plus a real ANDES baseline.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


def damping_ratio(s: complex) -> float:
    if abs(s) < 1e-12:
        return np.inf
    return float(-s.real / abs(s))


def qep_poles(M: np.ndarray, D: np.ndarray, L: np.ndarray) -> np.ndarray:
    """Solve M qddot + D qdot + L q = 0 via first-order companion form."""
    M = np.asarray(M, dtype=float)
    D = np.asarray(D, dtype=float)
    L = np.asarray(L, dtype=float)
    n = len(M)
    Minv = np.diag(1.0 / M)
    A = np.block(
        [
            [np.zeros((n, n)), np.eye(n)],
            [-Minv @ L, -Minv @ D],
        ]
    )
    return np.linalg.eigvals(A)


def oscillatory_margin(poles: Iterable[complex], imag_tol: float = 1e-7) -> tuple[float, complex]:
    candidates = [s for s in poles if abs(s.imag) > imag_tol and abs(s) > imag_tol and s.real < 1e-9]
    if not candidates:
        return np.inf, complex(np.nan, np.nan)
    ratios = np.array([damping_ratio(s) for s in candidates])
    idx = int(np.argmin(ratios))
    return float(ratios[idx]), complex(candidates[idx])


def spectral_estimates(
    M: np.ndarray,
    D: np.ndarray,
    L: np.ndarray,
    r: int = 6,
    eta_threshold: float = 0.25,
) -> dict:
    """Diagonal, second-order, reduced-QEP, and full-QEP estimates."""
    M = np.asarray(M, dtype=float)
    D = np.asarray(D, dtype=float)
    L = np.asarray(L, dtype=float)
    n = len(M)

    Minvh = np.diag(1.0 / np.sqrt(M))
    Lt = Minvh @ L @ Minvh
    Dt = Minvh @ D @ Minvh
    nu, Q = np.linalg.eigh((Lt + Lt.T) / 2)
    valid = np.where(nu > 1e-9)[0]
    Gamma = Q.T @ Dt @ Q
    delta = np.diag(Gamma)
    E = Gamma - np.diag(delta)

    zeta_modes = np.full(n, np.inf)
    zeta_modes[valid] = delta[valid] / (2 * np.sqrt(nu[valid]))
    c = int(np.argmin(zeta_modes))
    zeta_diag = float(zeta_modes[c])

    disc = delta[c] ** 2 - 4 * nu[c]
    if disc >= 0:
        s0 = complex((-delta[c] - np.sqrt(disc)) / 2, 0)
    else:
        s0 = complex(-delta[c] / 2, np.sqrt(4 * nu[c] - delta[c] ** 2) / 2)

    corr = 0.0 + 0.0j
    eta = 0.0
    singular_neighbor = False
    for ell in valid:
        ell = int(ell)
        if ell == c:
            continue
        pell = s0 * s0 + delta[ell] * s0 + nu[ell]
        if abs(pell) < 1e-10:
            singular_neighbor = True
            continue
        corr += E[c, ell] * E[ell, c] / pell
        eta = max(eta, abs(E[c, ell]) / (abs(np.sqrt(nu[ell]) - np.sqrt(nu[c])) + 1e-9))

    ds = (s0 * s0 / (2 * s0 + delta[c])) * corr
    s2 = s0 + ds
    zeta_second = damping_ratio(s2)

    # Damping-selected reduced QEP: critical mode + modes with lowest diagonal
    # damping and strongest coupling to critical mode.
    low = list(np.argsort(zeta_modes)[: min(r, n)])
    coupling = list(np.argsort(-np.abs(E[c, :]))[: min(r, n)])
    keep = []
    for k in [c] + low + coupling:
        if k in valid and k not in keep:
            keep.append(int(k))
        if len(keep) >= min(r, len(valid)):
            break
    Vr = Q[:, keep]
    Mr = np.eye(len(keep))
    Dr = Vr.T @ Dt @ Vr
    Lr = Vr.T @ Lt @ Vr
    rqep_poles = qep_poles(np.ones(len(keep)), Dr, Lr)
    zeta_rqep, pole_rqep = oscillatory_margin(rqep_poles)

    full_p = qep_poles(M, D, L)
    zeta_full, pole_full = oscillatory_margin(full_p)

    trigger_reject = bool(singular_neighbor or eta > eta_threshold or abs(ds) > 0.25 * max(abs(s0), 1e-12))
    adaptive = zeta_rqep if trigger_reject else zeta_second

    return {
        "zeta_full": zeta_full,
        "pole_full_real": pole_full.real,
        "pole_full_imag": pole_full.imag,
        "zeta_diag": zeta_diag,
        "zeta_second": float(zeta_second),
        "zeta_rqep": float(zeta_rqep),
        "pole_rqep_real": pole_rqep.real,
        "pole_rqep_imag": pole_rqep.imag,
        "zeta_adaptive": float(adaptive),
        "critical_mode": c,
        "nu_critical": float(nu[c]),
        "delta_critical": float(delta[c]),
        "eta": float(eta),
        "trigger_reject": trigger_reject,
        "kept_modes": keep,
    }


def load_andes_ieee39():
    import andes  # type: ignore

    case = andes.get_case("ieee39/ieee39_full.xlsx")
    ss = andes.load(case, setup=False, no_output=True)
    ss.setup()
    ss.PFlow.run()
    ss.EIG.run()
    return ss, case


def andes_inventory(ss) -> dict:
    names = [
        "Bus",
        "Line",
        "GENROU",
        "GENCLS",
        "TGOV1",
        "IEEEST",
        "PSS",
        "REGCA1",
        "REGCP1",
        "REGF1",
        "REGF2",
        "PLL1",
    ]
    return {name: int(getattr(getattr(ss, name), "n", 0)) for name in names if hasattr(ss, name)}


def andes_pole_summary(ss) -> dict:
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    stable_osc = [s for s in mu if abs(s.imag) > 1e-7 and s.real < 0]
    ratios = np.array([damping_ratio(s) for s in stable_osc])
    order = np.argsort(ratios)
    critical = stable_osc[int(order[0])] if len(order) else complex(np.nan, np.nan)
    top = []
    for idx in order[:10]:
        s = stable_osc[int(idx)]
        top.append(
            {
                "real": float(s.real),
                "imag": float(s.imag),
                "zeta": damping_ratio(s),
                "freq_hz": float(abs(s.imag) / (2 * np.pi)),
            }
        )
    return {
        "n_eigs": int(len(mu)),
        "n_positive_real": int(np.sum(mu.real > 1e-7)),
        "n_zero": int(np.sum(np.abs(mu) < 1e-7)),
        "n_oscillatory": int(np.sum(np.abs(mu.imag) > 1e-7)),
        "max_real": float(np.max(mu.real)) if len(mu) else np.nan,
        "critical_stable_osc": {
            "real": float(critical.real),
            "imag": float(critical.imag),
            "zeta": damping_ratio(critical),
            "freq_hz": float(abs(critical.imag) / (2 * np.pi)),
        },
        "least_damped_stable_oscillatory": top,
    }


def build_laplacian_from_andes(ss) -> tuple[np.ndarray, list[int], dict]:
    """Build lossless active-power stiffness Laplacian from solved ANDES lines."""
    bus_ids = [int(x) for x in ss.Bus.idx.v]
    bus_pos = {bus: i for i, bus in enumerate(bus_ids)}
    V = np.asarray(ss.Bus.v.v, dtype=float)
    theta = np.asarray(ss.Bus.a.v, dtype=float)
    n = len(bus_ids)
    L = np.zeros((n, n), dtype=float)

    line_rows = []
    for idx, (b1, b2, ghk, bhk) in enumerate(zip(ss.Line.bus1.v, ss.Line.bus2.v, ss.Line.ghk.v, ss.Line.bhk.v)):
        i = bus_pos[int(b1)]
        j = bus_pos[int(b2)]
        # Active-power Jacobian weight around solved operating point.
        # ANDES Line.b is shunt charging; ghk/bhk are the series branch
        # admittance terms. For mostly inductive lines, bhk is negative and
        # -bhk*cos(delta) gives the familiar positive synchronizing stiffness.
        angle = theta[i] - theta[j]
        w = float(V[i] * V[j] * ((-bhk) * np.cos(angle) + ghk * np.sin(angle)))
        if w <= 0:
            continue
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
        line_rows.append(
            {
                "line_index": idx,
                "from_bus": int(b1),
                "to_bus": int(b2),
                "weight": w,
            }
        )
    meta = {"bus_ids": bus_ids, "lines": line_rows}
    return L, bus_ids, meta


def kron_reduce(L: np.ndarray, keep: list[int]) -> np.ndarray:
    keep = list(keep)
    elim = [i for i in range(L.shape[0]) if i not in keep]
    Lkk = L[np.ix_(keep, keep)]
    if not elim:
        return Lkk
    Lke = L[np.ix_(keep, elim)]
    Lek = L[np.ix_(elim, keep)]
    Lee = L[np.ix_(elim, elim)]
    return Lkk - Lke @ np.linalg.pinv(Lee) @ Lek


def make_scenario(base_M: np.ndarray, rho: float, trial: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, dict]:
    """Create SG/IBR-like M,D vectors for the Kron-reduced generator network."""
    n = len(base_M)
    m = int(round(rho * n))
    converted = np.array([], dtype=int)
    if m > 0:
        converted = rng.choice(n, size=m, replace=False)

    M = base_M.copy()
    # SG baseline damping is mildly proportional. Converted IBRs have lower
    # effective inertia and heterogeneous control damping.
    damping_ratio_base = 0.45
    ratio = np.full(n, damping_ratio_base)
    if m > 0:
        M[converted] *= rng.uniform(0.25, 0.55, size=m)
        ratio[converted] = rng.lognormal(mean=np.log(0.55), sigma=0.55, size=m)
    D = np.diag(M * ratio)
    return M, D, {"rho": rho, "trial": trial, "converted": converted.tolist(), "damping_ratio": ratio.tolist()}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_figures(out: Path, estimator_rows: list[dict], line_rows: list[dict]) -> None:
    """Create lightweight reproducible figures when matplotlib is available."""
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    rhos = sorted({float(r["rho"]) for r in estimator_rows})
    methods = [
        ("err_diag", "Diagonal"),
        ("err_second", "Second order"),
        ("err_rqep", "Reduced QEP"),
        ("err_adaptive", "Adaptive"),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for key, label in methods:
        med = []
        p95 = []
        for rho in rhos:
            vals = np.array([float(r[key]) for r in estimator_rows if abs(float(r["rho"]) - rho) < 1e-12])
            med.append(np.median(vals) * 100)
            p95.append(np.percentile(vals, 95) * 100)
        ax.plot(rhos, med, marker="o", label=f"{label} median")
        ax.plot(rhos, p95, linestyle="--", alpha=0.55, label=f"{label} p95")
    ax.set_xlabel("IBR-like penetration in surrogate")
    ax.set_ylabel("Relative damping-margin error (%)")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "fig_estimator_errors_by_rho.pdf")
    fig.savefig(out / "fig_estimator_errors_by_rho.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    x = np.array([float(r["delta_nu_critical"]) for r in line_rows])
    y = np.array([float(r["delta_zeta"]) for r in line_rows])
    colors = ["#b22222" if bool(r["braess_like_reversal"]) else "#1f77b4" for r in line_rows]
    ax.scatter(x, y, c=colors, s=28, alpha=0.85)
    ax.axhline(0, color="0.2", linewidth=0.8)
    ax.axvline(0, color="0.2", linewidth=0.8)
    ax.set_xlabel(r"$\Delta \nu_c$ under 10% line reinforcement")
    ax.set_ylabel(r"$\Delta \zeta_{\min}$")
    ax.set_title("IEEE39-derived surrogate line reinforcements")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig_line_frequency_vs_damping.pdf")
    fig.savefig(out / "fig_line_frequency_vs_damping.png", dpi=200)
    plt.close(fig)


def run_surrogate_suite(ss, out: Path, trials: int = 60, seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    Lbus, bus_ids, net_meta = build_laplacian_from_andes(ss)

    gen_buses = [int(x) for x in ss.GENROU.bus.v]
    bus_pos = {bus: i for i, bus in enumerate(bus_ids)}
    keep = [bus_pos[b] for b in gen_buses]
    Lred = kron_reduce(Lbus, keep)

    # Remove numerical asymmetry and tiny row-sum drift.
    Lred = (Lred + Lred.T) / 2
    Lred -= np.diag(np.sum(Lred, axis=1))
    Lred = (Lred + Lred.T) / 2

    base_M = np.asarray(ss.GENROU.M.v, dtype=float)

    estimator_rows: list[dict] = []
    node_rows: list[dict] = []
    line_rows: list[dict] = []

    rhos = [0.0, 0.2, 0.4, 0.6, 0.8]
    for rho in rhos:
        for trial in range(trials):
            M, D, scenario = make_scenario(base_M, rho, trial, rng)
            est = spectral_estimates(M, D, Lred, r=6)
            row = {
                **{k: v for k, v in scenario.items() if k != "converted" and k != "damping_ratio"},
                "n_converted": len(scenario["converted"]),
                **{k: v for k, v in est.items() if k != "kept_modes"},
                "err_diag": abs(est["zeta_diag"] - est["zeta_full"]) / max(est["zeta_full"], 1e-12),
                "err_second": abs(est["zeta_second"] - est["zeta_full"]) / max(est["zeta_full"], 1e-12),
                "err_rqep": abs(est["zeta_rqep"] - est["zeta_full"]) / max(est["zeta_full"], 1e-12),
                "err_adaptive": abs(est["zeta_adaptive"] - est["zeta_full"]) / max(est["zeta_full"], 1e-12),
            }
            estimator_rows.append(row)

    # Weak-node and line/Braess diagnostics on one representative 60% scenario.
    M, D, scenario = make_scenario(base_M, 0.6, 0, np.random.default_rng(seed + 1000))
    base = spectral_estimates(M, D, Lred, r=6)
    zbase = base["zeta_full"]

    for i, bus in enumerate(gen_buses):
        D2 = D.copy()
        D2[i, i] *= 1.10
        e2 = spectral_estimates(M, D2, Lred, r=6)
        node_rows.append(
            {
                "bus": bus,
                "action": "increase_local_damping_10pct",
                "base_zeta": zbase,
                "post_zeta": e2["zeta_full"],
                "delta_zeta": e2["zeta_full"] - zbase,
                "sensitivity_fd": (e2["zeta_full"] - zbase) / max(D[i, i] * 0.10, 1e-12),
            }
        )

    # Rebuild every original line, perturb its susceptance weight, Kron reduce,
    # and evaluate damping vs frequency changes in the same surrogate scenario.
    for line in net_meta["lines"]:
        Lp = Lbus.copy()
        i = bus_pos[line["from_bus"]]
        j = bus_pos[line["to_bus"]]
        dw = 0.10 * line["weight"]
        Lp[i, i] += dw
        Lp[j, j] += dw
        Lp[i, j] -= dw
        Lp[j, i] -= dw
        Lpred = kron_reduce(Lp, keep)
        Lpred = (Lpred + Lpred.T) / 2
        est2 = spectral_estimates(M, D, Lpred, r=6)
        line_rows.append(
            {
                "line_index": line["line_index"],
                "from_bus": line["from_bus"],
                "to_bus": line["to_bus"],
                "base_weight": line["weight"],
                "base_zeta": zbase,
                "post_zeta": est2["zeta_full"],
                "delta_zeta": est2["zeta_full"] - zbase,
                "base_nu_critical": base["nu_critical"],
                "post_nu_critical": est2["nu_critical"],
                "delta_nu_critical": est2["nu_critical"] - base["nu_critical"],
                "braess_like_reversal": bool((est2["nu_critical"] - base["nu_critical"]) > 0 and (est2["zeta_full"] - zbase) < 0),
            }
        )

    write_csv(out / "estimator_table.csv", estimator_rows)
    write_csv(out / "weak_node_table.csv", node_rows)
    write_csv(out / "weak_link_table.csv", line_rows)
    with (out / "network_metadata.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "bus_ids": bus_ids,
                "generator_buses": gen_buses,
                "kron_reduced_laplacian": Lred.tolist(),
                "base_M": base_M.tolist(),
                "lines": net_meta["lines"],
            },
            f,
            indent=2,
        )

    write_figures(out, estimator_rows, line_rows)

    def summarize(rows: list[dict], key: str) -> dict:
        values = np.array([float(r[key]) for r in rows if np.isfinite(float(r[key]))])
        if len(values) == 0:
            return {"median": np.nan, "p95": np.nan, "mean": np.nan}
        return {
            "median": float(np.median(values)),
            "p95": float(np.percentile(values, 95)),
            "mean": float(np.mean(values)),
        }

    summary = {
        "n_trials": len(estimator_rows),
        "errors": {
            "diag": summarize(estimator_rows, "err_diag"),
            "second": summarize(estimator_rows, "err_second"),
            "rqep": summarize(estimator_rows, "err_rqep"),
            "adaptive": summarize(estimator_rows, "err_adaptive"),
        },
        "by_rho": {},
        "weak_node_best": sorted(node_rows, key=lambda r: r["delta_zeta"], reverse=True)[:5],
        "weak_link_best": sorted(line_rows, key=lambda r: r["delta_zeta"], reverse=True)[:5],
        "weak_link_worst": sorted(line_rows, key=lambda r: r["delta_zeta"])[:5],
        "braess_like_count": int(sum(1 for r in line_rows if r["braess_like_reversal"])),
        "n_lines_tested": len(line_rows),
        "representative_scenario": scenario,
    }
    for rho in rhos:
        subset = [r for r in estimator_rows if abs(float(r["rho"]) - rho) < 1e-12]
        summary["by_rho"][str(rho)] = {
            "diag": summarize(subset, "err_diag"),
            "second": summarize(subset, "err_second"),
            "rqep": summarize(subset, "err_rqep"),
            "adaptive": summarize(subset, "err_adaptive"),
            "trigger_reject_rate": float(np.mean([bool(r["trigger_reject"]) for r in subset])),
        }

    with (out / "surrogate_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("outputs") / "ieee39_rational_filter")
    parser.add_argument("--trials", type=int, default=60)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    ss, case = load_andes_ieee39()

    baseline = {
        "case": str(case),
        "inventory": andes_inventory(ss),
        "eigs": andes_pole_summary(ss),
        "limitation": (
            "The packaged ieee39_full.xlsx case has synchronous GENROU dynamics "
            "and no active REGCA1/REGCP1/REGF devices. IBR experiments below are "
            "IEEE39-derived surrogate experiments, not full ANDES IBR simulations."
        ),
    }
    with (args.out / "andes_baseline_summary.json").open("w", encoding="utf-8") as f:
        json.dump(baseline, f, indent=2)

    summary = run_surrogate_suite(ss, args.out, trials=args.trials)
    print(json.dumps({"baseline": baseline, "surrogate": summary}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
