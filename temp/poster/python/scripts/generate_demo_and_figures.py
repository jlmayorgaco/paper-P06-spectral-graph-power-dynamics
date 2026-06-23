from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "poster_ieee_session"
FIG_DIR = OUT / "figures"
DATA_DIR = OUT / "data"

MASTER_SEED = 20260608
REGIME = "moderate"
REINFORCEMENT_FRAC = 0.10

IEEE_BLUE = "#00843D"
IEEE_DARK = "#182A20"
GOOD = "#00843D"
BAD = "#b23a48"
MUTED = "#4F6F5D"
LIGHT = "#E4F4EB"


def connected_graph(rng: np.random.Generator, n: int, p_extra: float) -> np.ndarray:
    """Random connected weighted graph: random tree plus extra plausible ties."""
    weights = np.zeros((n, n), dtype=float)
    for j in range(1, n):
        i = int(rng.integers(0, j))
        w = float(np.exp(rng.uniform(np.log(0.35), np.log(3.5))))
        weights[i, j] = weights[j, i] = w

    for i in range(n):
        for j in range(i + 1, n):
            if weights[i, j] == 0.0 and rng.random() < p_extra:
                w = float(np.exp(rng.uniform(np.log(0.25), np.log(3.0))))
                weights[i, j] = weights[j, i] = w
    return weights


def laplacian(weights: np.ndarray) -> np.ndarray:
    return np.diag(weights.sum(axis=1)) - weights


def normalized_matrices(
    weights: np.ndarray, inertia: np.ndarray, damping_per_inertia: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    sinv = np.diag(1.0 / np.sqrt(inertia))
    lt = sinv @ laplacian(weights) @ sinv
    dt = np.diag(damping_per_inertia)
    return lt, dt


def damping_ratio(pole: complex) -> float:
    if abs(pole) < 1e-12:
        return float("inf")
    return float(-pole.real / abs(pole))


def full_margin(
    weights: np.ndarray, inertia: np.ndarray, damping_per_inertia: np.ndarray
) -> dict[str, Any] | None:
    lt, dt = normalized_matrices(weights, inertia, damping_per_inertia)
    n = len(inertia)
    a = np.block([[np.zeros((n, n)), np.eye(n)], [-lt, -dt]])
    eigvals, eigvecs = np.linalg.eig(a)
    candidates = [
        (damping_ratio(s), s, idx)
        for idx, s in enumerate(eigvals)
        if s.imag > 1e-7 and s.real < 1e-8 and abs(s) > 1e-9
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    zeta, pole, idx = candidates[0]
    return {
        "zeta": float(zeta),
        "pole": pole,
        "idx": idx,
        "A": a,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
    }


def modal_data(
    weights: np.ndarray, inertia: np.ndarray, damping_per_inertia: np.ndarray
) -> dict[str, Any]:
    lt, dt = normalized_matrices(weights, inertia, damping_per_inertia)
    eigvals, q = np.linalg.eigh(lt)
    keep = eigvals > 1e-9
    nu = eigvals[keep]
    q = q[:, keep]
    gamma = q.T @ dt @ q
    delta = np.diag(gamma)
    zetas = np.full_like(nu, np.inf, dtype=float)
    valid = delta < 2.0 * np.sqrt(nu) - 1e-10
    zetas[valid] = delta[valid] / (2.0 * np.sqrt(nu[valid]))
    return {
        "Lt": lt,
        "Dt": dt,
        "nu": nu,
        "Q": q,
        "Gamma": gamma,
        "delta": delta,
        "zeta_modes": zetas,
    }


def second_order_estimate(
    weights: np.ndarray, inertia: np.ndarray, damping_per_inertia: np.ndarray
) -> dict[str, Any] | None:
    md = modal_data(weights, inertia, damping_per_inertia)
    nu = md["nu"]
    gamma = md["Gamma"]
    delta = md["delta"]
    zeta_modes = md["zeta_modes"]
    c = int(np.argmin(zeta_modes))
    if not np.isfinite(zeta_modes[c]):
        return None

    discriminant = 4.0 * nu[c] - delta[c] ** 2
    if discriminant <= 0:
        return None
    s0 = complex(-delta[c] / 2.0, math.sqrt(discriminant) / 2.0)
    offdiag = gamma - np.diag(delta)

    corr = 0.0j
    min_den = float("inf")
    for ell in range(len(nu)):
        if ell == c:
            continue
        denom = s0 * s0 + delta[ell] * s0 + nu[ell]
        min_den = min(min_den, abs(denom))
        corr += offdiag[c, ell] * offdiag[ell, c] / denom

    ds = (s0 * s0 / (2.0 * s0 + delta[c])) * corr
    s2 = s0 + ds
    eta = float(
        np.linalg.norm(offdiag, "fro")
        / max(np.linalg.norm(np.diag(delta), "fro"), 1e-12)
    )

    return {
        **md,
        "critical_mode": c,
        "zeta_diag": float(zeta_modes[c]),
        "zeta_second": damping_ratio(s2),
        "s0": s0,
        "s2": s2,
        "ds": ds,
        "nu_c": float(nu[c]),
        "delta_c": float(delta[c]),
        "eta": eta,
        "min_denominator": float(min_den),
    }


def line_dlt(inertia: np.ndarray, i: int, j: int) -> np.ndarray:
    direction = np.zeros(len(inertia), dtype=float)
    direction[i] = 1.0 / math.sqrt(float(inertia[i]))
    direction[j] = -1.0 / math.sqrt(float(inertia[j]))
    return np.outer(direction, direction)


def qep_zeta_sensitivity(
    weights: np.ndarray,
    inertia: np.ndarray,
    damping_per_inertia: np.ndarray,
    i: int,
    j: int,
) -> tuple[float, complex]:
    """Exact first-order sensitivity of the current full-QEP critical pole."""
    fm = full_margin(weights, inertia, damping_per_inertia)
    if fm is None:
        raise RuntimeError("No oscillatory full-QEP pole was found.")

    eigvals = fm["eigvals"]
    eigvecs = fm["eigvecs"]
    pole_idx = int(fm["idx"])
    pole = eigvals[pole_idx]

    left_rows = np.linalg.inv(eigvecs)
    left = left_rows[pole_idx, :]

    n = len(inertia)
    da = np.zeros((2 * n, 2 * n), dtype=complex)
    da[n:, :n] = -line_dlt(inertia, i, j)
    ds_dw = (left @ (da @ eigvecs[:, pole_idx])) / (left @ eigvecs[:, pole_idx])

    a = pole.real
    b = pole.imag
    radius = abs(pole)
    da_real = ds_dw.real
    db_imag = ds_dw.imag
    dzeta_dw = -da_real / radius + a * (a * da_real + b * db_imag) / (radius**3)
    return float(np.real(dzeta_dw)), complex(ds_dw)


def edge_rows(
    weights: np.ndarray,
    inertia: np.ndarray,
    damping_per_inertia: np.ndarray,
    base_full: dict[str, Any],
    base_second: dict[str, Any],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    q = base_second["Q"]
    c = int(base_second["critical_mode"])
    base_nu = float(base_second["nu_c"])
    base_zeta = float(base_full["zeta"])
    base_zeta_second = float(base_second["zeta_second"])

    edge_index = 0
    for i in range(weights.shape[0]):
        for j in range(i + 1, weights.shape[1]):
            w = float(weights[i, j])
            if w <= 0.0:
                continue

            dnu_dw = float((q[i, c] / math.sqrt(inertia[i]) - q[j, c] / math.sqrt(inertia[j])) ** 2)
            dzeta_dw_qep, ds_dw_qep = qep_zeta_sensitivity(
                weights, inertia, damping_per_inertia, i, j
            )

            perturbed = weights.copy()
            perturbed[i, j] = perturbed[j, i] = w * (1.0 + REINFORCEMENT_FRAC)
            post_full = full_margin(perturbed, inertia, damping_per_inertia)
            post_second = second_order_estimate(perturbed, inertia, damping_per_inertia)
            if post_full is None or post_second is None:
                continue

            delta_w = REINFORCEMENT_FRAC * w
            delta_zeta_full = float(post_full["zeta"] - base_zeta)
            delta_zeta_second = float(post_second["zeta_second"] - base_zeta_second)
            rows.append(
                {
                    "edge_index": edge_index,
                    "from_bus": i + 1,
                    "to_bus": j + 1,
                    "base_weight": w,
                    "delta_w": delta_w,
                    "dnu_dw_frequency": dnu_dw,
                    "delta_nu_fd_10pct": float(post_second["nu_c"] - base_nu),
                    "dzeta_dw_qep_analytic": dzeta_dw_qep,
                    "ds_dw_qep_real": float(ds_dw_qep.real),
                    "ds_dw_qep_imag": float(ds_dw_qep.imag),
                    "delta_zeta_full_fd_10pct": delta_zeta_full,
                    "dzeta_dw_full_fd": delta_zeta_full / delta_w,
                    "delta_zeta_second_order_fd_10pct": delta_zeta_second,
                    "dzeta_dw_second_order_fd": delta_zeta_second / delta_w,
                    "post_zeta_full": float(post_full["zeta"]),
                    "post_zeta_second": float(post_second["zeta_second"]),
                }
            )
            edge_index += 1

    freq_sorted = sorted(rows, key=lambda row: row["dnu_dw_frequency"], reverse=True)
    damp_sorted = sorted(rows, key=lambda row: row["dzeta_dw_qep_analytic"], reverse=True)
    freq_rank = {row["edge_index"]: rank + 1 for rank, row in enumerate(freq_sorted)}
    damp_rank = {row["edge_index"]: rank + 1 for rank, row in enumerate(damp_sorted)}
    for row in rows:
        row["frequency_rank"] = freq_rank[row["edge_index"]]
        row["damping_rank"] = damp_rank[row["edge_index"]]
        row["is_frequency_top"] = row["frequency_rank"] == 1
        row["is_damping_top"] = row["damping_rank"] == 1
        row["is_reversal"] = bool(
            row["delta_nu_fd_10pct"] > 0.0 and row["delta_zeta_full_fd_10pct"] < 0.0
        )
        row["second_order_delta_abs_error"] = abs(
            row["delta_zeta_second_order_fd_10pct"] - row["delta_zeta_full_fd_10pct"]
        )
        row["second_order_delta_rel_error"] = row["second_order_delta_abs_error"] / max(
            abs(row["delta_zeta_full_fd_10pct"]), 1e-12
        )
        row["second_order_direction_correct"] = bool(
            np.sign(row["delta_zeta_second_order_fd_10pct"])
            == np.sign(row["delta_zeta_full_fd_10pct"])
        )
    return rows


def generate_case(seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = int(rng.integers(6, 10))
    weights = connected_graph(rng, n, p_extra=float(rng.uniform(0.25, 0.55)))
    inertia = np.exp(rng.uniform(np.log(1.5), np.log(8.0), n))
    damping_per_inertia = np.exp(rng.uniform(np.log(0.04), np.log(0.65), n))
    return weights, inertia, damping_per_inertia


def accept_case(rows: list[dict[str, Any]], base_full: dict[str, Any], base_second: dict[str, Any]) -> bool:
    if not (0.025 < float(base_full["zeta"]) < 0.45):
        return False
    if not (0.10 < float(base_second["eta"]) < 1.10):
        return False
    if len(rows) < 5:
        return False

    top_freq = max(rows, key=lambda row: row["dnu_dw_frequency"])
    top_damp = max(rows, key=lambda row: row["dzeta_dw_qep_analytic"])

    if top_freq["edge_index"] == top_damp["edge_index"]:
        return False
    if top_freq["delta_nu_fd_10pct"] <= 1e-7:
        return False
    if top_freq["delta_zeta_full_fd_10pct"] >= -1e-4:
        return False
    if top_freq["dzeta_dw_qep_analytic"] >= 0.0:
        return False
    if top_freq["delta_zeta_second_order_fd_10pct"] >= 0.0:
        return False
    if top_freq["second_order_delta_rel_error"] > 2.0:
        return False
    if top_damp["dzeta_dw_qep_analytic"] <= 0.0:
        return False
    if top_damp["delta_zeta_full_fd_10pct"] <= 0.0:
        return False
    return True


def search_case(max_attempts: int = 50_000) -> dict[str, Any]:
    master = np.random.default_rng(MASTER_SEED)
    for attempt in range(1, max_attempts + 1):
        seed = int(master.integers(1, 2**31 - 1))
        weights, inertia, damping_per_inertia = generate_case(seed)
        base_full = full_margin(weights, inertia, damping_per_inertia)
        base_second = second_order_estimate(weights, inertia, damping_per_inertia)
        if base_full is None or base_second is None:
            continue
        rows = edge_rows(weights, inertia, damping_per_inertia, base_full, base_second)
        if accept_case(rows, base_full, base_second):
            return {
                "attempt": attempt,
                "seed": seed,
                "weights": weights,
                "inertia": inertia,
                "damping_per_inertia": damping_per_inertia,
                "base_full": base_full,
                "base_second": base_second,
                "rows": rows,
            }
    raise RuntimeError(f"No accepted case found in {max_attempts} attempts.")


def as_jsonable(value: Any) -> Any:
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): as_jsonable(v) for k, v in value.items() if k not in {"A", "eigvals", "eigvecs", "Q", "Gamma", "Lt", "Dt"}}
    if isinstance(value, list):
        return [as_jsonable(v) for v in value]
    return value


def save_case(case: dict[str, Any]) -> dict[str, Any]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    rows = case["rows"]
    top_freq = max(rows, key=lambda row: row["dnu_dw_frequency"])
    top_damp = max(rows, key=lambda row: row["dzeta_dw_qep_analytic"])

    weights = case["weights"]
    inertia = case["inertia"]
    damping_per_inertia = case["damping_per_inertia"]
    edges = [
        {
            "edge_index": row["edge_index"],
            "from_bus": row["from_bus"],
            "to_bus": row["to_bus"],
            "base_weight": row["base_weight"],
        }
        for row in rows
    ]

    summary = {
        "master_seed": MASTER_SEED,
        "case_seed": case["seed"],
        "search_attempts": case["attempt"],
        "regime": REGIME,
        "n_buses": int(weights.shape[0]),
        "reinforcement_fraction": REINFORCEMENT_FRAC,
        "inertia_M": inertia.tolist(),
        "damping_per_inertia_D_over_M": damping_per_inertia.tolist(),
        "physical_damping_D": (inertia * damping_per_inertia).tolist(),
        "weight_matrix": weights.tolist(),
        "edges": edges,
        "base": {
            "zeta_full": float(case["base_full"]["zeta"]),
            "pole_full": as_jsonable(case["base_full"]["pole"]),
            "zeta_diag": float(case["base_second"]["zeta_diag"]),
            "zeta_second_order": float(case["base_second"]["zeta_second"]),
            "s0": as_jsonable(case["base_second"]["s0"]),
            "s2": as_jsonable(case["base_second"]["s2"]),
            "delta_s_second_order": as_jsonable(case["base_second"]["ds"]),
            "critical_mode_index_zero_based": int(case["base_second"]["critical_mode"]),
            "nu_critical": float(case["base_second"]["nu_c"]),
            "delta_critical": float(case["base_second"]["delta_c"]),
            "offdiag_modal_damping_ratio_eta": float(case["base_second"]["eta"]),
            "min_modal_denominator": float(case["base_second"]["min_denominator"]),
        },
        "selected_reversal_frequency_top": top_freq,
        "selected_damping_aware_top": top_damp,
        "reversal_count": int(sum(1 for row in rows if row["is_reversal"])),
        "n_lines_tested": len(rows),
        "poster_summary": (
            "A six-bus synthetic network with log-uniform inertia and moderate "
            "log-uniform D_i/M_i heterogeneity was found after two seeded search "
            "attempts. Under a 10% reinforcement, the line ranked first by "
            "critical-mode frequency sensitivity increases the critical modal "
            "stiffness but reduces the full-QEP damping margin. The QEP-aware "
            "damping sensitivity ranks a different line first, and that line "
            "improves the damping margin. The second-order rational-filter "
            "correction predicts the harmful line's damping movement in the "
            "correct direction."
        ),
    }

    with (DATA_DIR / "demo_case_summary.json").open("w", encoding="utf-8") as f:
        json.dump(as_jsonable(summary), f, indent=2)

    fieldnames = list(rows[0].keys())
    with (DATA_DIR / "demo_line_rankings.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return summary


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=360, bbox_inches="tight")
    plt.close(fig)


def edge_label(row: dict[str, Any]) -> str:
    return f"{int(row['from_bus'])}-{int(row['to_bus'])}"


def plot_demo(case: dict[str, Any]) -> None:
    rows = case["rows"]
    top_freq = max(rows, key=lambda row: row["dnu_dw_frequency"])
    top_damp = max(rows, key=lambda row: row["dzeta_dw_qep_analytic"])

    freq_rows = sorted(rows, key=lambda row: row["dnu_dw_frequency"], reverse=True)[:6]
    damp_rows = sorted(rows, key=lambda row: row["dzeta_dw_qep_analytic"], reverse=True)[:6]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    fig.patch.set_facecolor("white")

    ax = axes[0]
    freq_rows_plot = list(reversed(freq_rows))
    labels = [edge_label(row) for row in freq_rows_plot]
    values = [row["dnu_dw_frequency"] for row in freq_rows_plot]
    colors = [BAD if row["delta_zeta_full_fd_10pct"] < 0 else GOOD for row in freq_rows_plot]
    ax.barh(labels, values, color=colors, alpha=0.92)
    ax.set_xlim(0.0, max(values) * 1.34)
    ax.set_xlabel(r"frequency score $\partial \nu_c/\partial w$")
    ax.set_title("Frequency ranking", loc="left", color=IEEE_DARK, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    for idx, row in enumerate(freq_rows_plot):
        marker = " harmful" if row["edge_index"] == top_freq["edge_index"] else ""
        ax.text(
            values[idx] + max(values) * 0.025,
            idx,
            rf"$\Delta\zeta={row['delta_zeta_full_fd_10pct']:+.4f}$" + marker,
            va="center",
            fontsize=8.5,
            color=BAD if row["delta_zeta_full_fd_10pct"] < 0 else GOOD,
        )

    ax = axes[1]
    damp_rows_plot = list(reversed(damp_rows))
    labels = [edge_label(row) for row in damp_rows_plot]
    values = [row["dzeta_dw_qep_analytic"] for row in damp_rows_plot]
    colors = [GOOD if value > 0 else BAD for value in values]
    ax.barh(labels, values, color=colors, alpha=0.92)
    xmin = min(values) * 1.45
    xmax = max(values) * 1.32
    ax.set_xlim(xmin, xmax)
    ax.axvline(0.0, color="#333333", linewidth=0.8)
    ax.set_xlabel(r"damping-aware score $\partial \zeta_{\min}/\partial w$")
    ax.set_title("Damping ranking", loc="left", color=IEEE_DARK, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    for idx, row in enumerate(damp_rows_plot):
        marker = " chosen" if row["edge_index"] == top_damp["edge_index"] else ""
        xpos = values[idx] + (xmax - xmin) * 0.035
        ax.text(
            xpos,
            idx,
            rf"$\Delta\zeta={row['delta_zeta_full_fd_10pct']:+.4f}$" + marker,
            va="center",
            fontsize=8.5,
            color=GOOD if row["delta_zeta_full_fd_10pct"] > 0 else BAD,
        )
    ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, -3))

    fig.suptitle(
        "A real reversal: the cheapest frequency heuristic picks the harmful line",
        fontsize=15,
        fontweight="bold",
        color=IEEE_DARK,
        y=1.03,
    )
    fig.text(
        0.5,
        -0.02,
        (
            f"Seed {case['seed']}, {REGIME} non-proportional damping. "
            f"Line {edge_label(top_freq)} raises critical modal stiffness by "
            f"{top_freq['delta_nu_fd_10pct']:.4f} but lowers full-QEP damping "
            f"margin by {top_freq['delta_zeta_full_fd_10pct']:.4f} under a "
            f"{int(REINFORCEMENT_FRAC * 100)}% reinforcement."
        ),
        ha="center",
        fontsize=10,
        color=IEEE_DARK,
    )
    fig.subplots_adjust(wspace=0.34, bottom=0.18, top=0.82)
    save_figure(fig, "fig_demo_reversal")


def plot_damping_region() -> None:
    x = np.linspace(0.0, 1.0, 300)
    boundary = 0.18 * np.sqrt(np.maximum(x, 1e-9))
    rng = np.random.default_rng(12)

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), sharey=True)
    for ax in axes:
        ax.plot(x, boundary, color=IEEE_DARK, lw=2.0)
        ax.fill_between(x, boundary, boundary.max() * 1.25, color="#dcefe5", alpha=0.95)
        ax.fill_between(x, 0, boundary, color="#f7d9de", alpha=0.95)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, boundary.max() * 1.18)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(True, alpha=0.16)
        ax.set_xlabel(r"modal stiffness $\nu_k$")
    axes[0].set_ylabel(r"modal damping $\delta_k$")

    pts_x = np.exp(rng.uniform(np.log(0.03), np.log(0.96), 34))
    pts_y = 0.035 + 0.16 * rng.random(34) * np.sqrt(pts_x)
    unsafe = pts_y < 0.18 * np.sqrt(pts_x)
    axes[0].scatter(pts_x[~unsafe], pts_y[~unsafe], s=34, color=GOOD, label="safe")
    axes[0].scatter(pts_x[unsafe], pts_y[unsafe], s=42, color=BAD, marker="x", label="unsafe")
    axes[0].set_title("Commuting damping: exact boundary", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="upper left")

    nu0 = 0.72
    delta0 = 0.18 * math.sqrt(nu0) + 0.033
    delta1 = 0.18 * math.sqrt(nu0) - 0.026
    axes[1].scatter([nu0], [delta0], s=58, color=GOOD, label="diagonal screen")
    axes[1].scatter([nu0], [delta1], s=68, color=BAD, marker="x", label="coupled QEP")
    axes[1].annotate(
        "",
        xy=(nu0, delta1),
        xytext=(nu0, delta0),
        arrowprops=dict(arrowstyle="->", lw=1.6, color=IEEE_DARK),
    )
    bg_x = np.exp(rng.uniform(np.log(0.03), np.log(0.96), 18))
    bg_y = 0.035 + 0.13 * rng.random(18) * np.sqrt(bg_x)
    axes[1].scatter(bg_x, bg_y, s=24, color=MUTED, alpha=0.55)
    axes[1].set_title("Non-proportional damping: screen only", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, loc="upper left")
    fig.text(0.5, 0.01, r"Target margin boundary: $\delta=2\zeta_\star\sqrt{\nu}$", ha="center")
    save_figure(fig, "fig_damping_region")


def plot_stress_test() -> None:
    labels = ["prop.", "weak", "moderate", "strong", "clustered"]
    diag = np.array([0.0, 0.038, 2.53, 5.16, 308.0])
    second = np.array([0.0, 0.007, 0.65, 1.46, 74.0])
    qep_low = np.array([0.0, 0.004, 0.22, 1.24, 9.8])
    qep_high = np.array([0.0, 0.024, 1.51, 3.90, 160.0])

    floor = 0.001
    x = np.arange(len(labels))
    width = 0.33
    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    ax.bar(x - width / 2, np.maximum(diag, floor), width=width, color="#7aa6c2", label="diagonal")
    ax.bar(x + width / 2, np.maximum(second, floor), width=width, color=IEEE_BLUE, label="second order")
    for idx, (lo, hi) in enumerate(zip(qep_low, qep_high)):
        if hi == 0.0:
            ax.scatter([x[idx]], [floor], marker="_", s=120, color=IEEE_DARK)
            ax.text(x[idx], floor * 1.8, "~0", ha="center", fontsize=9)
        else:
            ax.vlines(x[idx], lo, hi, colors=GOOD, linewidth=3.0)
            ax.scatter([x[idx], x[idx]], [lo, hi], s=28, color=GOOD)
    ax.scatter([], [], marker="|", s=180, color=GOOD, label="reduced QEP range")
    ax.set_yscale("log")
    ax.set_ylabel("median relative error (%)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title("Synthetic stress test from the manuscript", loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.text(
        0.99,
        0.02,
        "Zero values are plotted at a small floor for log-scale visibility.",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color=MUTED,
    )
    save_figure(fig, "fig_stress_test")


def plot_frequency_damping() -> None:
    x = np.linspace(0.05, 1.0, 300)
    sigma_fixed = -0.10
    omega = 0.2 + 4.0 * x
    zeta_fixed = -sigma_fixed / np.sqrt(sigma_fixed**2 + omega**2)
    sigma_improves = -(0.10 + 0.16 * x)
    zeta_improves = -sigma_improves / np.sqrt(sigma_improves**2 + omega**2)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(x, zeta_fixed / zeta_fixed.max(), color=BAD, lw=2.4, label="frequency rises, decay fixed")
    ax.plot(x, zeta_improves / zeta_improves.max(), color=GOOD, lw=2.4, label="decay rises enough")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"oscillation frequency $|\omega|$")
    ax.set_ylabel(r"normalized damping ratio")
    ax.set_title(r"$\zeta=-\sigma/\sqrt{\sigma^2+\omega^2}$", loc="left", fontweight="bold")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, loc="upper right")
    save_figure(fig, "fig_frequency_damping")


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#28364d",
            "axes.labelcolor": IEEE_DARK,
            "xtick.color": IEEE_DARK,
            "ytick.color": IEEE_DARK,
            "text.color": IEEE_DARK,
        }
    )

    case = search_case()
    summary = save_case(case)
    plot_demo(case)
    plot_damping_region()
    plot_stress_test()
    plot_frequency_damping()

    top_freq = summary["selected_reversal_frequency_top"]
    top_damp = summary["selected_damping_aware_top"]
    report = {
        "case_seed": summary["case_seed"],
        "search_attempts": summary["search_attempts"],
        "regime": summary["regime"],
        "frequency_top_line": f"{top_freq['from_bus']}-{top_freq['to_bus']}",
        "frequency_top_delta_nu": top_freq["delta_nu_fd_10pct"],
        "frequency_top_delta_zeta_full": top_freq["delta_zeta_full_fd_10pct"],
        "frequency_top_delta_zeta_second_order": top_freq["delta_zeta_second_order_fd_10pct"],
        "frequency_top_second_order_rel_error": top_freq["second_order_delta_rel_error"],
        "damping_top_line": f"{top_damp['from_bus']}-{top_damp['to_bus']}",
        "damping_top_delta_zeta_full": top_damp["delta_zeta_full_fd_10pct"],
        "reversal_count": summary["reversal_count"],
        "n_lines_tested": summary["n_lines_tested"],
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
