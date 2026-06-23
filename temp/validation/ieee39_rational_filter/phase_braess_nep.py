"""Phase Braess-NEP: test for resonant network-control Braess mechanisms.

This script uses the Phase-0D controller-aware Schur NEP bridge as the
validation engine.  It perturbs physical ANDES line parameters, reruns power
flow/eigenanalysis, rebuilds the Schur NEP, and compares the contour-extracted
zeros before and after reinforcement.

The test is deliberately conservative:

* ANDES/NEP poles are the reference for finite differences.
* No fixed-point solver is used.
* A "resonant" classification requires all of:
  (i) damping decreases under reinforcement,
  (ii) the tracked network-family mode moves closer to an A_cc control pole,
  (iii) the condensed self-energy contribution dominates the first-order
       Schur-NEP sensitivity.

If these conditions do not occur in the tested calibrated cases, the report
keeps the resonant Braess mechanism as a theoretical prediction not observed in
this benchmark.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from phase0d_nep_contour_solver import (
    RNG_SEED,
    RectangleContour,
    annotate_zeros,
    beyn_contour,
    block_partition,
    classify_partition,
    cluster_zeros,
    control_participation,
    damping_ratio,
    default_contours,
    freq_hz,
    load_andes_case,
    null_vector,
    paired_delta_omega_positions,
    schur_nep,
    schur_nep_derivative,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase_braess_nep"
TMP = OUT / "perturbed_cases"
EPS_DEFAULT = 0.01


def write_json(path: Path, data: Any) -> None:
    def default(obj: Any) -> Any:
        if isinstance(obj, complex):
            return {"real": float(obj.real), "imag": float(obj.imag)}
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return str(obj)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=default)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def git_commit_id() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unavailable"


def case_catalog() -> dict[str, Path]:
    return {
        "mix60": ROOT / "validation" / "ieee39_rational_filter" / "cases" / "ieee39_ibr_mix60.xlsx",
        "mix60_soft_gfl": ROOT
        / "validation"
        / "ieee39_rational_filter"
        / "cases"
        / "phase1_calibrated"
        / "soft_gfl"
        / "ieee39_ibr_mix60.xlsx",
        "mix60_no_pss": ROOT
        / "validation"
        / "ieee39_rational_filter"
        / "cases"
        / "phase1_calibrated"
        / "no_pss"
        / "ieee39_ibr_mix60.xlsx",
        "mix40_no_exciter_pss": ROOT
        / "validation"
        / "ieee39_rational_filter"
        / "cases"
        / "phase1_calibrated"
        / "no_exciter_pss"
        / "ieee39_ibr_mix40.xlsx",
        "mix20_no_exciter_pss": ROOT
        / "validation"
        / "ieee39_rational_filter"
        / "cases"
        / "phase1_calibrated"
        / "no_exciter_pss"
        / "ieee39_ibr_mix20.xlsx",
    }


@dataclass
class CaseNEP:
    case: str
    path: Path
    names: list[str]
    e_idx: list[int]
    c_idx: list[int]
    e_names: list[str]
    c_names: list[str]
    blocks: dict[str, np.ndarray]
    zeros: list[dict[str, Any]]
    acc_eigs: np.ndarray
    max_real: float


def load_case_nep(case: str, path: Path, contours: list[RectangleContour]) -> CaseNEP:
    ss = load_andes_case(str(path))
    names = [str(x) for x in ss.EIG.x_name]
    A = np.asarray(ss.EIG.As, dtype=float)
    e_idx, c_idx, _ = classify_partition(names)
    blocks = block_partition(A, e_idx, c_idx)
    c_names = [names[i] for i in c_idx]
    e_names = [names[i] for i in e_idx]

    raw: list[dict[str, Any]] = []
    for ci, contour in enumerate(contours):
        result = beyn_contour(
            blocks,
            contour,
            m_probe=min(blocks["Aee"].shape[0], max(12, blocks["Aee"].shape[0])),
            rank_tol=1e-10,
            residual_tol=1e-5,
            seed=RNG_SEED + ci,
        )
        for z in result["zeros"]:
            z["contour"] = contour.name
        raw.extend(result["zeros"])
    zeros = annotate_zeros(cluster_zeros(raw), blocks, c_names, tau=0.70)
    return CaseNEP(
        case=case,
        path=path,
        names=names,
        e_idx=e_idx,
        c_idx=c_idx,
        e_names=e_names,
        c_names=c_names,
        blocks=blocks,
        zeros=zeros,
        acc_eigs=np.linalg.eigvals(blocks["Acc"]),
        max_real=float(np.max(np.asarray(ss.EIG.mu, dtype=complex).real)),
    )


def active_lines(path: Path) -> list[dict[str, Any]]:
    lines = pd.read_excel(path, sheet_name="Line")
    out: list[dict[str, Any]] = []
    for pos, row in lines.iterrows():
        if int(row.get("u", 1)) != 1:
            continue
        out.append(
            {
                "line_pos": int(pos),
                "line_idx": str(row["idx"]),
                "from_bus": int(row["bus1"]),
                "to_bus": int(row["bus2"]),
                "r": float(row["r"]),
                "x": float(row["x"]),
                "trans": int(row.get("trans", 0)),
            }
        )
    return out


def perturb_line_case(src: Path, dst: Path, line_pos: int, eps: float) -> None:
    sheets = pd.read_excel(src, sheet_name=None)
    line = sheets["Line"].copy()
    scale = 1.0 / (1.0 + eps)
    # Series reinforcement: reduce r and x.  Shunt charging is left unchanged
    # to avoid conflating reinforcement with a separate shunt compensation.
    for col in ("r", "x"):
        line.loc[line.index[line_pos], col] = float(line.loc[line.index[line_pos], col]) * scale
    sheets["Line"] = line
    dst.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(dst, engine="openpyxl") as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)


def selected_zero(zeros: list[dict[str, Any]], family: str | None = None) -> dict[str, Any] | None:
    filt = [z for z in zeros if z["zeta"] > 0 and np.isfinite(float(z["zeta"]))]
    if family is not None:
        filt = [z for z in filt if z["family"] == family]
    if not filt:
        return None
    return min(filt, key=lambda z: (float(z["zeta"]), float(z["freq_hz"])))


def nearest_zero(zeros: list[dict[str, Any]], s: complex) -> dict[str, Any] | None:
    if not zeros:
        return None
    return min(zeros, key=lambda z: abs(complex(z["real"], z["imag"]) - s))


def acc_pole_info(Acc: np.ndarray, names: list[str], target: complex) -> dict[str, Any]:
    vals, vecs = np.linalg.eig(Acc)
    k = int(np.argmin(np.abs(vals - target)))
    pole = complex(vals[k])
    v = vecs[:, k]
    top = np.argsort(np.abs(v))[-6:][::-1] if len(v) else []
    top_states = "; ".join(f"{names[i]}:{abs(v[i]):.4g}" for i in top)
    return {
        "pole": pole,
        "top_states": top_states,
        "is_pll": "PLL" in top_states.upper(),
    }


def matched_acc_info(base: CaseNEP, post: CaseNEP, s: complex) -> tuple[dict[str, Any], dict[str, Any]]:
    base_info = acc_pole_info(base.blocks["Acc"], base.c_names, s)
    post_info = acc_pole_info(post.blocks["Acc"], post.c_names, complex(base_info["pole"]))
    return base_info, post_info


def network_candidates(zeros: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fam = [z for z in zeros if z["family"] == "I_network" and np.isfinite(float(z["zeta"]))]
    if not fam:
        return []
    fam = sorted(fam, key=lambda z: (float(z["zeta"]), float(z["freq_hz"])))
    z0 = float(fam[0]["zeta"])
    # Include the network-margin mode plus nearby low-damping family-I modes.
    # This catches the predicted regime where a not-quite-critical network mode
    # is close to a controller pole and could become critical after reinforcement.
    out = [z for z in fam if float(z["zeta"]) <= max(0.20, 1.30 * z0)]
    return out[:8]


def left_right_null(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    U, _, Vh = np.linalg.svd(T)
    y = U[:, -1]
    x = Vh.conj().T[:, -1]
    x = x / max(np.linalg.norm(x), 1e-15)
    y = y / max(np.linalg.norm(y), 1e-15)
    return y, x


def d_zeta_from_ds(s: complex, ds: complex) -> float:
    a, b = float(s.real), float(s.imag)
    da, db = float(ds.real), float(ds.imag)
    r = max(abs(s), 1e-15)
    return float(-da / r + a * (a * da + b * db) / (r**3))


def dS_components(
    base: CaseNEP,
    post: CaseNEP,
    s: complex,
    eps: float,
) -> dict[str, np.ndarray]:
    b, p = base.blocks, post.blocks
    Aee, Aec, Ace, Acc = b["Aee"], b["Aec"], b["Ace"], b["Acc"]
    dAee = (p["Aee"] - Aee) / eps
    dAec = (p["Aec"] - Aec) / eps
    dAce = (p["Ace"] - Ace) / eps
    dAcc = (p["Acc"] - Acc) / eps

    nc = Acc.shape[0]
    R = np.linalg.solve(s * np.eye(nc, dtype=complex) - Acc, np.eye(nc, dtype=complex))
    dS_retained = -dAee
    dS_self = -(dAec @ R @ Ace + Aec @ R @ dAce + Aec @ R @ dAcc @ R @ Ace)

    dpos, wpos, _ = paired_delta_omega_positions(base.e_names)
    mask_A = np.zeros_like(dS_retained)
    if dpos and wpos:
        # Stiffness/frequency proxy: omega-row, delta-column retained block.
        mask_A[np.ix_(wpos, dpos)] = dS_retained[np.ix_(wpos, dpos)]
    dS_A = mask_A
    dS_B = dS_retained - dS_A
    return {"A_stiffness": dS_A, "B_retained_rotation": dS_B, "C_self_energy": dS_self}


def sensitivity_terms(base: CaseNEP, post: CaseNEP, mode: dict[str, Any], eps: float) -> dict[str, Any]:
    s = complex(mode["real"], mode["imag"])
    T = schur_nep(base.blocks, s)
    y, x = left_right_null(T)
    Ts = schur_nep_derivative(base.blocks, s)
    denom = y.conj().T @ Ts @ x
    comps = dS_components(base, post, s, eps)

    rows: dict[str, Any] = {}
    total_ds = 0.0 + 0.0j
    total_dz = 0.0
    for name, dS in comps.items():
        ds = -(y.conj().T @ dS @ x) / denom
        dz = d_zeta_from_ds(s, complex(ds))
        rows[f"{name}_ds_real"] = float(np.real(ds))
        rows[f"{name}_ds_imag"] = float(np.imag(ds))
        rows[f"{name}_dzeta"] = float(dz)
        total_ds += ds
        total_dz += dz

    # Signature of frequency-dependent control in T'(s): I + Aec R^2 Ace.
    Aec, Ace, Acc = base.blocks["Aec"], base.blocks["Ace"], base.blocks["Acc"]
    nc = Acc.shape[0]
    eye_c = np.eye(nc, dtype=complex)
    R_Ace = np.linalg.solve(s * eye_c - Acc, Ace)
    R2_Ace = np.linalg.solve(s * eye_c - Acc, R_Ace)
    Ts_control = Aec @ R2_Ace
    denom_control = y.conj().T @ Ts_control @ x
    denom_identity = y.conj().T @ x

    rows.update(
        {
            "total_analytic_ds_real": float(np.real(total_ds)),
            "total_analytic_ds_imag": float(np.imag(total_ds)),
            "total_analytic_dzeta": float(total_dz),
            "denominator_abs": float(abs(denom)),
            "control_denominator_abs": float(abs(denom_control)),
            "identity_denominator_abs": float(abs(denom_identity)),
            "control_denominator_fraction": float(abs(denom_control) / max(abs(denom), 1e-15)),
        }
    )
    return rows


def classify_result(row: dict[str, Any]) -> str:
    if not row["margin_reversal"] and not row["network_mode_reversal"]:
        return "no_reversal"
    if row["nonsmooth"]:
        return "nonsmooth_unclassified"
    terms = {
        "A": abs(row["A_stiffness_dzeta"]),
        "B": abs(row["B_retained_rotation_dzeta"]),
        "C": abs(row["C_self_energy_dzeta"]),
    }
    dominant = max(terms, key=terms.get)
    row["dominant_term"] = dominant
    near_control = bool(
        row["nearest_acc_is_pll"]
        and min(row["relative_distance_to_acc_base"], row["relative_distance_to_acc_post"]) <= 0.10
    )
    resonant = bool(
        row["network_mode_reversal"]
        and row["distance_to_acc_delta"] < 0
        and dominant == "C"
        and row["C_self_energy_dzeta"] < 0
        and near_control
    )
    if resonant:
        return "resonant"
    if dominant == "A":
        return "geometric_or_stiffness"
    if dominant == "C":
        return "self_energy_nonresonant"
    return "retained_rotation"


def analyze_line(
    base: CaseNEP,
    case_path: Path,
    line: dict[str, Any],
    eps: float,
    contours: list[RectangleContour],
) -> list[dict[str, Any]]:
    pert_path = TMP / base.case / f"{line['line_idx']}_eps{eps:.4f}.xlsx"
    perturb_line_case(case_path, pert_path, line["line_pos"], eps)
    post = load_case_nep(base.case + "__" + line["line_idx"], pert_path, contours)

    base_margin = selected_zero(base.zeros)
    post_margin = selected_zero(post.zeros)
    candidates = network_candidates(base.zeros)
    if base_margin is None or post_margin is None or not candidates:
        return [{
            "case": base.case,
            **line,
            "ok": False,
            "reason": "missing base/post margin or base network-family zero",
        }]

    rows: list[dict[str, Any]] = []
    post_family = [z for z in post.zeros if z["family"] == "I_network"]
    for mode_rank, base_net in enumerate(candidates, start=1):
        s_net = complex(base_net["real"], base_net["imag"])
        post_net = nearest_zero(post_family, s_net) or nearest_zero(post.zeros, s_net)
        assert post_net is not None

        base_acc, post_acc = matched_acc_info(base, post, s_net)
        s_post_net = complex(post_net["real"], post_net["imag"])
        base_acc_pole = complex(base_acc["pole"])
        post_acc_pole = complex(post_acc["pole"])
        base_distance = abs(s_net - base_acc_pole)
        post_distance = abs(s_post_net - post_acc_pole)

        terms = sensitivity_terms(base, post, base_net, eps)
        fd_network_dz = (float(post_net["zeta"]) - float(base_net["zeta"])) / eps
        fd_margin_dz = (float(post_margin["zeta"]) - float(base_margin["zeta"])) / eps
        analytic = float(terms["total_analytic_dzeta"])
        rel_error = abs(analytic - fd_network_dz) / max(abs(fd_network_dz), 1e-12)
        nonsmooth = bool(abs(s_post_net - s_net) > 0.40 * max(abs(s_net), 1.0))

        row: dict[str, Any] = {
            "case": base.case,
            **line,
            "ok": True,
            "eps": eps,
            "mode_rank": mode_rank,
            "base_max_real": base.max_real,
            "post_max_real": post.max_real,
            "base_margin_real": float(base_margin["real"]),
            "base_margin_imag": float(base_margin["imag"]),
            "base_margin_zeta": float(base_margin["zeta"]),
            "base_margin_family": base_margin["family"],
            "post_margin_real": float(post_margin["real"]),
            "post_margin_imag": float(post_margin["imag"]),
            "post_margin_zeta": float(post_margin["zeta"]),
            "post_margin_family": post_margin["family"],
            "delta_margin_zeta": float(post_margin["zeta"] - base_margin["zeta"]),
            "fd_margin_dzeta": float(fd_margin_dz),
            "base_network_real": float(base_net["real"]),
            "base_network_imag": float(base_net["imag"]),
            "base_network_freq_hz": float(base_net["freq_hz"]),
            "base_network_zeta": float(base_net["zeta"]),
            "base_network_pi_c": float(base_net["pi_c"]),
            "post_network_real": float(post_net["real"]),
            "post_network_imag": float(post_net["imag"]),
            "post_network_freq_hz": float(post_net["freq_hz"]),
            "post_network_zeta": float(post_net["zeta"]),
            "post_network_pi_c": float(post_net["pi_c"]),
            "delta_network_zeta": float(post_net["zeta"] - base_net["zeta"]),
            "fd_network_dzeta": float(fd_network_dz),
            "analytic_vs_fd_rel_error": float(rel_error),
            "nonsmooth": nonsmooth,
            "base_nearest_acc_real": float(base_acc_pole.real),
            "base_nearest_acc_imag": float(base_acc_pole.imag),
            "base_nearest_acc_top_states": base_acc["top_states"],
            "post_matched_acc_real": float(post_acc_pole.real),
            "post_matched_acc_imag": float(post_acc_pole.imag),
            "post_matched_acc_top_states": post_acc["top_states"],
            "nearest_acc_is_pll": bool(base_acc["is_pll"] or post_acc["is_pll"]),
            "distance_to_acc_base": float(base_distance),
            "distance_to_acc_post": float(post_distance),
            "relative_distance_to_acc_base": float(base_distance / max(abs(s_net), 1e-15)),
            "relative_distance_to_acc_post": float(post_distance / max(abs(s_post_net), 1e-15)),
            "distance_to_acc_delta": float(post_distance - base_distance),
            "margin_reversal": bool(post_margin["zeta"] < base_margin["zeta"]),
            "network_mode_reversal": bool(post_net["zeta"] < base_net["zeta"]),
            "dominant_term": "",
        }
        row.update(terms)
        row["classification"] = classify_result(row)
        rows.append(row)
    return rows


def plot_resonant_candidate(rows: list[dict[str, Any]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    candidates = [r for r in rows if r.get("classification") == "resonant"]
    if not candidates:
        candidates = [r for r in rows if r.get("network_mode_reversal")]
    if not candidates:
        return
    r = min(candidates, key=lambda x: x.get("fd_network_dzeta", 0.0))
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.scatter([r["base_network_real"]], [r["base_network_imag"]], label="base network mode", color="#1f77b4")
    ax.scatter([r["post_network_real"]], [r["post_network_imag"]], label="reinforced network mode", color="#ff7f0e")
    ax.scatter([r["base_nearest_acc_real"]], [r["base_nearest_acc_imag"]], label="base nearest Acc pole", color="#d62728", marker="x")
    ax.scatter([r["post_matched_acc_real"]], [r["post_matched_acc_imag"]], label="post matched Acc pole", color="#9467bd", marker="x")
    ax.annotate("", xy=(r["post_network_real"], r["post_network_imag"]), xytext=(r["base_network_real"], r["base_network_imag"]), arrowprops={"arrowstyle": "->"})
    ax.axvline(0, color="#888888", linewidth=0.8)
    ax.grid(True, alpha=0.25)
    ax.set_xlabel("Re(s)")
    ax.set_ylabel("Im(s) [rad/s]")
    ax.set_title(f"{r['case']} {r['line_idx']} classification={r['classification']}")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig_braess_nep_s_plane.png", dpi=220)
    fig.savefig(OUT / "fig_braess_nep_s_plane.pdf")
    plt.close(fig)


def write_report(status: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    resonant = [r for r in rows if r.get("classification") == "resonant"]
    reversals = [r for r in rows if r.get("margin_reversal")]
    net_reversals = [r for r in rows if r.get("network_mode_reversal")]
    smooth = [r for r in rows if r.get("ok") and not r.get("nonsmooth")]
    best = None
    if resonant:
        best = min(resonant, key=lambda r: r["fd_network_dzeta"])
    elif net_reversals:
        best = min(net_reversals, key=lambda r: r["fd_network_dzeta"])
    elif rows:
        best = min([r for r in rows if r.get("ok")], key=lambda r: r.get("fd_network_dzeta", math.inf), default=None)

    lines = [
        "# Phase Braess-NEP Report",
        "",
        f"Gate/verdict: **{status['verdict']}**",
        f"Random seed: `{status['random_seed']}`",
        f"Git commit: `{status['git_commit']}`",
        "",
        "## Method",
        "Each line reinforcement reduces the physical ANDES `r` and `x` fields by `1/(1+eps)` and reruns power flow, full eigensolve, and the Phase-0D Beyn Schur-NEP contour solver.  The finite-difference reference is the NEP/ANDES pole movement, not a scalar fixed-point model.",
        "",
        "The first-order Schur-NEP sensitivity is decomposed as:",
        "- A: retained omega-row / delta-column stiffness proxy;",
        "- B: residual retained-network/modal-rotation proxy;",
        "- C: condensed-control self-energy contribution from the Schur term;",
        "with the control-frequency signature reported as the fraction of the denominator coming from `A_ec (sI-A_cc)^(-2) A_ce`.",
        "",
        "A resonant Braess classification requires a network-family damping decrease, decreasing distance to a matched `A_cc` pole, dominant negative C contribution, relative pole distance below 10%, and a nearest `A_cc` pole whose eigenvector is PLL-dominated.  Otherwise the reversal is not claimed resonant.",
        "",
        "## Summary",
        f"- Cases tested: {', '.join(status['cases_tested'])}",
        f"- Line-mode perturbations completed: {len([r for r in rows if r.get('ok')])}",
        f"- Margin reversals: {len(reversals)}",
        f"- Network-family mode reversals: {len(net_reversals)}",
        f"- Resonant reversals: {len(resonant)}",
        f"- Smooth analytic-vs-NEP median relative error: {status.get('median_smooth_rel_error', 'nan')}",
        "",
    ]
    if best:
        lines.extend(
            [
                "## Clearest Case",
                f"- Case/line/mode-rank: `{best['case']}` `{best['line_idx']}` ({best['from_bus']}-{best['to_bus']}), mode {best.get('mode_rank')}",
                f"- Classification: **{best['classification']}**",
                f"- Margin zeta: {best['base_margin_zeta']:.6g} -> {best['post_margin_zeta']:.6g} (Delta={best['delta_margin_zeta']:.3e})",
                f"- Network zeta: {best['base_network_zeta']:.6g} -> {best['post_network_zeta']:.6g} (FD derivative={best['fd_network_dzeta']:.3e})",
                f"- Distance to matched Acc pole: {best['distance_to_acc_base']:.6g} -> {best['distance_to_acc_post']:.6g} (Delta={best['distance_to_acc_delta']:.3e})",
                f"- Relative Acc distance: {best['relative_distance_to_acc_base']:.3e} -> {best['relative_distance_to_acc_post']:.3e}; nearest Acc PLL-dominated={best['nearest_acc_is_pll']}",
                f"- Base nearest Acc top states: {best.get('base_nearest_acc_top_states', '')}",
                f"- A dzeta={best['A_stiffness_dzeta']:.3e}, B dzeta={best['B_retained_rotation_dzeta']:.3e}, C dzeta={best['C_self_energy_dzeta']:.3e}",
                f"- Analytic total dzeta={best['total_analytic_dzeta']:.3e}, FD network dzeta={best['fd_network_dzeta']:.3e}, relative error={best['analytic_vs_fd_rel_error']:.3e}",
                f"- Control denominator fraction={best['control_denominator_fraction']:.3e}",
                "",
            ]
        )
    lines.extend(
        [
            "## Verdict",
            status["interpretation"],
            "",
            "Artifacts:",
            f"- CSV: `{status['outputs']['line_csv']}`",
            f"- JSON: `{status['outputs']['status_json']}`",
            f"- Figure: `{status['outputs'].get('figure_png', 'not generated')}`",
        ]
    )
    (OUT / "phase_braess_nep_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="mix60,mix60_soft_gfl", help="Comma-separated case keys.")
    parser.add_argument("--eps", type=float, default=EPS_DEFAULT)
    parser.add_argument("--max-lines", type=int, default=0, help="0 means all active lines.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    np.random.seed(RNG_SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)

    catalog = case_catalog()
    contours = default_contours()
    selected = [x.strip() for x in args.cases.split(",") if x.strip()]
    rows: list[dict[str, Any]] = []
    case_summaries: list[dict[str, Any]] = []

    for case in selected:
        path = catalog.get(case)
        if path is None or not path.exists():
            case_summaries.append({"case": case, "ok": False, "reason": "case path missing", "path": str(path)})
            continue
        print(f"[braess-nep] base case {case}: {path}")
        try:
            base = load_case_nep(case, path, contours)
        except Exception as exc:  # noqa: BLE001 - validation report must capture blockers.
            case_summaries.append({"case": case, "ok": False, "reason": f"base load failed: {exc}", "path": str(path)})
            continue
        lines = active_lines(path)
        if args.max_lines and args.max_lines > 0:
            lines = lines[: args.max_lines]
        case_summaries.append(
            {
                "case": case,
                "ok": True,
                "path": str(path),
                "n_lines": len(lines),
                "base_max_real": base.max_real,
                "base_margin_zeta": selected_zero(base.zeros)["zeta"] if selected_zero(base.zeros) else None,
                "base_network_zeta": selected_zero(base.zeros, "I_network")["zeta"] if selected_zero(base.zeros, "I_network") else None,
            }
        )
        for k, line in enumerate(lines, start=1):
            print(f"[braess-nep] {case} line {k}/{len(lines)} {line['line_idx']} {line['from_bus']}-{line['to_bus']}")
            try:
                rows.extend(analyze_line(base, path, line, args.eps, contours))
            except Exception as exc:  # noqa: BLE001
                rows.append({"case": case, **line, "ok": False, "reason": str(exc), "eps": args.eps})

    ok_rows = [r for r in rows if r.get("ok")]
    smooth_errors = [r["analytic_vs_fd_rel_error"] for r in ok_rows if not r.get("nonsmooth") and np.isfinite(r.get("analytic_vs_fd_rel_error", np.nan))]
    resonant = [r for r in ok_rows if r.get("classification") == "resonant"]
    margin_reversals = [r for r in ok_rows if r.get("margin_reversal")]
    net_reversals = [r for r in ok_rows if r.get("network_mode_reversal")]
    if resonant:
        verdict = "PASS_RESONANT_FOUND"
        interp = (
            "At least one line reinforcement produced a network-family damping decrease, moved the tracked mode closer "
            "to a matched condensed-control pole, and was dominated by the self-energy term.  The resonant Braess "
            "mechanism is observed in this benchmark."
        )
    elif net_reversals or margin_reversals:
        verdict = "NO_RESONANT_OBSERVED"
        interp = (
            "Line reinforcement can reduce damping in the tested cases, but no completed line perturbation satisfied "
            "the full resonant criteria.  The observed reversals should be reported as non-resonant/geometric or "
            "retained-network effects; resonant Braess remains a theoretical prediction pending a case in the correct "
            "near-pole regime."
        )
    else:
        verdict = "NO_REVERSAL_OBSERVED"
        interp = (
            "No damping-margin or tracked network-family reversal was observed under the tested physical line "
            "reinforcements.  The resonant Braess claim is not supported by this benchmark."
        )

    status = {
        "gate": "Braess resonant network-control mechanism over Phase-0D NEP bridge",
        "verdict": verdict,
        "random_seed": RNG_SEED,
        "git_commit": git_commit_id(),
        "eps": args.eps,
        "cases_requested": selected,
        "cases_tested": [c["case"] for c in case_summaries if c.get("ok")],
        "case_summaries": case_summaries,
        "n_rows": len(rows),
        "n_ok": len(ok_rows),
        "n_margin_reversals": len(margin_reversals),
        "n_network_mode_reversals": len(net_reversals),
        "n_resonant": len(resonant),
        "median_smooth_rel_error": float(np.median(smooth_errors)) if smooth_errors else None,
        "interpretation": interp,
        "outputs": {
            "line_csv": str((OUT / "phase_braess_nep_line_decomposition.csv").relative_to(ROOT)),
            "status_json": str((OUT / "phase_braess_nep_status.json").relative_to(ROOT)),
            "report": str((OUT / "phase_braess_nep_report.md").relative_to(ROOT)),
            "figure_png": str((OUT / "fig_braess_nep_s_plane.png").relative_to(ROOT)),
        },
    }
    write_csv(OUT / "phase_braess_nep_line_decomposition.csv", rows)
    write_json(OUT / "phase_braess_nep_status.json", status)
    plot_resonant_candidate(ok_rows)
    write_report(status, ok_rows)
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
