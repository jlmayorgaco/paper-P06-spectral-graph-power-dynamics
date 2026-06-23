"""Gate: does inertia placement survive the controller-aware Schur bridge?

This script tests the boundary between two parts of the project:

1. The controller-aware Schur NEP bridge, already validated against full ANDES
   poles in diagnostic cases.
2. The network-inertia structural theory, which was derived without the
   frequency-dependent controller self-energy.

The gate is deliberately conservative.  It perturbs only documented inertia
parameters present in the ANDES Excel cases.  In the current IEEE39 IBR files,
REGF1 does not expose a clear virtual-inertia column, so this script tests
active GENROU.M perturbations and reports virtual-inertia validation as blocked.

No manuscript files are edited by this script.
"""

from __future__ import annotations

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
    one_sided_modes,
    schur_nep,
    schur_nep_derivative,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase_inertia_bridge_gate"
TMP = OUT / "perturbed_cases"
EPS_DEFAULT = 0.01
TAU_CONTROL = 0.70
SIGN_TOL = 1e-8


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
    import andes  # type: ignore

    return {
        "no_pss": ROOT / "outputs" / "ieee39_mode_diagnostics" / "ieee39_no_ieeest.xlsx",
        "base": Path(andes.get_case("ieee39/ieee39_full.xlsx")),
        "mix60": ROOT / "validation" / "ieee39_rational_filter" / "cases" / "ieee39_ibr_mix60.xlsx",
    }


def read_case(path: Path) -> dict[str, pd.DataFrame]:
    return pd.read_excel(path, sheet_name=None)


def write_case(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)


def active_mask(df: pd.DataFrame) -> pd.Series:
    if "u" not in df.columns:
        return pd.Series([True] * len(df), index=df.index)
    return pd.to_numeric(df["u"], errors="coerce").fillna(1).astype(float) > 0


def active_genrou_inertia(path: Path) -> list[dict[str, Any]]:
    sheets = read_case(path)
    if "GENROU" not in sheets:
        return []
    df = sheets["GENROU"]
    required = {"idx", "M"}
    if not required.issubset(df.columns):
        return []
    rows: list[dict[str, Any]] = []
    for pos, row in df.loc[active_mask(df)].iterrows():
        rows.append(
            {
                "row_pos": int(pos),
                "idx": str(row["idx"]),
                "bus": str(row.get("bus", "")),
                "M": float(row["M"]),
                "D": float(row["D"]) if "D" in df.columns else float("nan"),
            }
        )
    return rows


def virtual_inertia_columns(path: Path) -> dict[str, list[str]]:
    sheets = read_case(path)
    candidates: dict[str, list[str]] = {}
    names = {"M", "H", "Tj", "J", "Ta", "Tpm", "wdrp", "gammap"}
    for sheet_name, df in sheets.items():
        if sheet_name.upper().startswith(("REGF", "VSG", "GFM")):
            cols = [str(c) for c in df.columns if str(c) in names or "inert" in str(c).lower()]
            if cols:
                candidates[sheet_name] = cols
    return candidates


def perturb_genrou_m(src: Path, dst: Path, genrou_idx: str, eps: float) -> dict[str, Any]:
    sheets = read_case(src)
    if "GENROU" not in sheets:
        raise RuntimeError(f"No GENROU sheet in {src}")
    df = sheets["GENROU"].copy()
    if "idx" not in df.columns or "M" not in df.columns:
        raise RuntimeError(f"GENROU sheet lacks idx/M columns in {src}")
    mask = active_mask(df) & (df["idx"].astype(str) == genrou_idx)
    if int(mask.sum()) != 1:
        raise RuntimeError(f"Expected exactly one active GENROU row for {genrou_idx}; found {int(mask.sum())}")
    old_m = float(df.loc[mask, "M"].iloc[0])
    df["M"] = pd.to_numeric(df["M"], errors="coerce").astype(float)
    df.loc[mask, "M"] = old_m * (1.0 + eps)
    sheets["GENROU"] = df
    write_case(dst, sheets)
    return {"idx": genrou_idx, "old_M": old_m, "new_M": old_m * (1.0 + eps), "eps": eps}


@dataclass
class LoadedCase:
    label: str
    path: Path
    ss: Any
    names: list[str]
    blocks: dict[str, np.ndarray]
    c_names: list[str]
    full_modes: list[dict[str, Any]]
    bridge_zeros: list[dict[str, Any]]
    max_real: float


def full_modes_with_family(ss: Any, blocks: dict[str, np.ndarray], max_modes: int = 24) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    for rank, s in enumerate(one_sided_modes(mu)[:max_modes], start=1):
        x_e = null_vector(blocks, s)
        pi_c, pi_c_raw, _ = control_participation(blocks, s, x_e)
        rows.append(
            {
                "rank": rank,
                "pole": complex(s),
                "real": float(s.real),
                "imag": float(s.imag),
                "freq_hz": freq_hz(s),
                "zeta": damping_ratio(s),
                "pi_c": float(pi_c),
                "pi_c_raw": float(pi_c_raw),
                "family": "II_control" if pi_c >= TAU_CONTROL else "I_network",
            }
        )
    return rows


def bridge_zeros(blocks: dict[str, np.ndarray], c_names: list[str]) -> list[dict[str, Any]]:
    raw: list[dict[str, Any]] = []
    for ci, contour in enumerate(default_contours()):
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
    return annotate_zeros(cluster_zeros(raw), blocks, c_names, tau=TAU_CONTROL)


def load_case(label: str, path: Path, compute_bridge: bool = True) -> LoadedCase:
    ss = load_andes_case(str(path))
    names = [str(x) for x in ss.EIG.x_name]
    A = np.asarray(ss.EIG.As, dtype=float)
    e_idx, c_idx, _ = classify_partition(names)
    blocks = block_partition(A, e_idx, c_idx)
    c_names = [names[i] for i in c_idx]
    full = full_modes_with_family(ss, blocks)
    zeros = bridge_zeros(blocks, c_names) if compute_bridge else []
    return LoadedCase(
        label=label,
        path=path,
        ss=ss,
        names=names,
        blocks=blocks,
        c_names=c_names,
        full_modes=full,
        bridge_zeros=zeros,
        max_real=float(np.max(np.asarray(ss.EIG.mu, dtype=complex).real)),
    )


def first_mode(rows: list[dict[str, Any]]) -> dict[str, Any]:
    if not rows:
        raise RuntimeError("No oscillatory modes found")
    return rows[0]


def match_mode(rows: list[dict[str, Any]], s_ref: complex) -> dict[str, Any]:
    if not rows:
        raise RuntimeError("No oscillatory modes found")
    return min(rows, key=lambda row: abs(complex(row["pole"]) - s_ref))


def left_null_vector(blocks: dict[str, np.ndarray], s: complex) -> np.ndarray:
    """Left null vector y satisfying y^H T(s) ~= 0."""
    _, _, vh = np.linalg.svd(schur_nep(blocks, s).conj().T)
    return vh.conj().T[:, -1]


def dzeta_from_ds(s: complex, ds: complex) -> float:
    """Directional derivative of zeta(s)=-Re(s)/|s|."""
    r = abs(s)
    if r < 1e-14:
        return float("nan")
    return float(-ds.real / r + s.real * np.real(np.conj(s) * ds) / (r**3))


def local_nep_zeta_derivative(
    base_blocks: dict[str, np.ndarray],
    pert_blocks: dict[str, np.ndarray],
    s0: complex,
    eps: float,
) -> tuple[complex, float, float]:
    """Local NEP sensitivity using T_p from a finite difference of Schur operators.

    The available ANDES sheets do not expose a first-principles symbolic mapping
    from each Excel inertia parameter to every affected Schur block.  Instead of
    inventing one, this computes the local operator perturbation directly:

        T_p(s0) ~= (T_pert(s0) - T_base(s0)) / eps
        ds/deps = - y^H T_p(s0) x / y^H T_s(s0) x.

    This is still a NEP sensitivity at the base pole.  It does not use the
    perturbed eigensolve to choose or fit the answer.
    """
    x = null_vector(base_blocks, s0)
    y = left_null_vector(base_blocks, s0)
    t_p = (schur_nep(pert_blocks, s0) - schur_nep(base_blocks, s0)) / eps
    t_s = schur_nep_derivative(base_blocks, s0)
    den = y.conj() @ t_s @ x
    if abs(den) < 1e-14:
        return complex(np.nan, np.nan), float("nan"), float("nan")
    ds_deps = -complex(y.conj() @ t_p @ x) / complex(den)
    dzeta_deps = dzeta_from_ds(s0, ds_deps)
    residual = float(np.linalg.norm(schur_nep(base_blocks, s0) @ x))
    return ds_deps, dzeta_deps, residual


def sign_value(x: float, tol: float = SIGN_TOL) -> int:
    if x > tol:
        return 1
    if x < -tol:
        return -1
    return 0


def top_set(rows: list[dict[str, Any]], key: str, n: int = 3) -> set[str]:
    ranked = sorted(rows, key=lambda r: r[key], reverse=True)
    return {r["idx"] for r in ranked[: min(n, len(ranked))]}


def run_case_gate(case_name: str, case_path: Path, eps: float) -> dict[str, Any]:
    case_out = OUT / case_name
    pert_dir = TMP / case_name
    case_out.mkdir(parents=True, exist_ok=True)
    pert_dir.mkdir(parents=True, exist_ok=True)

    base = load_case(case_name, case_path, compute_bridge=True)
    genrous = active_genrou_inertia(case_path)
    virtual_cols = virtual_inertia_columns(case_path)
    if not genrous:
        return {
            "case": case_name,
            "case_path": str(case_path),
            "status": "BLOCKED_NO_INERTIA_PARAMETER",
            "reason": "No active GENROU.M parameter found.",
            "virtual_inertia_columns": virtual_cols,
        }

    base_full_crit = first_mode(base.full_modes)
    base_bridge_crit = first_mode(base.bridge_zeros)
    base_full_s = complex(base_full_crit["pole"])
    base_bridge_s = complex(base_bridge_crit["pole"])

    rows: list[dict[str, Any]] = []
    for gen in genrous:
        idx = gen["idx"]
        dst = pert_dir / f"{case_path.stem}_M_{idx}_eps{eps:.4g}.xlsx"
        perturb_info = perturb_genrou_m(case_path, dst, idx, eps)
        pert = load_case(f"{case_name}_{idx}", dst, compute_bridge=True)

        full_crit = first_mode(pert.full_modes)
        bridge_crit = first_mode(pert.bridge_zeros)
        full_matched = match_mode(pert.full_modes, base_full_s)
        bridge_matched = match_mode(pert.bridge_zeros, base_bridge_s)

        full_delta_global = float(full_crit["zeta"] - base_full_crit["zeta"])
        bridge_delta_global = float(bridge_crit["zeta"] - base_bridge_crit["zeta"])
        full_delta_matched = float(full_matched["zeta"] - base_full_crit["zeta"])
        bridge_delta_matched = float(bridge_matched["zeta"] - base_bridge_crit["zeta"])
        full_deriv = full_delta_global / eps
        bridge_deriv = bridge_delta_global / eps
        full_deriv_matched = full_delta_matched / eps
        bridge_deriv_matched = bridge_delta_matched / eps
        nep_ds_deps, nep_local_deriv, nep_local_residual = local_nep_zeta_derivative(
            base.blocks, pert.blocks, base_bridge_s, eps
        )
        nep_local_delta = nep_local_deriv * eps

        sign_match_global = sign_value(full_delta_global) == sign_value(bridge_delta_global)
        sign_match_matched = sign_value(full_delta_matched) == sign_value(bridge_delta_matched)
        sign_match_nep_local_full = sign_value(full_delta_matched) == sign_value(nep_local_delta)
        sign_match_nep_local_bridge = sign_value(bridge_delta_matched) == sign_value(nep_local_delta)
        mag_rel_error = abs(bridge_delta_global - full_delta_global) / max(abs(full_delta_global), 1e-12)
        mag_rel_error_matched = abs(bridge_delta_matched - full_delta_matched) / max(
            abs(full_delta_matched), 1e-12
        )
        mag_rel_error_nep_local_full = abs(nep_local_delta - full_delta_matched) / max(
            abs(full_delta_matched), 1e-12
        )
        mag_rel_error_nep_local_bridge = abs(nep_local_delta - bridge_delta_matched) / max(
            abs(bridge_delta_matched), 1e-12
        )
        rows.append(
            {
                "case": case_name,
                "idx": idx,
                "bus": gen["bus"],
                "old_M": perturb_info["old_M"],
                "new_M": perturb_info["new_M"],
                "eps": eps,
                "base_full_zeta": base_full_crit["zeta"],
                "base_full_family": base_full_crit["family"],
                "base_full_pi_c": base_full_crit["pi_c"],
                "base_bridge_zeta": base_bridge_crit["zeta"],
                "base_bridge_family": base_bridge_crit["family"],
                "base_bridge_pi_c": base_bridge_crit["pi_c"],
                "full_global_zeta": full_crit["zeta"],
                "bridge_global_zeta": bridge_crit["zeta"],
                "full_delta_global": full_delta_global,
                "bridge_delta_global": bridge_delta_global,
                "full_deriv_global": full_deriv,
                "bridge_deriv_global": bridge_deriv,
                "sign_match_global": bool(sign_match_global),
                "mag_rel_error_global": mag_rel_error,
                "full_matched_zeta": full_matched["zeta"],
                "bridge_matched_zeta": bridge_matched["zeta"],
                "full_delta_matched": full_delta_matched,
                "bridge_delta_matched": bridge_delta_matched,
                "full_deriv_matched": full_deriv_matched,
                "bridge_deriv_matched": bridge_deriv_matched,
                "nep_local_ds_deps_real": float(nep_ds_deps.real),
                "nep_local_ds_deps_imag": float(nep_ds_deps.imag),
                "nep_local_deriv_matched": nep_local_deriv,
                "nep_local_delta_matched": nep_local_delta,
                "nep_local_residual": nep_local_residual,
                "sign_match_matched": bool(sign_match_matched),
                "sign_match_nep_local_full": bool(sign_match_nep_local_full),
                "sign_match_nep_local_bridge": bool(sign_match_nep_local_bridge),
                "mag_rel_error_matched": mag_rel_error_matched,
                "mag_rel_error_nep_local_full": mag_rel_error_nep_local_full,
                "mag_rel_error_nep_local_bridge": mag_rel_error_nep_local_bridge,
                "full_global_family": full_crit["family"],
                "bridge_global_family": bridge_crit["family"],
                "global_mode_switched_full": bool(abs(complex(full_crit["pole"]) - complex(full_matched["pole"])) > 1e-5),
                "global_mode_switched_bridge": bool(
                    abs(complex(bridge_crit["pole"]) - complex(bridge_matched["pole"])) > 1e-5
                ),
                "perturbed_case": str(dst),
            }
        )

    n = len(rows)
    sign_rate_global = sum(1 for r in rows if r["sign_match_global"]) / max(n, 1)
    sign_rate_matched = sum(1 for r in rows if r["sign_match_matched"]) / max(n, 1)
    top_full = top_set(rows, "full_delta_global", n=3)
    top_bridge = top_set(rows, "bridge_delta_global", n=3)
    top_nep_local = top_set(rows, "nep_local_delta_matched", n=3)
    top3_overlap = len(top_full & top_bridge) / max(min(3, n), 1)
    top3_overlap_nep_local = len(top_full & top_nep_local) / max(min(3, n), 1)
    negative_full = [r for r in rows if r["full_delta_global"] < -SIGN_TOL]
    negative_bridge = [r for r in rows if r["bridge_delta_global"] < -SIGN_TOL]
    finite_signal = [r for r in rows if abs(r["full_delta_global"]) >= SIGN_TOL]
    sign_rate_nep_local_full = sum(1 for r in rows if r["sign_match_nep_local_full"]) / max(n, 1)
    sign_rate_nep_local_bridge = sum(1 for r in rows if r["sign_match_nep_local_bridge"]) / max(n, 1)

    status_gate1 = (
        "PASS"
        if sign_rate_global >= 0.80 and top3_overlap >= 1.0 and n >= 3
        else "WEAK_OR_BLOCKED"
    )
    if len(finite_signal) < max(2, math.ceil(0.5 * n)):
        status_gate1 = "WEAK_SIGNAL"

    gate2_status = "NETWORK_MODE_SCOPE"
    gate2_reason = "Critical mode is network-family; structural substitutability can be tested on retained network modes."
    if str(base_bridge_crit["family"]) == "II_control":
        gate2_status = "NOT_APPLICABLE_TO_CONTROL_CRITICAL"
        gate2_reason = (
            "Critical bridge mode is control-family. Projecting the network-inertia shunt residual "
            "onto this control-family mode would not test the original red-inertia classifier."
        )

    gate3_status = "NO_SIGN_REVERSAL_OBSERVED"
    if negative_full:
        gate3_status = "SIGN_REVERSAL_OBSERVED_ON_FULL_ANDES"
    elif negative_bridge:
        gate3_status = "BRIDGE_ONLY_NEGATIVE_NOT_CONFIRMED_BY_ANDES"

    summary = {
        "case": case_name,
        "case_path": str(case_path),
        "eps": eps,
        "n_candidates": n,
        "base_full_critical": {
            k: base_full_crit[k]
            for k in ["real", "imag", "freq_hz", "zeta", "pi_c", "family"]
        },
        "base_bridge_critical": {
            k: base_bridge_crit[k]
            for k in ["real", "imag", "freq_hz", "zeta", "pi_c", "family"]
        },
        "virtual_inertia_columns": virtual_cols,
        "virtual_inertia_validation_status": (
            "BLOCKED_NO_EXPLICIT_GFM_VIRTUAL_INERTIA_COLUMN"
            if not any("M" in cols or "H" in cols or "Tj" in cols for cols in virtual_cols.values())
            else "CANDIDATE_COLUMNS_PRESENT_REQUIRES_PHYSICAL_MAPPING"
        ),
        "gate1_bridge_vs_full": {
            "status": status_gate1,
            "sign_rate_global": sign_rate_global,
            "sign_rate_matched": sign_rate_matched,
            "top3_overlap": top3_overlap,
            "median_mag_rel_error_global": float(np.median([r["mag_rel_error_global"] for r in rows])),
            "median_mag_rel_error_matched": float(np.median([r["mag_rel_error_matched"] for r in rows])),
            "finite_signal_count": len(finite_signal),
            "acceptance_rule": "PASS if global sign rate >= 80%, top-3 overlap is complete, and signal is not numerically weak.",
        },
        "gate1_local_nep_sensitivity": {
            "status": status_gate1,
            "sign_rate_vs_full_matched": sign_rate_nep_local_full,
            "sign_rate_vs_bridge_matched": sign_rate_nep_local_bridge,
            "top3_overlap_vs_full_global": top3_overlap_nep_local,
            "median_mag_rel_error_vs_full_matched": float(
                np.median([r["mag_rel_error_nep_local_full"] for r in rows])
            ),
            "median_mag_rel_error_vs_bridge_matched": float(
                np.median([r["mag_rel_error_nep_local_bridge"] for r in rows])
            ),
            "note": (
                "Local NEP derivative uses ds/deps=-y^H T_p x/(y^H T_s x), "
                "with T_p estimated from the Schur-operator perturbation at the base pole."
            ),
        },
        "gate2_substitutability": {"status": gate2_status, "reason": gate2_reason},
        "gate3_inertia_reversal": {
            "status": gate3_status,
            "n_negative_full_global": len(negative_full),
            "n_negative_bridge_global": len(negative_bridge),
            "negative_full_indices": [r["idx"] for r in negative_full],
        },
    }
    write_csv(case_out / "candidate_results.csv", rows)
    write_json(case_out / "summary.json", summary)
    return {"summary": summary, "rows": rows}


def write_report(status: dict[str, Any], case_results: list[dict[str, Any]]) -> None:
    lines = [
        "# Phase Gate: Inertia Assignment over Controller-Aware Bridge",
        "",
        "This gate asks whether inertia-placement diagnostics survive when evaluated over the",
        "controller-aware Schur NEP bridge.  The script perturbs only documented inertia",
        "parameters present in the case files.  In the current IEEE39 IBR cases, REGF1 does",
        "not expose a clear `M`, `H`, or `Tj` virtual-inertia parameter, so GFM virtual-inertia",
        "assignment remains blocked unless a physical parameter mapping is supplied.",
        "",
        f"- Commit: `{status['git_commit']}`",
        f"- Random seed: `{status['random_seed']}`",
        f"- eps: `{status['eps']}`",
        "",
        "## Case Results",
        "",
    ]
    for item in case_results:
        s = item["summary"]
        g1 = s["gate1_bridge_vs_full"]
        g1loc = s["gate1_local_nep_sensitivity"]
        g2 = s["gate2_substitutability"]
        g3 = s["gate3_inertia_reversal"]
        lines.extend(
            [
                f"### {s['case']}",
                "",
                (
                    f"- Critical full mode: zeta={s['base_full_critical']['zeta']:.6g}, "
                    f"freq={s['base_full_critical']['freq_hz']:.6g} Hz, "
                    f"family={s['base_full_critical']['family']}, "
                    f"pi_c={s['base_full_critical']['pi_c']:.4g}."
                ),
                (
                    f"- Critical bridge mode: zeta={s['base_bridge_critical']['zeta']:.6g}, "
                    f"freq={s['base_bridge_critical']['freq_hz']:.6g} Hz, "
                    f"family={s['base_bridge_critical']['family']}, "
                    f"pi_c={s['base_bridge_critical']['pi_c']:.4g}."
                ),
                (
                    f"- Gate-1 bridge-vs-ANDES: **{g1['status']}**; "
                    f"sign_rate={g1['sign_rate_global']:.3f}, "
                    f"top3_overlap={g1['top3_overlap']:.3f}, "
                    f"finite_signal={g1['finite_signal_count']}/{s['n_candidates']}."
                ),
                (
                    f"- Gate-1 local NEP sensitivity: sign_rate_vs_full_matched="
                    f"{g1loc['sign_rate_vs_full_matched']:.3f}, "
                    f"top3_overlap_vs_full_global={g1loc['top3_overlap_vs_full_global']:.3f}, "
                    f"median_rel_error_vs_full_matched="
                    f"{g1loc['median_mag_rel_error_vs_full_matched']:.3g}."
                ),
                f"- Gate-2 substitutability: **{g2['status']}**. {g2['reason']}",
                (
                    f"- Gate-3 inertia sign reversal: **{g3['status']}**; "
                    f"negative_full={g3['n_negative_full_global']}/{s['n_candidates']}."
                ),
                f"- Virtual-inertia validation: **{s['virtual_inertia_validation_status']}**.",
                "",
            ]
        )

    lines.extend(
        [
            "## Verdict",
            "",
            f"- Overall status: **{status['overall_status']}**",
            f"- Abstract implication: {status['abstract_implication']}",
            "",
            "## What this gate did not test",
            "",
            "- It did not validate a real GFM virtual-inertia knob because the tested REGF1 sheets do not expose one.",
            "- It did not prove that the red-inertia substitutability classifier applies to control-family critical modes.",
            "- It did not replace the registered estimator/planning validation over the NEP bridge.",
        ]
    )
    (OUT / "phase_inertia_bridge_gate_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)

    eps = EPS_DEFAULT
    catalog = case_catalog()
    case_results: list[dict[str, Any]] = []
    for case_name in ["no_pss", "base", "mix60"]:
        path = catalog[case_name]
        if not path.exists():
            case_results.append({"summary": {"case": case_name, "status": "MISSING_CASE", "case_path": str(path)}})
            continue
        print(f"[gate] running {case_name}: {path}")
        case_results.append(run_case_gate(case_name, path, eps))

    summaries = [item["summary"] for item in case_results]
    g1_passes = [s for s in summaries if s.get("gate1_bridge_vs_full", {}).get("status") == "PASS"]
    g2_control_blocked = [
        s for s in summaries if s.get("gate2_substitutability", {}).get("status") == "NOT_APPLICABLE_TO_CONTROL_CRITICAL"
    ]
    virtual_blocked = [
        s
        for s in summaries
        if s.get("virtual_inertia_validation_status") == "BLOCKED_NO_EXPLICIT_GFM_VIRTUAL_INERTIA_COLUMN"
    ]

    if virtual_blocked:
        overall = "CONSERVATIVE_SCOPE"
        abstract = (
            "Do not claim validated virtual-inertia assignment over the controller-aware bridge. "
            "Report GENROU.M sensitivity checks and keep GFM virtual inertia as blocked pending a physical parameter mapping."
        )
    elif len(g1_passes) == len(summaries) and not g2_control_blocked:
        overall = "STRONG_SCOPE_CANDIDATE"
        abstract = "Bridge-based inertia assignment appears compatible with the tested controller-aware cases."
    else:
        overall = "CONSERVATIVE_SCOPE"
        abstract = (
            "Use the conservative abstract: the bridge is validated, but inertia substitutability remains validated only in the network-inertia regime or in limited GENROU.M checks."
        )

    status = {
        "overall_status": overall,
        "abstract_implication": abstract,
        "random_seed": RNG_SEED,
        "git_commit": git_commit_id(),
        "eps": eps,
        "cases": summaries,
    }
    write_json(OUT / "phase_inertia_bridge_gate_status.json", status)
    write_report(status, case_results)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
