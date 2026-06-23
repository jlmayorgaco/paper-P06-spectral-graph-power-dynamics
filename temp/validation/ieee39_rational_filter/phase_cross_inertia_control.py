"""Phase E0: inertia-control cross-term gate over the Schur/ANDES bridge.

This script tests the empirical question raised by the spectral-attribution
framing: does a joint inertia/control perturbation move a mode in a way that is
not explained by adding the inertia-only and control-only shifts?

The test is deliberately local and conservative.

Given a base ANDES case, it builds perturbed Excel cases:

* M-only: scale selected active GENROU.M rows;
* C-only: scale selected active PLL1 gains Kp/Ki;
* M+C: apply both perturbations.

It then uses the full ANDES eigensolve as the reference.  The Schur partition is
used only to report control participation pi_c for the tracked mode.

The non-additive cross derivative is estimated at a half step h:

    s_MC = (s(M+h,C+h) - s(M+h,C) - s(M,C+h) + s0) / h^2

and then tested at a larger step eps:

    additive prediction: s0 + eps (s_M + s_C)
    cross prediction:    additive + eps^2 s_MC

If the cross prediction improves the pole and damping-ratio prediction over the
additive one, the cross term is measurable for that case/mode.  If not, the
framework remains a useful attribution organizer, but the cross term is not an
empirical headline for the tested benchmark.

No manuscript files are edited by this script.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from phase0d_nep_contour_solver import (
    block_partition,
    classify_partition,
    control_participation,
    damping_ratio,
    freq_hz,
    load_andes_case,
    null_vector,
    one_sided_modes,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase_cross_inertia_control"
TMP = OUT / "perturbed_cases"
RNG_SEED = 20260608


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


def parse_targets(value: str | None) -> set[str] | None:
    if value is None or value.strip().lower() in {"", "all", "*"}:
        return None
    return {x.strip() for x in value.split(",") if x.strip()}


def scale_selected(
    df: pd.DataFrame,
    column: str,
    scale: float,
    targets: set[str] | None,
    target_column: str = "idx",
) -> int:
    if column not in df.columns:
        return 0
    mask = active_mask(df)
    if targets is not None and target_column in df.columns:
        mask &= df[target_column].astype(str).isin(targets)
    count = int(mask.sum())
    if count:
        df[column] = pd.to_numeric(df[column], errors="coerce").astype(float)
        df.loc[mask, column] = df.loc[mask, column] * scale
    return count


@dataclass(frozen=True)
class PerturbationSpec:
    inertia_scale: float
    control_scale: float
    genrou_targets: set[str] | None
    pll_targets: set[str] | None


def make_perturbed_case(base: Path, out: Path, spec: PerturbationSpec) -> dict[str, Any]:
    sheets = read_case(base)
    changed: dict[str, Any] = {
        "base_case": str(base),
        "out_case": str(out),
        "inertia_scale": spec.inertia_scale,
        "control_scale": spec.control_scale,
        "genrou_targets": sorted(spec.genrou_targets) if spec.genrou_targets else "all_active",
        "pll_targets": sorted(spec.pll_targets) if spec.pll_targets else "all_active",
        "n_genrou_m_scaled": 0,
        "n_pll_kp_scaled": 0,
        "n_pll_ki_scaled": 0,
    }

    if spec.inertia_scale != 1.0 and "GENROU" in sheets:
        changed["n_genrou_m_scaled"] = scale_selected(
            sheets["GENROU"], "M", spec.inertia_scale, spec.genrou_targets
        )

    if spec.control_scale != 1.0 and "PLL1" in sheets:
        changed["n_pll_kp_scaled"] = scale_selected(
            sheets["PLL1"], "Kp", spec.control_scale, spec.pll_targets
        )
        changed["n_pll_ki_scaled"] = scale_selected(
            sheets["PLL1"], "Ki", spec.control_scale, spec.pll_targets
        )

    write_case(out, sheets)
    return changed


def mode_table(case_path: Path, max_modes: int = 24, tau: float = 0.70) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ss = load_andes_case(str(case_path))
    names = [str(x) for x in ss.EIG.x_name]
    A = np.asarray(ss.EIG.As, dtype=float)
    e_idx, c_idx, rationale = classify_partition(names)
    blocks = block_partition(A, e_idx, c_idx)

    rows: list[dict[str, Any]] = []
    for rank, s in enumerate(one_sided_modes(np.asarray(ss.EIG.mu, dtype=complex))[:max_modes], start=1):
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
                "family": "II_control" if pi_c >= tau else "I_network",
            }
        )
    meta = {
        "case_path": str(case_path),
        "n_states": len(names),
        "n_retained": len(e_idx),
        "n_condensed": len(c_idx),
        "partition_rationale": rationale,
    }
    return rows, meta


def choose_mode(rows: list[dict[str, Any]], selector: str) -> dict[str, Any]:
    if not rows:
        raise RuntimeError("No oscillatory modes found in base case.")
    if selector == "critical":
        return rows[0]
    if selector == "network":
        for row in rows:
            if row["family"] == "I_network":
                return row
        raise RuntimeError("No network-family oscillatory mode found.")
    if selector == "control":
        for row in rows:
            if row["family"] == "II_control":
                return row
        raise RuntimeError("No control-family oscillatory mode found.")
    raise ValueError(f"Unknown selector: {selector}")


def match_mode(rows: list[dict[str, Any]], s_ref: complex) -> dict[str, Any]:
    if not rows:
        raise RuntimeError("No oscillatory modes found in perturbed case.")
    return min(rows, key=lambda row: abs(row["pole"] - s_ref))


def zeta_error(pred: complex, true: complex) -> float:
    return abs(damping_ratio(pred) - damping_ratio(true))


def complex_error(pred: complex, true: complex) -> float:
    return abs(pred - true)


def run_gate(
    base_case: Path,
    mode_selector: str,
    eps: float,
    genrou_targets: set[str] | None,
    pll_targets: set[str] | None,
    tau: float,
) -> dict[str, Any]:
    h = eps / 2.0
    TMP.mkdir(parents=True, exist_ok=True)

    cases = {
        "base": (1.0, 1.0),
        "M_h": (1.0 + h, 1.0),
        "C_h": (1.0, 1.0 + h),
        "MC_h": (1.0 + h, 1.0 + h),
        "MC_eps": (1.0 + eps, 1.0 + eps),
    }

    perturbation_rows = []
    case_paths: dict[str, Path] = {"base": base_case}
    stem = base_case.stem
    for label, (m_scale, c_scale) in cases.items():
        if label == "base":
            continue
        out = TMP / f"{stem}_{mode_selector}_{label}_eps{eps:.4g}.xlsx"
        info = make_perturbed_case(
            base_case,
            out,
            PerturbationSpec(m_scale, c_scale, genrou_targets, pll_targets),
        )
        info["label"] = label
        perturbation_rows.append(info)
        case_paths[label] = out

    tables: dict[str, list[dict[str, Any]]] = {}
    metas: dict[str, dict[str, Any]] = {}
    for label, path in case_paths.items():
        rows, meta = mode_table(path, tau=tau)
        tables[label] = rows
        metas[label] = meta

    base_mode = choose_mode(tables["base"], mode_selector)
    s0 = complex(base_mode["pole"])
    matched = {label: match_mode(rows, s0) for label, rows in tables.items()}

    s_m_h = complex(matched["M_h"]["pole"])
    s_c_h = complex(matched["C_h"]["pole"])
    s_mc_h = complex(matched["MC_h"]["pole"])
    s_mc_eps = complex(matched["MC_eps"]["pole"])

    d_m = (s_m_h - s0) / h
    d_c = (s_c_h - s0) / h
    d_cross = (s_mc_h - s_m_h - s_c_h + s0) / (h * h)

    pred_add = s0 + eps * (d_m + d_c)
    pred_cross = pred_add + eps * eps * d_cross

    add_complex_error = complex_error(pred_add, s_mc_eps)
    cross_complex_error = complex_error(pred_cross, s_mc_eps)
    add_zeta_error = zeta_error(pred_add, s_mc_eps)
    cross_zeta_error = zeta_error(pred_cross, s_mc_eps)
    true_shift = s_mc_eps - s0
    cross_shift = eps * eps * d_cross

    improvement_complex = (
        (add_complex_error - cross_complex_error) / add_complex_error
        if add_complex_error > 1e-15
        else 0.0
    )
    improvement_zeta = (
        (add_zeta_error - cross_zeta_error) / add_zeta_error if add_zeta_error > 1e-15 else 0.0
    )
    cross_fraction = abs(cross_shift) / max(abs(true_shift), 1e-15)
    measurable = bool(
        cross_fraction >= 0.05
        and (improvement_complex >= 0.10 or improvement_zeta >= 0.10)
        and cross_complex_error < add_complex_error
    )

    mode_rows = []
    for label, row in matched.items():
        mode_rows.append(
            {
                "label": label,
                "real": row["real"],
                "imag": row["imag"],
                "freq_hz": row["freq_hz"],
                "zeta": row["zeta"],
                "pi_c": row["pi_c"],
                "family": row["family"],
                "distance_to_base_pole": abs(complex(row["pole"]) - s0),
            }
        )

    prediction_rows = [
        {
            "prediction": "additive",
            "real": float(pred_add.real),
            "imag": float(pred_add.imag),
            "freq_hz": freq_hz(pred_add),
            "zeta": damping_ratio(pred_add),
            "complex_error": add_complex_error,
            "zeta_abs_error": add_zeta_error,
        },
        {
            "prediction": "additive_plus_cross",
            "real": float(pred_cross.real),
            "imag": float(pred_cross.imag),
            "freq_hz": freq_hz(pred_cross),
            "zeta": damping_ratio(pred_cross),
            "complex_error": cross_complex_error,
            "zeta_abs_error": cross_zeta_error,
        },
    ]

    return {
        "status": "PASS_MEASURABLE" if measurable else "NOT_MEASURABLE",
        "random_seed": RNG_SEED,
        "git_commit": git_commit_id(),
        "base_case": str(base_case),
        "mode_selector": mode_selector,
        "eps": eps,
        "half_step": h,
        "tau": tau,
        "base_mode": {
            "real": base_mode["real"],
            "imag": base_mode["imag"],
            "freq_hz": base_mode["freq_hz"],
            "zeta": base_mode["zeta"],
            "pi_c": base_mode["pi_c"],
            "family": base_mode["family"],
        },
        "perturbations": perturbation_rows,
        "mode_rows": mode_rows,
        "prediction_rows": prediction_rows,
        "cross_metrics": {
            "true_shift_abs": abs(true_shift),
            "cross_shift_abs": abs(cross_shift),
            "cross_fraction_of_true_shift": cross_fraction,
            "add_complex_error": add_complex_error,
            "cross_complex_error": cross_complex_error,
            "complex_error_improvement_fraction": improvement_complex,
            "add_zeta_abs_error": add_zeta_error,
            "cross_zeta_abs_error": cross_zeta_error,
            "zeta_error_improvement_fraction": improvement_zeta,
            "d_m_real": float(d_m.real),
            "d_m_imag": float(d_m.imag),
            "d_c_real": float(d_c.real),
            "d_c_imag": float(d_c.imag),
            "d_cross_real": float(d_cross.real),
            "d_cross_imag": float(d_cross.imag),
        },
        "partition_meta": metas["base"],
    }


def active_ids(case_path: Path, sheet_name: str) -> list[str]:
    sheets = read_case(case_path)
    if sheet_name not in sheets or "idx" not in sheets[sheet_name].columns:
        return []
    df = sheets[sheet_name]
    ids = df.loc[active_mask(df), "idx"].astype(str).tolist()
    return ids


def summarize_result(label: str, genrou_target: str, pll_target: str, result: dict[str, Any]) -> dict[str, Any]:
    cm = result["cross_metrics"]
    return {
        "label": label,
        "genrou_target": genrou_target,
        "pll_target": pll_target,
        "status": result["status"],
        "mode_selector": result["mode_selector"],
        "base_real": result["base_mode"]["real"],
        "base_imag": result["base_mode"]["imag"],
        "base_freq_hz": result["base_mode"]["freq_hz"],
        "base_zeta": result["base_mode"]["zeta"],
        "base_pi_c": result["base_mode"]["pi_c"],
        "base_family": result["base_mode"]["family"],
        "true_shift_abs": cm["true_shift_abs"],
        "cross_shift_abs": cm["cross_shift_abs"],
        "cross_fraction_of_true_shift": cm["cross_fraction_of_true_shift"],
        "add_complex_error": cm["add_complex_error"],
        "cross_complex_error": cm["cross_complex_error"],
        "complex_error_improvement_fraction": cm["complex_error_improvement_fraction"],
        "add_zeta_abs_error": cm["add_zeta_abs_error"],
        "cross_zeta_abs_error": cm["cross_zeta_abs_error"],
        "zeta_error_improvement_fraction": cm["zeta_error_improvement_fraction"],
    }


def write_report(result: dict[str, Any], out: Path) -> None:
    cm = result["cross_metrics"]
    lines = [
        "# Phase E0 Inertia-Control Cross-Term Gate",
        "",
        "This report tests whether the joint inertia/control perturbation is measurably non-additive.",
        "The full ANDES eigensolve is the pole reference.  The Schur partition is used only to",
        "report control participation of the tracked mode.",
        "",
        f"- Status: **{result['status']}**",
        f"- Base case: `{result['base_case']}`",
        f"- Mode selector: `{result['mode_selector']}`",
        f"- eps: `{result['eps']}`, half-step h: `{result['half_step']}`",
        f"- Commit: `{result['git_commit']}`",
        "",
        "## Base Tracked Mode",
        "",
        (
            f"- pole = {result['base_mode']['real']:.6g} + j{result['base_mode']['imag']:.6g}, "
            f"freq = {result['base_mode']['freq_hz']:.6g} Hz, "
            f"zeta = {result['base_mode']['zeta']:.6g}, "
            f"pi_c = {result['base_mode']['pi_c']:.4g}, "
            f"family = {result['base_mode']['family']}"
        ),
        "",
        "## Cross-Term Metrics",
        "",
        f"- |true joint shift| = `{cm['true_shift_abs']:.6e}`",
        f"- |estimated cross shift at eps| = `{cm['cross_shift_abs']:.6e}`",
        f"- cross/true shift fraction = `{cm['cross_fraction_of_true_shift']:.6g}`",
        f"- additive complex error = `{cm['add_complex_error']:.6e}`",
        f"- cross complex error = `{cm['cross_complex_error']:.6e}`",
        f"- complex-error improvement = `{cm['complex_error_improvement_fraction']:.6g}`",
        f"- additive zeta error = `{cm['add_zeta_abs_error']:.6e}`",
        f"- cross zeta error = `{cm['cross_zeta_abs_error']:.6e}`",
        f"- zeta-error improvement = `{cm['zeta_error_improvement_fraction']:.6g}`",
        "",
        "## Interpretation",
        "",
    ]
    if result["status"] == "PASS_MEASURABLE":
        lines.append(
            "The cross term is measurable under the registered local test: it is a nontrivial "
            "fraction of the true shift and improves the additive prediction."
        )
    else:
        lines.append(
            "The cross term is not measurable under this registered local test.  This does not "
            "invalidate the attribution framework, but it means the inertia-control cross term "
            "should not be used as a headline empirical claim for this benchmark without a "
            "stronger case."
        )

    out.mkdir(parents=True, exist_ok=True)
    (out / "phase_cross_inertia_control_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        type=Path,
        default=ROOT
        / "validation"
        / "ieee39_rational_filter"
        / "cases"
        / "phase1_calibrated"
        / "no_pss"
        / "ieee39_ibr_mix60.xlsx",
    )
    parser.add_argument("--mode", choices=["critical", "network", "control"], default="network")
    parser.add_argument("--eps", type=float, default=0.01)
    parser.add_argument("--genrou-targets", default="all", help="Comma-separated GENROU idx values, or all.")
    parser.add_argument("--pll-targets", default="all", help="Comma-separated PLL1 idx values, or all.")
    parser.add_argument("--tau", type=float, default=0.70)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument(
        "--sweep-singletons",
        action="store_true",
        help="Run all active single GENROU.M x single PLL1.Kp/Ki target pairs plus the all-active case.",
    )
    args = parser.parse_args()

    if not args.case.exists():
        raise FileNotFoundError(args.case)
    if args.eps <= 0:
        raise ValueError("--eps must be positive.")

    if args.out.exists():
        # Keep the destructive action tightly scoped to this script's output.
        shutil.rmtree(args.out)
    args.out.mkdir(parents=True, exist_ok=True)

    if args.sweep_singletons:
        rows = []
        genrou_ids = active_ids(args.case, "GENROU")
        pll_ids = active_ids(args.case, "PLL1")
        sweep_items: list[tuple[str, set[str] | None, str, set[str] | None]] = [
            ("all_active", None, "all_active", None)
        ]
        for g in genrou_ids:
            for p in pll_ids:
                sweep_items.append((g, {g}, p, {p}))

        results = []
        for i, (g_label, g_targets, p_label, p_targets) in enumerate(sweep_items, start=1):
            result = run_gate(
                base_case=args.case,
                mode_selector=args.mode,
                eps=args.eps,
                genrou_targets=g_targets,
                pll_targets=p_targets,
                tau=args.tau,
            )
            result["sweep_index"] = i
            result["genrou_target_label"] = g_label
            result["pll_target_label"] = p_label
            results.append(result)
            rows.append(summarize_result(f"{g_label}__{p_label}", g_label, p_label, result))

        write_json(args.out / "phase_cross_inertia_control_sweep.json", results)
        write_csv(args.out / "sweep_summary.csv", rows)
        best = max(rows, key=lambda r: (r["cross_fraction_of_true_shift"], r["complex_error_improvement_fraction"]))
        status = {
            "status": "PASS_MEASURABLE" if any(r["status"] == "PASS_MEASURABLE" for r in rows) else "NOT_MEASURABLE",
            "n_cases": len(rows),
            "n_measurable": sum(1 for r in rows if r["status"] == "PASS_MEASURABLE"),
            "best_by_cross_fraction": best,
        }
        write_json(args.out / "phase_cross_inertia_control_status.json", status)
        print(json.dumps(status, indent=2))
        return

    result = run_gate(
        base_case=args.case,
        mode_selector=args.mode,
        eps=args.eps,
        genrou_targets=parse_targets(args.genrou_targets),
        pll_targets=parse_targets(args.pll_targets),
        tau=args.tau,
    )

    write_json(args.out / "phase_cross_inertia_control_status.json", result)
    write_csv(args.out / "tracked_modes.csv", result["mode_rows"])
    write_csv(args.out / "prediction_errors.csv", result["prediction_rows"])
    write_csv(args.out / "perturbations.csv", result["perturbations"])
    write_report(result, args.out)

    print(json.dumps({"status": result["status"], "cross_metrics": result["cross_metrics"]}, indent=2))


if __name__ == "__main__":
    main()
