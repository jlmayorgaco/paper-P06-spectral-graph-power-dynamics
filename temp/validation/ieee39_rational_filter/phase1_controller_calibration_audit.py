"""
Phase-1 controller calibration audit for the derived IEEE39 IBR cases.

The goal is conservative: sweep a small set of documented GFL/GFM controller
profiles, run full ANDES PFlow/EIG, record stability and participation of the
dominant/critical modes, and make clear whether calibration is sufficient to
unlock estimator validation. This is not an optimizer and it does not tune
parameters to favor the proposed estimator.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE_CASE_DIR = ROOT / "validation" / "ieee39_rational_filter" / "cases"
OUT = ROOT / "outputs" / "ieee39_phase1_calibration"
OUT_CASES = BASE_CASE_DIR / "phase1_calibrated"


@dataclass(frozen=True)
class Profile:
    name: str
    description: str
    apply: Callable[[dict[str, pd.DataFrame]], None]


def damping_ratio(s: complex) -> float:
    return float("inf") if abs(s) < 1e-12 else float(-s.real / abs(s))


def freq_hz(s: complex) -> float:
    return float(abs(s.imag) / (2 * np.pi)) if np.isfinite(s.imag) else float("nan")


def set_cols(df: pd.DataFrame, values: dict[str, float]) -> None:
    for col, val in values.items():
        if col in df.columns:
            df[col] = val


def scale_cols(df: pd.DataFrame, values: dict[str, float]) -> None:
    for col, scale in values.items():
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce") * scale


def profile_baseline(sheets: dict[str, pd.DataFrame]) -> None:
    return None


def profile_slow_pll(sheets: dict[str, pd.DataFrame]) -> None:
    if "PLL1" in sheets:
        set_cols(sheets["PLL1"], {"Kp": 0.35, "Ki": 0.035, "Tf": 0.10, "Tp": 0.10})


def profile_soft_gfl(sheets: dict[str, pd.DataFrame]) -> None:
    if "PLL1" in sheets:
        set_cols(sheets["PLL1"], {"Kp": 0.35, "Ki": 0.035, "Tf": 0.10, "Tp": 0.10})
    if "REECA1" in sheets:
        set_cols(sheets["REECA1"], {"Kqv": 2.0, "Kqp": 0.2, "Kqi": 0.02, "Kvp": 0.2, "Kvi": 0.02})
        set_cols(sheets["REECA1"], {"Trv": 0.05, "Tp": 0.05, "Tpord": 0.05, "Tiq": 0.05})
    if "REPCA1" in sheets:
        set_cols(sheets["REPCA1"], {"Kp": 0.2, "Ki": 0.02, "Kpg": 0.2, "Kig": 0.02, "Tfltr": 0.05})


def profile_soft_gfm(sheets: dict[str, pd.DataFrame]) -> None:
    if "REGF1" in sheets:
        set_cols(
            sheets["REGF1"],
            {
                "KPi": 0.15,
                "KIi": 4.0,
                "KPv": 0.8,
                "KIv": 3.0,
                "KPplim": 1.0,
                "KIplim": 8.0,
                "KPqlim": 0.05,
                "KIqlim": 0.6,
                "Tpm": 0.08,
                "Tr": 0.02,
                "Te": 0.02,
                "wdrp": 0.05,
                "Qdrp": 0.06,
            },
        )


def profile_soft_all(sheets: dict[str, pd.DataFrame]) -> None:
    profile_soft_gfl(sheets)
    profile_soft_gfm(sheets)


def profile_very_soft_all(sheets: dict[str, pd.DataFrame]) -> None:
    profile_soft_all(sheets)
    if "PLL1" in sheets:
        set_cols(sheets["PLL1"], {"Kp": 0.15, "Ki": 0.01, "Tf": 0.15, "Tp": 0.15})
    if "REECA1" in sheets:
        set_cols(sheets["REECA1"], {"Kqv": 0.8, "Kqp": 0.08, "Kqi": 0.008, "Kvp": 0.08, "Kvi": 0.008})
    if "REPCA1" in sheets:
        set_cols(sheets["REPCA1"], {"Kp": 0.08, "Ki": 0.008, "Kpg": 0.08, "Kig": 0.008})
    if "REGF1" in sheets:
        set_cols(sheets["REGF1"], {"KPi": 0.08, "KIi": 1.5, "KPv": 0.4, "KIv": 1.2, "Tpm": 0.12})


def profile_droop_damped_gfm(sheets: dict[str, pd.DataFrame]) -> None:
    if "REGF1" in sheets:
        set_cols(sheets["REGF1"], {"wdrp": 0.08, "Qdrp": 0.09, "Tpm": 0.12, "Tr": 0.03, "Te": 0.03})
        scale_cols(sheets["REGF1"], {"KPi": 0.4, "KIi": 0.4, "KPv": 0.4, "KIv": 0.4})


def profile_no_pss(sheets: dict[str, pd.DataFrame]) -> None:
    sheets.pop("IEEEST", None)


def profile_no_exciter_pss(sheets: dict[str, pd.DataFrame]) -> None:
    sheets.pop("IEEEST", None)
    sheets.pop("IEEEX1", None)


def profile_no_pss_soft_all(sheets: dict[str, pd.DataFrame]) -> None:
    profile_no_pss(sheets)
    profile_soft_all(sheets)


PROFILES = [
    Profile("baseline", "Generated generic parameters from Phase 0.", profile_baseline),
    Profile("slow_pll", "GFL PLL gains reduced and PLL filters slowed.", profile_slow_pll),
    Profile("soft_gfl", "Slow PLL plus softer REECA1/REPCA1 voltage/reactive controls.", profile_soft_gfl),
    Profile("soft_gfm", "Softer REGF1 inner controls and slower droop response.", profile_soft_gfm),
    Profile("soft_all", "soft_gfl and soft_gfm combined.", profile_soft_all),
    Profile("very_soft_all", "More conservative version of soft_all.", profile_very_soft_all),
    Profile("droop_damped_gfm", "Softer GFM droop and scaled REGF1 gains only.", profile_droop_damped_gfm),
    Profile("no_pss", "Diagnostic removal of retained synchronous-machine IEEEST/PSS rows.", profile_no_pss),
    Profile("no_exciter_pss", "Diagnostic removal of retained IEEEX1 and IEEEST rows.", profile_no_exciter_pss),
    Profile("no_pss_soft_all", "Diagnostic no_pss plus soft_all inverter controls.", profile_no_pss_soft_all),
]


def read_case(path: Path) -> dict[str, pd.DataFrame]:
    return pd.read_excel(path, sheet_name=None)


def write_case(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)


def dominant_family(state_name: str) -> str:
    families = [
        "PLL1",
        "REGCP1",
        "REECA1",
        "REPCA1",
        "REGF1",
        "GENROU",
        "TGOV1N",
        "IEEEX1",
        "IEEEST",
        "BusFreq",
    ]
    for fam in families:
        if fam in state_name:
            return fam
    return "unknown"


def mode_kind(freq: float, family: str, state_name: str) -> str:
    if family in {"PLL1", "REGCP1", "REECA1", "REPCA1", "REGF1"}:
        return "ibr-control"
    if family in {"IEEEX1", "IEEEST", "TGOV1N"}:
        return "sg-control"
    if family == "GENROU" and ("omega" in state_name or "delta" in state_name):
        return "electromechanical"
    if freq >= 5.0:
        return "high-frequency-control"
    if 0.1 <= freq <= 3.0:
        return "electromechanical-or-control"
    return "unclassified"


def top_participation(ss, mode_idx: int, top_n: int = 8) -> list[dict[str, object]]:
    pf = np.asarray(ss.EIG.pfactors)
    names = [str(x) for x in ss.EIG.x_name]
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    vals = np.abs(pf[:, mode_idx]) if pf.shape[1] == len(mu) else np.abs(pf[mode_idx, :])
    rows = []
    for rank, state_idx in enumerate(np.argsort(vals)[-top_n:][::-1], start=1):
        state_name = names[state_idx] if state_idx < len(names) else str(state_idx)
        rows.append(
            {
                "rank": rank,
                "state_index": int(state_idx),
                "state_name": state_name,
                "family": dominant_family(state_name),
                "participation": float(vals[state_idx]),
            }
        )
    return rows


def summarize_case(case_path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    import andes  # type: ignore

    result: dict[str, object] = {"case_path": str(case_path), "ok": False}
    part_rows: list[dict[str, object]] = []
    try:
        ss = andes.load(str(case_path), setup=False, no_output=True)
        ss.setup()
        pflow_ok = bool(ss.PFlow.run())
        eig_ok = bool(ss.EIG.run()) if pflow_ok else False
        result.update({"pflow_ok": pflow_ok, "eig_ok": eig_ok, "ok": bool(pflow_ok and eig_ok)})
        if not eig_ok:
            return result, part_rows

        mu = np.asarray(ss.EIG.mu, dtype=complex)
        pos_idx = np.where(mu.real > 1e-7)[0]
        stable_osc_idx = np.where((mu.real < 0) & (np.abs(mu.imag) > 0.1))[0]
        dom_pos_idx = int(pos_idx[np.argmax(mu[pos_idx].real)]) if len(pos_idx) else None
        if len(stable_osc_idx):
            zetas = np.array([damping_ratio(mu[i]) for i in stable_osc_idx])
            crit_idx = int(stable_osc_idx[int(np.argmin(zetas))])
        else:
            crit_idx = None

        def mode_record(idx: int | None, prefix: str) -> dict[str, object]:
            if idx is None:
                return {
                    f"{prefix}_index": math.nan,
                    f"{prefix}_real": math.nan,
                    f"{prefix}_imag": math.nan,
                    f"{prefix}_zeta": math.nan,
                    f"{prefix}_freq_hz": math.nan,
                    f"{prefix}_top_state": "",
                    f"{prefix}_top_family": "",
                    f"{prefix}_kind": "",
                }
            s = mu[idx]
            top = top_participation(ss, idx, top_n=1)[0]
            f_hz = freq_hz(s)
            return {
                f"{prefix}_index": idx,
                f"{prefix}_real": float(s.real),
                f"{prefix}_imag": float(s.imag),
                f"{prefix}_zeta": damping_ratio(s),
                f"{prefix}_freq_hz": f_hz,
                f"{prefix}_top_state": top["state_name"],
                f"{prefix}_top_family": top["family"],
                f"{prefix}_kind": mode_kind(f_hz, str(top["family"]), str(top["state_name"])),
            }

        result.update(
            {
                "n_eigs": int(len(mu)),
                "n_positive_real": int(len(pos_idx)),
                "n_positive_oscillatory": int(np.sum(np.abs(mu[pos_idx].imag) > 0.1)) if len(pos_idx) else 0,
                "n_zero": int(np.sum(np.abs(mu) < 1e-7)),
                "max_real": float(np.max(mu.real)) if len(mu) else math.nan,
                "small_signal_stable": bool(len(pos_idx) == 0),
            }
        )
        result.update(mode_record(dom_pos_idx, "dominant_positive"))
        result.update(mode_record(crit_idx, "critical_stable_osc"))

        for label, idx in (("dominant_positive", dom_pos_idx), ("critical_stable_osc", crit_idx)):
            if idx is None:
                continue
            for row in top_participation(ss, idx, top_n=8):
                row.update(
                    {
                        "mode_label": label,
                        "mode_index": idx,
                        "mode_real": float(mu[idx].real),
                        "mode_imag": float(mu[idx].imag),
                        "mode_freq_hz": freq_hz(mu[idx]),
                    }
                )
                part_rows.append(row)
    except Exception as exc:
        result.update({"error_type": type(exc).__name__, "error": str(exc)})
    return result, part_rows


def run_profile(case_path: Path, profile: Profile) -> tuple[Path, dict[str, object], list[dict[str, object]]]:
    sheets = read_case(case_path)
    profile.apply(sheets)
    out_case = OUT_CASES / profile.name / case_path.name
    write_case(out_case, sheets)
    summary, parts = summarize_case(out_case)
    summary.update({"base_case": case_path.stem, "profile": profile.name, "profile_description": profile.description})
    for row in parts:
        row.update({"base_case": case_path.stem, "profile": profile.name})
    return out_case, summary, parts


def write_report(summary_rows: list[dict[str, object]], path: Path) -> None:
    df = pd.DataFrame(summary_rows)
    stable = df[df["small_signal_stable"] == True]
    best = best_by_case(df)

    lines = [
        "# Phase-1 Controller Calibration Audit",
        "",
        "This audit sweeps documented controller profiles for the derived IEEE39 IBR cases. "
        "The reference is the full ANDES eigensolve for each modified case.",
        "",
        f"- Profiles tested: {df['profile'].nunique()}",
        f"- Case/profile runs: {len(df)}",
        f"- Stable case/profile runs: {len(stable)}",
        "",
        "## Best profile by case",
    ]
    for _, row in best.iterrows():
        lines.append(
            "- "
            f"{row['base_case']}: controller-only {row['best_controller_profile']} "
            f"(pos={row['controller_pos']}, max_real={row['controller_max_real']:.4g}); "
            f"diagnostic {row['best_diagnostic_profile']} "
            f"(pos={row['diagnostic_pos']}, max_real={row['diagnostic_max_real']:.4g}); "
            f"critical top={row['diagnostic_critical_top_state']}"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "- A profile that makes a case stable is a calibration candidate, not final tuning.",
            "- High-frequency critical modes are classified by participation, not by frequency alone.",
            "- Estimator validation remains blocked until a calibrated full-ANDES IBR case is matched by a reduced `(L,M,D)`/control-coupling surrogate.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def best_by_case(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for case, group in df.groupby("base_case"):
        controller = group[~group["profile"].astype(str).str.startswith("no_")]
        diagnostic = group[group["profile"].astype(str).str.startswith("no_")]
        controller_best = controller.sort_values(["n_positive_real", "max_real"], ascending=[True, True]).iloc[0]
        diagnostic_best = diagnostic.sort_values(["n_positive_real", "max_real"], ascending=[True, True]).iloc[0]
        rows.append(
            {
                "base_case": case,
                "best_controller_profile": controller_best["profile"],
                "controller_stable": bool(controller_best["small_signal_stable"]),
                "controller_pos": int(controller_best["n_positive_real"]),
                "controller_max_real": float(controller_best["max_real"]),
                "controller_critical_zeta": float(controller_best["critical_stable_osc_zeta"]),
                "controller_critical_freq_hz": float(controller_best["critical_stable_osc_freq_hz"]),
                "controller_critical_top_state": controller_best["critical_stable_osc_top_state"],
                "controller_critical_kind": controller_best["critical_stable_osc_kind"],
                "best_diagnostic_profile": diagnostic_best["profile"],
                "diagnostic_stable": bool(diagnostic_best["small_signal_stable"]),
                "diagnostic_pos": int(diagnostic_best["n_positive_real"]),
                "diagnostic_max_real": float(diagnostic_best["max_real"]),
                "diagnostic_critical_zeta": float(diagnostic_best["critical_stable_osc_zeta"]),
                "diagnostic_critical_freq_hz": float(diagnostic_best["critical_stable_osc_freq_hz"]),
                "diagnostic_critical_top_state": diagnostic_best["critical_stable_osc_top_state"],
                "diagnostic_critical_kind": diagnostic_best["critical_stable_osc_kind"],
            }
        )
    return pd.DataFrame(rows)


def profile_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for profile, group in df.groupby("profile"):
        rows.append(
            {
                "profile": profile,
                "n_runs": int(len(group)),
                "n_stable": int((group["small_signal_stable"] == True).sum()),
                "median_positive_modes": float(group["n_positive_real"].median()),
                "max_positive_modes": int(group["n_positive_real"].max()),
                "min_max_real": float(group["max_real"].min()),
                "max_max_real": float(group["max_real"].max()),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_stable", "min_max_real"], ascending=[False, True])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    OUT_CASES.mkdir(parents=True, exist_ok=True)
    case_paths = sorted(BASE_CASE_DIR.glob("ieee39_ibr_*.xlsx"))

    summary_rows: list[dict[str, object]] = []
    participation_rows: list[dict[str, object]] = []
    for case_path in case_paths:
        for profile in PROFILES:
            _, summary, parts = run_profile(case_path, profile)
            summary_rows.append(summary)
            participation_rows.extend(parts)
            print(
                f"{case_path.stem:18s} {profile.name:16s} "
                f"ok={summary.get('ok')} stable={summary.get('small_signal_stable')} "
                f"pos={summary.get('n_positive_real')} max={summary.get('max_real')}"
            )

    pd.DataFrame(summary_rows).to_csv(OUT / "phase1_calibration_sweep.csv", index=False)
    pd.DataFrame(participation_rows).to_csv(OUT / "phase1_mode_participation.csv", index=False)
    best_by_case(pd.DataFrame(summary_rows)).to_csv(OUT / "phase1_best_by_case.csv", index=False)
    profile_summary(pd.DataFrame(summary_rows)).to_csv(OUT / "phase1_profile_summary.csv", index=False)
    with (OUT / "phase1_calibration_sweep.json").open("w", encoding="utf-8") as f:
        json.dump(summary_rows, f, indent=2)
    write_report(summary_rows, OUT / "phase1_calibration_report.md")
    print(f"Wrote {OUT / 'phase1_calibration_sweep.csv'}")
    print(f"Wrote {OUT / 'phase1_mode_participation.csv'}")
    print(f"Wrote {OUT / 'phase1_calibration_report.md'}")


if __name__ == "__main__":
    main()
