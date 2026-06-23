"""Diagnose positive-real modes in the packaged ANDES IEEE 39 dynamic case.

The ANDES ``ieee39_full.xlsx`` case has positive-real non-oscillatory modes in
the full eigenanalysis. This script isolates whether those modes are associated
with network electromechanical dynamics or with specific dynamic model families.
It also records dominant participation factors for the positive modes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def damping_ratio(s: complex) -> float:
    return float("inf") if abs(s) < 1e-12 else float(-s.real / abs(s))


def summarize_system(ss, imag_tol: float = 0.1) -> dict[str, object]:
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    pos_idx = np.where(mu.real > 1e-7)[0]
    pos = mu[pos_idx]
    pos_osc = pos[np.abs(pos.imag) > imag_tol]
    stable_osc = [z for z in mu if z.real < 0 and abs(z.imag) > imag_tol]
    if stable_osc:
        ratios = np.array([damping_ratio(z) for z in stable_osc])
        crit = stable_osc[int(np.argmin(ratios))]
    else:
        crit = complex(np.nan, np.nan)
    return {
        "n_eigs": int(len(mu)),
        "n_positive": int(len(pos)),
        "n_positive_osc": int(len(pos_osc)),
        "max_real": float(mu.real.max()) if len(mu) else float("nan"),
        "critical_real": float(crit.real),
        "critical_imag": float(crit.imag),
        "critical_zeta": damping_ratio(crit),
        "critical_freq_hz": float(abs(crit.imag) / (2 * np.pi)) if np.isfinite(crit.imag) else float("nan"),
    }


def participation_report(ss, top_n: int = 10) -> list[dict[str, object]]:
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    pf = np.asarray(ss.EIG.pfactors)
    names = [str(x) for x in ss.EIG.x_name]
    rows: list[dict[str, object]] = []
    for mode_idx in np.where(mu.real > 1e-7)[0]:
        vals = np.abs(pf[:, mode_idx]) if pf.shape[1] == len(mu) else np.abs(pf[mode_idx, :])
        for rank, state_idx in enumerate(np.argsort(vals)[-top_n:][::-1], start=1):
            rows.append(
                {
                    "mode_index": int(mode_idx),
                    "mode_real": float(mu[mode_idx].real),
                    "mode_imag": float(mu[mode_idx].imag),
                    "rank": rank,
                    "state_index": int(state_idx),
                    "state_name": names[state_idx] if state_idx < len(names) else str(state_idx),
                    "participation": float(vals[state_idx]),
                }
            )
    return rows


def run_variant(base_sheets: dict[str, pd.DataFrame], base_path: str, remove: list[str], out: Path):
    import andes  # type: ignore

    if remove:
        path = out / ("ieee39_" + "_".join(["no"] + [x.lower() for x in remove]) + ".xlsx")
        sheets = {k: v.copy() for k, v in base_sheets.items() if k not in remove}
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            for sheet_name, df in sheets.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        path = Path(base_path)
    ss = andes.load(str(path), setup=False, no_output=True)
    ss.setup()
    ss.PFlow.run()
    ss.EIG.run()
    return path, ss


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("outputs") / "ieee39_mode_diagnostics")
    return parser.parse_args()


def main() -> None:
    import andes  # type: ignore

    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    base = andes.get_case("ieee39/ieee39_full.xlsx")
    base_sheets = pd.read_excel(base, sheet_name=None)
    variants = {
        "base": [],
        "no_tgov": ["TGOV1N"],
        "no_pss": ["IEEEST"],
        "no_exciter_pss": ["IEEEX1", "IEEEST"],
        "genrou_only": ["TGOV1N", "IEEEX1", "IEEEST", "ACEc", "Toggler"],
    }

    summary_rows = []
    participation_rows = []
    for name, remove in variants.items():
        path, ss = run_variant(base_sheets, base, remove, args.out)
        rec = {"variant": name, "removed": ",".join(remove), "path": str(path)}
        rec.update(summarize_system(ss))
        summary_rows.append(rec)
        if name == "base":
            participation_rows.extend(participation_report(ss))

    pd.DataFrame(summary_rows).to_csv(args.out / "baseline_positive_mode_diagnostics.csv", index=False)
    pd.DataFrame(participation_rows).to_csv(args.out / "positive_mode_participation.csv", index=False)
    with (args.out / "baseline_positive_mode_diagnostics.json").open("w", encoding="utf-8") as f:
        json.dump(summary_rows, f, indent=2)
    print(pd.DataFrame(summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()
