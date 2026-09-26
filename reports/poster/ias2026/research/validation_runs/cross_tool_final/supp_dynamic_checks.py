# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""Supplementary Phase D checks the pack does not perform.

The pack's run_andes_dynamic.py re-runs ANDES but compares it with the FROZEN internal
columns of results/F1_eigenvalue_reconciliation.csv. This script:
- re-solves the internal model live for every F1 stage and subset (same code path as
  experiments/F1_reconciliation.py; nothing is tuned);
- compares with the FRESH ANDES spectra written by the pack run
  (pack/results/F1_andes_spectra.npz):
    * whole dynamic spectrum, nearest-neighbour worst distance (as in F1);
    * electromechanical band (0.3-1.5 Hz): every internal band mode to its nearest
      ANDES eigenvalue AND every ANDES band mode to its nearest internal eigenvalue;
      the worst of both directions is the band mismatch;
- checks repeatability: fresh ANDES spectra vs the frozen results/F1 spectra, and live
  internal band alpha/frequency vs the frozen internal columns.

Writes supp/dynamic_supp.csv and supp/dynamic_supp.json.
Usage (tx3-analysis venv, from the research root): python supp_dynamic_checks.py .
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "src"))
sys.path.insert(0, str(repo / "experiments"))
import F1_reconciliation as f1  # noqa: E402  (module import only; main() is not called)
from ibr_cycles.dynamics.modal_family import band_candidates  # noqa: E402
from ibr_cycles.dynamics.modes import eigen_analysis  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402

HERE = Path(__file__).resolve().parent
out = HERE / "supp"
out.mkdir(exist_ok=True)
fresh = np.load(HERE / "pack/results/F1_andes_spectra.npz")
frozen = np.load(repo / "results/F1/F1_andes_spectra.npz")
frozen_table = pd.read_csv(repo / "results/F1_eigenvalue_reconciliation.csv")
fresh_table = pd.read_csv(HERE / "pack/results/andes_equivalent.csv")
LO, HI = f1.BAND_HZ


def band(values: np.ndarray) -> np.ndarray:
    f = values.imag / (2 * np.pi)
    return values[(f >= LO) & (f <= HI) & (np.abs(values) > f1.ZERO)]


def nearest(src: np.ndarray, dst: np.ndarray) -> float:
    if src.size == 0:
        return 0.0
    if dst.size == 0:
        return float("inf")
    return float(np.abs(src[:, None] - dst[None, :]).min(axis=1).max())


rows = []
for stage, services in f1.STAGES.items():
    for members in f1.SUBSETS:
        label = "+".join(map(str, members)) or "BASE"
        plan = (
            ReplacementPlan.of({b: 1.0 for b in members}, device="static_power")
            if members
            else ReplacementPlan.of({})
        )
        case = solve_case(plan, machine_services=services)
        internal = f1.dynamic_spectrum(case.system.A)
        spectrum = eigen_analysis(case.system.A)
        bc = band_candidates(spectrum, f1.BAND_HZ)
        worst = max(bc, key=lambda m: m.real) if bc else None
        key = f"{stage}|{label}"
        an = fresh[key]
        an = np.sort_complex(an[np.abs(an) > f1.ZERO])
        whole_abs, whole_rel = f1.match(internal, an)
        b_int, b_an = band(internal), band(an)
        fz = (
            np.sort_complex(frozen[key])
            if key in frozen.files
            else np.zeros(0, complex)
        )
        fr = np.sort_complex(fresh[key])
        frozen_row = frozen_table[
            (frozen_table.stage == stage) & (frozen_table.members == label)
        ].iloc[0]
        fresh_row = fresh_table[
            (fresh_table.stage == stage) & (fresh_table.members == label)
        ].iloc[0]
        rows.append(
            {
                "stage": stage,
                "members": label,
                "internal_alpha_band_live": worst.real,
                "internal_freq_band_hz_live": worst.frequency_hz,
                "internal_band_rhp_live": sum(1 for m in bc if m.real > 0.0),
                "andes_alpha_band_fresh": float(fresh_row.alpha_band),
                "andes_freq_band_hz_fresh": float(fresh_row.freq_band_hz),
                "andes_band_rhp_fresh": int(fresh_row.band_rhp_count),
                "alpha_error_live_vs_fresh": abs(
                    worst.real - float(fresh_row.alpha_band)
                ),
                "freq_error_hz_live_vs_fresh": abs(
                    worst.frequency_hz - float(fresh_row.freq_band_hz)
                ),
                "rhp_agrees_live_vs_fresh": sum(1 for m in bc if m.real > 0.0)
                == int(fresh_row.band_rhp_count),
                "n_band_internal": int(b_int.size),
                "n_band_andes": int(b_an.size),
                "band_mismatch_internal_to_andes": nearest(b_int, an),
                "band_mismatch_andes_to_internal": nearest(b_an, internal),
                "whole_spectrum_worst_match": whole_abs,
                "whole_spectrum_worst_relative": whole_rel,
                "fresh_vs_frozen_andes_spectrum_maxabs": float(np.max(np.abs(fr - fz)))
                if fz.size == fr.size
                else float("nan"),
                "internal_live_vs_frozen_alpha": abs(
                    worst.real - float(frozen_row.internal_alpha_band)
                ),
                "internal_live_vs_frozen_freq_hz": abs(
                    worst.frequency_hz - float(frozen_row.internal_freq_band_hz)
                ),
            }
        )
        print(
            f"{stage:22s} {label:12s} a_int={worst.real:+.6f} a_andes={float(fresh_row.alpha_band):+.6f} "
            f"band_mis={max(rows[-1]['band_mismatch_internal_to_andes'], rows[-1]['band_mismatch_andes_to_internal']):.2e} "
            f"whole={whole_abs:.2e}",
            flush=True,
        )

df = pd.DataFrame(rows)
df["band_mismatch"] = df[
    ["band_mismatch_internal_to_andes", "band_mismatch_andes_to_internal"]
].max(axis=1)
df.to_csv(out / "dynamic_supp.csv", index=False)
summary = {}
for stage in f1.STAGES:
    s = df[df.stage == stage]
    summary[stage] = {
        "n_subsets": int(len(s)),
        "max_alpha_error": float(s.alpha_error_live_vs_fresh.max()),
        "max_freq_error_hz": float(s.freq_error_hz_live_vs_fresh.max()),
        "rhp_agree_all": bool(s.rhp_agrees_live_vs_fresh.all()),
        "max_band_mismatch": float(s.band_mismatch.max()),
        "max_whole_spectrum_match": float(s.whole_spectrum_worst_match.max()),
        "max_fresh_vs_frozen_andes": float(
            s.fresh_vs_frozen_andes_spectrum_maxabs.max()
        ),
        "max_internal_live_vs_frozen_alpha": float(
            s.internal_live_vs_frozen_alpha.max()
        ),
    }
(out / "dynamic_supp.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
