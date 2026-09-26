"""E31 - overnight 2. Independent ANDES validation of the flagship lattice.

The ANDES side is produced by ``E31_andes_worker.py`` under its own interpreter
and imports none of this project's code. This driver runs the matching lattice
here and compares, against bands declared in ``BANDS`` below BEFORE any lattice
comparison was inspected.

The comparison is deliberately like for like. ANDES has no model of this
project's grid-following converter, so the shared configuration is the one both
tools can express: the machine at a replaced bus is removed and its power keeps
being delivered as a constant-power injection. That is exactly the negative
control this project calls ``static_power``, and it is the configuration in which
CLAIMS.md N4c asserts the failure persists without any converter.
"""

from __future__ import annotations

import subprocess
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, pin_blas_threads
from _v2c_common import (
    BAND_HZ,
    CORE,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)

NAME = "E31_ANDES_validation"
ANDES_PYTHON = Path(
    "C:/Users/walla/Documents/Github/paper-P06-spectral-graph-power-dynamics"
    "/.venv/tx3-andes/Scripts/python.exe"
)

#: Declared before any lattice comparison was looked at.
BANDS = {
    "base_frequency_relative": 0.10,
    "base_damping_sign": "must agree",
    "ordering_full_portfolio_is_worst": "must agree",
    "ordering_spearman_alpha_vs_size": 0.5,
    "per_subset_shift_sign_agreement": 0.80,
    "note": (
        "device models differ by construction, so magnitudes are not compared; "
        "only the frequency, the sign of the damping and the ordering are"
    ),
}


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question=(
            "Does an independent dynamic implementation reproduce the qualitative "
            "flagship ordering?"
        ),
        config={
            "andes_version": "2.0.0",
            "andes_case": "andes/cases/ieee39/ieee39_full.xlsx",
            "shared_configuration": "machine removed, constant-power injection kept",
            "bands_declared_before_comparison": BANDS,
            "band_hz": list(BAND_HZ),
        },
        workers=1,
    )
    started = time.time()
    raw = experiment.path("E31_andes_raw.csv")
    if not raw.exists():
        experiment.note("running the ANDES worker under its own interpreter")
        subprocess.run(
            [
                str(ANDES_PYTHON),
                str(Path(__file__).with_name("E31_andes_worker.py")),
                str(raw),
                str(experiment.path("cases")),
            ],
            check=True,
            capture_output=True,
        )
    andes_table = pd.read_csv(raw)
    andes_ok = andes_table[andes_table.status == "OK"]
    experiment.note(
        f"ANDES: {len(andes_ok)} of {len(andes_table)} cases solved"
    )

    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    pinned = sorted(machine_labels(flagship.system.labels))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)

    rows = []
    for size in range(len(CORE) + 1):
        for members in combinations(CORE, size):
            label = "+".join(map(str, members)) or "BASE"
            entry = {"members": label, "size": size}
            for device, tag in (("static_power", "static"), ("gfl", "gfl")):
                try:
                    case = (
                        base
                        if not members
                        else solve_case(
                            ReplacementPlan.of(
                                {b: 1.0 for b in members}, device=device
                            )
                        )
                    )
                    spectrum = eigen_analysis(case.system.A)
                    reading = read_family(
                        anchor, base, case, spectrum, common_labels=pinned
                    )
                    entry[f"alpha_{tag}"] = (
                        reading.family.alpha if reading else float("nan")
                    )
                    entry[f"freq_{tag}_hz"] = (
                        reading.family.frequency_worst_hz if reading else float("nan")
                    )
                    entry[f"band_rhp_{tag}"] = int(
                        sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0)
                    )
                except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
                    entry[f"alpha_{tag}"] = float("nan")
                    entry[f"freq_{tag}_hz"] = float("nan")
                    entry[f"band_rhp_{tag}"] = -1
            rows.append(entry)
    internal = pd.DataFrame(rows)

    merged = internal.merge(
        andes_ok[andes_ok.governors == False][  # noqa: E712 - explicit column filter
            ["members", "alpha_band", "freq_band_hz", "band_rhp_count", "rhp_total"]
        ].rename(
            columns={
                "alpha_band": "alpha_andes",
                "freq_band_hz": "freq_andes_hz",
                "band_rhp_count": "band_rhp_andes",
                "rhp_total": "rhp_total_andes",
            }
        ),
        on="members",
        how="left",
    )
    base_row = merged[merged.members == "BASE"].iloc[0]
    merged["shift_static"] = merged.alpha_static - base_row.alpha_static
    merged["shift_andes"] = merged.alpha_andes - base_row.alpha_andes
    merged["shift_sign_agrees"] = np.sign(merged.shift_static) == np.sign(
        merged.shift_andes
    )
    merged["abs_freq_error"] = (merged.freq_andes_hz - merged.freq_static_hz).abs()
    merged["rel_freq_error"] = merged.abs_freq_error / merged.freq_static_hz
    base_row = merged[merged.members == "BASE"].iloc[0]  # now with derived columns
    experiment.save_table(merged, "E31_ANDES_validation.csv")

    non_empty = merged[merged["size"] > 0]
    checks = {
        "base_frequency_relative": float(base_row.rel_freq_error)
        if np.isfinite(base_row.rel_freq_error)
        else float("nan"),
        "base_damping_sign_agrees": bool(
            np.sign(base_row.alpha_static) == np.sign(base_row.alpha_andes)
        ),
        "internal_full_portfolio_is_worst": bool(
            merged.alpha_static.idxmax() == merged.index[-1]
        ),
        "andes_full_portfolio_is_worst": bool(
            merged.alpha_andes.idxmax() == merged.index[-1]
        ),
        "spearman_internal_alpha_vs_size": float(
            pd.Series(merged.alpha_static).corr(pd.Series(merged["size"]), method="spearman")
        ),
        "spearman_andes_alpha_vs_size": float(
            pd.Series(merged.alpha_andes).corr(pd.Series(merged["size"]), method="spearman")
        ),
        "shift_sign_agreement": float(non_empty.shift_sign_agrees.mean()),
        "internal_static_flagship_unstable": bool(
            merged.iloc[-1].alpha_static > 0
        ),
        "andes_flagship_unstable": bool(merged.iloc[-1].alpha_andes > 0),
        "andes_spurious_rhp_total_base": int(base_row.rhp_total_andes),
        "median_rel_freq_error": float(merged.rel_freq_error.median()),
    }

    frequency_ok = checks["base_frequency_relative"] <= BANDS["base_frequency_relative"]
    ordering_ok = (
        checks["internal_full_portfolio_is_worst"]
        and checks["andes_full_portfolio_is_worst"]
        and checks["shift_sign_agreement"] >= BANDS["per_subset_shift_sign_agreement"]
    )
    verdict = "PASS" if (frequency_ok and ordering_ok) else (
        "G2_MIXED" if frequency_ok else "G2_FAIL"
    )

    print()
    print(
        merged[
            [
                "members",
                "size",
                "alpha_static",
                "alpha_andes",
                "shift_static",
                "shift_andes",
                "shift_sign_agrees",
                "freq_static_hz",
                "freq_andes_hz",
                "rel_freq_error",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:9.5f}")
    )
    print()
    for key, value in checks.items():
        print(f"  {key:42s} {value}")

    experiment.finish(verdict, checks=checks, bands=BANDS, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
