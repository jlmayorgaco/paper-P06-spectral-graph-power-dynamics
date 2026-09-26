"""F1 - equation-equivalent reconciliation of the internal model against ANDES.

The ANDES side is produced by ``F1_andes_equivalent_worker.py`` under its own
interpreter and imports none of this project's code.

Two configurations both tools can express exactly:

    R2  manual excitation, no stabilizer, no governor
    R3  first-order AVR, no stabilizer, no governor

ANDES realises the first-order AVR as SEXS with TA/TB = 1, which collapses its
lead-lag to unity and leaves ``K / (1 + s TE)`` — the internal AVR term for term
— with limits opened so that neither model saturates.

The comparison is at the level of the WHOLE dynamic spectrum, not one tracked
mode. The internal reduced matrix and the ANDES state matrix have different
dimensions and different state orderings, so the invariant compared is the
multiset of eigenvalues: every internal dynamic eigenvalue is matched to its
nearest ANDES eigenvalue and the worst matching distance is reported. A model
difference of any consequence shows up there.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
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
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

NAME = "F1_reconciliation"
ZERO = 1e-3
SUBSETS = [()] + [(b,) for b in CORE] + [tuple(CORE)]
#: internal service settings matching each ANDES stage
STAGES = {
    "R2 manual excitation": {"pss": 0.0, "avr_manual": 1.0},
    "R3 first-order AVR": {"pss": 0.0},
}


def dynamic_spectrum(matrix):
    values = np.linalg.eigvals(np.asarray(matrix, dtype=np.complex128))
    return np.sort_complex(values[np.abs(values) > ZERO])


def match(internal, andes):
    """Nearest-neighbour matching distance from each internal eigenvalue."""

    if internal.size == 0 or andes.size == 0:
        return float("nan"), float("nan")
    distances = np.abs(internal[:, None] - andes[None, :]).min(axis=1)
    scale = max(float(np.abs(internal).max()), 1.0)
    return float(distances.max()), float(distances.max() / scale)


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Are ANDES and the internal simulator solving the same dynamic model?",
        config={
            "stages": list(STAGES),
            "andes_avr": "SEXS with TA/TB = 1, K = case KA, TE = case TE, limits opened",
            "andes_machine": "GENROU with xd2 = xd1, xq2 = xq1, subtransient constants 1e-4",
            "loads": "constant power on both sides, pq2z disabled",
            "governor": "removed on both sides",
            "comparison": "whole dynamic spectrum, nearest-neighbour matching",
        },
        workers=1,
    )
    started = time.time()

    andes_table = pd.read_csv(RESULTS / "F1" / "F1_andes_equivalent.csv")
    spectra = np.load(RESULTS / "F1" / "F1_andes_spectra.npz")

    rows = []
    for stage, services in STAGES.items():
        reference = solve_case(ReplacementPlan.of({}), machine_services=services)
        flagship = solve_case(
            ReplacementPlan.of({b: 1.0 for b in CORE}, device="static_power"),
            machine_services=services,
        )
        pinned = sorted(machine_labels(flagship.system.labels))
        anchor = base_anchor(
            eigen_analysis(reference.system.A), reference, nominal_reference()
        )
        for members in SUBSETS:
            label = "+".join(map(str, members)) or "BASE"
            case = (
                reference
                if not members
                else solve_case(
                    ReplacementPlan.of({b: 1.0 for b in members}, device="static_power"),
                    machine_services=services,
                )
            )
            spectrum = eigen_analysis(case.system.A)
            reading = read_family(
                anchor, reference, case, spectrum, common_labels=pinned
            )
            band = band_candidates(spectrum, BAND_HZ)
            worst_band = max(band, key=lambda m: m.real) if band else None
            internal = dynamic_spectrum(case.system.A)

            key = f"{stage}|{label}"
            reference_row = andes_table[
                (andes_table.stage == stage) & (andes_table.members == label)
            ]
            andes_values = spectra[key] if key in spectra.files else np.zeros(0, complex)
            andes_dynamic = np.sort_complex(
                andes_values[np.abs(andes_values) > ZERO]
            )
            absolute, relative = match(internal, andes_dynamic)

            row = {
                "stage": stage,
                "members": label,
                "size": len(members),
                "internal_alpha_band": worst_band.real if worst_band else np.nan,
                "internal_freq_band_hz": worst_band.frequency_hz if worst_band else np.nan,
                "internal_band_rhp": sum(1 for m in band if m.real > 0.0),
                "internal_alpha_family": reading.family.alpha if reading else np.nan,
                "internal_n_dynamic": int(internal.size),
                "andes_n_dynamic": int(andes_dynamic.size),
                "worst_eigenvalue_match": absolute,
                "worst_relative_match": relative,
            }
            if not reference_row.empty:
                entry = reference_row.iloc[0]
                row.update(
                    {
                        "andes_alpha_band": float(entry.alpha_band),
                        "andes_freq_band_hz": float(entry.freq_band_hz),
                        "andes_band_rhp": int(entry.band_rhp_count),
                        "alpha_band_error": abs(
                            float(entry.alpha_band) - row["internal_alpha_band"]
                        ),
                        "freq_band_error_hz": abs(
                            float(entry.freq_band_hz) - row["internal_freq_band_hz"]
                        ),
                        "rhp_agrees": bool(
                            int(entry.band_rhp_count) == row["internal_band_rhp"]
                        ),
                    }
                )
            rows.append(row)
            print(
                "  %-22s %-12s internal a=%+.5f f=%.4f rhp=%d | andes a=%+.5f f=%.4f rhp=%d"
                " | worst spectrum match %.2e"
                % (
                    stage, label,
                    row["internal_alpha_band"], row["internal_freq_band_hz"],
                    row["internal_band_rhp"],
                    row.get("andes_alpha_band", np.nan),
                    row.get("andes_freq_band_hz", np.nan),
                    row.get("andes_band_rhp", -1),
                    absolute,
                ),
                flush=True,
            )

    table = pd.DataFrame(rows)
    target = RESULTS / "F1_eigenvalue_reconciliation.csv"
    table.to_csv(target, index=False)
    experiment.save_table(table, "F1_eigenvalue_reconciliation.csv")

    r3 = table[table.stage == "R3 first-order AVR"]
    r2 = table[table.stage == "R2 manual excitation"]
    checks = {
        "R3_max_alpha_band_error": float(r3.alpha_band_error.max()),
        "R3_max_freq_band_error_hz": float(r3.freq_band_error_hz.max()),
        "R3_rhp_agrees_everywhere": bool(r3.rhp_agrees.all()),
        "R3_worst_spectrum_match": float(r3.worst_eigenvalue_match.max()),
        "R2_max_alpha_band_error": float(r2.alpha_band_error.max()),
        "R2_rhp_agrees_everywhere": bool(r2.rhp_agrees.all()),
        "R2_worst_spectrum_match": float(r2.worst_eigenvalue_match.max()),
        "internal_flagship_R3": float(
            r3[r3.members == "30+33+35+37"].internal_alpha_band.iloc[0]
        ),
        "andes_flagship_R3": float(
            r3[r3.members == "30+33+35+37"].andes_alpha_band.iloc[0]
        ),
        "internal_flagship_R2": float(
            r2[r2.members == "30+33+35+37"].internal_alpha_band.iloc[0]
        ),
        "andes_flagship_R2": float(
            r2[r2.members == "30+33+35+37"].andes_alpha_band.iloc[0]
        ),
    }
    agree = (
        checks["R3_rhp_agrees_everywhere"]
        and checks["R3_max_alpha_band_error"] < 1e-4
        and checks["R3_max_freq_band_error_hz"] < 1e-3
    )
    verdict = "HARMONIZED_AGREE" if agree else "HARMONIZED_DISAGREE"

    print()
    for key, value in checks.items():
        print(f"  {key:32s} {value}")
    print()
    print("F1 GATE:", verdict)
    experiment.finish(verdict, checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
