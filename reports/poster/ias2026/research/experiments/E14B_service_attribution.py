"""E14B - PHASE E0. Which synchronous service does the Track-A failure need?

Question
    E14 showed the order-4 inter-area failure is not caused by converter control.
    It did not say what was lost. Which synchronous service, restored at the four
    retired buses, removes the failure?

Method
    A synchronous condenser is retained at each retired bus: the full machine
    apparatus at a declared rating, producing no active power, while the
    converter carries every displaced megawatt. Individual services are then
    scaled off one at a time, and in a small factorial.

Tracked quantity
    The flagship inter-area branch is followed by MAC on the machine states
    common to both configurations. Reporting the spectral abscissa alone would
    be wrong: restoring the condenser moves the binding constraint to a
    different mode, and the branch of interest must still be located.

Status of the result
    NUMERICAL OBSERVATION at one operating point, one core, rho = 1.

Usage
    python experiments/E14B_service_attribution.py
"""

from __future__ import annotations

import argparse
import time

import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modal_tracking import modal_assurance
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

EXPERIMENT = "E14B_service_attribution"
CORE = (30, 33, 35, 37)
#: Inertia is scaled rather than removed: a condenser with zero inertia is a
#: singular machine, not a physical service ablation.
INERTIA_GRID = (0.10, 0.25, 0.50, 0.75, 1.00)
PSS_GRID = (0.00, 0.25, 0.50, 0.75, 1.00)
CONDENSER_GRID = (0.25, 0.50, 1.00)


def machine_labels(labels: tuple[str, ...]) -> set[str]:
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def solve(**kwargs) -> tuple[object, object]:
    case = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}, **kwargs))
    return case, eigen_analysis(case.system.A)


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E0 service attribution")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    started = time.time()

    # Reference branch: the critical mode of the unrepaired flagship.
    flagship, flagship_spectrum = solve()
    dynamic = [m for m in flagship_spectrum.modes if abs(m.value) > 1e-3]
    reference = max(dynamic, key=lambda m: m.real)
    ref_labels = machine_labels(flagship.system.labels)
    ref_index = {n: i for i, n in enumerate(flagship.system.labels)}

    def track(case, spectrum) -> dict[str, float]:
        common = sorted(ref_labels & machine_labels(case.system.labels))
        left = reference.right[[ref_index[n] for n in common]]
        index = {n: i for i, n in enumerate(case.system.labels)}
        select = [index[n] for n in common]
        candidates = [m for m in spectrum.modes if abs(m.value) > 1e-3]
        best = max(candidates, key=lambda m: modal_assurance(left, m.right[select]))
        return {
            "tracked_real": best.real,
            "tracked_frequency_hz": best.frequency_hz,
            "tracked_damping": best.damping,
            "tracked_mac": modal_assurance(left, best.right[select]),
            "common_machines": len(common) // 2,
        }

    rows: list[dict[str, object]] = []

    def record(group: str, label: str, case, spectrum, **extra):
        verdict = assess(spectrum, case.system.labels)
        rows.append(
            {
                "group": group,
                "configuration": label,
                "n_states": case.n_states,
                "spectral_abscissa": verdict.spectral_abscissa,
                "critical_frequency_hz": verdict.critical_frequency_hz,
                "zeta_min": verdict.zeta_min,
                "unstable": verdict.unstable,
                "gz_condition": case.gz_condition,
                **track(case, spectrum),
                **extra,
            }
        )

    base_case = solve_case(ReplacementPlan.of({}))
    base_spectrum = eigen_analysis(base_case.system.A)
    record("A0", "original synchronous generators", base_case, base_spectrum)
    record("A1", "full PV-GFL replacement", flagship, flagship_spectrum)

    full = {b: 1.0 for b in CORE}
    for rating in CONDENSER_GRID:
        case, spectrum = solve(condenser={b: rating for b in CORE})
        record("A2", f"condenser rating {rating:g}", case, spectrum, rating=rating)

    for scale in INERTIA_GRID:
        case, spectrum = solve(condenser=full, condenser_services={"inertia": scale})
        record("A3", f"inertia scale {scale:g}", case, spectrum, inertia=scale)

    for scale in PSS_GRID:
        case, spectrum = solve(condenser=full, condenser_services={"pss": scale})
        record("A4", f"PSS scale {scale:g}", case, spectrum, pss=scale)

    case, spectrum = solve(condenser=full, condenser_services={"avr_manual": 1.0})
    record("A5", "manual excitation (AVR frozen)", case, spectrum)

    case, spectrum = solve(
        condenser=full, condenser_services={"avr_manual": 1.0, "pss": 0.0}
    )
    record("A6", "manual excitation and no PSS", case, spectrum)

    for inertia in (0.10, 1.00):
        for avr_manual in (0.0, 1.0):
            for pss in (0.0, 1.0):
                case, spectrum = solve(
                    condenser=full,
                    condenser_services={
                        "inertia": inertia,
                        "avr_manual": avr_manual,
                        "pss": pss,
                    },
                )
                record(
                    "A7",
                    f"H x{inertia:g}, AVR {'manual' if avr_manual else 'auto'}, "
                    f"PSS {'off' if pss == 0 else 'on'}",
                    case,
                    spectrum,
                    inertia=inertia,
                    avr_manual=avr_manual,
                    pss=pss,
                )

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table.to_csv(tables / f"{EXPERIMENT}_services.csv", index=False)

    a1 = table[table.group == "A1"].iloc[0]
    a7 = table[table.group == "A7"]
    single_service = {
        "inertia_only_matters": bool(
            (table[table.group == "A3"].tracked_real.max() > 0)
            and (table[table.group == "A3"].tracked_real.min() < 0)
        ),
        "pss_only_matters": bool(
            (table[table.group == "A4"].tracked_real.max() > 0)
            and (table[table.group == "A4"].tracked_real.min() < 0)
        ),
    }

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={
            "core": list(CORE),
            "inertia_grid": list(INERTIA_GRID),
            "pss_grid": list(PSS_GRID),
            "condenser_grid": list(CONDENSER_GRID),
        },
    )
    manifest.finish(
        "SUCCESS",
        cases=len(table),
        reference_branch={
            "real": reference.real,
            "frequency_hz": reference.frequency_hz,
        },
        unrepaired_tracked_real=float(a1.tracked_real),
        condenser_restores_stability=bool(
            table[table.group == "A2"].tracked_real.max() < 0
        ),
        minimum_condenser_rating_that_works=float(
            table[(table.group == "A2") & (table.tracked_real < 0)].rating.min()
        )
        if (table[table.group == "A2"].tracked_real < 0).any()
        else None,
        factorial=a7[
            ["configuration", "tracked_real", "tracked_frequency_hz", "unstable"]
        ].to_dict("records"),
        single_service=single_service,
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    columns = [
        "configuration",
        "n_states",
        "spectral_abscissa",
        "tracked_real",
        "tracked_frequency_hz",
        "tracked_damping",
        "tracked_mac",
    ]
    print(
        "reference branch (unrepaired flagship): %+.5f at %.4f Hz"
        % (reference.real, reference.frequency_hz)
    )
    print()
    for group in ("A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7"):
        block = table[table.group == group]
        if block.empty:
            continue
        print(f"--- {group}")
        print(
            block[columns].to_string(index=False, float_format=lambda v: f"{v:10.5f}")
        )
        print()
    print("%d cases in %.0f s" % (len(table), time.time() - started))
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
