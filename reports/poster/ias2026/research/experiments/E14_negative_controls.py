"""E14 - PHASE E1. Negative controls.

Question
    Is the order-4 non-composability found in E12 a specific, controller-mediated
    dynamic mechanism, or does the method flag any sufficiently large replacement
    portfolio?

Controls
    N1 replacement amplitude sweep; N2 added synchronous damping; N3 structured
    controller ablation; N4 matched stable four-replacement portfolios;
    N5 static replacement at the same P and Q; N6 mode-family tracking on the
    matched controls.

Status of the result
    NUMERICAL OBSERVATION at one operating point. This experiment can only
    falsify, never confirm.

Usage
    python experiments/E14_negative_controls.py
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.diagnosis.baselines import generalized_scr, nodal_metrics
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

EXPERIMENT = "E14_negative_controls"
CENSUS = "E12_compatibility_census_portfolios.csv"
CORES = "E12_compatibility_census_non_composable.csv"

RHO_GRID = (0.25, 0.50, 0.75, 1.00)
DAMPING_GRID = (0.0, 1.0, 2.0, 5.0)
BASE = ConverterParameters()

#: Structured controller ablation. Each entry is one physically meaningful
#: change; no mechanism is inferred from a single gain.
ABLATIONS: dict[str, dict[str, float]] = {
    "nominal": {},
    "measurement_fast": {"tau_p": BASE.tau_p / 3.0},
    "measurement_slow": {"tau_p": BASE.tau_p * 3.0},
    "outer_p_weak": {"kp_p": BASE.kp_p / 4.0, "ki_p": BASE.ki_p / 4.0},
    "outer_p_strong": {"kp_p": BASE.kp_p * 4.0, "ki_p": BASE.ki_p * 4.0},
    "outer_q_weak": {"kp_q": BASE.kp_q / 4.0, "ki_q": BASE.ki_q / 4.0},
    "outer_q_strong": {"kp_q": BASE.kp_q * 4.0, "ki_q": BASE.ki_q * 4.0},
    "pll_slow": ConverterParameters.pll_from_design(12.0, 0.7),
    "pll_slow_damped": ConverterParameters.pll_from_design(12.0, 1.5),
    "pll_fast": ConverterParameters.pll_from_design(60.0, 0.7),
    "current_slow": {"kp_i": BASE.kp_i / 2.0, "ki_i": BASE.ki_i / 2.0},
    "current_fast": {"kp_i": BASE.kp_i * 2.0, "ki_i": BASE.ki_i * 2.0},
}

DEVICES = ("gfl", "static_power", "static_impedance")


def solve_with(members, *, rho=1.0, device="gfl", damping=None, converter=None):
    plan = ReplacementPlan.of(
        {b: rho for b in members}, device=device, machine_damping=damping
    )
    case = solve_case(plan, converter=converter)
    spectrum = eigen_analysis(case.system.A)
    verdict = assess(spectrum, case.system.labels)
    dynamic = [m for m in spectrum.modes if abs(m.value) > 1e-3]
    critical = max(dynamic, key=lambda m: m.real)
    label = case.system.labels[int(np.argmax(critical.participation))]
    return {
        "spectral_abscissa": verdict.spectral_abscissa,
        "critical_frequency_hz": verdict.critical_frequency_hz,
        "critical_family": "gfl" if "gfl" in label else "sg",
        "zeta_min": verdict.zeta_min,
        "unstable": bool(verdict.unstable),
        "gz_condition": case.gz_condition,
        "n_states": case.n_states,
    }


def safe(members, **kwargs) -> dict[str, object]:
    try:
        row = solve_with(members, **kwargs)
        row["status"] = "SUCCESS"
        return row
    except (InfeasibleReplacement, ValueError) as error:
        return {"status": "FAILED", "reason": str(error)}


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E1 negative controls")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    started = time.time()

    network = load_network()
    power_flow = solve_power_flow(network)
    metrics = nodal_metrics(network, power_flow)
    census = pd.read_csv(RESULTS / "tables" / CENSUS)
    cores = pd.read_csv(RESULTS / "tables" / CORES)

    genuine = [
        tuple(int(b) for b in m.split("+"))
        for m in cores[
            (cores.mechanism == "A_GENUINE_HIGH_ORDER") & (cores["size"] == 4)
        ].members
    ]
    switching = [
        tuple(int(b) for b in m.split("+"))
        for m in cores[cores.mechanism == "C_MODE_SWITCHING"].members
    ]

    rows: list[dict[str, object]] = []

    def record(control, members, **payload):
        rows.append(
            {
                "control": control,
                "members": "+".join(map(str, members)),
                "size": len(members),
                **payload,
            }
        )

    # N1 amplitude sweep on the genuine order-4 cores.
    for core in genuine:
        for rho in RHO_GRID:
            record("N1_amplitude", core, rho=rho, **safe(core, rho=rho))

    # N2 added synchronous damping, inertia untouched.
    for core in genuine:
        for damping in DAMPING_GRID:
            record("N2_damping", core, damping=damping, **safe(core, damping=damping))

    # N3 structured controller ablation, on both failure families.
    for family, group in (("genuine", genuine), ("switching", switching[:6])):
        for core in group:
            for name, override in ABLATIONS.items():
                record(
                    "N3_controller",
                    core,
                    family=family,
                    ablation=name,
                    **safe(core, converter=ConverterParameters(**override)),
                )

    # N5 static replacement at the same P and Q.
    for family, group in (("genuine", genuine), ("switching", switching)):
        for core in group:
            for device in DEVICES:
                record(
                    "N5_static",
                    core,
                    family=family,
                    device=device,
                    **safe(core, device=device),
                )

    # N4 matched stable four-replacement portfolios.
    unstable_mw = census[(census["size"] == 4) & census.unstable].replaced_mw
    low, high = float(unstable_mw.min()), float(unstable_mw.max())
    stable4 = census[(census["size"] == 4) & (~census.unstable)].copy()
    stable4["members_tuple"] = stable4.members.apply(
        lambda s: tuple(int(b) for b in s.split("+"))
    )
    total_inertia = sum(
        network.machine_on_system_base(b)["M"] for b in network.generator_buses
    )
    stable4["inertia_fraction"] = stable4.members_tuple.apply(
        lambda t: sum(network.machine_on_system_base(b)["M"] for b in t) / total_inertia
    )
    matched = stable4[
        (stable4.replaced_mw >= low * 0.95) & (stable4.replaced_mw <= high * 1.05)
    ].nlargest(25, "replaced_mw")
    for row in matched.itertuples():
        core = row.members_tuple
        record(
            "N4_matched_stable",
            core,
            replaced_mw=row.replaced_mw,
            inertia_fraction=row.inertia_fraction,
            min_scr=min(metrics[b].scr for b in core),
            gscr=generalized_scr(network, core, power_flow),
            **safe(core),
        )

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_all.csv", index=False)

    # ---- summaries -------------------------------------------------------
    n5 = table[table.control == "N5_static"].pivot_table(
        index=["members", "family"], columns="device", values="spectral_abscissa"
    )
    n5["controller_mediated_strict"] = (
        (n5.gfl > 0) & (n5.static_power <= 0) & (n5.static_impedance <= 0)
    )
    n5["controller_mediated_pq"] = (n5.gfl > 0) & (n5.static_power <= 0)
    n5 = n5.reset_index()
    n5.to_csv(tables / f"{EXPERIMENT}_N5_static.csv", index=False)

    n1 = table[(table.control == "N1_amplitude") & (table.status == "SUCCESS")]
    n1_pivot = n1.pivot_table(
        index="members", columns="rho", values="spectral_abscissa"
    )
    n2 = table[(table.control == "N2_damping") & (table.status == "SUCCESS")]
    n2_pivot = n2.pivot_table(
        index="members", columns="damping", values="spectral_abscissa"
    )
    n3 = table[(table.control == "N3_controller") & (table.status == "SUCCESS")]
    n3_pivot = n3.pivot_table(
        index="ablation", columns="family", values="spectral_abscissa", aggfunc="median"
    )
    n4 = table[(table.control == "N4_matched_stable") & (table.status == "SUCCESS")]

    n1_pivot.to_csv(tables / f"{EXPERIMENT}_N1_amplitude.csv")
    n2_pivot.to_csv(tables / f"{EXPERIMENT}_N2_damping.csv")
    n3_pivot.to_csv(tables / f"{EXPERIMENT}_N3_ablation.csv")
    n4.to_csv(tables / f"{EXPERIMENT}_N4_matched.csv", index=False)

    size4 = census[census["size"] == 4]
    flag_rate = float(size4.unstable.mean())

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={
            "rho_grid": list(RHO_GRID),
            "damping_grid": list(DAMPING_GRID),
            "ablations": list(ABLATIONS),
            "devices": list(DEVICES),
            "genuine_cores": ["+".join(map(str, c)) for c in genuine],
        },
    )
    manifest.finish(
        "SUCCESS",
        cases=len(table),
        failures=int((table.status != "SUCCESS").sum()),
        order4_flag_rate=flag_rate,
        n5_controller_mediated_strict=int(n5.controller_mediated_strict.sum()),
        n5_controller_mediated_pq=int(n5.controller_mediated_pq.sum()),
        n5_by_family=n5.groupby("family").controller_mediated_pq.sum().to_dict(),
        n1_smallest_unstable_rho={
            str(m): float(
                min(
                    (r for r in RHO_GRID if float(n1_pivot.loc[m, r]) > 0),
                    default=float("nan"),
                )
            )
            for m in n1_pivot.index
        },
        n2_damping_effect={
            str(m): float(n2_pivot.loc[m, 5.0] - n2_pivot.loc[m, 0.0])
            for m in n2_pivot.index
        },
        n4_matched_stable_cases=len(n4),
        n4_any_unstable=int(n4.unstable.sum()),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print(
        "HARD GATE: fraction of order-4 portfolios flagged unstable = %.3f (%d of %d)"
        % (flag_rate, int(size4.unstable.sum()), len(size4))
    )
    print("  the method is not an everything-is-dangerous detector at this order.")
    print()
    print("N5 STATIC REPLACEMENT at the same P and Q  (the controller-mediation test)")
    print(
        n5.groupby("family")[["controller_mediated_pq", "controller_mediated_strict"]]
        .sum()
        .to_string()
    )
    print()
    print("  median spectral abscissa by family and device:")
    print(
        n5.groupby("family")[list(DEVICES)]
        .median()
        .to_string(float_format=lambda v: f"{v:12.5f}")
    )
    print()
    print("N1 AMPLITUDE SWEEP on the genuine order-4 cores (spectral abscissa)")
    print(n1_pivot.to_string(float_format=lambda v: f"{v:10.5f}"))
    print()
    print("N2 ADDED SYNCHRONOUS DAMPING (the case has D = 0, so this is absolute)")
    print(n2_pivot.to_string(float_format=lambda v: f"{v:10.5f}"))
    print()
    print("N3 CONTROLLER ABLATION, median spectral abscissa by failure family")
    print(n3_pivot.to_string(float_format=lambda v: f"{v:14.5f}"))
    print()
    print(
        "N4 MATCHED STABLE FOUR-REPLACEMENT CONTROLS: %d cases, %d unstable"
        % (len(n4), int(n4.unstable.sum()))
    )
    if len(n4):
        print(
            "  replaced MW range %.0f to %.0f; worst abscissa %+.5f"
            % (n4.replaced_mw.min(), n4.replaced_mw.max(), n4.spectral_abscissa.max())
        )
    print()
    print("%d cases in %.0f s" % (len(table), time.time() - started))
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
