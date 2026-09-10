"""E10 - PHASE A. IEEE-39 base model and numerical validation.

Question
    Is the IEEE-39 replacement benchmark numerically sound and dynamically
    eligible as a starting point, before any replacement is evaluated?

Controls
    Ybus and AC power flow against an independent ANDES solution; equilibrium
    residuals; index-1 conditioning; identification of the reference modes left
    free by the absence of a governor; eligibility of the base damping.

Status of the result
    NUMERICAL OBSERVATION. This experiment produces no scientific claim; it
    decides whether the benchmark may be used at all.

Usage
    python experiments/E10_ieee39_baseline.py
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from _bootstrap import RESULTS, ROOT
from ibr_cycles.diagnosis.baselines import nodal_metrics
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_network import load_network, solve_power_flow

EXPERIMENT = "E10_ieee39_baseline"
REFERENCE = ROOT / "configs" / "ias2026" / "ieee39_andes_powerflow_reference.json"
YBUS_REFERENCE = ROOT / "configs" / "ias2026" / "ieee39_ybus_andes.npy"


def network_validation(network, power_flow) -> dict[str, float | bool]:
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    lines_only = np.load(YBUS_REFERENCE)
    shunts = np.zeros_like(lines_only)
    for bus, value in ((4, 1.0), (5, 2.0)):
        position = network.position(bus)
        shunts[position, position] += 1j * value
    magnitudes = np.abs(power_flow.voltages)
    angles = np.angle(power_flow.voltages)
    return {
        "ybus_max_deviation": float(np.abs(network.ybus - (lines_only + shunts)).max()),
        "voltage_max_deviation": float(
            np.abs(magnitudes - np.array(reference["v"])).max()
        ),
        "angle_max_deviation": float(np.abs(angles - np.array(reference["a"])).max()),
        "power_flow_mismatch": power_flow.max_mismatch,
        "min_voltage": power_flow.min_voltage,
        "max_voltage": power_flow.max_voltage,
        "reference_tool": reference["tool"],
    }


def generator_table(network, power_flow, metrics) -> pd.DataFrame:
    rows = []
    for bus in network.generator_buses:
        generation = power_flow.injection(bus, network.ybus) + network.loads.get(
            bus, 0j
        )
        rating = network.machines[bus]["Sn"] / 100.0
        machine = network.machine_on_system_base(bus)
        row = {
            "bus": bus,
            "rating_mva": network.machines[bus]["Sn"],
            "p_pu": generation.real,
            "q_pu": generation.imag,
            "apparent_over_rating": abs(generation) / rating,
            "reactive_over_rating": generation.imag / rating,
            "voltage_pu": abs(power_flow.at(bus)),
            "inertia_m_system_base": machine["M"],
            "xd1_system_base": machine["xd1"],
            "replacement_candidate": bus in network.replacement_candidates,
        }
        if bus in metrics:
            row.update(metrics[bus].as_row())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE A IEEE-39 base model")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()

    network = load_network()
    power_flow = solve_power_flow(network)
    validation = network_validation(network, power_flow)
    metrics = nodal_metrics(network, power_flow)

    case = solve_case(ReplacementPlan.of({}))
    spectrum = eigen_analysis(case.system.A)
    verdict = assess(spectrum, case.system.labels)

    table = generator_table(network, power_flow, metrics)
    modes = pd.DataFrame(
        [
            {
                "real": m.real,
                "imag": m.imag,
                "frequency_hz": m.frequency_hz,
                "damping": m.damping,
                "condition": m.condition,
                "dominant_state": case.system.labels[int(np.argmax(m.participation))],
                "participation": float(m.participation.max()),
            }
            for m in spectrum.modes
        ]
    ).sort_values("real", ascending=False)

    tables = RESULTS / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table.to_csv(tables / f"{EXPERIMENT}_generators.csv", index=False)
    modes.to_csv(tables / f"{EXPERIMENT}_base_modes.csv", index=False)

    eligible = (
        validation["ybus_max_deviation"] < 1e-9
        and validation["voltage_max_deviation"] < 1e-5
        and case.equilibrium.norm_f < 1e-9
        and case.equilibrium.norm_g < 1e-9
        and case.gz_condition < 1e8
        and verdict.reference_exclusion_is_clean
        and not verdict.unstable
    )

    config = {
        "network_sha256": network.source_sha256,
        "fidelity": "two_axis_4th_order + first_order_avr + ieeest_2_state_pss",
        "governor": "not_retained",
        "candidates": list(network.replacement_candidates),
    }
    manifest = Manifest(experiment=EXPERIMENT, seed=args.seed, config=config)
    manifest.finish(
        "SUCCESS" if eligible else "BASE_INELIGIBLE",
        validation=validation,
        n_differential_states=case.n_states,
        n_algebraic_variables=case.dae.n_z,
        gz_condition=case.gz_condition,
        equilibrium={
            "norm_f": case.equilibrium.norm_f,
            "norm_g": case.equilibrium.norm_g,
        },
        base_assessment=verdict.__dict__,
        max_device_loading=case.max_loading,
        nodal_metrics=[m.as_row() for m in metrics.values()],
    )
    path = manifest.write(RESULTS / "manifests")

    print(
        "NETWORK VALIDATION (independent reference: %s)" % validation["reference_tool"]
    )
    print("  Ybus max deviation      %.3e" % validation["ybus_max_deviation"])
    print("  |V| max deviation       %.3e" % validation["voltage_max_deviation"])
    print("  angle max deviation     %.3e rad" % validation["angle_max_deviation"])
    print("  power-flow mismatch     %.3e" % validation["power_flow_mismatch"])
    print(
        "  voltage range           %.4f to %.4f pu"
        % (validation["min_voltage"], validation["max_voltage"])
    )
    print()
    print("MODEL")
    print("  differential states     %d" % case.n_states)
    print("  algebraic variables     %d" % case.dae.n_z)
    print(
        "  equilibrium |f|,|g|     %.2e , %.2e"
        % (case.equilibrium.norm_f, case.equilibrium.norm_g)
    )
    print("  cond(gz)                %.3e" % case.gz_condition)
    print()
    print("REFERENCE MODES (no governor retained)")
    print("  excluded                %d  (expected 2)" % verdict.reference_modes)
    print("  largest excluded |lam|  %.3e" % verdict.reference_max_abs)
    print("  mean delta/omega share  %.3f" % verdict.reference_participation)
    print()
    print("BASE ELIGIBILITY")
    print(
        "  spectral abscissa       %+.6f  (non-reference)" % verdict.spectral_abscissa
    )
    print("  RHP modes               %d" % verdict.rhp_modes)
    print(
        "  zeta_min (0.1-5 Hz)     %+.6f at %.4f Hz"
        % (verdict.zeta_min, verdict.zeta_min_frequency_hz)
    )
    print("  oscillatory modes       %d" % verdict.oscillatory_modes)
    print()
    print(
        table[
            [
                "bus",
                "rating_mva",
                "p_pu",
                "q_pu",
                "apparent_over_rating",
                "scr",
                "thevenin_magnitude_pu",
                "replacement_candidate",
            ]
        ].to_string(index=False, na_rep="-")
    )
    print()
    print(f"manifest -> {path}")
    print("PHASE A ELIGIBLE" if eligible else "PHASE A BASE INELIGIBLE")
    return 0 if eligible else 1


if __name__ == "__main__":
    raise SystemExit(main())
