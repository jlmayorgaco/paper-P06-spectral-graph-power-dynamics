"""E22 - Track-A validation gate 1. Does the phenomenon survive re-equilibration?

Question
    The Track-A result was found at a frozen operating point, where a matched
    reactive policy leaves the AC solution bit-identical. Does the same inter-area
    branch still cross when the operating point is genuinely allowed to move, and
    what happens under a realistic plant voltage regulator?

Policies
    matched          converter supplies its rating share of the machine's Q.
                     The AC solution is unchanged; this is the discovery case.
    voltage_control  converter regulates its terminal voltage through a PI on
                     |V| that generates the reactive reference. The AC solution
                     is again unchanged, but the dynamics are not: this is what a
                     real plant controller does, and it is the physically
                     realistic policy.
    unity_pf         converter operates at unity power factor. The bus is retyped
                     to PQ, the power flow genuinely moves, and every device is
                     reinitialized on the new operating point.

Reported
    feasibility, voltage excursion, required reactive capability, spectral
    abscissa, the tracked inter-area branch, and the gauge-invariant closure
    distance of E17.

Status of the result
    NUMERICAL OBSERVATION at one loading level, rho = 1.

Usage
    python experiments/E22_reequilibrated_policies.py
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modal_tracking import modal_assurance
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.port_admittance import build_action_space

EXPERIMENT = "E22_reequilibrated_policies"
#: Reactive capability of the converter as a fraction of its rating, used only to
#: report whether the matched policy is physically askable of a real plant.
CAPABILITY = 0.33

POLICIES = {
    "matched": ({}, None),
    "voltage_control": ({}, ConverterParameters(voltage_control=True)),
    "unity_pf": ({"q_policy": "unity_pf"}, None),
}


def machine_labels(labels):
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def main() -> int:
    parser = argparse.ArgumentParser(description="Track-A gate 1")
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    base_voltages = np.abs(base.dae.power_flow.voltages)
    network = base.dae.network
    cores = pd.read_csv(
        RESULTS / "tables" / "E12_compatibility_census_non_composable.csv"
    )
    genuine = [
        tuple(int(b) for b in m.split("+"))
        for m in cores[
            (cores.mechanism == "A_GENUINE_HIGH_ORDER") & (cores["size"] == 4)
        ].members
    ]

    rows = []
    for core in genuine:
        reference_case = solve_case(ReplacementPlan.of({b: 1.0 for b in core}))
        reference_spectrum = eigen_analysis(reference_case.system.A)
        reference = max(
            (m for m in reference_spectrum.modes if abs(m.value) > 1e-3),
            key=lambda m: m.real,
        )
        ref_labels = machine_labels(reference_case.system.labels)
        ref_index = {n: i for i, n in enumerate(reference_case.system.labels)}
        required_q = {
            b: abs(
                (
                    base.dae.power_flow.injection(b, network.ybus)
                    + network.loads.get(b, 0j)
                ).imag
            )
            / (network.machines[b]["Sn"] / 100.0)
            for b in core
        }

        for name, (kwargs, converter) in POLICIES.items():
            entry = {
                "core": "+".join(map(str, core)),
                "policy": name,
                "max_reactive_over_rating": max(required_q.values()),
                "within_capability": max(required_q.values()) <= CAPABILITY,
            }
            try:
                case = solve_case(
                    ReplacementPlan.of({b: 1.0 for b in core}, **kwargs),
                    converter=converter,
                )
            except InfeasibleReplacement as error:
                entry.update({"status": "INFEASIBLE", "reason": str(error)})
                rows.append(entry)
                continue
            spectrum = eigen_analysis(case.system.A)
            verdict = assess(spectrum, case.system.labels)
            common = sorted(ref_labels & machine_labels(case.system.labels))
            left = reference.right[[ref_index[n] for n in common]]
            index = {n: i for i, n in enumerate(case.system.labels)}
            select = [index[n] for n in common]
            best = max(
                (m for m in spectrum.modes if abs(m.value) > 1e-3),
                key=lambda m: modal_assurance(left, m.right[select]),
            )
            try:
                space = build_action_space(base, case, core)
                closure = float(
                    abs(space.split(reference.value)["closest_to_minus_one"] + 1.0)
                )
            except (np.linalg.LinAlgError, ValueError):
                closure = float("nan")
            entry.update(
                {
                    "status": "SUCCESS",
                    "n_states": case.n_states,
                    "spectral_abscissa": verdict.spectral_abscissa,
                    "zeta_min": verdict.zeta_min,
                    "unstable": verdict.unstable,
                    "tracked_real": best.real,
                    "tracked_frequency_hz": best.frequency_hz,
                    "tracked_mac": modal_assurance(left, best.right[select]),
                    "max_voltage_change": float(
                        np.abs(
                            np.abs(case.dae.power_flow.voltages) - base_voltages
                        ).max()
                    ),
                    "min_voltage": case.dae.power_flow.min_voltage,
                    "closure_distance": closure,
                    "equilibrium_norm_g": case.equilibrium.norm_g,
                    # The MAC-best match sometimes lands on a fast real mode:
                    # the branch has merged and its real part is not comparable.
                    "branch_resolved": bool(best.frequency_hz > 0.1),
                }
            )
            rows.append(entry)

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_policies.csv", index=False)

    ok = table[table.status == "SUCCESS"]

    def stats(name: str) -> dict[str, float]:
        block = ok[ok.policy == name]
        resolved = block[block.branch_resolved]
        return {
            "cases": int(len(block)),
            "system_unstable": int(block.unstable.sum()),
            "branch_resolved": int(len(resolved)),
            "tracked_branch_unstable": int((resolved.tracked_real > 0).sum()),
            "median_abscissa": float(block.spectral_abscissa.median()),
            "min_mac": float(block.tracked_mac.min()),
            "min_voltage": float(block.min_voltage.min()),
        }

    survival = {name: stats(name) for name in POLICIES}

    # The gate is about the TRACKED inter-area branch, not about the system being
    # unstable for any reason: under unity power factor the operating point gets
    # weak enough that a converter-driven divergence can dominate, which is the
    # Track-B mechanism and must not be counted as survival of Track A.
    flagship_row = ok[(ok.core == "30+33+35+37")].set_index("policy")
    flagship_survives = bool(
        flagship_row.loc["unity_pf", "tracked_real"] > 0
        and flagship_row.loc["unity_pf", "tracked_mac"] > 0.9
    )
    survives = flagship_survives

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={
            "cores": ["+".join(map(str, c)) for c in genuine],
            "capability": CAPABILITY,
        },
    )
    manifest.finish(
        "SUCCESS" if survives else "PHENOMENON_IS_OPERATING_POINT_SPECIFIC",
        survival=survival,
        flagship_survives_reequilibration=flagship_survives,
        flagship_rows=flagship_row.reset_index()[
            [
                "policy",
                "spectral_abscissa",
                "tracked_real",
                "tracked_frequency_hz",
                "tracked_mac",
                "min_voltage",
            ]
        ].to_dict("records"),
        matched_reactive_within_capability=bool(ok.within_capability.all()),
        max_reactive_over_rating=float(ok.max_reactive_over_rating.max()),
        unity_pf_min_voltage=float(ok[ok.policy == "unity_pf"].min_voltage.min()),
        unity_pf_max_voltage_change=float(
            ok[ok.policy == "unity_pf"].max_voltage_change.max()
        ),
        voltage_control_stabilizes=int(
            (~ok[ok.policy == "voltage_control"].unstable).sum()
        ),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("nine genuine order-4 cores under three reactive policies")
    print()
    for name in POLICIES:
        block = ok[ok.policy == name]
        print("--- %s" % name)
        print(
            block[
                [
                    "core",
                    "spectral_abscissa",
                    "tracked_real",
                    "tracked_frequency_hz",
                    "tracked_mac",
                    "max_voltage_change",
                    "min_voltage",
                    "closure_distance",
                ]
            ].to_string(index=False, float_format=lambda v: f"{v:10.5f}")
        )
        s = survival[name]
        print(
            "    system unstable %d/%d | tracked branch resolved %d/%d and unstable"
            " %d | min MAC %.3f | min voltage %.4f"
            % (
                s["system_unstable"],
                s["cases"],
                s["branch_resolved"],
                s["cases"],
                s["tracked_branch_unstable"],
                s["min_mac"],
                s["min_voltage"],
            )
        )
        print()
    print(
        "matched policy reactive requirement, worst over the nine cores: %.3f of rating"
        % ok.max_reactive_over_rating.max()
    )
    print(
        "  within a %.0f%% capability: %s"
        % (CAPABILITY * 100, ok.within_capability.all())
    )
    print()
    print("flagship 30+33+35+37 across the three policies:")
    print(
        flagship_row.reset_index()[
            [
                "policy",
                "spectral_abscissa",
                "tracked_real",
                "tracked_frequency_hz",
                "tracked_mac",
                "min_voltage",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:10.5f}")
    )
    print()
    print(
        "GATE 1 PASSED, the flagship branch still crosses after re-equilibration"
        if survives
        else "GATE 1 FAILED, the flagship branch is operating-point specific"
    )
    print(f"manifest -> {path}")
    return 0 if survives else 1


if __name__ == "__main__":
    raise SystemExit(main())
