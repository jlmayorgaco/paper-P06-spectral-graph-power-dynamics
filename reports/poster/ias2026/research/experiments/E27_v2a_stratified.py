"""E27 - Track-A v2A. Stratified load by availability interaction.

Protocol: configs/ias2026/trackA_v2.yaml. v1 remains frozen and refuted; no v1
sample is reused and the seed is deliberately different.

What v2A asks, all declared before the run:

    H1   delta_alpha_Q depends on the operating condition
    H2   PRIMARY, conditional: alpha_Q > 0 implies alpha_V < 0
    H2b  SECONDARY, global: D_V < D_Q on the COMPLEX branch position
    H4   there is a frontier in (load, availability)
    H5   d_4 tightens as the branch approaches instability

Notation. ``lambda`` is a system mode; ``mu`` is an eigenvalue of the port
interaction operator; ``d_4 = min |mu + 1|`` is evaluated at the tracked
inter-area branch of that sample, never at a generic critical mode.

Usage
    python experiments/E27_v2a_stratified.py
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _bootstrap import RESULTS
from ibr_cycles.dynamics.modal_tracking import modal_assurance
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space
from ibr_cycles.uncertainty.sampling import (
    audit_limits,
    sample_operating_point,
    stratified_grid,
)

EXPERIMENT = "E27_v2a_stratified"
CORE = (30, 33, 35, 37)
SEED = 20260910
LOAD_STRATA = np.array([0.85, 0.95, 1.00, 1.05])
AVAILABILITY_STRATA = np.array([0.70, 0.80, 0.90, 0.95, 1.00])
PER_CELL = 30
MAC_FLOOR = 0.80
MIN_VOLTAGE, MAX_VOLTAGE = 0.85, 1.15
BOOTSTRAP = 5000


def machine_labels(labels):
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def pick_anchor(spectrum, case, nominal):
    reference_vector, reference_labels = nominal
    common = sorted(set(reference_labels) & machine_labels(case.system.labels))
    if len(common) < 4:
        return None
    left = np.array([reference_vector[n] for n in common])
    index = {n: i for i, n in enumerate(case.system.labels)}
    select = [index[n] for n in common]
    candidates = [
        m
        for m in spectrum.modes
        if abs(m.value) > 1e-3 and 0.3 <= m.frequency_hz <= 1.5
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda m: modal_assurance(left, m.right[select]))


def track(anchor, ref_labels, ref_index, case, spectrum):
    common = sorted(ref_labels & machine_labels(case.system.labels))
    if len(common) < 4:
        return None
    left = anchor.right[[ref_index[n] for n in common]]
    index = {n: i for i, n in enumerate(case.system.labels)}
    select = [index[n] for n in common]
    best = max(
        (m for m in spectrum.modes if abs(m.value) > 1e-3),
        key=lambda m: modal_assurance(left, m.right[select]),
    )
    return best, modal_assurance(left, best.right[select])


def bootstrap_ci(flags, rng, draws=BOOTSTRAP):
    flags = np.asarray(flags, dtype=float)
    if flags.size == 0:
        return float("nan"), float("nan"), float("nan")
    means = flags[rng.integers(0, flags.size, size=(draws, flags.size))].mean(axis=1)
    return (
        float(flags.mean()),
        float(np.percentile(means, 2.5)),
        float(np.percentile(means, 97.5)),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Track-A v2A")
    parser.add_argument("--per-cell", type=int, default=PER_CELL)
    args = parser.parse_args()
    started = time.time()

    base_network = load_network()
    rng = np.random.default_rng(SEED)
    points = stratified_grid(LOAD_STRATA, AVAILABILITY_STRATA, args.per_cell, rng)

    nominal_q = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    nominal_mode = max(
        (m for m in eigen_analysis(nominal_q.system.A).modes if abs(m.value) > 1e-3),
        key=lambda m: m.real,
    )
    nominal = (
        {
            n: nominal_mode.right[i]
            for i, n in enumerate(nominal_q.system.labels)
            if n.startswith(("delta_sg", "omega_sg"))
        },
        machine_labels(nominal_q.system.labels),
    )

    plan_q = ReplacementPlan.of({b: 1.0 for b in CORE})
    plan_pf = ReplacementPlan.of({b: 1.0 for b in CORE}, q_policy="unity_pf")
    voltage = ConverterParameters(voltage_control=True)

    rows = []
    for index, (load, reactive, availability) in enumerate(points):
        point = sample_operating_point(
            base_network,
            active_load=load,
            reactive_load=reactive,
            availability=availability,
            availability_buses=CORE,
            rng=rng,
        )
        entry = {
            "sample": index,
            "active_load": load,
            "reactive_load": reactive,
            "availability": availability,
            "load_cell": int(np.searchsorted(LOAD_STRATA, load, side="right") - 1),
            "availability_cell": int(
                np.searchsorted(AVAILABILITY_STRATA, availability, side="right") - 1
            ),
            "slack_estimate_pu": point.slack_estimate_pu,
            "unserved_pu": point.unserved_pu,
            "dispatch_saturated": ";".join(map(str, point.saturated)),
            "curtailed": ";".join(map(str, point.curtailed)),
        }
        if not point.dispatchable:
            entry.update({"status": "REJECTED", "reason": "not dispatchable"})
            rows.append(entry)
            continue
        try:
            cases = {
                "base": solve_case(ReplacementPlan.of({}), network=point.network),
                "q": solve_case(plan_q, network=point.network),
                "v": solve_case(plan_q, network=point.network, converter=voltage),
                "pf": solve_case(plan_pf, network=point.network),
            }
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            entry.update({"status": "REJECTED", "reason": str(error)[:90]})
            rows.append(entry)
            continue

        voltages = np.abs(cases["base"].dae.power_flow.voltages)
        if voltages.min() < MIN_VOLTAGE or voltages.max() > MAX_VOLTAGE:
            entry.update({"status": "REJECTED", "reason": "voltage limits"})
            rows.append(entry)
            continue

        spectra = {k: eigen_analysis(c.system.A) for k, c in cases.items()}
        anchor = pick_anchor(spectra["base"], cases["base"], nominal)
        if anchor is None:
            entry.update({"status": "REJECTED", "reason": "no anchor in the band"})
            rows.append(entry)
            continue
        ref_labels = machine_labels(cases["base"].system.labels)
        ref_index = {n: i for i, n in enumerate(cases["base"].system.labels)}

        ok = True
        for name in ("base", "q", "v", "pf"):
            found = track(anchor, ref_labels, ref_index, cases[name], spectra[name])
            if found is None or found[1] < MAC_FLOOR:
                ok = False
                break
            entry[f"alpha_{name}"] = found[0].real
            entry[f"imag_{name}"] = found[0].imag
            entry[f"freq_{name}"] = found[0].frequency_hz
            entry[f"mac_{name}"] = found[1]
        if not ok:
            entry.update({"status": "REJECTED", "reason": "mode tracking below MAC"})
            rows.append(entry)
            continue

        # d_4 is evaluated AT the tracked inter-area branch of this sample.
        lambda_q = complex(entry["alpha_q"], entry["imag_q"])
        lambda_base = complex(entry["alpha_base"], entry["imag_base"])
        lambda_v = complex(entry["alpha_v"], entry["imag_v"])
        try:
            space = build_action_space(cases["base"], cases["q"], CORE)
            entry["d_4"] = float(
                abs(space.split(lambda_q)["closest_to_minus_one"] + 1.0)
            )
        except (np.linalg.LinAlgError, ValueError):
            entry["d_4"] = float("nan")

        audit = audit_limits(point.network, cases["base"].dae.power_flow)
        entry["slack_active_solved_pu"] = audit.slack_active_pu
        entry["slack_within_limits"] = audit.slack_within_limits
        entry["reactive_binding"] = ";".join(
            f"{bus}:{kind}" for bus, kind in audit.reactive_binding
        )

        entry["delta_alpha_q"] = entry["alpha_q"] - entry["alpha_base"]
        entry["delta_alpha_v"] = entry["alpha_v"] - entry["alpha_q"]
        # H2b uses the COMPLEX branch position: voltage control moves frequency
        # as well as damping, and a real-part-only measure would hide that.
        entry["D_Q"] = float(abs(lambda_q - lambda_base))
        entry["D_V"] = float(abs(lambda_v - lambda_base))
        entry["damping_distance_q"] = abs(entry["alpha_q"] - entry["alpha_base"])
        entry["damping_distance_v"] = abs(entry["alpha_v"] - entry["alpha_base"])
        entry["degrades"] = bool(entry["delta_alpha_q"] > 0.0)
        entry["q_unstable"] = bool(entry["alpha_q"] > 0.0)
        entry["v_stable"] = bool(entry["alpha_v"] < 0.0)
        entry["modal_restoration"] = bool(entry["D_V"] < entry["D_Q"])
        entry["damping_restoration"] = bool(
            entry["damping_distance_v"] < entry["damping_distance_q"]
        )
        entry["status"] = "ACCEPTED"
        rows.append(entry)

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_samples.csv", index=False)
    accepted = table[table.status == "ACCEPTED"]

    # ---- cell map ---------------------------------------------------------
    cells = []
    for i in range(len(LOAD_STRATA) - 1):
        for j in range(len(AVAILABILITY_STRATA) - 1):
            block = accepted[
                (accepted.load_cell == i) & (accepted.availability_cell == j)
            ]
            point, low, high = bootstrap_ci(block.degrades.to_numpy(), rng)
            cells.append(
                {
                    "load_cell": i,
                    "availability_cell": j,
                    "load_range": f"{LOAD_STRATA[i]:.2f}-{LOAD_STRATA[i + 1]:.2f}",
                    "availability_range": (
                        f"{AVAILABILITY_STRATA[j]:.2f}-{AVAILABILITY_STRATA[j + 1]:.2f}"
                    ),
                    "n": len(block),
                    "fraction_degrades": point,
                    "ci_low": low,
                    "ci_high": high,
                    "median_delta_alpha_q": float(block.delta_alpha_q.median())
                    if len(block)
                    else np.nan,
                    "q_unstable_fraction": float(block.q_unstable.mean())
                    if len(block)
                    else np.nan,
                }
            )
    cell_table = pd.DataFrame(cells)
    cell_table.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_cells.csv", index=False)

    feasible = cell_table[cell_table.n >= 5].sort_values(
        ["load_cell", "availability_cell"]
    )
    top = feasible.iloc[-1] if len(feasible) else cell_table.iloc[-1]
    bottom = feasible.iloc[0] if len(feasible) else cell_table.iloc[0]

    # H1  the effect depends on the operating condition
    h1 = bool(
        top.fraction_degrades > bottom.fraction_degrades and top.ci_low > bottom.ci_high
    )
    # H4  a frontier exists in (load, availability)
    above = cell_table[(cell_table.n >= 5) & (cell_table.ci_low > 0.5)]
    below = cell_table[(cell_table.n >= 5) & (cell_table.ci_high < 0.5)]
    h4 = bool(len(above) > 0 and len(below) > 0)
    # H2  PRIMARY, conditional mitigation
    unstable = accepted[accepted.q_unstable]
    h2_point, h2_low, h2_high = bootstrap_ci(unstable.v_stable.to_numpy(), rng)
    h2 = bool(h2_low > 0.90)
    # H2b SECONDARY, global modal restoration on the complex branch
    h2b_point, h2b_low, h2b_high = bootstrap_ci(
        accepted.modal_restoration.to_numpy(), rng
    )
    h2b = bool(h2b_low > 0.80)
    aux_point, aux_low, aux_high = bootstrap_ci(
        accepted.damping_restoration.to_numpy(), rng
    )
    # H5  the closure descriptor tightens as the branch approaches instability
    valid = accepted[np.isfinite(accepted.d_4)]
    rho, pvalue = (
        spearmanr(valid.d_4, valid.alpha_q) if len(valid) > 10 else (np.nan, np.nan)
    )
    h5 = bool(rho < 0 and pvalue < 0.01)

    passes = bool(h1 and h2 and h4 and h5)

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=SEED,
        config={
            "protocol": "configs/ias2026/trackA_v2.yaml",
            "load_strata": LOAD_STRATA.tolist(),
            "availability_strata": AVAILABILITY_STRATA.tolist(),
            "per_cell": args.per_cell,
        },
    )
    manifest.finish(
        "SUCCESS" if passes else "V2A_INCOMPLETE",
        drawn=len(table),
        accepted=len(accepted),
        rejection_rate=float(1.0 - len(accepted) / max(len(table), 1)),
        rejection_reasons=table[table.status == "REJECTED"]
        .reason.value_counts()
        .to_dict()
        if (table.status == "REJECTED").any()
        else {},
        cells=cell_table.to_dict("records"),
        H1_operating_dependence={
            "passed": h1,
            "top_cell": top.to_dict(),
            "bottom_cell": bottom.to_dict(),
        },
        H2_conditional_mitigation={
            "passed": h2,
            "fraction": h2_point,
            "ci": [h2_low, h2_high],
            "n_unstable": len(unstable),
        },
        H2b_modal_restoration={
            "passed": h2b,
            "fraction": h2b_point,
            "ci": [h2b_low, h2b_high],
            "auxiliary_damping_only": {
                "fraction": aux_point,
                "ci": [aux_low, aux_high],
            },
        },
        H4_frontier={
            "passed": h4,
            "cells_above_half": len(above),
            "cells_below_half": len(below),
        },
        H5_closure_tracks={
            "passed": h5,
            "spearman": float(rho),
            "p_value": float(pvalue),
        },
        reactive_binding_samples=int((accepted.reactive_binding != "").sum()),
        slack_within_limits=int(accepted.slack_within_limits.sum()),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("v2A stratified campaign, protocol trackA_v2.yaml, seed %d" % SEED)
    print(
        "  drawn %d, accepted %d, rejection rate %.1f%%"
        % (len(table), len(accepted), 100 * (1 - len(accepted) / max(len(table), 1)))
    )
    if (table.status == "REJECTED").any():
        print(
            "  rejection reasons: %s"
            % table[table.status == "REJECTED"].reason.value_counts().head(4).to_dict()
        )
    print()
    print("fraction of samples where the replacement DEGRADES the branch")
    grid = cell_table.pivot(
        index="load_range", columns="availability_range", values="fraction_degrades"
    )
    print(grid.to_string(float_format=lambda v: f"{v:8.3f}", na_rep="     -"))
    print()
    print("median d_alpha_Q by cell")
    print(
        cell_table.pivot(
            index="load_range",
            columns="availability_range",
            values="median_delta_alpha_q",
        ).to_string(float_format=lambda v: f"{v:8.4f}", na_rep="     -")
    )
    print()
    print("samples per cell")
    print(
        cell_table.pivot(
            index="load_range", columns="availability_range", values="n"
        ).to_string(na_rep="  -")
    )
    print()
    print("PREREGISTERED HYPOTHESES")
    print(
        "  H1  operating dependence  %-5s  top %s %.3f [%.3f, %.3f]"
        "  vs bottom %s %.3f [%.3f, %.3f]"
        % (
            h1,
            top.load_range,
            top.fraction_degrades,
            top.ci_low,
            top.ci_high,
            bottom.load_range,
            bottom.fraction_degrades,
            bottom.ci_low,
            bottom.ci_high,
        )
    )
    print(
        "  H2  conditional repair    %-5s  %.3f [%.3f, %.3f] over n=%d unstable samples"
        % (h2, h2_point, h2_low, h2_high, len(unstable))
    )
    print(
        "  H2b modal restoration     %-5s  D_V < D_Q in %.3f [%.3f, %.3f]"
        % (h2b, h2b_point, h2b_low, h2b_high)
    )
    print(
        "      auxiliary, damping only        %.3f [%.3f, %.3f]"
        % (aux_point, aux_low, aux_high)
    )
    print(
        "  H4  frontier              %-5s  %d cells above one half, %d below"
        % (h4, len(above), len(below))
    )
    print(
        "  H5  d_4 tracks alpha_Q    %-5s  Spearman %+.3f, p = %.2e" % (h5, rho, pvalue)
    )
    print()
    print(
        "  limit audit: %d of %d accepted samples have a generator on a reactive"
        " limit; slack within limits in %d"
        % (
            int((accepted.reactive_binding != "").sum()),
            len(accepted),
            int(accepted.slack_within_limits.sum()),
        )
    )
    print()
    print("V2A PASSED" if passes else "V2A DID NOT MEET ITS PREREGISTERED CRITERIA")
    print(f"manifest -> {path}")
    return 0 if passes else 1


if __name__ == "__main__":
    raise SystemExit(main())
