"""E26 - Track-A validation gate 5. Held-out Monte Carlo of the DIRECTION.

Frozen protocol: configs/ias2026/trackA_final_v1.yaml. Nothing in this script
selects, tunes or filters on anything seen after that file was written.

Hypotheses, declared in the frozen file before any sample was drawn, all on the
SAME tracked branch of each individual sample:

    H1  reactive-setpoint replacement degrades it        d_alpha_Q > 0
    H2  voltage regulation mitigates that degradation    d_alpha_V < 0
    H3  ordering  alpha_unityPF > alpha_Qset > alpha_Vctrl

The question is not how many operating points are unstable. It is whether the
direction of the effect and the direction of the repair generalize away from the
operating point where the mechanism was found.

Usage
    python experiments/E26_holdout_monte_carlo.py [--samples 500]
"""

from __future__ import annotations

import argparse
import time
from dataclasses import replace as dataclass_replace

import numpy as np
import pandas as pd
from scipy.stats import qmc

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
from ibr_cycles.models.ieee39_network import load_network

EXPERIMENT = "E26_holdout_monte_carlo"
CORE = (30, 33, 35, 37)
SEED = 20260909
MAC_FLOOR = 0.80
MIN_VOLTAGE = 0.85
MAX_VOLTAGE = 1.15
BOOTSTRAP = 5000

BOUNDS = {
    "active_load": (0.85, 1.15),
    "reactive_load": (0.90, 1.10),
    "pv_availability": (0.70, 1.00),
}
SCATTER = 0.05
DISPATCH = 0.10


def machine_labels(labels):
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def perturbed_network(base_network, sample, rng):
    """One held-out operating point: loads, dispatch and PV availability."""

    active, reactive, availability = sample
    loads = {}
    for bus, load in base_network.loads.items():
        jitter = 1.0 + rng.uniform(-SCATTER, SCATTER)
        loads[bus] = complex(load.real * active * jitter, load.imag * reactive * jitter)
    # Generation follows load. Without this the slack machine absorbs the whole
    # load increase and is rejected for exceeding its own rating, which is an
    # artefact of the dispatch model rather than an infeasible operating point.
    pv = {}
    for bus, spec in base_network.pv.items():
        factor = active * (1.0 + rng.uniform(-DISPATCH, DISPATCH))
        if bus in CORE:
            factor *= availability
        pv[bus] = {**spec, "p": spec["p"] * factor}
    return dataclass_replace(base_network, loads=loads, pv=pv)


def pick_anchor(spectrum, case, nominal):
    """Base-case mode closest to the nominal inter-area branch, 0.3 to 1.5 Hz."""

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


def track(reference, ref_labels, ref_index, case, spectrum):
    common = sorted(ref_labels & machine_labels(case.system.labels))
    if len(common) < 4:
        return None
    left = reference.right[[ref_index[n] for n in common]]
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
    parser = argparse.ArgumentParser(description="Track-A gate 5")
    parser.add_argument("--samples", type=int, default=500)
    args = parser.parse_args()
    started = time.time()

    base_network = load_network()
    engine = qmc.LatinHypercube(d=3, seed=SEED)
    unit = engine.random(args.samples)
    lower = np.array(
        [BOUNDS[k][0] for k in ("active_load", "reactive_load", "pv_availability")]
    )
    upper = np.array(
        [BOUNDS[k][1] for k in ("active_load", "reactive_load", "pv_availability")]
    )
    samples = qmc.scale(unit, lower, upper)
    rng = np.random.default_rng(SEED)

    # The nominal inter-area eigenvector, expressed by state name so that it can
    # be compared across configurations with different state vectors.
    nominal_base = solve_case(ReplacementPlan.of({}))
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
    del nominal_base

    plan_q = ReplacementPlan.of({b: 1.0 for b in CORE})
    plan_pf = ReplacementPlan.of({b: 1.0 for b in CORE}, q_policy="unity_pf")
    voltage = ConverterParameters(voltage_control=True)

    rows = []
    for index, sample in enumerate(samples):
        network = perturbed_network(base_network, sample, rng)
        entry = {
            "sample": index,
            "active_load": float(sample[0]),
            "reactive_load": float(sample[1]),
            "pv_availability": float(sample[2]),
        }
        try:
            base_case = solve_case(ReplacementPlan.of({}), network=network)
            q_case = solve_case(plan_q, network=network)
            v_case = solve_case(plan_q, network=network, converter=voltage)
            pf_case = solve_case(plan_pf, network=network)
        except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
            entry.update({"status": "REJECTED", "reason": str(error)[:120]})
            rows.append(entry)
            continue

        voltages = np.abs(base_case.dae.power_flow.voltages)
        if voltages.min() < MIN_VOLTAGE or voltages.max() > MAX_VOLTAGE:
            entry.update({"status": "REJECTED", "reason": "voltage limits"})
            rows.append(entry)
            continue

        # Anchor: the mode of THIS sample's base case closest to the nominal
        # inter-area eigenvector. Anchoring on the replacement case's critical
        # mode is wrong when that case is stable, because its critical mode is
        # then not necessarily the branch under study.
        base_spectrum = eigen_analysis(base_case.system.A)
        anchor = pick_anchor(base_spectrum, base_case, nominal)
        if anchor is None:
            entry.update({"status": "REJECTED", "reason": "no anchor in the band"})
            rows.append(entry)
            continue
        ref_labels = machine_labels(base_case.system.labels)
        ref_index = {n: i for i, n in enumerate(base_case.system.labels)}
        q_critical = anchor

        tracked = {}
        ok = True
        for name, case in (
            ("base", base_case),
            ("q", q_case),
            ("v", v_case),
            ("pf", pf_case),
        ):
            spectrum = eigen_analysis(case.system.A)
            found = track(q_critical, ref_labels, ref_index, case, spectrum)
            if found is None or found[1] < MAC_FLOOR:
                ok = False
                break
            tracked[name] = found
            entry[f"alpha_{name}"] = found[0].real
            entry[f"freq_{name}"] = found[0].frequency_hz
            entry[f"mac_{name}"] = found[1]
            entry[f"abscissa_{name}"] = assess(
                spectrum, case.system.labels
            ).spectral_abscissa
        if not ok:
            entry.update(
                {"status": "REJECTED", "reason": "mode tracking below MAC floor"}
            )
            rows.append(entry)
            continue

        entry["delta_alpha_q"] = entry["alpha_q"] - entry["alpha_base"]
        entry["delta_alpha_v"] = entry["alpha_v"] - entry["alpha_q"]
        entry["h1"] = bool(entry["delta_alpha_q"] > 0.0)
        entry["h2"] = bool(entry["delta_alpha_v"] < 0.0)
        entry["h3"] = bool(entry["alpha_pf"] > entry["alpha_q"] > entry["alpha_v"])
        entry["status"] = "ACCEPTED"
        rows.append(entry)

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_samples.csv", index=False)
    accepted = table[table.status == "ACCEPTED"]

    statistics = {}
    for name in ("h1", "h2", "h3"):
        point, low, high = bootstrap_ci(accepted[name].to_numpy(), rng)
        statistics[name] = {"fraction": point, "ci_low": low, "ci_high": high}

    passes = bool(
        statistics["h1"]["ci_low"] > 0.80 and statistics["h2"]["ci_low"] > 0.80
    )

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=SEED,
        config={
            "protocol": "configs/ias2026/trackA_final_v1.yaml",
            "samples": args.samples,
            "bounds": BOUNDS,
            "mac_floor": MAC_FLOOR,
        },
    )
    manifest.finish(
        "SUCCESS" if passes else "DIRECTION_DOES_NOT_GENERALIZE",
        drawn=len(table),
        accepted=len(accepted),
        rejection_rate=float(1.0 - len(accepted) / max(len(table), 1)),
        rejection_reasons=table[table.status == "REJECTED"]
        .reason.value_counts()
        .to_dict()
        if (table.status == "REJECTED").any()
        else {},
        statistics=statistics,
        median_delta_alpha_q=float(accepted.delta_alpha_q.median()),
        median_delta_alpha_v=float(accepted.delta_alpha_v.median()),
        q_unstable_fraction=float((accepted.alpha_q > 0).mean()),
        v_unstable_fraction=float((accepted.alpha_v > 0).mean()),
        min_mac=float(
            accepted[["mac_base", "mac_q", "mac_v", "mac_pf"]].to_numpy().min()
        ),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("held-out Monte Carlo, protocol frozen in trackA_final_v1.yaml")
    print(
        "  drawn %d, accepted %d, rejection rate %.1f%%"
        % (len(table), len(accepted), 100 * (1 - len(accepted) / max(len(table), 1)))
    )
    if (table.status == "REJECTED").any():
        print(
            "  rejection reasons: %s"
            % table[table.status == "REJECTED"].reason.value_counts().to_dict()
        )
    print()
    print(
        "  minimum MAC over every accepted sample and configuration: %.4f"
        % accepted[["mac_base", "mac_q", "mac_v", "mac_pf"]].to_numpy().min()
    )
    print()
    print("  tracked branch real part, median over accepted samples")
    for name, label in (
        ("base", "reference"),
        ("q", "Q setpoint"),
        ("pf", "unity power factor"),
        ("v", "voltage control"),
    ):
        block = accepted[f"alpha_{name}"]
        print(
            "    %-22s %+9.5f   (5th %+9.5f, 95th %+9.5f)  unstable %5.1f%%"
            % (
                label,
                block.median(),
                block.quantile(0.05),
                block.quantile(0.95),
                100 * (block > 0).mean(),
            )
        )
    print()
    print("  FROZEN DIRECTIONAL HYPOTHESES")
    names = {
        "h1": "H1  Q setpoint degrades the branch      d_alpha_Q > 0",
        "h2": "H2  voltage control mitigates it        d_alpha_V < 0",
        "h3": "H3  ordering  unityPF > Qset > Vctrl",
    }
    for key, label in names.items():
        s = statistics[key]
        print(
            "    %-46s %.3f  95%% CI [%.3f, %.3f]"
            % (label, s["fraction"], s["ci_low"], s["ci_high"])
        )
    print()
    print(
        "    median d_alpha_Q = %+.5f, median d_alpha_V = %+.5f"
        % (accepted.delta_alpha_q.median(), accepted.delta_alpha_v.median())
    )
    print()
    print(
        "GATE 5 PASSED, the direction of the effect and of the repair generalize"
        if passes
        else "GATE 5 FAILED, the direction does not generalize; downgrade"
        " the claim to the discovery operating point"
    )
    print(f"manifest -> {path}")
    return 0 if passes else 1


if __name__ == "__main__":
    raise SystemExit(main())
