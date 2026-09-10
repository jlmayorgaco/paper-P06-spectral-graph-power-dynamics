"""E30 - overnight 1. Does the order-4 mechanism survive realistic Q dispatch?

The frozen campaigns hold the operating point fixed by giving the converter the
reactive output the machine was producing. That isolates the mechanism and it is
not how a plant is dispatched. This runs the complete 16-subset lattice of the
flagship under four reactive policies, frozen in
``configs/ias2026/overnight_policies.yaml`` before any case was solved.

Observable: the frozen v2C inter-area modal family, unchanged, with the
comparison basis pinned to the flagship survivors so the overlap threshold means
the same thing at every node of every lattice.

Classification per policy:
    A SAME MECHANISM, B SAME FAMILY DIFFERENT SEVERITY, C DIFFERENT FAMILY,
    D MECHANISM DISAPPEARS, E INFEASIBLE OPERATING POINT
"""

from __future__ import annotations

import time
from itertools import combinations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS  # noqa: F401  (path bootstrap)
from _overnight import Experiment, pin_blas_threads
from _v2c_common import (
    BAND_HZ,
    CORE,
    FAMILY_THRESHOLD,
    base_anchor,
    machine_labels,
    nominal_reference,
    read_family,
)
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.dynamics.modal_family import band_candidates
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.qpolicy import POLICIES, solve_policy

NAME = "E30_rebalanced_policy_lattice"
ZERO_MODE = 1e-3


def describe(case, anchor, base_case, pinned, base_voltages, network):
    spectrum = eigen_analysis(case.system.A)
    reading = read_family(anchor, base_case, case, spectrum, common_labels=pinned)
    voltages = np.abs(case.dae.power_flow.voltages)
    row = {
        "spectral_abscissa": float(
            max(m.real for m in spectrum.modes if abs(m.value) > ZERO_MODE)
        ),
        "rhp_count_total": int(
            sum(
                1
                for m in spectrum.modes
                if abs(m.value) > ZERO_MODE and m.real > 0.0 and m.imag >= 0.0
            )
        ),
        "band_rhp_count": int(
            sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0.0)
        ),
        "min_voltage": float(voltages.min()),
        "max_voltage": float(voltages.max()),
        "operating_point_drift": float(
            np.abs(np.abs(case.dae.power_flow.voltages) - base_voltages).max()
        ),
        "max_device_loading": float(case.max_loading),
        "gz_condition": float(case.gz_condition),
    }
    if reading is None:
        row["tracking"] = "TRACKING_FAILURE"
        return row, None
    family = reading.family
    row.update(
        {
            "tracking": "ok",
            "alpha_IA": family.alpha,
            "freq_IA_worst_hz": family.frequency_worst_hz,
            "damping_IA_worst": float(family.worst.damping),
            "family_size": family.size,
            "overlap_min": min(family.overlaps),
            "independence": family.independence,
            "descendant_reals": ";".join(f"{v:.6f}" for v in family.reals),
            "descendant_freqs_hz": ";".join(f"{v:.4f}" for v in family.frequencies_hz),
        }
    )
    return row, family.alpha


def effective_order(truncations):
    """Lowest interaction order whose reconstruction already predicts instability.

    4 means the effect is irreducibly fourth order. A smaller number means a
    lower-order model already sees it, so the fourth-order description is not
    what makes the portfolio unstable.
    """

    for k in range(5):
        if truncations[k] > 0.0:
            return k
    return None


def classify(truncations, mu4, exact, flagship_band_rhp, families_ok, all_feasible):
    """The five declared categories.

    D is decided on whether an instability exists at all, not on whether damping
    moved. A flagship that stays in the left half plane with no right-half-plane
    mode in the band is a portfolio where the mechanism does not occur, however
    its damping compares with the base case.
    """

    if not all_feasible:
        return "E"
    if not families_ok:
        return "C"
    unstable = exact > 0.0 or flagship_band_rhp > 0
    if not unstable:
        return "D"
    if (
        truncations[1] < 0
        and truncations[2] < 0
        and truncations[3] < 0
        and truncations[4] > 0
        and mu4 > 0
    ):
        return "A"
    return "B"


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Does the Track-A order-4 mechanism survive realistic Q dispatch?",
        seed=None,
        config={
            "policies": list(POLICIES),
            "core": list(CORE),
            "capability_config": "configs/ias2026/overnight_policies.yaml",
            "observable": "frozen v2C inter-area modal family envelope",
            "family_threshold": FAMILY_THRESHOLD,
            "band_hz": list(BAND_HZ),
        },
        workers=1,
    )
    started = time.time()

    network = load_network()
    base_case = solve_case(ReplacementPlan.of({}))
    base_voltages = np.abs(base_case.dae.power_flow.voltages)
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    pinned = sorted(machine_labels(flagship.system.labels))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base_case.system.A), base_case, nominal)
    if anchor is None:
        experiment.finish("ERROR", reason="no base anchor")
        return 1
    experiment.note(
        f"base anchor {anchor.real:+.6f} at {anchor.frequency_hz:.4f} Hz; "
        f"pinned basis {len(pinned)} machine states"
    )

    subsets = [
        members
        for size in range(len(CORE) + 1)
        for members in combinations(CORE, size)
    ]

    rows = []
    per_policy: dict[str, dict[frozenset, float]] = {p: {} for p in POLICIES}
    feasible: dict[str, bool] = {p: True for p in POLICIES}
    families_ok: dict[str, bool] = {p: True for p in POLICIES}

    for policy in POLICIES:
        for members in subsets:
            label = "+".join(map(str, members)) or "BASE"
            row = {
                "policy": policy,
                "members": label,
                "size": len(members),
            }
            if not members:
                described, alpha = describe(
                    base_case, anchor, base_case, pinned, base_voltages, network
                )
                row.update(described)
                row.update(
                    {
                        "status": "ACCEPTED",
                        "saturated": "",
                        "voltage_regulating": "",
                        "worst_q_over_capability": 0.0,
                        "worst_utilisation": 0.0,
                        "q_required_pu": "",
                    }
                )
                per_policy[policy][frozenset()] = alpha
                rows.append(row)
                continue
            try:
                outcome = solve_policy(members, policy, network=network)
            except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError) as error:
                row.update(
                    {"status": "REJECTED", "reason": str(error)[:120], "tracking": "-"}
                )
                feasible[policy] = False
                rows.append(row)
                continue
            described, alpha = describe(
                outcome.case, anchor, base_case, pinned, base_voltages, network
            )
            row.update(described)
            row.update(
                {
                    "status": "ACCEPTED",
                    "saturated": ";".join(map(str, outcome.saturated)),
                    "voltage_regulating": ";".join(map(str, outcome.voltage_regulating)),
                    "worst_q_over_capability": outcome.worst_q_ratio,
                    "worst_utilisation": outcome.worst_utilisation,
                    "q_required_pu": ";".join(
                        f"{b}:{outcome.q_required_pu[b]:+.4f}"
                        for b in sorted(outcome.q_required_pu)
                    ),
                    "q_capability_pu": ";".join(
                        f"{b}:{outcome.q_capability_pu[b]:.4f}"
                        for b in sorted(outcome.q_capability_pu)
                    ),
                    "clamp_passes": outcome.clamp_passes,
                }
            )
            if alpha is None:
                families_ok[policy] = False
            else:
                per_policy[policy][frozenset(map(str, members))] = alpha
            rows.append(row)
        done = len(per_policy[policy])
        experiment.note(f"{policy}: {done} of 16 subsets tracked")

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E30_rebalanced_policy_lattice.csv")

    names = tuple(map(str, CORE))
    summary_rows = []
    for policy in POLICIES:
        values = per_policy[policy]
        if len(values) != 16:
            summary_rows.append(
                {
                    "policy": policy,
                    "classification": "E",
                    "tracked_subsets": len(values),
                    "note": "lattice incomplete; Moebius reconstruction not legitimate",
                }
            )
            continue
        decomposition = decompose(lambda s: values[frozenset(s)], names)
        truncations = {k: float(decomposition.truncated(names, k)) for k in range(5)}
        mu4 = float(decomposition.terms[frozenset(names)])
        exact = values[frozenset(names)]
        base_alpha = values[frozenset()]
        block = table[(table.policy == policy) & (table.status == "ACCEPTED")]
        proper = block[block["size"] < 4]
        flagship_rhp = int(block[block["size"] == 4].band_rhp_count.iloc[0])
        verdict = classify(
            truncations,
            mu4,
            exact,
            flagship_rhp,
            families_ok[policy],
            feasible[policy],
        )
        order = effective_order(truncations)
        summary_rows.append(
            {
                "policy": policy,
                "classification": verdict,
                "tracked_subsets": 16,
                "alpha_base": base_alpha,
                "alpha_exact": exact,
                "alpha_le_1": truncations[1],
                "alpha_le_2": truncations[2],
                "alpha_le_3": truncations[3],
                "alpha_le_4": truncations[4],
                "mu4": mu4,
                "effective_order": order if order is not None else -1,
                "kappa": 4 if verdict == "A" else np.nan,
                "reconstruction_residual": abs(truncations[4] - exact),
                "proper_subsets_stable": bool((proper.alpha_IA < 0).all()),
                "proper_subsets_band_rhp_zero": bool((proper.band_rhp_count == 0).all()),
                "flagship_band_rhp": flagship_rhp,
                "worst_proper_subset_alpha": float(proper.alpha_IA.max()),
                "unstable_proper_subsets": ";".join(
                    proper[proper.alpha_IA > 0].members.tolist()
                ),
                "family_size_set": sorted(set(block.family_size.dropna().astype(int))),
                "overlap_min": float(block.overlap_min.min()),
                "worst_drift": float(block.operating_point_drift.max()),
                "min_voltage": float(block.min_voltage.min()),
                "max_voltage": float(block.max_voltage.max()),
                "worst_loading": float(block.max_device_loading.max()),
                "any_saturated": bool((block.saturated.fillna("") != "").any()),
                "flagship_freq_hz": float(
                    block[block["size"] == 4].freq_IA_worst_hz.iloc[0]
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    experiment.save_table(summary, "E30_policy_summary.csv")

    print()
    print(summary.to_string(index=False, float_format=lambda v: f"{v:9.5f}"))

    classifications = dict(zip(summary.policy, summary.classification, strict=True))
    survives = [p for p, c in classifications.items() if c == "A"]
    engineering = {p: c for p, c in classifications.items() if p != "matched"}
    # G1 asks whether re-equilibrated engineering policies destroy the mechanism
    # ENTIRELY. A policy that still produces an unstable portfolio has not
    # destroyed it, even when the interaction order collapses.
    any_unstable = bool(
        (summary[summary.policy != "matched"].flagship_band_rhp > 0).any()
    )
    if not any_unstable:
        gate = "G1_FAIL"
    elif len(survives) == 1 and "matched" in survives:
        gate = "G1_QUALIFIED"
    else:
        gate = "PASS"

    experiment.finish(
        gate,
        classifications=classifications,
        survives_order4=survives,
        summary={r["policy"]: r for r in summary_rows},
        elapsed_s=time.time() - started,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
