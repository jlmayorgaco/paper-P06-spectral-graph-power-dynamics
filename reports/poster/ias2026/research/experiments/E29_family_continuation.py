"""E29 - v2C step 1. Deterministic continuation check of the inter-area FAMILY.

This runs BEFORE the v2C protocol is frozen and BEFORE any v2C seed is drawn. It
is deterministic, uses no random sample, and its only purpose is to establish
whether the descendant SET is a stable object under a subspace metric where the
individual members are not.

Three questions:

  1. how many descendants does the frozen base inter-area reference acquire, and
     does that number stay constant along a load continuation?
  2. is the subspace the family spans continuous along the continuation, measured
     by the worst principal-angle cosine between consecutive steps, while the
     individual assurance values chatter?
  3. is the region abs(alpha_IA) <= 0.05 actually populated under the frozen
     sampling domain? This decides between the two candidate near-boundary
     strata of the v2C preregistration.

Usage
    python experiments/E29_family_continuation.py
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from _v2c_common import (
    BAND_HZ,
    CORE,
    FAMILY_THRESHOLD,
    base_anchor,
    nominal_reference,
    read_family,
)
from ibr_cycles.dynamics.modal_family import subspace_similarity
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.uncertainty.sampling import sample_operating_point

EXPERIMENT = "E29_family_continuation"
THRESHOLD = FAMILY_THRESHOLD
AVAILABILITIES = (0.75, 0.85, 0.92, 1.00)
LOADS = np.round(np.arange(0.900, 1.0501, 0.0025), 6)


def main() -> int:
    started = time.time()
    base_network = load_network()
    plan_q = ReplacementPlan.of({b: 1.0 for b in CORE})
    voltage = ConverterParameters(voltage_control=True)
    nominal = nominal_reference()

    rows = []
    for availability in AVAILABILITIES:
        previous_basis = None
        previous_load = None
        for load in LOADS:
            rng = np.random.default_rng(0)  # deterministic: no jitter, no scatter
            point = sample_operating_point(
                base_network,
                active_load=float(load),
                reactive_load=1.0,
                availability=float(availability),
                availability_buses=CORE,
                rng=rng,
                dispatch_jitter=0.0,
                load_scatter=0.0,
            )
            row = {
                "availability": availability,
                "load": load,
                "dispatchable": bool(point.dispatchable),
            }
            if not point.dispatchable:
                rows.append(row)
                previous_basis, previous_load = None, None
                continue
            try:
                base = solve_case(ReplacementPlan.of({}), network=point.network)
                replaced = solve_case(plan_q, network=point.network)
                regulated = solve_case(
                    plan_q, network=point.network, converter=voltage
                )
            except Exception as error:  # noqa: BLE001 - recorded, not swallowed
                row["reason"] = str(error)[:80]
                rows.append(row)
                previous_basis, previous_load = None, None
                continue

            base_spectrum = eigen_analysis(base.system.A)
            anchor = base_anchor(base_spectrum, base, nominal)
            if anchor is None:
                row["reason"] = "no base anchor in the band"
                rows.append(row)
                continue

            for name, case, spectrum in (
                ("q", replaced, eigen_analysis(replaced.system.A)),
                ("v", regulated, eigen_analysis(regulated.system.A)),
            ):
                reading = read_family(anchor, base, case, spectrum)
                if reading is None:
                    row[f"{name}_size"] = 0
                    continue
                family = reading.family
                row[f"{name}_size"] = family.size
                row[f"alpha_{name}_family"] = family.alpha
                row[f"freq_{name}_worst"] = family.frequency_worst_hz
                row[f"{name}_independence"] = family.independence
                row[f"{name}_overlap_min"] = min(family.overlaps)
                row[f"{name}_overlap_max"] = max(family.overlaps)
                # the single-branch observable that v2B used, for contrast
                row[f"alpha_{name}_argmax"] = reading.argmax_real
                row[f"freq_{name}_argmax"] = reading.argmax_frequency_hz
                if name == "q":
                    if previous_basis is not None and previous_load is not None:
                        row["subspace_similarity_to_previous"] = subspace_similarity(
                            previous_basis, family.basis
                        )
                        row["previous_load"] = previous_load
                    previous_basis, previous_load = family.basis, load
            row["alpha_base"] = anchor.real
            row["freq_base"] = anchor.frequency_hz
            rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_continuation.csv", index=False)
    served = table[table.dispatchable & table.q_size.notna()]

    sizes = served.q_size.value_counts().sort_index()
    similarity = served.subspace_similarity_to_previous.dropna()
    swaps = served[
        np.sign(served.alpha_q_argmax) != np.sign(served.alpha_q_family)
    ]
    near = served[served.alpha_q_family.abs() <= 0.05]
    near_argmax = served[served.alpha_q_argmax.abs() <= 0.05]

    print("E29 deterministic family continuation, %d dispatchable points" % len(served))
    print("  loads %.3f to %.3f in steps of 0.0025, availabilities %s"
          % (LOADS.min(), LOADS.max(), list(AVAILABILITIES)))
    print()
    print("1. family size after replacement")
    for size, count in sizes.items():
        print("     %d descendants : %4d points" % (int(size), int(count)))
    print("     base case family size: %s"
          % sorted(set(served.get("base_size", pd.Series(dtype=float)).dropna())))
    print()
    print("2. subspace continuity along the continuation (consecutive steps)")
    print("     worst principal-angle cosine: min %.6f, 1st pct %.6f, median %.6f"
          % (similarity.min(), similarity.quantile(0.01), similarity.median()))
    print("     member independence (smallest singular value): min %.4f, median %.4f"
          % (served.q_independence.min(), served.q_independence.median()))
    print("     per-member overlap to the base reference: min %.4f, max %.4f"
          % (served.q_overlap_min.min(), served.q_overlap_max.max()))
    print("     points where the argmax selector disagrees in SIGN with the"
          " envelope: %d of %d" % (len(swaps), len(served)))
    print()
    print("3. is abs(alpha_IA) <= 0.05 populated?")
    print("     family envelope : %d of %d points" % (len(near), len(served)))
    print("     argmax selector : %d of %d points" % (len(near_argmax), len(served)))
    print()
    for availability in AVAILABILITIES:
        block = served[served.availability == availability]
        if block.empty:
            continue
        sign_change = block[block.alpha_q_family > 0]
        crossing = sign_change.load.min() if not sign_change.empty else np.nan
        print("     availability %.2f : envelope crosses zero near load %.4f,"
              " family size %s"
              % (availability, crossing, sorted(set(block.q_size.astype(int)))))

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=0,
        config={
            "band_hz": list(BAND_HZ),
            "threshold": THRESHOLD,
            "deterministic": True,
            "purpose": "pre-registration evidence for v2C; no random sample",
        },
    )
    manifest.finish(
        "REPORTED",
        dispatchable_points=len(served),
        family_sizes={int(k): int(v) for k, v in sizes.items()},
        subspace_similarity_min=float(similarity.min()),
        subspace_similarity_median=float(similarity.median()),
        independence_min=float(served.q_independence.min()),
        argmax_sign_disagreements=int(len(swaps)),
        near_boundary_points_envelope=int(len(near)),
        near_boundary_points_argmax=int(len(near_argmax)),
        elapsed_s=time.time() - started,
    )
    print(f"manifest -> {manifest.write(RESULTS / 'manifests')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
