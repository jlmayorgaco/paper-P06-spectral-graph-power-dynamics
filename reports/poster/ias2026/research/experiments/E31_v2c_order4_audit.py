"""E31 - v2C item 10. Does the Track-A order-4 claim survive the family definition?

E15 established the order-4 claim on a SINGLE argmax-MAC branch. v2B showed that
selector is defective once the inter-area reference acquires two descendants, so
every earlier result that depended on it has to be re-audited rather than
assumed. This re-runs the whole 16-subset lattice with the frozen v2C observable:

    alpha_IA(S) = max over the descendant family of the base inter-area mode

The comparison basis is pinned to the flagship survivors for every subset, so the
overlap threshold means the same thing at every node of the lattice.

The claim under audit, from E15:

    reconstruction through order 1 stable
    reconstruction through order 2 stable
    reconstruction through order 3 stable
    exact, i.e. through order 4, unstable

If it survives, the Track-A higher-order result holds under a strictly stronger
modal definition. If it does not, the scalar Moebius order claim is downgraded
immediately. The branch-independent spectral facts are reported beside it either
way.

Usage
    python experiments/E31_v2c_order4_audit.py
"""

from __future__ import annotations

import time
from itertools import combinations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
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
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

EXPERIMENT = "E31_v2c_order4_audit"
ZERO_MODE = 1e-3


def main() -> int:
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in CORE}))
    pinned = sorted(machine_labels(flagship.system.labels))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)
    if anchor is None:
        print("no base anchor in the band")
        return 1

    rows = []
    envelope: dict[frozenset[str], float] = {}
    single: dict[frozenset[str], float] = {}
    for size in range(len(CORE) + 1):
        for members in combinations(CORE, size):
            case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
            spectrum = eigen_analysis(case.system.A)
            reading = read_family(
                anchor, base, case, spectrum, common_labels=pinned
            )
            key = frozenset(map(str, members))
            row = {
                "members": "+".join(map(str, members)) or "BASE",
                "size": size,
                "spectral_abscissa": float(
                    max(m.real for m in spectrum.modes if abs(m.value) > ZERO_MODE)
                ),
                "rhp_count": int(
                    sum(
                        1
                        for m in spectrum.modes
                        if abs(m.value) > ZERO_MODE and m.real > 0.0 and m.imag >= 0.0
                    )
                ),
                "band_rhp_count": int(
                    sum(1 for m in band_candidates(spectrum, BAND_HZ) if m.real > 0.0)
                ),
            }
            if reading is None:
                row["family_size"] = 0
                row["tracking"] = "TRACKING_FAILURE"
                rows.append(row)
                continue
            family = reading.family
            row.update(
                {
                    "tracking": "ok",
                    "family_size": family.size,
                    "alpha_IA": family.alpha,
                    "freq_worst_hz": family.frequency_worst_hz,
                    "overlap_min": min(family.overlaps),
                    "overlap_max": max(family.overlaps),
                    "independence": family.independence,
                    "max_mode_condition": family.max_condition,
                    "descendant_reals": ";".join(f"{v:.6f}" for v in family.reals),
                    "descendant_freqs_hz": ";".join(
                        f"{v:.4f}" for v in family.frequencies_hz
                    ),
                    "argmax_alpha": reading.argmax_real,
                    "argmax_freq_hz": reading.argmax_frequency_hz,
                }
            )
            envelope[key] = family.alpha
            single[key] = reading.argmax_real
            rows.append(row)

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_lattice.csv", index=False)

    names = tuple(map(str, CORE))
    if len(envelope) != 16:
        print("lattice incomplete: %d of 16 subsets tracked" % len(envelope))
        return 1

    family_decomposition = decompose(lambda s: envelope[frozenset(s)], names)
    single_decomposition = decompose(lambda s: single[frozenset(s)], names)
    truncations = {
        k: family_decomposition.truncated(names, k) for k in range(5)
    }
    single_truncations = {
        k: single_decomposition.truncated(names, k) for k in range(5)
    }
    mu4 = float(family_decomposition.terms[frozenset(names)])
    single_mu4 = float(single_decomposition.terms[frozenset(names)])

    terms = pd.DataFrame(
        [
            {
                "order": len(k),
                "members": "+".join(sorted(k)) or "BASE",
                "term_family": v,
                "term_single_branch": single_decomposition.terms[k],
            }
            for k, v in sorted(
                family_decomposition.terms.items(),
                key=lambda kv: (len(kv[0]), sorted(kv[0])),
            )
        ]
    )
    terms.to_csv(tables / f"{EXPERIMENT}_mobius_terms.csv", index=False)

    tracked = table[table.tracking == "ok"]
    overlap_min = float(tracked.overlap_min.min())
    survives = bool(
        truncations[1] < 0
        and truncations[2] < 0
        and truncations[3] < 0
        and truncations[4] > 0
        and mu4 > 0
        and overlap_min >= FAMILY_THRESHOLD
    )
    residual = abs(truncations[4] - envelope[frozenset(names)])

    print("E31 order-4 audit under the inter-area FAMILY definition")
    print("  base anchor %+.6f at %.4f Hz; comparison basis pinned to the %d"
          " flagship survivors" % (anchor.real, anchor.frequency_hz, len(pinned)))
    print()
    print(table[[
        "members", "size", "family_size", "alpha_IA", "freq_worst_hz",
        "overlap_min", "argmax_alpha", "spectral_abscissa", "band_rhp_count",
    ]].to_string(index=False, float_format=lambda v: f"{v:10.5f}"))
    print()
    print("  family size across the lattice: %s"
          % dict(sorted(tracked.family_size.value_counts().items())))
    print("  minimum member overlap: %.4f (threshold %.2f)" % (overlap_min, FAMILY_THRESHOLD))
    print("  tracking failures: %d" % int((table.tracking != "ok").sum()))
    print()
    print("Moebius reconstruction of alpha_IA, truncated at each interaction order")
    print("   order   family envelope        single branch (v2B observable)")
    for k in range(5):
        print("   <= %d   %+12.6f  %-8s  %+12.6f  %s"
              % (k, truncations[k], "UNSTABLE" if truncations[k] > 0 else "stable",
                 single_truncations[k],
                 "UNSTABLE" if single_truncations[k] > 0 else "stable"))
    print("   irreducible fourth-order term  mu4 = %+.6f (single branch %+.6f)"
          % (mu4, single_mu4))
    print("   reconstruction residual %.2e" % residual)
    print()
    print("ORDER-4 CLAIM SURVIVES THE FAMILY DEFINITION" if survives
          else "ORDER-4 CLAIM DOES NOT SURVIVE - DOWNGRADE IT")

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=0,
        config={
            "core": list(CORE),
            "observable": "alpha_IA family envelope",
            "family_threshold": FAMILY_THRESHOLD,
            "pinned_basis": len(pinned),
        },
    )
    manifest.finish(
        "SUCCESS" if survives else "DOWNGRADE",
        survives=survives,
        truncations={f"alpha_le_{k}": float(v) for k, v in truncations.items()},
        truncations_single_branch={
            f"alpha_le_{k}": float(v) for k, v in single_truncations.items()
        },
        mu4=mu4,
        mu4_single_branch=single_mu4,
        overlap_min=overlap_min,
        family_sizes={
            int(k): int(v) for k, v in tracked.family_size.value_counts().items()
        },
        reconstruction_residual=float(residual),
        band_rhp_by_size={
            int(k): sorted(map(int, v))
            for k, v in table.groupby("size").band_rhp_count.apply(set).items()
        },
        elapsed_s=time.time() - started,
    )
    print(f"manifest -> {manifest.write(RESULTS / 'manifests')}")
    return 0 if survives else 1


if __name__ == "__main__":
    raise SystemExit(main())
