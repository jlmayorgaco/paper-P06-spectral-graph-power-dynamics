"""E15 - PHASE E2. Exact order-4 audit of the Track-A flagship.

Question
    Does the flagship satisfy the order-4 claim on ONE tracked mode branch, and
    is the irreducible fourth-order term genuinely required for the crossing?

Method
    All 16 subsets of the core. The flagship critical branch is followed by MAC
    on the machine states that survive in the flagship, which are a subset of the
    survivors of every proper subset, so the comparison basis is well defined
    throughout. The Moebius decomposition is then taken of the TRACKED
    eigenvalue real part, not of the raw spectral abscissa: the abscissa is a
    maximum over modes and mixes branches.

Acceptance
    alpha(<=1) < 0, alpha(<=2) < 0, alpha(<=3) < 0, alpha(exact) > 0, with the
    same tracked branch and MAC above a declared floor. Otherwise AMBIGUOUS.

Status of the result
    NUMERICAL OBSERVATION at one operating point, rho = 1.

Usage
    python experiments/E15_order4_audit.py [--core 30 33 35 37]
"""

from __future__ import annotations

import argparse
import time
from itertools import combinations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.cycles.mobius import decompose
from ibr_cycles.diagnosis.screening import assess
from ibr_cycles.dynamics.modal_tracking import modal_assurance
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case

EXPERIMENT = "E15_order4_audit"
MAC_FLOOR = 0.80


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E2 order-4 audit")
    parser.add_argument("--core", type=int, nargs=4, default=[30, 33, 35, 37])
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    core = tuple(args.core)
    started = time.time()

    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in core}))
    spectrum = eigen_analysis(flagship.system.A)
    dynamic = [m for m in spectrum.modes if abs(m.value) > 1e-3]
    reference = max(dynamic, key=lambda m: m.real)
    labels = [
        n for n in flagship.system.labels if n.startswith(("delta_sg", "omega_sg"))
    ]
    index = {n: i for i, n in enumerate(flagship.system.labels)}
    left = reference.right[[index[n] for n in labels]]

    rows = []
    tracked: dict[frozenset[str], float] = {}
    for size in range(len(core) + 1):
        for members in combinations(core, size):
            case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
            local = eigen_analysis(case.system.A)
            verdict = assess(local, case.system.labels)
            position = {n: i for i, n in enumerate(case.system.labels)}
            select = [position[n] for n in labels]
            candidates = [m for m in local.modes if abs(m.value) > 1e-3]
            best = max(candidates, key=lambda m: modal_assurance(left, m.right[select]))
            mac = modal_assurance(left, best.right[select])
            dominant = case.system.labels[int(np.argmax(best.participation))]
            rows.append(
                {
                    "members": "+".join(map(str, members)) or "BASE",
                    "size": size,
                    "spectral_abscissa": verdict.spectral_abscissa,
                    "tracked_real": best.real,
                    "tracked_imag": best.imag,
                    "tracked_frequency_hz": best.frequency_hz,
                    "tracked_damping": best.damping,
                    "tracked_mac": mac,
                    "tracked_dominant_state": dominant,
                    "tracked_family": "gfl" if "gfl" in dominant else "sg",
                    "mode_condition": best.condition,
                    "unstable": verdict.unstable,
                    "gz_condition": case.gz_condition,
                }
            )
            tracked[frozenset(map(str, members))] = best.real

    table = pd.DataFrame(rows)
    names = tuple(map(str, core))
    decomposition = decompose(lambda s: tracked[frozenset(s)], names)
    truncations = {f"alpha_le_{k}": decomposition.truncated(names, k) for k in range(5)}
    mu4 = decomposition.terms[frozenset(names)]

    mac_min = float(table.tracked_mac.min())
    families = set(table.tracked_family)
    accepted = bool(
        truncations["alpha_le_1"] < 0
        and truncations["alpha_le_2"] < 0
        and truncations["alpha_le_3"] < 0
        and truncations["alpha_le_4"] > 0
        and mu4 > 0
        and mac_min >= MAC_FLOOR
        and len(families) == 1
    )

    terms = pd.DataFrame(
        [
            {
                "order": len(k),
                "members": "+".join(sorted(k)) or "BASE",
                "term": v,
            }
            for k, v in sorted(
                decomposition.terms.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))
            )
        ]
    )

    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_lattice.csv", index=False)
    terms.to_csv(tables / f"{EXPERIMENT}_mobius_terms.csv", index=False)

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={"core": list(core), "mac_floor": MAC_FLOOR},
    )
    manifest.finish(
        "SUCCESS" if accepted else "AMBIGUOUS",
        reference_branch={
            "real": reference.real,
            "frequency_hz": reference.frequency_hz,
        },
        truncations=truncations,
        mu4=float(mu4),
        mac_min=mac_min,
        tracked_families=sorted(families),
        reconstruction_residual=float(
            abs(truncations["alpha_le_4"] - tracked[frozenset(names)])
        ),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print(
        "core %s   reference branch %+.6f at %.4f Hz"
        % ("+".join(map(str, core)), reference.real, reference.frequency_hz)
    )
    print()
    print(
        table[
            [
                "members",
                "size",
                "spectral_abscissa",
                "tracked_real",
                "tracked_frequency_hz",
                "tracked_damping",
                "tracked_mac",
                "tracked_family",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:10.5f}")
    )
    print()
    print("Moebius terms of the TRACKED branch real part")
    print(terms.to_string(index=False, float_format=lambda v: f"{v:12.6f}"))
    print()
    print("reconstruction truncated at each interaction order:")
    for k in range(5):
        value = truncations[f"alpha_le_{k}"]
        print(
            "  order <= %d : %+.6f  %s"
            % (k, value, "UNSTABLE" if value > 0 else "stable")
        )
    print("  irreducible fourth-order term mu4 = %+.6f" % mu4)
    print("  minimum MAC across the lattice = %.4f (floor %.2f)" % (mac_min, MAC_FLOOR))
    print("  tracked families = %s" % sorted(families))
    print()
    print("ORDER-4 CLAIM ACCEPTED" if accepted else "ORDER-4 CLAIM AMBIGUOUS")
    print(f"manifest -> {path}")
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
