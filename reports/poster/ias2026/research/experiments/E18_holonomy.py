"""E18 - PHASE E2. Closed-cycle holonomies of the Track-A replacement core.

Question
    E17 established that the closure condition at the created mode is
    representation invariant. Can it be attributed to a closed structure over the
    four replacement ports, and does that structure distinguish the flagship from
    stable four-replacement portfolios matched on size and megawatts?

Objects
    With ``M_ab(s) = dY_a(s) K_ab(s)`` the directed cycle
    ``gamma: a1 -> a2 -> ... -> ap -> a1`` carries

        H_gamma(s) = M_a1a2 M_a2a3 ... M_apa1

    Under an admissible port basis change ``M -> S^-1 M S``, so ``H_gamma``
    transforms by similarity within its base port and its trace, determinant and
    eigenvalues are invariant. Entries and norms are not, and are not reported.

Discriminating test
    If a cycle score is equally large for stable matched portfolios it explains
    nothing. The comparison against matched controls is the point of this
    experiment, not the flagship numbers on their own.

Status of the result
    NUMERICAL OBSERVATION at one operating point, rho = 1.

Usage
    python experiments/E18_holonomy.py [--core 30 33 35 37]
"""

from __future__ import annotations

import argparse
import time
from itertools import combinations, permutations

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.port_admittance import build_action_space

EXPERIMENT = "E18_holonomy"
BAND_HZ = (0.1, 5.0)
N_MATCHED = 12


def directed_cycles(order: int) -> list[tuple[int, ...]]:
    """Every simple directed cycle of length 2 to ``order``, once per orbit."""

    cycles: list[tuple[int, ...]] = []
    for length in range(2, order + 1):
        for members in combinations(range(order), length):
            head, rest = members[0], members[1:]
            for tail in permutations(rest):
                cycle = (head, *tail)
                if length > 2 and cycle[::-1][-1:] + cycle[::-1][:-1] in cycles:
                    continue
                cycles.append(cycle)
    return cycles


def holonomy(space, m: np.ndarray, cycle: tuple[int, ...]) -> np.ndarray:
    product = None
    for position, node in enumerate(cycle):
        following = cycle[(position + 1) % len(cycle)]
        block = space.block(m, node, following)
        product = block if product is None else product @ block
    return product


def scores(space, m: np.ndarray, cycle: tuple[int, ...]) -> dict[str, float]:
    h = holonomy(space, m, cycle)
    eigenvalues = np.linalg.eigvals(h)
    return {
        "abs_trace": float(abs(np.trace(h))),
        "abs_det": float(abs(np.linalg.det(h))),
        "spectral_radius": float(np.max(np.abs(eigenvalues))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E2 holonomies")
    parser.add_argument("--core", type=int, nargs="+", default=[30, 33, 35, 37])
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    core = tuple(args.core)
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in core}))
    space = build_action_space(base, flagship, core)
    created = max(
        (m for m in eigen_analysis(flagship.system.A).modes if abs(m.value) > 1e-3),
        key=lambda m: m.real,
    ).value
    cycles = directed_cycles(len(core))
    names = {i: core[i] for i in range(len(core))}

    # --- cycle spectrum at the created mode --------------------------------
    m_created = space.m(created)
    rows = []
    for cycle in cycles:
        rows.append(
            {
                "cycle": "->".join(str(names[i]) for i in cycle)
                + f"->{names[cycle[0]]}",
                "length": len(cycle),
                **scores(space, m_created, cycle),
            }
        )
    at_mode = pd.DataFrame(rows).sort_values("spectral_radius", ascending=False)

    # --- gauge invariance of the cycle invariants ---------------------------
    rng = np.random.default_rng(args.seed)
    bases = []
    for _ in core:
        candidate = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        bases.append(candidate / np.linalg.norm(candidate))
    moved = space.with_bases(bases)
    m_moved = moved.m(created)
    drift = []
    for cycle in cycles:
        a = scores(space, m_created, cycle)
        b = scores(moved, m_moved, cycle)
        drift.append(
            max(
                abs(b[k] - a[k]) / max(a[k], 1e-300)
                for k in ("abs_trace", "abs_det", "spectral_radius")
            )
        )
    gauge_drift = float(max(drift))

    # --- frequency sweep of the longest cycle -------------------------------
    sigma = created.real
    frequencies = np.unique(
        np.concatenate(
            [
                np.linspace(BAND_HZ[0], BAND_HZ[1], 60),
                np.linspace(
                    max(BAND_HZ[0], abs(created.imag) / (2 * np.pi) - 0.15),
                    abs(created.imag) / (2 * np.pi) + 0.15,
                    40,
                ),
            ]
        )
    )
    sweep_rows = []
    for f in frequencies:
        s = complex(sigma, 2.0 * np.pi * f)
        m = space.m(s)
        entry = {"frequency_hz": float(f)}
        for length in range(2, len(core) + 1):
            same = [c for c in cycles if len(c) == length]
            entry[f"max_rho_len{length}"] = max(
                scores(space, m, c)["spectral_radius"] for c in same
            )
        split = space.split(s)
        entry["distance_of_q_eig_to_minus_one"] = float(
            abs(split["closest_to_minus_one"] + 1.0)
        )
        sweep_rows.append(entry)
    sweep = pd.DataFrame(sweep_rows)

    # --- discrimination against matched stable portfolios -------------------
    census = pd.read_csv(RESULTS / "tables" / "E12_compatibility_census_portfolios.csv")
    target_mw = float(
        census[census.members == "+".join(map(str, core))].replaced_mw.iloc[0]
    )
    stable = census[(census["size"] == 4) & (~census.unstable)].copy()
    stable["gap"] = (stable.replaced_mw - target_mw).abs()
    matched = stable.nsmallest(N_MATCHED, "gap")

    comparison = []
    for row in matched.itertuples():
        members = tuple(int(b) for b in row.members.split("+"))
        case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
        other = build_action_space(base, case, members)
        m_other = other.m(created)
        record = {
            "portfolio": row.members,
            "replaced_mw": row.replaced_mw,
            "unstable": False,
        }
        for length in range(2, len(members) + 1):
            same = [c for c in cycles if len(c) == length]
            record[f"max_rho_len{length}"] = max(
                scores(other, m_other, c)["spectral_radius"] for c in same
            )
        split = other.split(created)
        record["distance_of_q_eig_to_minus_one"] = float(
            abs(split["closest_to_minus_one"] + 1.0)
        )
        comparison.append(record)
    flagship_record = {
        "portfolio": "+".join(map(str, core)),
        "replaced_mw": target_mw,
        "unstable": True,
        "distance_of_q_eig_to_minus_one": float(
            abs(space.split(created)["closest_to_minus_one"] + 1.0)
        ),
    }
    for length in range(2, len(core) + 1):
        same = [c for c in cycles if len(c) == length]
        flagship_record[f"max_rho_len{length}"] = max(
            scores(space, m_created, c)["spectral_radius"] for c in same
        )
    comparison.insert(0, flagship_record)
    contrast = pd.DataFrame(comparison)

    tables = RESULTS / "tables"
    at_mode.to_csv(tables / f"{EXPERIMENT}_cycles_at_mode.csv", index=False)
    sweep.to_csv(tables / f"{EXPERIMENT}_frequency_sweep.csv", index=False)
    contrast.to_csv(tables / f"{EXPERIMENT}_matched_contrast.csv", index=False)

    length_column = f"max_rho_len{len(core)}"
    flagship_value = float(contrast.iloc[0][length_column])
    controls = contrast.iloc[1:]
    separation = flagship_value / max(float(controls[length_column].max()), 1e-300)
    closure_separation = float(controls.distance_of_q_eig_to_minus_one.min()) / max(
        flagship_record["distance_of_q_eig_to_minus_one"], 1e-300
    )
    discriminative = bool(separation > 2.0)

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={"core": list(core), "band_hz": list(BAND_HZ), "matched": N_MATCHED},
    )
    manifest.finish(
        "SUCCESS",
        created_mode={"real": created.real, "imag": created.imag},
        cycles=len(cycles),
        gauge_drift_of_cycle_invariants=gauge_drift,
        strongest_cycle_at_mode=at_mode.iloc[0].to_dict(),
        flagship_len4_spectral_radius=flagship_value,
        matched_len4_max=float(controls[length_column].max()),
        cycle_separation_ratio=separation,
        closure_separation_ratio=closure_separation,
        cycle_score_is_discriminative=discriminative,
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print(
        "created mode %+.6f%+.6fj   %d directed cycles over %d ports"
        % (created.real, created.imag, len(cycles), len(core))
    )
    print("  gauge drift of trace, determinant and spectral radius: %.3e" % gauge_drift)
    print()
    print("cycle invariants at the created mode, strongest first")
    print(at_mode.to_string(index=False, float_format=lambda v: f"{v:12.5f}"))
    print()
    print(
        "flagship against %d matched stable four-replacement portfolios" % len(controls)
    )
    print(contrast.to_string(index=False, float_format=lambda v: f"{v:12.5f}"))
    print()
    print("  flagship max rho of a %d-cycle       %.5f" % (len(core), flagship_value))
    print("  largest among matched controls      %.5f" % controls[length_column].max())
    print("  separation ratio                    %.2f" % separation)
    print(
        "  flagship distance of a Q eigenvalue to -1   %.3e"
        % flagship_record["distance_of_q_eig_to_minus_one"]
    )
    print(
        "  smallest such distance among controls      %.3e"
        % controls.distance_of_q_eig_to_minus_one.min()
    )
    print()
    if discriminative:
        print("CYCLE SCORE DISCRIMINATES THE FLAGSHIP")
    else:
        print(
            "CYCLE SCORE DOES NOT DISCRIMINATE - do not use it as an explanatory metric"
        )
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
