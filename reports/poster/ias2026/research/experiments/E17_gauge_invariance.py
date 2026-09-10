"""E17 - PHASE E2. Is the Track-A collective result representation invariant?

Question
    E16 found that at the created mode the collective factor vanishes while the
    individual factor does not. Does that survive a change of port coordinates,
    or is it an artefact of writing the ports in rectangular components?

Admissible transformation
    An invertible 2x2 basis change at each replacement port, ``E_a -> E_a S_a``,
    with the update transforming as ``dY_a -> S_a^-1 dY_a S_a^-T`` so that the
    operator itself is untouched. Then ``M -> S^-1 M S`` block-conformally.

Predicted invariants
    det(I+M); every diagonal block determinant; their product, i.e. the
    individual factor; det(I+Q); the eigenvalues of Q; and therefore the
    closure condition "-1 is an eigenvalue of Q".

How invariance is measured, and why it matters here
    At the created mode ``det(I+M)`` and ``det(I+Q)`` are themselves about 1e-8,
    because that is what "the mode is a zero of the collective factor" means. A
    RELATIVE error on a quantity evaluated at its own zero is not a test of
    invariance, it is a test of conditioning: an ill-conditioned basis amplifies
    roundoff by the square of its condition number, so a 1e-14 perturbation shows
    up as a large relative change of a 1e-8 number. Determinant invariance is
    therefore checked at generic probe points away from the spectrum, where the
    determinants are order one. At the created mode the reported invariant is the
    one that is well scaled there: the distance from the nearest eigenvalue of Q
    to -1, measured absolutely against eigenvalues of order one.

Predicted non-invariants
    sigma_min(I+Q), the norm of Q, and the entries of Q. These are measured too,
    so that the report can state which quantities may be published.

Bases tested
    Identity; random invertible bases at condition numbers 1, 10, 100 and 1000;
    and the physically motivated polar basis, the Jacobian of (vx, vy) with
    respect to (|V|, angle) at the operating point.

Status of the result
    EXACT IDENTITY where invariance is predicted; any deviation beyond the
    finite-difference floor of the port model is an implementation error.

Usage
    python experiments/E17_gauge_invariance.py [--core 30 33 35 37]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.models.port_admittance import build_action_space

EXPERIMENT = "E17_gauge_invariance"
CONDITIONS = (1.0, 10.0, 100.0, 1000.0)
SEEDS = (0, 1, 2)


def random_basis(rng: np.random.Generator, condition: float) -> np.ndarray:
    """A 2x2 invertible basis with a prescribed condition number."""

    if condition == 1.0:
        angle = rng.uniform(0.0, 2.0 * np.pi)
        return np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]],
            dtype=np.complex128,
        )
    left = np.linalg.qr(rng.normal(size=(2, 2)))[0]
    right = np.linalg.qr(rng.normal(size=(2, 2)))[0]
    values = np.diag([condition, 1.0])
    return (left @ values @ right).astype(np.complex128)


def relative(a: complex, b: complex) -> float:
    return float(abs(a - b) / max(abs(b), 1e-300))


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E2 gauge invariance")
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

    reference_probes = (3.0 + 7.0j, -1.0 + 2.5j, 5.0 - 4.0j)
    reference = space.split(created)
    voltages = {
        b: complex(base.dae.voltages(base.equilibrium.z)[base.dae.network.position(b)])
        for b in core
    }

    generic = tuple(reference_probes)
    generic_reference = [space.split(s) for s in generic]

    def compare(label: str, transformed, condition: float) -> dict[str, object]:
        moved = transformed.split(created)
        eig_reference = np.sort_complex(reference["q_eigenvalues"])
        eig_moved = np.sort_complex(moved["q_eigenvalues"])
        eig_scale = max(float(np.abs(eig_reference).max()), 1e-300)
        away = [transformed.split(s) for s in generic]
        return {
            "basis": label,
            "condition": condition,
            # determinants, measured where they are order one
            "err_full_generic": max(
                relative(m["full"], r["full"])
                for m, r in zip(away, generic_reference, strict=True)
            ),
            "err_collective_generic": max(
                relative(m["collective"], r["collective"])
                for m, r in zip(away, generic_reference, strict=True)
            ),
            # well-scaled at the created mode
            "err_individual": relative(moved["individual"], reference["individual"]),
            "err_q_eigenvalues": float(
                np.abs(eig_moved - eig_reference).max() / eig_scale
            ),
            "closure_distance": float(abs(moved["closest_to_minus_one"] + 1.0)),
            # predicted NOT to be invariant
            "change_sigma_min": abs(
                moved["sigma_min_i_plus_q"] - reference["sigma_min_i_plus_q"]
            )
            / max(reference["sigma_min_i_plus_q"], 1e-300),
            "change_norm_q": abs(moved["norm_q"] - reference["norm_q"])
            / max(reference["norm_q"], 1e-300),
        }

    rows = [compare("identity", space, 1.0)]
    for condition in CONDITIONS:
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            bases = [random_basis(rng, condition) for _ in core]
            rows.append(
                compare(
                    f"random cond={condition:g} seed={seed}",
                    space.with_bases(bases),
                    condition,
                )
            )
    polar = space.with_bases(space.polar_bases(voltages))
    polar_condition = max(float(np.linalg.cond(b)) for b in space.polar_bases(voltages))
    rows.append(compare("polar (|V|, angle)", polar, polar_condition))

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_invariance.csv", index=False)

    invariant_columns = [
        "err_full_generic",
        "err_collective_generic",
        "err_individual",
        "err_q_eigenvalues",
    ]
    # Acceptance is judged on bases whose conditioning does not itself destroy
    # the arithmetic. A similarity at condition 1000 amplifies roundoff by about
    # 1e6, which is a numerical limit of the test and not gauge dependence; the
    # row is still reported.
    usable = table[table.condition <= 100.0]
    worst_invariant = float(usable[invariant_columns].to_numpy().max())
    tolerance = 1e-6
    closure_spread = float(table.closure_distance.max())
    accepted = bool(worst_invariant < tolerance and closure_spread < 1e-5)
    non_invariant_range = {
        "sigma_min_i_plus_q": float(table.change_sigma_min.max()),
        "norm_q": float(table.change_norm_q.max()),
    }

    closure = reference["closest_to_minus_one"]
    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={
            "core": list(core),
            "conditions": list(CONDITIONS),
            "seeds": list(SEEDS),
        },
    )
    manifest.finish(
        "SUCCESS" if accepted else "GAUGE_DEPENDENT",
        created_mode={"real": created.real, "imag": created.imag},
        worst_invariant_error=worst_invariant,
        worst_invariant_error_all_conditions=float(
            table[invariant_columns].to_numpy().max()
        ),
        closure_distance_max=closure_spread,
        tolerance=tolerance,
        non_invariant_relative_change=non_invariant_range,
        reference_split={
            "full": abs(reference["full"]),
            "individual": abs(reference["individual"]),
            "collective": abs(reference["collective"]),
            "spectral_radius_q": reference["spectral_radius_q"],
            "closest_q_eigenvalue_to_minus_one": [closure.real, closure.imag],
            "distance_to_minus_one": abs(closure + 1.0),
        },
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("created mode %+.6f%+.6fj" % (created.real, created.imag))
    print()
    print("reference split in rectangular port coordinates")
    print("  |full|                 %.6e" % abs(reference["full"]))
    print("  |individual|           %.6f" % abs(reference["individual"]))
    print("  |collective|           %.6e" % abs(reference["collective"]))
    print("  rho(Q)                 %.8f" % reference["spectral_radius_q"])
    print(
        "  Q eigenvalue nearest -1: %+.8f%+.8fj  (distance %.3e)"
        % (closure.real, closure.imag, abs(closure + 1.0))
    )
    print()
    print("relative change under an admissible port basis change")
    print(table.to_string(index=False, float_format=lambda v: f"{v:12.3e}"))
    print()
    print(
        "  worst error over the PREDICTED INVARIANTS   %.3e  (tolerance %.0e,"
        " bases at condition <= 100)" % (worst_invariant, tolerance)
    )
    print(
        "  largest distance of a Q eigenvalue to -1     %.3e  over every basis"
        % closure_spread
    )
    print(
        "  relative change of sigma_min(I+Q)  %.3e  <- not invariant, as predicted"
        % non_invariant_range["sigma_min_i_plus_q"]
    )
    print(
        "  relative change of ||Q||           %.3e  <- not invariant, as predicted"
        % non_invariant_range["norm_q"]
    )
    print()
    print("GAUGE INVARIANCE CONFIRMED" if accepted else "RESULT IS GAUGE DEPENDENT")
    print(f"manifest -> {path}")
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
