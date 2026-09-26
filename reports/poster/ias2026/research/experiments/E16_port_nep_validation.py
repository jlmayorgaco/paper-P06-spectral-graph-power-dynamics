"""E16 - PHASE E2. Port-space representation of the Track-A flagship.

Question
    Can the order-4 loss-of-synchronous-support failure be represented and
    localized in a low-dimensional replacement-port space, or does it only exist
    in the full state-space model?

Method
    Build the bus-space operator T(s) of the base case and of the flagship, check
    that it reproduces the full DAE spectrum, measure the rank of the replacement
    update, and evaluate the determinant lemma in the resulting action space.

Accuracy note
    The device port models are built by central differences, so the lemma is
    verified only to about 1e-6 relative. That is the finite-difference floor of
    this implementation, not a property of the identity; it is reported next to
    every residual rather than being described as machine precision.

Status of the result
    NUMERICAL OBSERVATION at one operating point, rho = 1.

Usage
    python experiments/E16_port_nep_validation.py [--core 30 33 35 37]
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
from ibr_cycles.models.port_admittance import build_port_operator

EXPERIMENT = "E16_port_nep_validation"
PROBES = (3.0 + 7.0j, -1.0 + 2.5j, 0.2 + 3.6j, 5.0 - 4.0j, -0.5 + 1.0j)


def port_visibility(operator, spectrum, generic=3.0 + 7.0j) -> pd.DataFrame:
    reference = operator.smallest_singular_value(generic)
    rows = []
    for mode in spectrum.modes:
        if abs(mode.value) <= 1e-3:
            continue
        sigma = operator.smallest_singular_value(mode.value)
        rows.append(
            {
                "real": mode.real,
                "imag": mode.imag,
                "frequency_hz": mode.frequency_hz,
                "sigma_min": sigma,
                "visible": sigma < 1e-6 * max(reference, 1.0),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="PHASE E2 port NEP validation")
    parser.add_argument("--core", type=int, nargs="+", default=[30, 33, 35, 37])
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    core = tuple(args.core)
    started = time.time()

    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in core}))
    t0 = build_port_operator(base.dae, base.equilibrium.x, base.equilibrium.z)
    ts = build_port_operator(
        flagship.dae, flagship.equilibrium.x, flagship.equilibrium.z
    )

    base_visibility = port_visibility(t0, eigen_analysis(base.system.A))
    flagship_spectrum = eigen_analysis(flagship.system.A)
    flagship_visibility = port_visibility(ts, flagship_spectrum)
    created = max(
        (m for m in flagship_spectrum.modes if abs(m.value) > 1e-3),
        key=lambda m: m.real,
    ).value

    position = {b: base.dae.network.position(b) for b in core}
    size = 2 * len(core)
    dimension = t0.dimension
    selector = np.zeros((dimension, size), dtype=np.complex128)
    for k, bus in enumerate(core):
        selector[2 * position[bus], 2 * k] = 1.0
        selector[2 * position[bus] + 1, 2 * k + 1] = 1.0

    def update(s: complex) -> np.ndarray:
        """T_S(s) - T_0(s) restricted to the replacement ports."""

        block = np.zeros((size, size), dtype=np.complex128)
        for k, bus in enumerate(core):
            block[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = t0.bus_admittance(
                s, bus
            ) - ts.bus_admittance(s, bus)
        return block

    def green(s: complex) -> np.ndarray:
        return selector.T @ np.linalg.solve(t0.evaluate(s), selector)

    def exact_ratio(s: complex) -> complex:
        sign_s, log_s = np.linalg.slogdet(ts.evaluate(s))
        sign_0, log_0 = np.linalg.slogdet(t0.evaluate(s))
        return complex((sign_s / sign_0) * np.exp(log_s - log_0))

    # structure of the update, measured against the finite-difference floor
    difference = ts.evaluate(PROBES[0]) - t0.evaluate(PROBES[0])
    off_core = np.ones((dimension, dimension), dtype=bool)
    for p in position.values():
        off_core[2 * p : 2 * p + 2, 2 * p : 2 * p + 2] = False
    fd_floor = float(np.abs(difference[off_core]).max())
    scale = float(np.abs(difference).max())
    rank = int(np.linalg.matrix_rank(difference, tol=5.0 * fd_floor))

    residuals = []
    for s in PROBES:
        lemma = complex(
            np.linalg.det(np.eye(size, dtype=np.complex128) + update(s) @ green(s))
        )
        exact = exact_ratio(s)
        residuals.append(
            {
                "s": str(s),
                "lemma": lemma,
                "exact": exact,
                "relative_error": abs(lemma - exact) / max(abs(exact), 1e-300),
            }
        )
    residual_table = pd.DataFrame(residuals)

    identity = np.eye(size, dtype=np.complex128)
    m = update(created) @ green(created)
    total = identity + m
    self_block = np.zeros_like(total)
    individual = 1.0 + 0.0j
    for k in range(len(core)):
        block = total[2 * k : 2 * k + 2, 2 * k : 2 * k + 2]
        self_block[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = block
        individual *= complex(np.linalg.det(block))
    q = np.linalg.solve(self_block, total) - identity

    split = {
        "abs_full": float(abs(np.linalg.det(total))),
        "sigma_min_full": float(np.linalg.svd(total, compute_uv=False)[-1]),
        "abs_individual": float(abs(individual)),
        "abs_collective": float(abs(np.linalg.det(identity + q))),
        "spectral_radius_q": float(np.max(np.abs(np.linalg.eigvals(q)))),
        "sigma_min_i_plus_q": float(np.linalg.svd(identity + q, compute_uv=False)[-1]),
    }

    tables = RESULTS / "tables"
    base_visibility.to_csv(tables / f"{EXPERIMENT}_base_modes.csv", index=False)
    flagship_visibility.to_csv(tables / f"{EXPERIMENT}_flagship_modes.csv", index=False)
    residual_table.to_csv(tables / f"{EXPERIMENT}_lemma.csv", index=False)

    accepted = bool(
        base_visibility.visible.all()
        and flagship_visibility.visible.all()
        and rank == size
        and residual_table.relative_error.max() < 1e-4
        and split["abs_collective"] < 1e-5
        and split["abs_individual"] > 1e-3
    )

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=args.seed,
        config={"core": list(core), "probes": [str(s) for s in PROBES]},
    )
    manifest.finish(
        "SUCCESS" if accepted else "PORT_REPRESENTATION_FAILED",
        port_dimension=dimension,
        action_dimension=size,
        state_space_dimension=flagship.n_states,
        base_modes_visible=int(base_visibility.visible.sum()),
        base_modes_total=len(base_visibility),
        flagship_modes_visible=int(flagship_visibility.visible.sum()),
        flagship_modes_total=len(flagship_visibility),
        update_rank=rank,
        finite_difference_floor_relative=fd_floor / scale,
        lemma_max_relative_error=float(residual_table.relative_error.max()),
        created_mode={"real": created.real, "imag": created.imag},
        split=split,
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print(
        "port operator %dx%d for a %d-state flagship"
        % (dimension, dimension, flagship.n_states)
    )
    print(
        "  base modes reproduced      %d of %d"
        % (base_visibility.visible.sum(), len(base_visibility))
    )
    print(
        "  flagship modes reproduced  %d of %d"
        % (flagship_visibility.visible.sum(), len(flagship_visibility))
    )
    print()
    print("replacement update structure")
    print("  finite-difference floor    %.2e relative" % (fd_floor / scale))
    print(
        "  numerical rank             %d  (2 per replaced bus, maximum %d)"
        % (rank, size)
    )
    print()
    print("determinant lemma in the %dx%d action space" % (size, size))
    print(residual_table.to_string(index=False))
    print()
    print("at the created mode %+.6f%+.6fj" % (created.real, created.imag))
    for key, value in split.items():
        print("  %-22s %.6g" % (key, value))
    print()
    print("PORT REPRESENTATION ACCEPTED" if accepted else "PORT REPRESENTATION FAILED")
    print(f"manifest -> {path}")
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
