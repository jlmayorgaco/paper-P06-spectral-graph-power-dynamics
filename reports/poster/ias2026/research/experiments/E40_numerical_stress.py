"""E40 - overnight 12. Numerical stress test of the closure computation.

Wider than E30: three perturbation amplitudes, three grid densities, three solve
algorithms, local minimisation, and the variation of the closure EIGENVALUE
itself rather than only its distance to -1.

Extended precision is not available on this platform. numpy's ``longdouble`` has
a 53-bit mantissa here, identical to ``float64``, and nothing may be installed
tonight, so there is no higher-precision arithmetic to fall back on. The
perturbation ensembles and the three-algorithm comparison are the substitute and
the limitation is stated rather than hidden.

Cohorts: 20 generic and 20 near-boundary points from the E33 deterministic
continuation, which is exactly reproducible, plus the flagship, the RC repair and
the 25 % condenser at the nominal point.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from _bootstrap import RESULTS  # noqa: F401
from _overnight import Experiment, pin_blas_threads, run_directory
from _v2c_common import BAND_HZ, CORE, base_anchor, nominal_reference, read_family
from E30_v2c_numerics import margin_from_solution, solve_lu, solve_qr, solve_svd
from E32_tds_validation import rc_converter
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import (
    build_action_space,
    closure_at,
    closure_profile,
    refine_closure,
)
from ibr_cycles.uncertainty.sampling import sample_operating_point

NAME = "E40_numerical_stress"
BAND_SAMPLES = 61
PERTURBATIONS = (1e-14, 1e-13, 1e-12)
ENSEMBLE = 30
LONGDOUBLE_BITS = int(np.finfo(np.longdouble).nmant) + 1


def closure_eigenvalue(space, omega):
    return complex(space.split(complex(0.0, omega))["closest_to_minus_one"])


def closure_from_solution(space, s, solution):
    """Margin AND the interaction eigenvalue nearest -1, from a given solve.

    The perturbation ensemble has to read the eigenvalue off the PERTURBED
    operator. Calling ``closure_eigenvalue(space, omega)`` inside the loop reads
    it off the unperturbed one, which makes the reported variation identically
    zero and says nothing at all.
    """

    u = space.selector()
    m = space.update(s) @ (u.T @ solution)
    size = m.shape[0]
    identity = np.eye(size, dtype=np.complex128)
    total = identity + m
    blocks = np.zeros_like(total)
    for a in range(space.order):
        rows = slice(2 * a, 2 * a + 2)
        blocks[rows, rows] = total[rows, rows]
    q = np.linalg.solve(blocks, total) - identity
    eigenvalues = np.linalg.eigvals(q)
    index = int(np.argmin(np.abs(eigenvalues + 1.0)))
    return float(np.abs(eigenvalues[index] + 1.0)), complex(eigenvalues[index])


def stress(space, label, extra, rng):
    profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
    coarse = min(profile, key=lambda s: s.margin)
    omega = coarse.omega
    s = complex(0.0, omega)
    operator = space.t0.evaluate(s)
    selector = space.selector()

    solutions = {
        "lu": solve_lu(operator, selector),
        "qr": solve_qr(operator, selector),
        "svd": solve_svd(operator, selector),
    }
    margins = {
        name: margin_from_solution(space, s, value)
        for name, value in solutions.items()
    }
    residuals = {
        name: float(np.linalg.norm(operator @ value - selector, 2))
        / (float(np.linalg.norm(operator, 2)) * float(np.linalg.norm(value, 2)))
        for name, value in solutions.items()
    }

    grids = {}
    for samples in (BAND_SAMPLES, 2 * BAND_SAMPLES - 1, 4 * BAND_SAMPLES - 3):
        best = min(
            closure_profile(space, band_hz=BAND_HZ, samples=samples),
            key=lambda x: x.margin,
        )
        grids[samples] = (best.margin, best.frequency_hz)
    half_width = 2.0 * np.pi * (BAND_HZ[1] - BAND_HZ[0]) / (BAND_SAMPLES - 1)
    local = refine_closure(space, omega, half_width=half_width)

    mu = closure_eigenvalue(space, omega)
    scale = float(np.linalg.norm(operator, "fro")) / np.sqrt(operator.size)
    ensembles = {}
    for amplitude in PERTURBATIONS:
        values, eigenvalue_shifts = [], []
        for _ in range(ENSEMBLE):
            noise = amplitude * scale * (
                rng.standard_normal(operator.shape)
                + 1j * rng.standard_normal(operator.shape)
            )
            try:
                perturbed = np.linalg.solve(operator + noise, selector)
            except np.linalg.LinAlgError:
                continue
            margin, perturbed_mu = closure_from_solution(space, s, perturbed)
            values.append(margin)
            eigenvalue_shifts.append(abs(perturbed_mu - mu))
        ensembles[amplitude] = (
            float(np.std(values)) if values else float("nan"),
            float(np.max(np.abs(np.array(values) - margins["lu"])))
            if values
            else float("nan"),
            float(np.max(eigenvalue_shifts)) if eigenvalue_shifts else float("nan"),
        )

    row = {
        "case": label,
        "m4": coarse.margin,
        "freq_port_hz": coarse.frequency_hz,
        "cond_T0": coarse.condition,
        "backward_residual": coarse.solve_residual,
        "mu_real": mu.real,
        "mu_imag": mu.imag,
        "m4_lu": margins["lu"],
        "m4_qr": margins["qr"],
        "m4_svd": margins["svd"],
        "algorithm_spread": float(max(margins.values()) - min(margins.values())),
        "residual_lu": residuals["lu"],
        "residual_qr": residuals["qr"],
        "residual_svd": residuals["svd"],
        "m4_grid_61": grids[BAND_SAMPLES][0],
        "m4_grid_121": grids[2 * BAND_SAMPLES - 1][0],
        "m4_grid_241": grids[4 * BAND_SAMPLES - 3][0],
        "freq_grid_61_hz": grids[BAND_SAMPLES][1],
        "freq_grid_241_hz": grids[4 * BAND_SAMPLES - 3][1],
        "grid_variation": float(
            max(v[0] for v in grids.values()) - min(v[0] for v in grids.values())
        ),
        "freq_variation_hz": float(
            max(v[1] for v in grids.values()) - min(v[1] for v in grids.values())
        ),
        "m4_local_min": local.margin,
        "freq_local_hz": local.frequency_hz,
        "local_reduction": float(coarse.margin - local.margin),
    }
    for amplitude, (std, worst, eigen_shift) in ensembles.items():
        tag = f"{amplitude:.0e}"
        row[f"ensemble_std_{tag}"] = std
        row[f"ensemble_worst_shift_{tag}"] = worst
        row[f"ensemble_worst_eigenvalue_shift_{tag}"] = eigen_shift
    row.update(extra)
    return row


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Is the closure computation numerically trustworthy?",
        config={
            "perturbations": list(PERTURBATIONS),
            "ensemble": ENSEMBLE,
            "grids": [BAND_SAMPLES, 2 * BAND_SAMPLES - 1, 4 * BAND_SAMPLES - 3],
            "algorithms": ["LU", "QR", "SVD"],
            "no_explicit_inverse": True,
            "longdouble_mantissa_bits": LONGDOUBLE_BITS,
            "extended_precision": (
                "unavailable: numpy longdouble is 53-bit on this platform and no "
                "package may be installed"
            ),
        },
        workers=1,
    )
    started = time.time()
    rng = np.random.default_rng(20260918)

    continuation = pd.read_csv(
        run_directory() / "E33_closure_continuation" / "E33_closure_continuation.csv"
    )
    good = continuation[
        (continuation.status == "OK") & continuation.m4_frozen.notna()
    ].copy()
    good["abs_alpha"] = good.alpha_IA.abs()
    near = good.nsmallest(20, "abs_alpha")
    generic = good.sort_values("abs_alpha").iloc[::max(len(good) // 20, 1)][:20]

    network = load_network()
    plan = ReplacementPlan.of({b: 1.0 for b in CORE})
    rows = []

    for cohort, block in (("near boundary", near), ("generic", generic)):
        for _, entry in block.iterrows():
            point = sample_operating_point(
                network,
                active_load=float(entry.load),
                reactive_load=1.0,
                availability=float(entry.availability),
                availability_buses=CORE,
                rng=np.random.default_rng(0),
                dispatch_jitter=0.0,
                load_scatter=0.0,
            )
            try:
                base = solve_case(ReplacementPlan.of({}), network=point.network)
                replaced = solve_case(plan, network=point.network)
                space = build_action_space(base, replaced, CORE)
            except (InfeasibleReplacement, ValueError, np.linalg.LinAlgError):
                continue
            rows.append(
                stress(
                    space,
                    cohort,
                    {
                        "load": float(entry.load),
                        "availability": float(entry.availability),
                        "alpha_IA": float(entry.alpha_IA),
                    },
                    rng,
                )
            )
        print(f"  {cohort}: {len(rows)} cumulative", flush=True)

    base = solve_case(ReplacementPlan.of({}))
    nominal = nominal_reference()
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal)
    for label, case_plan, converter in (
        ("flagship", plan, None),
        ("RC repair", plan, rc_converter()),
        (
            "condenser 25%",
            ReplacementPlan.of(
                {b: 1.0 for b in CORE}, condenser={b: 0.25 for b in CORE}
            ),
            None,
        ),
    ):
        case = solve_case(case_plan, converter=converter)
        reading = read_family(anchor, base, case, eigen_analysis(case.system.A))
        space = build_action_space(base, case, CORE)
        rows.append(
            stress(
                space,
                label,
                {
                    "load": 1.0,
                    "availability": 1.0,
                    "alpha_IA": reading.family.alpha if reading else float("nan"),
                },
                rng,
            )
        )
        print(f"  {label} done", flush=True)

    table = pd.DataFrame(rows)
    experiment.save_table(table, "E40_numerical_stress.csv")

    checks = {
        "cases": int(len(table)),
        "max_backward_residual": float(table.backward_residual.max()),
        "max_algorithm_spread": float(table.algorithm_spread.max()),
        "max_grid_variation": float(table.grid_variation.max()),
        "max_freq_variation_hz": float(table.freq_variation_hz.max()),
        "max_local_reduction": float(table.local_reduction.max()),
        "min_m4": float(table.m4.min()),
        "min_m4_local": float(table.m4_local_min.min()),
        "max_cond_T0": float(table.cond_T0.max()),
        "longdouble_mantissa_bits": LONGDOUBLE_BITS,
    }
    for amplitude in PERTURBATIONS:
        tag = f"{amplitude:.0e}"
        checks[f"max_ensemble_shift_{tag}"] = float(
            table[f"ensemble_worst_shift_{tag}"].max()
        )
        checks[f"max_ensemble_eigenvalue_shift_{tag}"] = float(
            table[f"ensemble_worst_eigenvalue_shift_{tag}"].max()
        )
    worst_numerical = max(
        checks["max_algorithm_spread"],
        checks["max_ensemble_shift_1e-12"],
        float(table.backward_residual.max()),
    )
    checks["headroom_factor"] = float(checks["min_m4"] / max(worst_numerical, 1e-300))
    checks["grid_is_the_dominant_error"] = bool(
        checks["max_grid_variation"] > 100.0 * worst_numerical
    )

    print()
    print(
        table[
            [
                "case",
                "alpha_IA",
                "m4",
                "m4_grid_241",
                "m4_local_min",
                "cond_T0",
                "backward_residual",
                "algorithm_spread",
                "ensemble_worst_shift_1e-12",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:11.3e}")
    )
    print()
    for key, value in checks.items():
        print(f"  {key:34s} {value}")
    experiment.finish("REPORTED", checks=checks, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
