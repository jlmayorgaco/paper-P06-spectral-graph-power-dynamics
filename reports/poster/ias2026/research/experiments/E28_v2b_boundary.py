"""E28 - Track-A v2B. Boundary-oriented campaign with a corrected port margin.

Protocol: configs/ias2026/trackA_v2B.yaml, seed 20260911. v1 and v2A stay frozen;
no earlier sample is reused. v2A's global criterion remains unsatisfied.

Three independent verdicts, no combined criterion:

    H4B  second-order response surface for delta_alpha_Q and for the indicator
         of alpha_Q > 0, with the two zero-level curves and their bands
    H2B  conditional repair: alpha_Q > 0 implies alpha_V < 0
    H5B  imaginary-axis port margin m_4 is smaller near the dynamic boundary,
         and omega_port agrees with the tracked frequency there

Three boundaries are kept apart: dispatchability, effect sign, dynamic
stability. A non-dispatchable point is never counted as unstable.

Usage
    python experiments/E28_v2b_boundary.py [--per-cell 20]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

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
from ibr_cycles.models.port_admittance import (
    build_action_space,
    imaginary_axis_closure,
)
from ibr_cycles.uncertainty.regression import least_squares, logistic
from ibr_cycles.uncertainty.sampling import (
    audit_limits,
    sample_operating_point,
    stratified_grid,
)

EXPERIMENT = "E28_v2b_boundary"
CORE = (30, 33, 35, 37)
SEED = 20260911
LOAD_STRATA = np.array([0.90, 0.95, 0.98, 1.00, 1.02, 1.05])
AVAILABILITY_STRATA = np.array([0.75, 0.85, 0.92, 0.96, 1.00])
PER_CELL = 20
BAND_HZ = (0.3, 1.5)
BAND_SAMPLES = 61
MAC_FLOOR = 0.80
MIN_VOLTAGE, MAX_VOLTAGE = 0.85, 1.15
NEAR_BOUNDARY = 0.05
FAR_FROM_BOUNDARY = 0.20
BOOTSTRAP = 5000


def machine_labels(labels):
    return {n for n in labels if n.startswith(("delta_sg", "omega_sg"))}


def pick_anchor(spectrum, case, nominal):
    vector, names = nominal
    common = sorted(set(names) & machine_labels(case.system.labels))
    if len(common) < 4:
        return None
    left = np.array([vector[n] for n in common])
    index = {n: i for i, n in enumerate(case.system.labels)}
    select = [index[n] for n in common]
    candidates = [
        m
        for m in spectrum.modes
        if abs(m.value) > 1e-3 and BAND_HZ[0] <= m.frequency_hz <= BAND_HZ[1]
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


def band_envelope(anchor, ref_labels, ref_index, case, spectrum, floor=MAC_FLOOR):
    """Worst damping among every band mode of MAC-comparable shape.

    DIAGNOSTIC, not a hypothesis variable. The preregistered observable is the
    single argmax-MAC branch; this records whether that branch is alone in the
    band or shares it with another of indistinguishable shape.
    """

    common = sorted(ref_labels & machine_labels(case.system.labels))
    if len(common) < 4:
        return None
    left = anchor.right[[ref_index[n] for n in common]]
    index = {n: i for i, n in enumerate(case.system.labels)}
    select = [index[n] for n in common]
    kept = [
        (m.real, m.frequency_hz, modal_assurance(left, m.right[select]))
        for m in spectrum.modes
        if abs(m.value) > 1e-3
        and m.value.imag >= 0.0
        and BAND_HZ[0] <= m.frequency_hz <= BAND_HZ[1]
    ]
    kept = [k for k in kept if k[2] >= floor]
    if not kept:
        return None
    worst = max(kept, key=lambda k: k[0])
    return worst[0], worst[1], len(kept)


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


def design_matrix(load, availability):
    return np.column_stack(
        [
            np.ones_like(load),
            load,
            availability,
            load * availability,
            load**2,
            availability**2,
        ]
    )


def zero_curve(coefficients, availability_grid, load_grid):
    """Load at which the fitted surface crosses zero UPWARD, per availability.

    The same routine serves both models: the logistic surface crosses one half
    exactly where its linear predictor crosses zero.

    The upcrossing is the one that is asked for, from improvement to degradation
    and from stable to unstable as load rises. A convex quadratic also has a
    second root on the far side of its vertex, and at high availability that root
    lands on the first grid point, which is the fit reaching the edge of its own
    data rather than a located boundary.
    """

    out = []
    for availability in availability_grid:
        column = design_matrix(load_grid, np.full_like(load_grid, availability))
        prediction = column @ coefficients
        upward = np.flatnonzero((prediction[:-1] < 0.0) & (prediction[1:] >= 0.0))
        out.append(float(load_grid[upward[0]]) if upward.size else np.nan)
    return np.array(out)


def zero_curve_band(
    fit_once, load, availability, response, availability_grid, load_grid, rng, draws=400
):
    """Percentile bootstrap band for a zero-level curve.

    Resamples operating points with replacement, refits, and re-locates the
    crossing. A draw whose surface has no crossing at a given availability
    contributes nothing there, so the band is reported with the count of draws
    that actually produced a curve.
    """

    n = load.size
    curves = np.full((draws, availability_grid.size), np.nan)
    for draw in range(draws):
        pick = rng.integers(0, n, size=n)
        try:
            coefficients = fit_once(load[pick], availability[pick], response[pick])
        except np.linalg.LinAlgError:
            continue
        curves[draw] = zero_curve(coefficients, availability_grid, load_grid)
    support = np.sum(np.isfinite(curves), axis=0)
    low = np.full(availability_grid.size, np.nan)
    high = np.full(availability_grid.size, np.nan)
    floor = max(draws // 10, 1)
    for column in range(availability_grid.size):
        if support[column] < floor:
            continue
        finite = curves[np.isfinite(curves[:, column]), column]
        low[column] = float(np.percentile(finite, 2.5))
        high[column] = float(np.percentile(finite, 97.5))
    return low, high, support


def main() -> int:
    parser = argparse.ArgumentParser(description="Track-A v2B")
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
        headroom = sum(
            base_network.pv[b]["pmax"] - point.network.pv[b]["p"]
            for b in base_network.pv
        )
        entry = {
            "sample": index,
            "load": load,
            "availability": availability,
            "reactive_load": reactive,
            "aggregate_headroom_pu": float(headroom),
            "dispatchable": bool(point.dispatchable),
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
        for name in ("base", "q", "v"):
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

        # Diagnostics. They consume no randomness and feed no verdict.
        for name in ("base", "q", "v"):
            envelope = band_envelope(
                anchor, ref_labels, ref_index, cases[name], spectra[name]
            )
            if envelope is None:
                entry[f"alpha_{name}_worst_band"] = float("nan")
                entry[f"freq_{name}_worst_band"] = float("nan")
                entry[f"n_band_{name}"] = 0
                continue
            entry[f"alpha_{name}_worst_band"] = envelope[0]
            entry[f"freq_{name}_worst_band"] = envelope[1]
            entry[f"n_band_{name}"] = envelope[2]

        audit = audit_limits(point.network, cases["base"].dae.power_flow)
        entry["binding_generators"] = ";".join(
            f"{bus}:{kind}" for bus, kind in audit.reactive_binding
        )
        try:
            space = build_action_space(cases["base"], cases["q"], CORE)
            margin = imaginary_axis_closure(
                space, band_hz=BAND_HZ, samples=BAND_SAMPLES
            )
            entry["m_4"] = margin.margin
            entry["omega_port"] = margin.omega_port
            entry["freq_port_hz"] = margin.frequency_port_hz
            entry["worst_condition_T0"] = margin.worst_condition
            descriptive = space.split(complex(0.0, entry["imag_q"]))
            entry["d4_jwQ"] = float(
                abs(descriptive["closest_to_minus_one"] + 1.0)
            )
        except (np.linalg.LinAlgError, ValueError):
            for key in ("m_4", "omega_port", "freq_port_hz", "worst_condition_T0", "d4_jwQ"):
                entry[key] = float("nan")

        entry["delta_alpha_q"] = entry["alpha_q"] - entry["alpha_base"]
        entry["degrades"] = bool(entry["delta_alpha_q"] > 0.0)
        entry["q_unstable"] = bool(entry["alpha_q"] > 0.0)
        entry["v_stable"] = bool(entry["alpha_v"] < 0.0)
        entry["status"] = "ACCEPTED"
        rows.append(entry)

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_samples.csv", index=False)
    accepted = table[table.status == "ACCEPTED"].copy()

    # ---- H4B response surface ---------------------------------------------
    load = accepted.load.to_numpy(float)
    availability = accepted.availability.to_numpy(float)
    design = design_matrix(load, availability)
    names = ("const", "L", "A", "L*A", "L^2", "A^2")
    effect_response = accepted.delta_alpha_q.to_numpy(float)
    unstable_response = accepted.q_unstable.to_numpy(float)
    linear = least_squares(design, effect_response, names)
    binary = logistic(design, unstable_response, names)
    interaction_p = linear.p_of("L*A")
    logistic_interaction_p = binary.p_of("L*A")
    interval = 1.959963985
    linear_ci = {
        name: [
            float(c - interval * s),
            float(c + interval * s),
        ]
        for name, c, s in zip(
            names, linear.coefficients, linear.standard_errors, strict=True
        )
    }

    load_grid = np.linspace(load.min(), load.max(), 400)
    availability_grid = np.linspace(availability.min(), availability.max(), 25)
    effect_curve = zero_curve(linear.coefficients, availability_grid, load_grid)
    stability_curve = zero_curve(binary.coefficients, availability_grid, load_grid)

    effect_low, effect_high, effect_support = zero_curve_band(
        lambda l, a, y: least_squares(design_matrix(l, a), y, names).coefficients,
        load,
        availability,
        effect_response,
        availability_grid,
        load_grid,
        rng,
    )
    stability_low, stability_high, stability_support = zero_curve_band(
        lambda l, a, y: logistic(design_matrix(l, a), y, names).coefficients,
        load,
        availability,
        unstable_response,
        availability_grid,
        load_grid,
        rng,
    )

    # The dispatchable support, so no curve is drawn where the fleet cannot serve.
    # Availability is continuous inside each stratum, so the edge is measured per
    # stratum: the largest load the fleet actually served there.
    served = table[table.dispatchable]
    centres, edges = [], []
    for low_edge, high_edge in zip(
        AVAILABILITY_STRATA[:-1], AVAILABILITY_STRATA[1:], strict=True
    ):
        inside = served[
            (served.availability >= low_edge) & (served.availability <= high_edge)
        ]
        if inside.empty:
            continue
        centres.append(0.5 * (low_edge + high_edge))
        edges.append(float(inside.load.max()))
    dispatch_edge = (
        np.interp(availability_grid, centres, edges)
        if len(centres) >= 2
        else np.full(availability_grid.size, edges[0] if edges else np.nan)
    )

    curves = pd.DataFrame(
        {
            "availability": availability_grid,
            "load_dispatchable_max": dispatch_edge,
            "load_effect_sign_zero": effect_curve,
            "load_effect_sign_zero_low": effect_low,
            "load_effect_sign_zero_high": effect_high,
            "effect_band_support": effect_support,
            "load_stability_zero": stability_curve,
            "load_stability_zero_low": stability_low,
            "load_stability_zero_high": stability_high,
            "stability_band_support": stability_support,
        }
    )
    curves.to_csv(RESULTS / "tables" / f"{EXPERIMENT}_boundaries.csv", index=False)

    # A curve counts only where it lies inside the dispatchable region.
    inside_effect = np.isfinite(effect_curve) & np.isfinite(dispatch_edge) & (
        effect_curve <= dispatch_edge
    )
    inside_stability = np.isfinite(stability_curve) & np.isfinite(dispatch_edge) & (
        stability_curve <= dispatch_edge
    )
    h4b_curves_exist = bool(inside_effect.any() and inside_stability.any())
    h4b = (
        "PASS"
        if (interaction_p < 0.05 and h4b_curves_exist)
        else "PARTIAL"
        if (interaction_p < 0.05 or h4b_curves_exist)
        else "FAIL"
    )

    # ---- H2B conditional repair -------------------------------------------
    unstable = accepted[accepted.q_unstable]
    h2_point, h2_low, h2_high = bootstrap_ci(unstable.v_stable.to_numpy(), rng)
    h2b = (
        "PASS"
        if h2_low > 0.90
        else "PARTIAL"
        if (h2_point > 0.90 and h2_low > 0.80)
        else "FAIL"
    )

    # ---- H5B port margin ---------------------------------------------------
    valid = accepted[np.isfinite(accepted.m_4)]
    near = valid[valid.alpha_q.abs() <= NEAR_BOUNDARY]
    far = valid[valid.alpha_q.abs() >= FAR_FROM_BOUNDARY]
    if len(near) >= 5 and len(far) >= 5:
        statistic, p_margin = mannwhitneyu(near.m_4, far.m_4, alternative="less")
    else:
        statistic, p_margin = np.nan, np.nan
    frequency_gap = (
        (near.freq_port_hz - near.imag_q.abs() / (2 * np.pi)).abs().median()
        if len(near)
        else np.nan
    )
    primary = bool(np.isfinite(p_margin) and p_margin < 0.01)
    secondary = bool(np.isfinite(frequency_gap) and frequency_gap < 0.05)
    h5b = "PASS" if (primary and secondary) else "PARTIAL" if primary else "FAIL"

    # ---- diagnostics, no verdict depends on them ---------------------------
    split = accepted[accepted.n_band_q >= 2]
    disagree = accepted[
        np.sign(accepted.alpha_q) != np.sign(accepted.alpha_q_worst_band)
    ]
    band_unstable = accepted[accepted.alpha_q_worst_band > 0.0]
    band_repaired = band_unstable[band_unstable.alpha_v_worst_band < 0.0]
    diagnostics = {
        "samples_with_a_split_band": int(len(split)),
        "mean_band_modes_above_mac_floor_q": float(accepted.n_band_q.mean()),
        "mean_band_modes_above_mac_floor_base": float(accepted.n_band_base.mean()),
        "tracked_and_band_envelope_disagree_in_sign": int(len(disagree)),
        "band_envelope_unstable": int(len(band_unstable)),
        "band_envelope_unstable_and_repaired": int(len(band_repaired)),
        "min_abs_alpha_q_tracked": float(accepted.alpha_q.abs().min()),
        "min_abs_alpha_q_band_envelope": float(
            accepted.alpha_q_worst_band.abs().min()
        ),
    }

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=SEED,
        config={
            "protocol": "configs/ias2026/trackA_v2B.yaml",
            "load_strata": LOAD_STRATA.tolist(),
            "availability_strata": AVAILABILITY_STRATA.tolist(),
            "per_cell": args.per_cell,
            "band_hz": list(BAND_HZ),
        },
    )
    manifest.finish(
        "REPORTED",
        drawn=len(table),
        accepted=len(accepted),
        rejection_rate=float(1 - len(accepted) / max(len(table), 1)),
        H4B={
            "verdict": h4b,
            "linear": linear.summary(),
            "linear_ci95": linear_ci,
            "logistic": binary.summary(),
            "interaction_p": interaction_p,
            "logistic_interaction_p": logistic_interaction_p,
            "r_squared": linear.goodness,
            "pseudo_r_squared": binary.goodness,
            "logistic_converged": binary.converged,
            "curves_exist_inside_dispatchable": h4b_curves_exist,
            "effect_curve_points_inside": int(inside_effect.sum()),
            "stability_curve_points_inside": int(inside_stability.sum()),
        },
        H2B={
            "verdict": h2b,
            "fraction": h2_point,
            "ci": [h2_low, h2_high],
            "n_unstable": len(unstable),
        },
        H5B={
            "verdict": h5b,
            "n_near": len(near),
            "n_far": len(far),
            "median_m4_near": float(near.m_4.median()) if len(near) else None,
            "median_m4_far": float(far.m_4.median()) if len(far) else None,
            "mannwhitney_p": float(p_margin),
            "median_frequency_gap_hz": float(frequency_gap),
        },
        diagnostics=diagnostics,
        worst_condition_T0=float(accepted.worst_condition_T0.max()),
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    print("v2B boundary campaign, seed %d" % SEED)
    print("  drawn %d, accepted %d, rejection %.1f%% (dispatchability is a boundary,"
          " not an instability)"
          % (len(table), len(accepted), 100 * (1 - len(accepted) / max(len(table), 1))))
    print("  worst cond(T0(jw)) over the band: %.2e" % accepted.worst_condition_T0.max())
    print()
    print("H4B  response surface for delta_alpha_Q      verdict %s" % h4b)
    for name in names:
        i = linear.index(name)
        print("       %-6s %+10.4f  [%+9.4f, %+9.4f]  p = %.3e"
              % (name, linear.coefficients[i], linear_ci[name][0],
                 linear_ci[name][1], linear.p_values[i]))
    print("       R^2 = %.3f ; logistic pseudo-R^2 = %.3f ; logistic interaction"
          " p = %.3e" % (linear.goodness, binary.goodness, logistic_interaction_p))
    print("       zero-level curves inside the dispatchable region: %s"
          " (effect %d/%d, stability %d/%d availability points)"
          % (h4b_curves_exist, inside_effect.sum(), availability_grid.size,
             inside_stability.sum(), availability_grid.size))
    print()
    print("H2B  conditional repair                      verdict %s" % h2b)
    print("       alpha_V < 0 in %.3f [%.3f, %.3f] over n=%d unstable samples"
          % (h2_point, h2_low, h2_high, len(unstable)))
    print()
    print("H5B  imaginary-axis port margin              verdict %s" % h5b)
    print("       near the boundary  n=%d  median m_4 = %s"
          % (len(near), f"{near.m_4.median():.4f}" if len(near) else "-"))
    print("       far from it        n=%d  median m_4 = %s"
          % (len(far), f"{far.m_4.median():.4f}" if len(far) else "-"))
    print("       one-sided Mann-Whitney p = %.3e" % p_margin)
    print("       median |f_port - f_mode| near the boundary = %.4f Hz" % frequency_gap)
    print()
    print("  the three verdicts are independent; no combined criterion is applied")
    print()
    print("DIAGNOSTICS  (no verdict depends on these)")
    print("       band modes of MAC-comparable shape: %.2f in the base case,"
          " %.2f after replacement"
          % (diagnostics["mean_band_modes_above_mac_floor_base"],
             diagnostics["mean_band_modes_above_mac_floor_q"]))
    print("       samples whose band carries two such branches: %d of %d"
          % (len(split), len(accepted)))
    print("       tracked branch and band envelope disagree in sign: %d samples"
          % len(disagree))
    print("       smallest |alpha_Q|: tracked %.4f, band envelope %.4f"
          % (diagnostics["min_abs_alpha_q_tracked"],
             diagnostics["min_abs_alpha_q_band_envelope"]))
    print("       band envelope unstable in %d samples, of which voltage control"
          " leaves the whole band stable in %d"
          % (len(band_unstable), len(band_repaired)))
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
