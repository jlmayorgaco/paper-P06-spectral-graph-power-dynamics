"""E28 - Track-A v2C. Independent replication on the inter-area FAMILY envelope.

Protocol: configs/ias2026/v2c_modal_family_protocol.yaml, seed 20260912, frozen
before execution. Exactly one methodological change from v2B: the stability
endpoint is the envelope of the descendant modal family, not one selected
descendant. Everything else is inherited unchanged.

v1, v2A and v2B stay frozen. H5B stays NON-EXECUTABLE. The v2B alpha_Q curves
stay invalid as physical boundaries and nothing here rehabilitates them.

Three independent verdicts:

    H4C  does the load-by-availability interaction survive on alpha_IA?
    H2C  does the frozen repair still work, on alpha_IA and on the whole band?
    H5C  does the port margin tighten as alpha_IA approaches zero, and does
         omega_port track the worst descendant frequency?

Usage
    python experiments/E28_v2c_modal_family.py [--per-cell 20]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import beta, mannwhitneyu, spearmanr

from _bootstrap import RESULTS
from _v2c_common import (
    BAND_HZ,
    CORE,
    FAMILY_THRESHOLD,
    band_worst,
    base_anchor,
    nominal_reference,
    read_family,
)
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.io.manifest import Manifest
from ibr_cycles.models.ieee39_case import (
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters
from ibr_cycles.models.ieee39_network import load_network
from ibr_cycles.models.port_admittance import build_action_space, closure_profile
from ibr_cycles.uncertainty.regression import least_squares, logistic
from ibr_cycles.uncertainty.sampling import (
    audit_limits,
    sample_operating_point,
    stratified_grid,
)

EXPERIMENT = "E28_v2c_modal_family"
SEED = 20260912
LOAD_STRATA = np.array([0.90, 0.95, 0.98, 1.00, 1.02, 1.05])
AVAILABILITY_STRATA = np.array([0.75, 0.85, 0.92, 0.96, 1.00])
PER_CELL = 20
BAND_SAMPLES = 61
MIN_VOLTAGE, MAX_VOLTAGE = 0.85, 1.15
NEAR_K = 50
FAR_K = 50
NEAR_ABS = 0.05
BOOTSTRAP = 5000
TERMS = ("const", "L", "A", "L*A", "L^2", "A^2")


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

    The corrected upcrossing convention, inherited from the v2B fix: a convex
    quadratic also has a root on the far side of its vertex, and at the edge of
    the data that root is the fit leaving its own support, not a boundary.
    """

    out = []
    for availability in availability_grid:
        column = design_matrix(load_grid, np.full_like(load_grid, availability))
        prediction = column @ coefficients
        upward = np.flatnonzero((prediction[:-1] < 0.0) & (prediction[1:] >= 0.0))
        out.append(float(load_grid[upward[0]]) if upward.size else np.nan)
    return np.array(out)


def zero_curve_band(fit_once, load, availability, response, ag, lg, rng, draws=400):
    n = load.size
    curves = np.full((draws, ag.size), np.nan)
    for draw in range(draws):
        pick = rng.integers(0, n, size=n)
        try:
            curves[draw] = zero_curve(
                fit_once(load[pick], availability[pick], response[pick]), ag, lg
            )
        except np.linalg.LinAlgError:
            continue
    support = np.sum(np.isfinite(curves), axis=0)
    low = np.full(ag.size, np.nan)
    high = np.full(ag.size, np.nan)
    floor = max(draws // 10, 1)
    for column in range(ag.size):
        if support[column] < floor:
            continue
        finite = curves[np.isfinite(curves[:, column]), column]
        low[column] = float(np.percentile(finite, 2.5))
        high[column] = float(np.percentile(finite, 97.5))
    return low, high, support


def clopper_pearson(successes: int, trials: int, level: float = 0.95):
    """Exact binomial interval. Never widens a claim, never narrows one."""

    if trials == 0:
        return float("nan"), 0.0, 1.0
    alpha = 1.0 - level
    low = beta.ppf(alpha / 2.0, successes, trials - successes + 1) if successes else 0.0
    high = (
        beta.ppf(1.0 - alpha / 2.0, successes + 1, trials - successes)
        if successes < trials
        else 1.0
    )
    return successes / trials, float(low), float(high)


def record_family(entry, tag, reading, spectrum):
    family = reading.family
    entry[f"alpha_{tag}"] = family.alpha
    entry[f"freq_{tag}_worst"] = family.frequency_worst_hz
    entry[f"n_family_{tag}"] = family.size
    entry[f"independence_{tag}"] = family.independence
    entry[f"overlap_min_{tag}"] = min(family.overlaps)
    entry[f"overlap_max_{tag}"] = max(family.overlaps)
    entry[f"max_condition_{tag}"] = family.max_condition
    entry[f"band_worst_{tag}"] = band_worst(spectrum)
    # every descendant individually, DIAGNOSTIC only
    for k, (real, frequency) in enumerate(
        zip(family.reals, family.frequencies_hz, strict=True)
    ):
        entry[f"descendant{k}_real_{tag}"] = real
        entry[f"descendant{k}_freq_{tag}"] = frequency
    # the v2B single-branch observable, for contrast only
    entry[f"argmax_alpha_{tag}"] = reading.argmax_real
    entry[f"argmax_freq_{tag}"] = reading.argmax_frequency_hz


def main() -> int:
    parser = argparse.ArgumentParser(description="Track-A v2C")
    parser.add_argument("--per-cell", type=int, default=PER_CELL)
    args = parser.parse_args()
    started = time.time()

    base_network = load_network()
    rng = np.random.default_rng(SEED)
    points = stratified_grid(LOAD_STRATA, AVAILABILITY_STRATA, args.per_cell, rng)
    nominal = nominal_reference()
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
        anchor = base_anchor(spectra["base"], cases["base"], nominal)
        if anchor is None:
            entry.update({"status": "REJECTED", "reason": "TRACKING_FAILURE: no anchor"})
            rows.append(entry)
            continue
        entry["alpha_base_anchor"] = anchor.real
        entry["freq_base_anchor"] = anchor.frequency_hz

        failed = None
        for tag in ("base", "q", "v"):
            reading = read_family(anchor, cases["base"], cases[tag], spectra[tag])
            if reading is None:
                failed = tag
                break
            record_family(entry, tag, reading, spectra[tag])
        if failed is not None:
            entry.update(
                {"status": "REJECTED", "reason": f"TRACKING_FAILURE: empty family ({failed})"}
            )
            rows.append(entry)
            continue

        audit = audit_limits(point.network, cases["base"].dae.power_flow)
        entry["binding_generators"] = ";".join(
            f"{bus}:{kind}" for bus, kind in audit.reactive_binding
        )

        try:
            space = build_action_space(cases["base"], cases["q"], CORE)
            profile = closure_profile(space, band_hz=BAND_HZ, samples=BAND_SAMPLES)
            best = min(profile, key=lambda sample: sample.margin)
            entry["m_4"] = best.margin
            entry["omega_port"] = best.omega
            entry["freq_port_hz"] = best.frequency_hz
            entry["cond_T0_at_omega_port"] = best.condition
            entry["solve_residual_at_omega_port"] = best.solve_residual
            entry["worst_condition_T0"] = max(s.condition for s in profile)
            entry["worst_solve_residual"] = max(s.solve_residual for s in profile)
        except (np.linalg.LinAlgError, ValueError):
            for key in (
                "m_4",
                "omega_port",
                "freq_port_hz",
                "cond_T0_at_omega_port",
                "solve_residual_at_omega_port",
                "worst_condition_T0",
                "worst_solve_residual",
            ):
                entry[key] = float("nan")

        entry["delta_alpha_ia"] = entry["alpha_q"] - entry["alpha_base"]
        entry["degrades"] = bool(entry["delta_alpha_ia"] > 0.0)
        entry["q_unstable"] = bool(entry["alpha_q"] > 0.0)
        entry["v_stable"] = bool(entry["alpha_v"] < 0.0)
        entry["v_band_stable"] = bool(entry["band_worst_v"] < 0.0)
        entry["delta_f_port_hz"] = entry["freq_port_hz"] - entry["freq_q_worst"]
        entry["status"] = "ACCEPTED"
        rows.append(entry)

    table = pd.DataFrame(rows)
    tables = RESULTS / "tables"
    table.to_csv(tables / f"{EXPERIMENT}_samples.csv", index=False)
    accepted = table[table.status == "ACCEPTED"].copy()

    tracking_failures = int(
        table.reason.fillna("").str.startswith("TRACKING_FAILURE").sum()
    )

    # ---- H4C response surface on alpha_IA ----------------------------------
    load = accepted.load.to_numpy(float)
    availability = accepted.availability.to_numpy(float)
    design = design_matrix(load, availability)
    effect = accepted.delta_alpha_ia.to_numpy(float)
    unstable_flag = accepted.q_unstable.to_numpy(float)
    linear = least_squares(design, effect, TERMS)
    binary = logistic(design, unstable_flag, TERMS)
    interval = 1.959963985
    linear_ci = {
        name: [float(c - interval * s), float(c + interval * s)]
        for name, c, s in zip(
            TERMS, linear.coefficients, linear.standard_errors, strict=True
        )
    }
    residual = effect - design @ linear.coefficients
    fitted = design @ linear.coefficients
    residual_diagnostics = {
        "mean": float(residual.mean()),
        "std": float(residual.std(ddof=len(TERMS))),
        "max_abs": float(np.abs(residual).max()),
        "skew": float(((residual - residual.mean()) ** 3).mean() / residual.std() ** 3),
        "kurtosis": float(
            ((residual - residual.mean()) ** 4).mean() / residual.std() ** 4
        ),
        "spearman_abs_residual_vs_fitted": float(
            spearmanr(np.abs(residual), fitted).statistic
        ),
    }

    load_grid = np.linspace(load.min(), load.max(), 400)
    availability_grid = np.linspace(availability.min(), availability.max(), 25)
    effect_curve = zero_curve(linear.coefficients, availability_grid, load_grid)
    stability_curve = zero_curve(binary.coefficients, availability_grid, load_grid)
    effect_low, effect_high, effect_support = zero_curve_band(
        lambda l, a, y: least_squares(design_matrix(l, a), y, TERMS).coefficients,
        load, availability, effect, availability_grid, load_grid, rng,
    )
    stability_low, stability_high, stability_support = zero_curve_band(
        lambda l, a, y: logistic(design_matrix(l, a), y, TERMS).coefficients,
        load, availability, unstable_flag, availability_grid, load_grid, rng,
    )

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

    # worst descendant frequency at the physical boundary, per availability
    boundary_frequency = []
    for column, availability_value in enumerate(availability_grid):
        target = stability_curve[column]
        if not np.isfinite(target):
            boundary_frequency.append(np.nan)
            continue
        distance = (accepted.availability - availability_value).abs() / max(
            availability.std(), 1e-9
        ) + (accepted.load - target).abs() / max(load.std(), 1e-9)
        nearest = accepted.loc[distance.nsmallest(10).index]
        boundary_frequency.append(float(nearest.freq_q_worst.median()))

    curves = pd.DataFrame(
        {
            "availability": availability_grid,
            "load_dispatchable_max": dispatch_edge,
            "load_effect_sign_zero": effect_curve,
            "load_effect_sign_zero_low": effect_low,
            "load_effect_sign_zero_high": effect_high,
            "effect_band_support": effect_support,
            "load_alpha_ia_zero": stability_curve,
            "load_alpha_ia_zero_low": stability_low,
            "load_alpha_ia_zero_high": stability_high,
            "stability_band_support": stability_support,
            "freq_worst_at_boundary_hz": boundary_frequency,
        }
    )
    curves.to_csv(tables / f"{EXPERIMENT}_boundaries.csv", index=False)

    inside_effect = (
        np.isfinite(effect_curve)
        & np.isfinite(dispatch_edge)
        & (effect_curve <= dispatch_edge)
    )
    inside_stability = (
        np.isfinite(stability_curve)
        & np.isfinite(dispatch_edge)
        & (stability_curve <= dispatch_edge)
    )
    interaction_p = linear.p_of("L*A")
    logistic_interaction_p = binary.p_of("L*A")
    same_sign = np.sign(linear.coefficients[linear.index("L*A")]) == np.sign(
        binary.coefficients[binary.index("L*A")]
    )
    h4c = (
        "PASS"
        if interaction_p < 0.05
        else "PARTIAL"
        if (interaction_p < 0.20 and same_sign)
        else "FAIL"
    )

    # ---- H2C repair --------------------------------------------------------
    unstable = accepted[accepted.q_unstable]
    repaired = int(unstable.v_stable.sum())
    h2_point, h2_low, h2_high = clopper_pearson(repaired, len(unstable))
    band_unstable = accepted[accepted.band_worst_q > 0.0]
    band_repaired = int(band_unstable.v_band_stable.sum())
    band_point, band_low, band_high = clopper_pearson(
        band_repaired, len(band_unstable)
    )
    h2c = (
        "PASS"
        if h2_low > 0.90
        else "PARTIAL"
        if (h2_point > 0.90 and h2_low > 0.80)
        else "FAIL"
    )

    # ---- H5C port margin ---------------------------------------------------
    valid = accepted[np.isfinite(accepted.m_4)].copy()
    valid["abs_alpha"] = valid.alpha_q.abs()
    ordered = valid.sort_values("abs_alpha")
    near = ordered.head(NEAR_K)
    far = ordered.tail(FAR_K)
    # A fixed-count stratum cannot be empty, but it CAN overlap its own control
    # group when too few samples survive. That would make the comparison vacuous
    # in exactly the way v2B's empty stratum was, so it is refused, not reported.
    groups_disjoint = bool(
        len(valid) >= NEAR_K + FAR_K
        and not set(near.index) & set(far.index)
    )
    statistic, p_margin = mannwhitneyu(near.m_4, far.m_4, alternative="less")
    correlation = spearmanr(valid.m_4, valid.abs_alpha)
    log_fit = least_squares(
        np.column_stack([np.ones(len(valid)), valid.abs_alpha.to_numpy(float)]),
        np.log(valid.m_4.to_numpy(float)),
        ("const", "abs_alpha_IA"),
    )
    gap = near.delta_f_port_hz.abs()
    median_gap = float(gap.median())
    p95_gap_near = float(np.percentile(gap, 95))
    p95_gap = float(np.percentile(valid.delta_f_port_hz.abs(), 95))
    robustness = valid[valid.abs_alpha <= NEAR_ABS]
    robustness_p = (
        float(mannwhitneyu(robustness.m_4, far.m_4, alternative="less").pvalue)
        if len(robustness) >= 5
        else float("nan")
    )
    primary_1 = bool(groups_disjoint and p_margin < 0.01)
    primary_2 = bool(correlation.statistic > 0 and correlation.pvalue < 0.01)
    secondary = bool(median_gap < 0.05)
    h5c = (
        "PASS"
        if (primary_1 and primary_2 and secondary)
        else "PARTIAL"
        if (primary_1 and primary_2)
        else "FAIL"
    )

    manifest = Manifest(
        experiment=EXPERIMENT,
        seed=SEED,
        config={
            "protocol": "configs/ias2026/v2c_modal_family_protocol.yaml",
            "observable": "alpha_IA, inter-area modal family envelope",
            "band_hz": list(BAND_HZ),
            "family_threshold": FAMILY_THRESHOLD,
            "per_cell": args.per_cell,
            "near_k": NEAR_K,
        },
    )
    manifest.finish(
        "REPORTED",
        drawn=len(table),
        accepted=len(accepted),
        tracking_failures=tracking_failures,
        family_sizes={
            int(k): int(v) for k, v in accepted.n_family_q.value_counts().items()
        },
        H4C={
            "verdict": h4c,
            "linear": linear.summary(),
            "linear_ci95": linear_ci,
            "logistic": binary.summary(),
            "interaction_p": interaction_p,
            "logistic_interaction_p": logistic_interaction_p,
            "interaction_signs_agree": bool(same_sign),
            "r_squared": linear.goodness,
            "pseudo_r_squared": binary.goodness,
            "residual_diagnostics": residual_diagnostics,
            "curves_inside_dispatchable": bool(
                inside_effect.any() and inside_stability.any()
            ),
        },
        H2C={
            "verdict": h2c,
            "alpha_ia": {
                "n_unstable": len(unstable),
                "repaired": repaired,
                "fraction": h2_point,
                "clopper_pearson95": [h2_low, h2_high],
                "min_repaired_alpha_ia": float(unstable.alpha_v.max())
                if len(unstable)
                else None,
            },
            "whole_band": {
                "n_unstable": len(band_unstable),
                "repaired": band_repaired,
                "fraction": band_point,
                "clopper_pearson95": [band_low, band_high],
                "worst_band_alpha_v_over_all": float(accepted.band_worst_v.max()),
            },
        },
        H5C={
            "verdict": h5c,
            "near_k": NEAR_K,
            "median_m4_near": float(near.m_4.median()),
            "median_m4_far": float(far.m_4.median()),
            "near_abs_alpha_max": float(near.abs_alpha.max()),
            "far_abs_alpha_min": float(far.abs_alpha.min()),
            "mannwhitney_p": float(p_margin),
            "spearman_m4_vs_abs_alpha": float(correlation.statistic),
            "spearman_p": float(correlation.pvalue),
            "log_m4_fit": log_fit.summary(),
            "median_abs_delta_f_near_hz": median_gap,
            "p95_abs_delta_f_near_hz": p95_gap_near,
            "p95_abs_delta_f_all_hz": p95_gap,
            "near_far_groups_disjoint": groups_disjoint,
            "robustness_option_A_n": int(len(robustness)),
            "robustness_option_A_p": robustness_p,
        },
        numerics={
            "worst_condition_T0": float(accepted.worst_condition_T0.max()),
            "worst_solve_residual": float(accepted.worst_solve_residual.max()),
        },
        elapsed_s=time.time() - started,
    )
    path = manifest.write(RESULTS / "manifests")

    rejection = 1 - len(accepted) / max(len(table), 1)
    print("v2C modal-family campaign, seed %d" % SEED)
    print("  drawn %d, accepted %d, rejection %.1f%%" % (len(table), len(accepted), 100 * rejection))
    print("  tracking failures: %d" % tracking_failures)
    print("  family size after replacement: %s"
          % dict(sorted(accepted.n_family_q.value_counts().items())))
    print("  worst cond(T0) %.2e ; worst backward residual %.2e"
          % (accepted.worst_condition_T0.max(), accepted.worst_solve_residual.max()))
    print()
    print("H4C  response surface on alpha_IA            verdict %s" % h4c)
    for name in TERMS:
        i = linear.index(name)
        print("       %-6s %+10.4f  [%+9.4f, %+9.4f]  p = %.3e"
              % (name, linear.coefficients[i], linear_ci[name][0],
                 linear_ci[name][1], linear.p_values[i]))
    print("       R^2 = %.3f ; logistic pseudo-R^2 = %.3f" % (linear.goodness, binary.goodness))
    print("       logistic interaction p = %.3e, signs agree: %s"
          % (logistic_interaction_p, same_sign))
    print("       residuals: mean %+.2e, sd %.4f, skew %+.2f, kurtosis %.2f,"
          " |res| vs fitted rho = %+.3f"
          % (residual_diagnostics["mean"], residual_diagnostics["std"],
             residual_diagnostics["skew"], residual_diagnostics["kurtosis"],
             residual_diagnostics["spearman_abs_residual_vs_fitted"]))
    print()
    print("H2C  repair                                  verdict %s" % h2c)
    print("       alpha_IA  : %d of %d, %.4f  Clopper-Pearson [%.4f, %.4f]"
          % (repaired, len(unstable), h2_point, h2_low, h2_high))
    print("       whole band: %d of %d, %.4f  Clopper-Pearson [%.4f, %.4f]"
          % (band_repaired, len(band_unstable), band_point, band_low, band_high))
    print("       worst full-band alpha under voltage control, all samples: %+.4f"
          % accepted.band_worst_v.max())
    if len(unstable):
        print("       least negative repaired alpha_IA: %+.4f" % unstable.alpha_v.max())
    print()
    print("H5C  port margin against alpha_IA            verdict %s" % h5c)
    print("       near k=%d (abs alpha <= %.4f) median m_4 = %.4f"
          % (NEAR_K, near.abs_alpha.max(), near.m_4.median()))
    print("       far  k=%d (abs alpha >= %.4f) median m_4 = %.4f"
          % (FAR_K, far.abs_alpha.min(), far.m_4.median()))
    print("       one-sided Mann-Whitney p = %.3e" % p_margin)
    print("       Spearman m_4 vs abs(alpha_IA) = %+.4f, p = %.3e"
          % (correlation.statistic, correlation.pvalue))
    print("       log(m_4) on abs(alpha_IA): slope %+.4f, p = %.3e (preregistered"
          " convenience form, NOT a claimed power law)"
          % (log_fit.coefficients[1], log_fit.p_of("abs_alpha_IA")))
    print("       median abs(delta f) near = %.4f Hz, 95th pct near = %.4f Hz,"
          " 95th pct over all = %.4f Hz" % (median_gap, p95_gap_near, p95_gap))
    if not groups_disjoint:
        print("       NEAR AND FAR GROUPS OVERLAP: the comparison is vacuous and"
              " primary_1 is refused, not reported as a result")
    print("       robustness, option A stratum abs(alpha_IA) <= %.2f: n = %d, p = %.3e"
          % (NEAR_ABS, len(robustness), robustness_p))
    print()
    print("  the three verdicts are independent; no combined criterion is applied")
    print(f"manifest -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
