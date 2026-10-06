"""Deterministic post-processing for frozen IAS26-060B outputs.

This script consumes the immutable Julia case table and frozen manifest. It
performs no power-flow, DAE, or eigenvalue solves and never edits input files.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TAU = 1e-8
Z95 = 1.959963984540054
ROOT_REPO = Path(__file__).resolve().parents[6]
CONTRACT = ROOT_REPO / "reports/poster/ias2026/research/bnd_h4_mechanism"
F1_NOMINAL = CONTRACT / "results/20260926T161116Z_d0fecb32_ias26_020_f1_audit_v1/tables/F1_FINAL_AUDIT.csv"
F050_NOMINAL = CONTRACT / "results/20260926T174848Z_d0fecb32_ias26_050_nominal_v1/tables/IAS26-050_NOMINAL_COMPARISON.csv"


def dump_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return math.nan, math.nan
    p = k / n
    den = 1.0 + Z95 * Z95 / n
    center = (p + Z95 * Z95 / (2 * n)) / den
    half = Z95 * math.sqrt(p * (1 - p) / n + Z95 * Z95 / (4 * n * n)) / den
    return center - half, center + half


def classify_alpha(value: float) -> str:
    if not np.isfinite(value):
        return "INVALID"
    if value < -TAU:
        return "STABLE"
    if value > TAU:
        return "UNSTABLE"
    return "INDETERMINATE"


def strict_subsets(mask: int) -> list[int]:
    return [m for m in range(16) if m != mask and (m & mask) == m]


def minimal_unstable_masks(alpha_by_mask: dict[int, float]) -> list[int]:
    out = []
    for mask, alpha in alpha_by_mask.items():
        if mask.bit_count() < 2 or classify_alpha(alpha) != "UNSTABLE":
            continue
        vals = [alpha_by_mask.get(sub, math.nan) for sub in strict_subsets(mask)]
        if vals and all(classify_alpha(v) == "STABLE" for v in vals):
            out.append(mask)
    return sorted(out, key=lambda m: (m.bit_count(), m))


def mask_name(mask: int) -> str:
    buses = [str(bus) for i, bus in enumerate((30, 33, 35, 37)) if mask & (1 << i)]
    return "+".join(buses) if buses else "BASE"


def nominal_point() -> tuple[float, float]:
    f1 = pd.read_csv(F1_NOMINAL)
    alpha = dict(zip(f1.portfolio_mask_CORE_order_30_33_35_37.astype(int), f1["alpha_perp_s-1"], strict=True))
    x = max(alpha[m] for m in range(16) if m != 15)
    return float(x), float(alpha[15])


def event_row(name: str, states: pd.Series, valid_ids: pd.Series, n_total: int) -> dict:
    known = states.isin(["TRUE", "FALSE"])
    k = int((states == "TRUE").sum())
    n = int(known.sum())
    lo, hi = wilson(k, n)
    unknown_all = n_total - n
    return {
        "event": name,
        "estimate": k / n if n else math.nan,
        "numerator": k,
        "denominator": n,
        "wilson95_low": lo,
        "wilson95_high": hi,
        "invalid_count": int((~valid_ids).sum()),
        "indeterminate_count": int((states == "INDETERMINATE").sum()),
        "all_id_conservative_low": k / n_total,
        "all_id_conservative_high": (k + unknown_all) / n_total,
        "ensemble_interpretation": "descriptive synthetic ensemble; not a real-world network probability",
    }


def stats_with_bootstrap(values: pd.Series, rng: np.random.Generator, reps: int) -> dict:
    a = values.to_numpy(dtype=float)
    a = a[np.isfinite(a)]
    if not a.size:
        return {"N": 0}
    quantiles = {q: float(np.quantile(a, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}
    boot = {q: np.empty(reps) for q in (0.05, 0.5, 0.95)}
    for start in range(0, reps, 250):
        stop = min(reps, start + 250)
        indices = rng.integers(0, a.size, size=(stop - start, a.size))
        sample = a[indices]
        for q in boot:
            boot[q][start:stop] = np.quantile(sample, q, axis=1)
    return {
        "N": int(a.size),
        "mean": float(np.mean(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "median": quantiles[0.5],
        "q05": quantiles[0.05],
        "q25": quantiles[0.25],
        "q75": quantiles[0.75],
        "q95": quantiles[0.95],
        "min": float(np.min(a)),
        "max": float(np.max(a)),
        "bootstrap_percentile_95": {
            str(q): [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))]
            for q, v in boot.items()
        },
        "bootstrap_replicates": reps,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_root.resolve()
    config_path = run / "config/IAS26-060_MC_OPERATING_V1.json"
    manifest_path = run / "inputs/IAS26-060_SCENARIOS_V1.csv"
    case_path = run / "derived/MC_CASES.csv"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    manifest = pd.read_csv(manifest_path)
    cases = pd.read_csv(case_path)
    if len(manifest) != 1000 or manifest.scenario_id.nunique() != 1000:
        raise RuntimeError("frozen manifest must have exactly 1,000 unique IDs")
    if len(cases) != 17000 or cases.scenario_id.nunique() != 1000:
        raise RuntimeError(f"full MC integrity failure: {len(cases)} rows")
    if cases.groupby("scenario_id").size().ne(17).any():
        raise RuntimeError("each scenario must have exactly 17 retained evaluations")
    treatment_counts = cases.groupby(["scenario_id", "treatment"]).size().unstack(fill_value=0)
    original_col = treatment_counts.get("ORIGINAL", pd.Series(0, index=treatment_counts.index))
    intervention_cols = treatment_counts.drop(columns=["ORIGINAL"], errors="ignore").sum(axis=1)
    if original_col.ne(16).any():
        raise RuntimeError("each scenario must contain all 16 original portfolios")
    if intervention_cols.ne(1).any():
        raise RuntimeError("each scenario must contain exactly one paired intervention")

    derived = run / "derived"
    tables = run / "tables"
    figures = run / "figures"
    claims = run / "claims"
    for p in (derived, tables, figures, claims):
        p.mkdir(parents=True, exist_ok=True)

    mask_col = "portfolio_mask"
    original = cases[(cases.treatment == "ORIGINAL")].copy()
    intervention = cases[(cases.treatment != "ORIGINAL")].copy()
    summaries = []
    for sid, group in cases.groupby("scenario_id", sort=False):
        orig = original[original.scenario_id == sid]
        ret = intervention[intervention.scenario_id == sid]
        all_complete = len(orig) == 16 and len(ret) == 1 and (cases.loc[cases.scenario_id == sid, "status"] == "PASS").all()
        bymask = {int(r[mask_col]): float(r.alpha_perp) for _, r in orig.iterrows() if pd.notna(r.alpha_perp)}
        finite_full = len(bymask) == 16 and len(ret) == 1 and pd.notna(ret.iloc[0].alpha_perp)
        if len(orig[orig[mask_col] == 15]) != 1:
            raise RuntimeError(f"missing H4 original case for {sid}")
        h4 = orig[orig[mask_col] == 15].iloc[0]
        y = float(h4.alpha_perp) if pd.notna(h4.alpha_perp) else math.nan
        strict = [bymask.get(m, math.nan) for m in range(15)]
        x = float(max(strict)) if all(np.isfinite(strict)) else math.nan
        mval = float(min(y, -x)) if np.isfinite(x) and np.isfinite(y) else math.nan
        alpha_ret = float(ret.iloc[0].alpha_perp) if len(ret) == 1 and pd.notna(ret.iloc[0].alpha_perp) else math.nan
        all_determinate = finite_full and all(classify_alpha(v) != "INDETERMINATE" for v in [*bymask.values(), alpha_ret])
        blockers = minimal_unstable_masks(bymask) if finite_full else []
        other = [m for m in blockers if m != 15]
        if not all_complete or not finite_full:
            h4_event = phenomenon = rescue = "INVALID"
            improve = deteriorate = "INVALID"
        elif not all_determinate:
            h4_event = phenomenon = rescue = improve = deteriorate = "INDETERMINATE"
        else:
            h4_event = "TRUE" if x < -TAU and y > TAU else "FALSE"
            phenomenon = "TRUE" if blockers else "FALSE"
            rescue = "TRUE" if y > TAU and alpha_ret < -TAU else "FALSE"
            delta = alpha_ret - y
            improve = "TRUE" if delta < -TAU else ("FALSE" if delta > TAU else "INDETERMINATE")
            deteriorate = "TRUE" if delta > TAU else ("FALSE" if delta < -TAU else "INDETERMINATE")
        delta = alpha_ret - y if np.isfinite(alpha_ret) and np.isfinite(y) else math.nan
        alphas = [bymask.get(m, math.nan) for m in range(16)]
        h4_class = classify_alpha(y)
        summaries.append({
            "scenario_id": sid,
            "scenario_valid": bool(all_complete and finite_full),
            "load_scale": float(h4.load_scale),
            "epsilon_l2": float(np.linalg.norm(list(json.loads(h4.epsilon_vector_json).values()))),
            "saturation_count": int(h4.saturation_count),
            "voltage_min": float(h4.voltage_min) if pd.notna(h4.voltage_min) else math.nan,
            "voltage_max": float(h4.voltage_max) if pd.notna(h4.voltage_max) else math.nan,
            "voltage_limit_margin": min(float(h4.voltage_min) - 0.9, 1.1 - float(h4.voltage_max)) if pd.notna(h4.voltage_min) and pd.notna(h4.voltage_max) else math.nan,
            "x_s_max_proper_alpha": x,
            "y_s_h4_alpha_original": y,
            "m_s_signed_blocker_margin": mval,
            "alpha_perp_original": y,
            "alpha_perp_intervention": alpha_ret,
            "delta_alpha": delta,
            "h4_class": h4_class,
            "h4_minimal_blocker_persists": h4_event,
            "phenomenon_persists": phenomenon,
            "intervention_rescues_h4": rescue,
            "intervention_improves_h4": improve,
            "intervention_deteriorates_h4": deteriorate,
            "minimal_unstable_blocker_masks": ";".join(map(str, blockers)) if blockers else ("NONE" if finite_full else "UNRESOLVED_INVALID"),
            "minimal_unstable_blocker_ids": ";".join(mask_name(m) for m in blockers) if blockers else ("NONE" if finite_full else "UNRESOLVED_INVALID"),
            "other_minimal_blocker": bool(other),
            "critical_mode_tracking": "IDENTITY_UNRESOLVED_NO_CROSS_SCENARIO_MAC_IN_RUNNER",
            "case_count": len(cases[cases.scenario_id == sid]),
        })
    scen = pd.DataFrame(summaries).merge(manifest, on="scenario_id", validate="one_to_one", suffixes=("", "_manifest"))
    scen.to_csv(derived / "SCENARIO_METRICS.csv", index=False)
    n_total = len(scen)
    valid = scen.scenario_valid.astype(bool)

    # Event estimates preserve invalid IDs and deadband cases separately.
    events = []
    event_map = {
        "E1_H4_MINIMAL_BLOCKER_PERSISTS": scen.h4_minimal_blocker_persists,
        "E2_COLLECTIVE_PHENOMENON_PERSISTS": scen.phenomenon_persists,
        "E3_INTERVENTION_IMPROVES": scen.intervention_improves_h4,
        "E4_INTERVENTION_RESCUES_H4": scen.intervention_rescues_h4,
        "E5_INTERVENTION_DETERIORATES": scen.intervention_deteriorates_h4,
    }
    for name, states in event_map.items():
        events.append(event_row(name, states, valid, n_total))
    pd.DataFrame(events).to_csv(tables / "MC_EVENT_RATES.csv", index=False)

    stat_seed = int(config["statistical_analysis"]["scenario_level_bootstrap"]["seed"])
    reps = int(config["statistical_analysis"]["scenario_level_bootstrap"]["replicates"])
    rng = np.random.Generator(np.random.PCG64(stat_seed))
    stat_rows = []
    for metric in ("x_s_max_proper_alpha", "y_s_h4_alpha_original", "m_s_signed_blocker_margin", "delta_alpha"):
        eligible = scen.loc[valid, metric]
        entry = stats_with_bootstrap(eligible, rng, reps)
        entry.update({"metric": metric, "missing_or_invalid_count": n_total - int(eligible.notna().sum()), "deadband_count": int(eligible.abs().le(TAU).sum())})
        stat_rows.append(entry)
    dump_json(derived / "MC_MARGIN_BOOTSTRAP_STATS.json", {"seed": stat_seed, "bit_generator": "NumPy PCG64", "unit": "scenario_id cluster", "metrics": stat_rows})

    # Descriptive strata using fixed deciles of frozen L and epsilon norm.
    scen["load_decile"] = pd.qcut(scen.load_scale, 10, labels=False, duplicates="drop") + 1
    scen["epsilon_decile"] = pd.qcut(scen.epsilon_l2, 10, labels=False, duplicates="drop") + 1
    strat_rows = []
    for variable, col in (("load_level_decile", "load_decile"), ("spatial_stress_decile", "epsilon_decile"), ("generator_saturation_count", "saturation_count")):
        for level, grp in scen.groupby(col, dropna=False, sort=True):
            gvalid = grp[grp.scenario_valid]
            strat_rows.append({
                "stratifier": variable, "level": level, "N_total": len(grp), "N_valid": len(gvalid),
                "H4_persistence_n": int((gvalid.h4_minimal_blocker_persists == "TRUE").sum()),
                "H4_persistence_rate_valid": float((gvalid.h4_minimal_blocker_persists == "TRUE").mean()) if len(gvalid) else math.nan,
                "phenomenon_persistence_n": int((gvalid.phenomenon_persists == "TRUE").sum()),
                "intervention_rescue_n": int((gvalid.intervention_rescues_h4 == "TRUE").sum()),
                "median_x": float(gvalid.x_s_max_proper_alpha.median()) if len(gvalid) else math.nan,
                "median_y": float(gvalid.y_s_h4_alpha_original.median()) if len(gvalid) else math.nan,
                "median_delta_alpha": float(gvalid.delta_alpha.median()) if len(gvalid) else math.nan,
            })
    pd.DataFrame(strat_rows).to_csv(tables / "MC_STRATIFIED_DESCRIPTIVE.csv", index=False)

    # Deterministic counterexample gallery. Each row is an existing scenario/case.
    gallery: list[dict] = []
    def add(label: str, subset: pd.DataFrame, sort_col: str | None = None, ascending: bool = True, note: str = "") -> None:
        pool = subset.copy()
        if sort_col and len(pool):
            pool = pool.sort_values([sort_col, "scenario_id"], ascending=[ascending, True], kind="mergesort")
        if len(pool):
            r = pool.iloc[0]
            gallery.append({"gallery_id": label, "scenario_id": r.scenario_id, "selection": note, "load_scale": r.load_scale, "x_s": r.x_s_max_proper_alpha, "y_s": r.y_s_h4_alpha_original, "m_s": r.m_s_signed_blocker_margin, "alpha_intervention": r.alpha_perp_intervention, "delta_alpha": r.delta_alpha, "blocker_ids": r.minimal_unstable_blocker_ids, "saturation_count": r.saturation_count, "voltage_limit_margin": r.voltage_limit_margin, "status": "SELECTED"})
        else:
            gallery.append({"gallery_id": label, "scenario_id": "NONE_OBSERVED", "selection": note, "status": "NONE_OBSERVED"})
    good = scen[valid]
    blockers = good[good.h4_minimal_blocker_persists == "TRUE"]
    nx, ny = nominal_point()
    gallery.append({"gallery_id": "C1_CANONICAL_H4", "scenario_id": "CANONICAL_NOMINAL", "selection": "Existing frozen nominal F1/F050 evidence; no rerun", "load_scale": 1.0, "x_s": nx, "y_s": ny, "m_s": min(ny, -nx), "alpha_intervention": -0.017384676061556643, "delta_alpha": -0.14439114388986633, "blocker_ids": "30+33+35+37", "saturation_count": 0, "voltage_limit_margin": math.nan, "status": "SELECTED"})
    add("C2_STRONG_BLOCKER", blockers, "m_s_signed_blocker_margin", False, "maximum positive m_s among H4 blockers")
    add("C3_WEAK_BLOCKER_NEAR_BOUNDARY", blockers[blockers.m_s_signed_blocker_margin > 0], "m_s_signed_blocker_margin", True, "minimum positive m_s among H4 blockers")
    add("C4_WITNESS_CHANGED", good[(good.h4_minimal_blocker_persists == "FALSE") & good.other_minimal_blocker], "scenario_id", True, "H4 is not minimal but another V4 minimal blocker exists")
    add("C5_H4_DISAPPEARS", good[good.h4_minimal_blocker_persists == "FALSE"], "y_s_h4_alpha_original", True, "valid scenario without H4 minimal-blocker persistence")
    rescued = good[good.intervention_rescues_h4 == "TRUE"]
    add("C6_SUCCESSFUL_RESCUE", rescued, "alpha_perp_intervention", True, "smallest (most stable) retuned alpha among rescues")
    add("C7_MARGINAL_RESCUE", rescued, "alpha_perp_intervention", False, "largest retuned alpha still below -tau_dec")
    add("C8_REPAIR_FAILS", good[(good.y_s_h4_alpha_original > TAU) & (good.alpha_perp_intervention > TAU)], "alpha_perp_intervention", False, "original unstable and retuned remains unstable; largest residual positive alpha after the fixed intervention")
    add("C9_REPAIR_DETERIORATES", good[good.delta_alpha > TAU], "delta_alpha", False, "maximum positive delta_alpha")
    gallery.append({"gallery_id": "C10_MODE_SWITCH", "scenario_id": "NONE_OBSERVED", "selection": "No cross-scenario MAC alignment was computed; family identity is unresolved, so no mode-switch case is asserted.", "status": "IDENTITY_UNRESOLVED"})
    add("C11_EXTREME_LOW_LOAD", good, "load_scale", True, "minimum L_s")
    add("C12_EXTREME_HIGH_LOAD", good, "load_scale", False, "maximum L_s")
    add("C13_HIGHEST_SATURATION", good, "saturation_count", False, "maximum saturation_count; stable tie by scenario ID")
    add("C14_CLOSEST_VOLTAGE_LIMIT", good, "voltage_limit_margin", True, "minimum distance to 0.9/1.1 pu")
    pd.DataFrame(gallery).to_csv(tables / "CASE_GALLERY.csv", index=False)

    # Every plot's CSV keeps all valid and invalid IDs, so no outlier is hidden.
    f1x, f1y = nominal_point()
    p3csv = scen[["scenario_id", "scenario_valid", "x_s_max_proper_alpha", "y_s_h4_alpha_original", "h4_minimal_blocker_persists", "phenomenon_persists"]].copy()
    p3csv.to_csv(figures / "P3_MC_BLOCKER_QUADRANT.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 5.6), constrained_layout=True)
    ax.fill_between([-0.35, 0], 0, 0.35, color="#e9f2ed", zorder=0)
    mask = p3csv.scenario_valid
    ax.scatter(p3csv.loc[mask, "x_s_max_proper_alpha"], p3csv.loc[mask, "y_s_h4_alpha_original"], s=18, alpha=0.42, linewidths=0, rasterized=True, c="#315c72", label="valid MC scenario")
    ax.scatter([f1x], [f1y], marker="*", s=180, color="#c43b3b", edgecolor="white", linewidth=0.7, zorder=5, label="canonical nominal")
    ax.axvline(0, color="#333333", lw=0.9); ax.axhline(0, color="#333333", lw=0.9)
    ax.set(xlabel=r"$x_s=\max_{R\subsetneq H_4}\alpha_{\perp,s}(R)$ (s$^{-1}$)", ylabel=r"$y_s=\alpha_{\perp,s}(H_4)$ (s$^{-1}$)", title="H4 blocker robustness across the frozen synthetic ensemble")
    ax.legend(loc="best", frameon=False)
    er = pd.DataFrame(events).set_index("event")
    e1, e2 = er.loc["E1_H4_MINIMAL_BLOCKER_PERSISTS"], er.loc["E2_COLLECTIVE_PHENOMENON_PERSISTS"]
    ax.text(0.02, 0.98, f"N valid={int(mask.sum())}/1000\nH4 persistence={e1.estimate:.1%} (Wilson 95% {e1.wilson95_low:.1%}–{e1.wilson95_high:.1%})\nPhenomenon persistence={e2.estimate:.1%} (Wilson 95% {e2.wilson95_low:.1%}–{e2.wilson95_high:.1%})\nInvalid={int((~mask).sum())}; indeterminate={int((scen.h4_minimal_blocker_persists=='INDETERMINATE').sum())}", transform=ax.transAxes, va="top", ha="left", fontsize=8, bbox={"boxstyle":"round,pad=0.4", "fc":"white", "ec":"#bbbbbb", "alpha":0.94})
    for ext in ("pdf", "svg", "png"):
        fig.savefig(figures / f"P3_MC_BLOCKER_QUADRANT.{ext}", dpi=240)
    plt.close(fig)

    p4csv = scen[["scenario_id", "scenario_valid", "alpha_perp_original", "alpha_perp_intervention", "delta_alpha", "intervention_improves_h4", "intervention_rescues_h4", "intervention_deteriorates_h4"]].copy()
    p4csv.to_csv(figures / "P4A_MC_CONTROL_RESCUE.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 5.6), constrained_layout=True)
    ax.scatter(p4csv.loc[mask, "alpha_perp_original"], p4csv.loc[mask, "alpha_perp_intervention"], s=18, alpha=0.42, linewidths=0, rasterized=True, c="#315c72")
    bounds = [float(min(p4csv.loc[mask, "alpha_perp_original"].min(), p4csv.loc[mask, "alpha_perp_intervention"].min())), float(max(p4csv.loc[mask, "alpha_perp_original"].max(), p4csv.loc[mask, "alpha_perp_intervention"].max()))]
    ax.plot(bounds, bounds, ls="--", lw=0.9, color="#555555"); ax.axvline(0, color="#333333", lw=0.9); ax.axhline(0, color="#333333", lw=0.9)
    ax.set(xlabel=r"Original H4 $\alpha_\perp$ (s$^{-1}$)", ylabel=r"Retuned H4 $\alpha_\perp$ (s$^{-1}$)", title="Paired fixed-retuning response; same scenario IDs")
    ev = er
    e3, e4, e5 = (ev.loc[k] for k in ("E3_INTERVENTION_IMPROVES", "E4_INTERVENTION_RESCUES_H4", "E5_INTERVENTION_DETERIORATES"))
    delta_stats = next(x for x in stat_rows if x["metric"] == "delta_alpha")
    ax.text(0.02, 0.98, f"Improvement={e3.estimate:.1%} (95% {e3.wilson95_low:.1%}–{e3.wilson95_high:.1%})\nRescue={e4.estimate:.1%} (95% {e4.wilson95_low:.1%}–{e4.wilson95_high:.1%})\nDeterioration={e5.estimate:.1%} (95% {e5.wilson95_low:.1%}–{e5.wilson95_high:.1%})\nMedian Δα={delta_stats['median']:.4g} s$^{{-1}}$", transform=ax.transAxes, va="top", ha="left", fontsize=8, bbox={"boxstyle":"round,pad=0.4", "fc":"white", "ec":"#bbbbbb", "alpha":0.94})
    for ext in ("pdf", "svg", "png"):
        fig.savefig(figures / f"P4A_MC_CONTROL_RESCUE.{ext}", dpi=240)
    plt.close(fig)

    # S1 is descriptive only: fixed decile means and 10–90% spread.
    valid_scen = scen[valid].copy()
    loadbins = valid_scen.groupby("load_decile", observed=True).agg(load_median=("load_scale", "median"), y_median=("y_s_h4_alpha_original", "median"), y_q10=("y_s_h4_alpha_original", lambda s: s.quantile(.1)), y_q90=("y_s_h4_alpha_original", lambda s: s.quantile(.9)), m_median=("m_s_signed_blocker_margin", "median"), n=("scenario_id", "size")).reset_index()
    loadbins.to_csv(figures / "S1_MC_LOAD_SENSITIVITY.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.8, 4.4), constrained_layout=True)
    ax.errorbar(loadbins.load_median, loadbins.y_median, yerr=[loadbins.y_median-loadbins.y_q10, loadbins.y_q90-loadbins.y_median], marker="o", capsize=2, color="#315c72")
    ax.axhline(0, color="#333333", lw=0.9)
    ax.set(xlabel=r"Load scale $L_s$ (decile median)", ylabel=r"H4 $\alpha_\perp$ (s$^{-1}$; median and 10–90% range)", title="Descriptive H4 sensitivity versus load level")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(figures / f"S1_MC_LOAD_SENSITIVITY.{ext}", dpi=240)
    plt.close(fig)

    blocker_rows = []
    for _, r in good.iterrows():
        for bid in str(r.minimal_unstable_blocker_ids).split(";"):
            blocker_rows.append({"scenario_id": r.scenario_id, "minimal_blocker_id": bid})
    freq = pd.DataFrame(blocker_rows).groupby("minimal_blocker_id").size().rename("scenario_count").reset_index() if blocker_rows else pd.DataFrame(columns=["minimal_blocker_id", "scenario_count"])
    freq.to_csv(tables / "MINIMAL_BLOCKER_IDENTITY_FREQUENCY.csv", index=False)
    if (freq.minimal_blocker_id.dropna() != "30+33+35+37").any() and len(freq):
        fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
        ax.bar(freq.minimal_blocker_id, freq.scenario_count, color="#315c72"); ax.set_ylabel("Scenario count"); ax.set_xlabel("Minimal unstable portfolio within frozen V4"); ax.tick_params(axis="x", rotation=30)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(figures / f"S2_MINIMAL_BLOCKER_IDENTITY.{ext}", dpi=240)
        plt.close(fig)

    # Automatic cross-code holdout selection: five disjoint ranked strata.
    available = good.copy(); selections: list[dict] = []; already: set[str] = set()
    def pick_stratum(name: str, frame: pd.DataFrame, col: str | None, ascending: bool = True, fallback: pd.DataFrame | None = None) -> None:
        frame = frame[~frame.scenario_id.isin(already)].copy()
        if col and len(frame):
            frame = frame.sort_values([col, "scenario_id"], ascending=[ascending, True], kind="mergesort")
        chosen = frame.head(10)
        if len(chosen) < 10 and fallback is not None:
            remainder = fallback[~fallback.scenario_id.isin(already | set(chosen.scenario_id))].copy()
            if col and len(remainder):
                remainder = remainder.sort_values([col, "scenario_id"], ascending=[ascending, True], kind="mergesort")
            chosen = pd.concat([chosen, remainder.head(10 - len(chosen))], ignore_index=True)
        for rank, (_, rr) in enumerate(chosen.iterrows(), 1):
            already.add(rr.scenario_id)
            actual = "SELECTED" if rr.scenario_id in set(frame.scenario_id) else "PREDECLARED_NEAREST_AVAILABLE_FALLBACK"
            selections.append({"stratum": name, "rank_within_stratum": rank, "scenario_id": rr.scenario_id, "x_s": rr.x_s_max_proper_alpha, "y_s": rr.y_s_h4_alpha_original, "m_s": rr.m_s_signed_blocker_margin, "alpha_intervention": rr.alpha_perp_intervention, "selection_status": actual})
        for rank in range(len(chosen) + 1, 11):
            selections.append({"stratum": name, "rank_within_stratum": rank, "scenario_id": "NONE_AVAILABLE", "selection_status": "EMPTY_STRATUM_NO_REPLACEMENT"})
    pick_stratum("strongest_blockers", available[available.h4_minimal_blocker_persists == "TRUE"], "m_s_signed_blocker_margin", False, available)
    pick_stratum("near_boundary_blockers", available[(available.h4_minimal_blocker_persists == "TRUE") & (available.m_s_signed_blocker_margin > 0)], "m_s_signed_blocker_margin", True, available[available.h4_minimal_blocker_persists == "TRUE"])
    pick_stratum("h4_not_persistent", available[available.h4_minimal_blocker_persists == "FALSE"], "y_s_h4_alpha_original", True, available)
    pick_stratum("successful_rescue", available[available.intervention_rescues_h4 == "TRUE"], "alpha_perp_intervention", True, available[available.h4_minimal_blocker_persists == "TRUE"])
    repair = available[(available.y_s_h4_alpha_original > TAU) & (available.alpha_perp_intervention > -TAU)]
    if not len(repair):
        repair = available.assign(_distance=(available.alpha_perp_intervention + TAU).abs()).sort_values(["_distance", "scenario_id"])
    pick_stratum("repair_failure_or_nearest_fallback", repair, "alpha_perp_intervention", False, available)
    pd.DataFrame(selections).to_csv(tables / "CROSSCODE_HOLDOUT_SELECTION.csv", index=False)

    # Initial falsification ledger from MC evidence; TDS-dependent rows are finalized later.
    h4_true = int((good.h4_minimal_blocker_persists == "TRUE").sum())
    h4_false = int((good.h4_minimal_blocker_persists == "FALSE").sum())
    e2_true = int((good.phenomenon_persists == "TRUE").sum())
    e3_true = int((good.intervention_improves_h4 == "TRUE").sum())
    e4_false = int((good.intervention_rescues_h4 == "FALSE").sum())
    e5_true = int((good.intervention_deteriorates_h4 == "TRUE").sum())
    valid_cases = cases[(cases.status == "PASS") & cases.physical_feasible.astype(bool)]
    outside_band_unstable = valid_cases[
        (valid_cases.alpha_perp > TAU)
        & (valid_cases.alpha_omega < -TAU)
        & ((valid_cases.critical_frequency_hz < 0.3) | (valid_cases.critical_frequency_hz > 1.5))
    ]
    n9_scenarios = int(outside_band_unstable.scenario_id.nunique())
    invalid_share = float((~valid).mean())
    falsification = [
        ["N1", "Do proper subsets become unstable before H4?", "SUPPORTED" if any(good.x_s_max_proper_alpha > TAU) else "REFUTED", f"{int((good.x_s_max_proper_alpha > TAU).sum())} valid scenarios have an unstable strict H4 subset", "tables/CASE_GALLERY.csv; derived/SCENARIO_METRICS.csv"],
        ["N2", "Is H4 stable in some operating scenarios?", "SUPPORTED" if any(good.y_s_h4_alpha_original < -TAU) else "REFUTED", f"{int((good.y_s_h4_alpha_original < -TAU).sum())} valid scenarios have stable H4", "derived/SCENARIO_METRICS.csv"],
        ["N3", "Does g=0.25 worsen alpha?", "SUPPORTED" if e5_true else "REFUTED", f"{e5_true} material deteriorations among valid paired scenarios", "tables/MC_EVENT_RATES.csv"],
        ["N4", "Does g=0.25 fail to rescue unstable H4?", "SUPPORTED" if e4_false else "REFUTED", f"{e4_false} valid determinate non-rescues; rescue={int((good.intervention_rescues_h4=='TRUE').sum())}", "tables/MC_EVENT_RATES.csv"],
        ["N5", "Do modal families switch?", "NOT_TESTABLE", "The primary runner retains eigenvectors, but no cross-scenario quotient-basis MAC alignment is computed in this postprocessor; no identity claim is made.", "derived/SCENARIO_METRICS.csv"],
        ["N6", "Do results change materially with load level?", "NOT_TESTABLE", "Load-decile results are descriptive; no material-effect threshold was preregistered. No causal or inferential load trend is claimed.", "tables/MC_STRATIFIED_DESCRIPTIVE.csv"],
        ["N7", "Does maximum generator saturation explain extreme outcomes?", "NOT_TESTABLE", "Saturation strata are descriptive and cannot establish explanation/causality without a preregistered attribution contrast.", "tables/MC_STRATIFIED_DESCRIPTIVE.csv"],
        ["N8", "Does nonlinear TDS contradict eigenanalysis?", "NOT_TESTABLE", "Awaiting completed IAS26-080 trajectory analysis.", "raw/tds"],
        ["N9", "Is there a rightmost full-spectrum instability outside 0.3–1.5 Hz while alpha_Omega is stable?", "SUPPORTED" if n9_scenarios else "REFUTED", f"{len(outside_band_unstable)} valid case cells in {n9_scenarios} scenarios satisfy alpha_perp>tau, alpha_Omega<-tau, and critical frequency outside [0.3,1.5] Hz", "derived/MC_CASES.csv"],
        ["N10", "Is a substantial share physically infeasible?", "NOT_TESTABLE", f"Observed {int((~valid).sum())}/{n_total} ({invalid_share:.1%}) physically infeasible IDs; 'substantial' has no preregistered cutoff, so no threshold is invented", "tables/MC_EVENT_RATES.csv; derived/SCENARIO_METRICS.csv"],
    ]
    pd.DataFrame(falsification, columns=["question_id", "question", "status", "finding", "evidence_path"]).to_csv(claims / "FALSIFICATION_LEDGER.csv", index=False)

    full_audit = {
        "run_id": run.name,
        "scenario_ids_in_manifest": int(manifest.scenario_id.nunique()),
        "case_rows": int(len(cases)),
        "scenario_ids_with_17_cases": int((cases.groupby("scenario_id").size() == 17).sum()),
        "status_counts": cases.status.value_counts().to_dict(),
        "valid_scenario_count": int(valid.sum()),
        "invalid_scenario_count": int((~valid).sum()),
        "scenario_id_replacements": 0,
        "outcome_dependent_filtering": False,
        "spectrum_files": len(list((run / "raw/julia/scenarios").glob("IAS26-060-S*.jld2"))),
        "expected_spectrum_files": 1000,
        "max_pf_residual": float(cases.pf_residual.max()),
        "max_dae_f_residual": float(cases.dae_f_residual.max()),
        "max_dae_g_residual": float(cases.dae_g_residual.max()),
        "max_eigenpair_residual": float(cases.eigenpair_residual.max()),
        "max_complement_condition": float(cases.complement_condition.max()),
        "bootstrap_seed": stat_seed,
        "bootstrap_replicates": reps,
        "TDS_model_scope": "nonlinear phasor-domain DAE, not EMT",
    }
    dump_json(derived / "MC_FULL_AUDIT.json", full_audit)
    print(json.dumps(full_audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
