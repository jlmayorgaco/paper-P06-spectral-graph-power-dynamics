"""UC02: the E13/E39 baseline audit and its matched analyses, with true active power.

The frozen "replaced MW" predictor (E13, E39) is the converter RATING in MVA.
This re-audit keeps it (relabelled replaced_sn_mva) and adds the measured
active dispatch replaced_pg_mw, the documented SG limit replaced_pmax_mw and the
reactive dispatch replaced_q_mvar (UC01). Nothing is tuned; no frozen file is
rewritten.

1. Pipeline check: the frozen E39 loop is re-run exactly (same predictors, sizes,
   permutations and seed 20260917) and compared row by row with the frozen CSV.
2. New predictors, sizes 4/5/6: ROC AUC, sign-flipped AUC, PR AUC (both
   orientations), precision/recall at k = 10 (E39) and at k = number unstable
   (E13 screening), permutation p (5000 draws, fresh seed 20260918).
3. Matched analyses redone on Pg:
   a. E14 N4: stable four-replacement portfolios inside the failing portfolios'
      replaced-quantity range (+-5 percent), top 25, on Sn (frozen rule, check)
      and on Pg (corrected).
   b. nearest-neighbour Pg matching of every failing size-4 portfolio to a
      distinct stable size-4 portfolio.
   c. E18: the 12 stable size-4 portfolios nearest to the flagship in Pg; the
      closure distance and the longest-cycle spectral radius at the flagship's
      created mode are RECOMPUTED for them (same functions as E18).
   d. E38: the order-4 logistic for the bus-30 indicator, adjusted for Sn (frozen,
      reproduced) and for Pg (corrected).
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from _uc import (
    FROZEN_TABLES,
    OVERNIGHT,
    RESULTS,
    UCExperiment,
    out_dir,
    parse,
    sha256,
    write_json,
)

from _overnight import pin_blas_threads  # noqa: E402
from ibr_cycles.uncertainty.classification import (  # noqa: E402
    average_precision,
    permutation_auc,
    precision_at_k,
    recall_at_k,
    roc_auc,
)
from ibr_cycles.uncertainty.regression import logistic  # noqa: E402

OUT = out_dir("UC02")
SIZES = (4, 5, 6)
PERMUTATIONS = 5000
FROZEN_E39 = OVERNIGHT / "E39_baseline_audit" / "E39_baseline_audit.csv"
FROZEN_E18 = FROZEN_TABLES / "E18_holonomy_matched_contrast.csv"
PREDICTORS = FROZEN_TABLES / "E13_baseline_challenge_predictors.csv"
CENSUS = FROZEN_TABLES / "E12_compatibility_census_portfolios.csv"
E39_EXPECTATION = {
    "additive_prediction": "higher is more unstable",
    "pairwise_prediction": "higher is more unstable",
    "replaced_mw": "higher is more unstable",
    "replaced_inertia_fraction": "higher is more unstable",
    "min_scr": "LOWER is more unstable",
    "mean_scr": "LOWER is more unstable",
    "gscr": "LOWER is more unstable",
    "max_miif": "higher is more unstable",
}
NEW = {
    "replaced_sn_mva": "higher is more unstable",
    "replaced_pg_mw": "higher is more unstable",
    "replaced_pmax_mw": "higher is more unstable",
    "replaced_q_mvar": "higher is more unstable",
}


def score_rows(data, predictors, rng, tag):
    rows = []
    for size in SIZES:
        block = data[data["size"] == size]
        labels = block.unstable.to_numpy(bool)
        if labels.sum() in (0, labels.size):
            continue
        k_star = int(labels.sum())
        for predictor, expectation in predictors.items():
            scores = block[predictor].to_numpy(float)
            auc = roc_auc(labels, scores)
            perm = permutation_auc(labels, scores, draws=PERMUTATIONS, rng=rng)
            rows.append(
                {
                    "pass": tag,
                    "size": size,
                    "n": int(labels.size),
                    "unstable": k_star,
                    "base_rate": float(labels.mean()),
                    "predictor": predictor,
                    "expectation": expectation,
                    "roc_auc": auc,
                    "roc_auc_sign_flipped": 1.0 - auc,
                    "pr_auc": average_precision(labels, scores),
                    "pr_auc_sign_flipped": average_precision(labels, -scores),
                    "precision_at_10": precision_at_k(labels, scores, 10),
                    "recall_at_10": recall_at_k(labels, scores, 10),
                    "precision_at_k_unstable": precision_at_k(labels, scores, k_star),
                    "recall_at_k_unstable": recall_at_k(labels, scores, k_star),
                    "permutation_null_mean": perm["null_mean"],
                    "permutation_p": perm["p_two_sided"],
                }
            )
    return rows


def matched_e14(table, column):
    size4 = table[table["size"] == 4]
    fail = size4[size4.unstable]
    low, high = float(fail[column].min()), float(fail[column].max())
    stable = size4[~size4.unstable]
    window = stable[(stable[column] >= low * 0.95) & (stable[column] <= high * 1.05)]
    top = window.nlargest(25, column)
    scr_lo, scr_hi = float(fail.min_scr.min()), float(fail.min_scr.max())
    joint = window[
        (window.min_scr >= scr_lo - 1e-9) & (window.min_scr <= scr_hi + 1e-9)
    ]
    return {
        "failing_min_scr_range": [scr_lo, scr_hi],
        "stable_in_window_and_failing_min_scr_range": int(len(joint)),
        "joint_worst_abscissa": float(joint.spectral_abscissa.max())
        if len(joint)
        else None,
        "match_on": column,
        "failing_range": [low, high],
        "stable_in_window": int(len(window)),
        "selected": int(len(top)),
        "selected_range": [float(top[column].min()), float(top[column].max())],
        "selected_min_scr_range": [float(top.min_scr.min()), float(top.min_scr.max())],
        "selected_gscr_range": [float(top.gscr.min()), float(top.gscr.max())],
        "selected_worst_abscissa": float(top.spectral_abscissa.max()),
        "all_selected_stable": bool((~top.unstable).all()),
        "members": top.members.tolist(),
    }


def nearest_pairs(table, column):
    size4 = table[table["size"] == 4]
    fail = size4[size4.unstable].sort_values(column)
    pool = size4[~size4.unstable].copy()
    rows = []
    for r in fail.itertuples():
        gap = (pool[column] - getattr(r, column)).abs()
        j = gap.idxmin()
        rows.append(
            {
                "failing": r.members,
                "failing_value": getattr(r, column),
                "matched_stable": pool.loc[j, "members"],
                "matched_value": pool.loc[j, column],
                "relative_gap": float(gap.loc[j] / getattr(r, column)),
                "failing_abscissa": r.spectral_abscissa,
                "matched_abscissa": pool.loc[j, "spectral_abscissa"],
            }
        )
        pool = pool.drop(index=j)
    return pd.DataFrame(rows)


def e18_rematch(table):
    """Recompute E18's matched contrast with Pg matching (same functions as E18)."""

    from E18_holonomy import N_MATCHED, directed_cycles, scores
    from ibr_cycles.dynamics.modes import eigen_analysis
    from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
    from ibr_cycles.models.port_admittance import build_action_space

    core = (30, 33, 35, 37)
    base = solve_case(ReplacementPlan.of({}))
    flagship = solve_case(ReplacementPlan.of({b: 1.0 for b in core}))
    created = max(
        (m for m in eigen_analysis(flagship.system.A).modes if abs(m.value) > 1e-3),
        key=lambda m: m.real,
    ).value
    cycles = directed_cycles(len(core))
    flag_row = table[table.members == "+".join(map(str, core))].iloc[0]
    size4 = table[(table["size"] == 4) & (~table.unstable)].copy()
    size4["gap"] = (size4.replaced_pg_mw - flag_row.replaced_pg_mw).abs()
    matched = size4.nsmallest(N_MATCHED, "gap")
    rows = []
    space = build_action_space(base, flagship, core)
    for label, members, case, sp in [("flagship", core, flagship, space)] + [
        ("control", parse(m), None, None) for m in matched.members
    ]:
        if case is None:
            case = solve_case(ReplacementPlan.of({b: 1.0 for b in members}))
            sp = build_action_space(base, case, members)
        m = sp.m(created)
        rec = {
            "role": label,
            "portfolio": "+".join(map(str, members)),
            "replaced_pg_mw": float(
                table.loc[
                    table.members == "+".join(map(str, members)), "replaced_pg_mw"
                ].iloc[0]
            ),
            "replaced_sn_mva": float(
                table.loc[
                    table.members == "+".join(map(str, members)), "replaced_sn_mva"
                ].iloc[0]
            ),
            "distance_of_q_eig_to_minus_one": float(
                abs(sp.split(created)["closest_to_minus_one"] + 1.0)
            ),
        }
        for length in range(2, len(members) + 1):
            same = [c for c in cycles if len(c) == length]
            rec[f"max_rho_len{length}"] = max(
                scores(sp, m, c)["spectral_radius"] for c in same
            )
        rows.append(rec)
    frame = pd.DataFrame(rows)
    controls = frame[frame.role == "control"]
    flag = frame.iloc[0]
    return frame, {
        "controls": int(len(controls)),
        "controls_pg_range": [
            float(controls.replaced_pg_mw.min()),
            float(controls.replaced_pg_mw.max()),
        ],
        "flagship_pg": float(flag.replaced_pg_mw),
        "flagship_closure_distance": float(flag.distance_of_q_eig_to_minus_one),
        "controls_closure_distance_range": [
            float(controls.distance_of_q_eig_to_minus_one.min()),
            float(controls.distance_of_q_eig_to_minus_one.max()),
        ],
        "closure_separation_ratio": float(
            controls.distance_of_q_eig_to_minus_one.min()
            / flag.distance_of_q_eig_to_minus_one
        ),
        "cycle_separation_ratio_len4": float(
            flag.max_rho_len4 / controls.max_rho_len4.max()
        ),
        "overlap_with_frozen_sn_matched": sorted(
            set(controls.portfolio) & set(pd.read_csv(FROZEN_E18).portfolio.iloc[1:])
        ),
    }


def e38_logistic(merged, column):
    order4 = merged[merged["size"] == 4].copy()
    order4["bus30"] = order4.members.map(lambda m: float(30 in parse(m)))
    terms = ("const", "bus30", column, "inertia", "min_scr", "gscr")
    design = np.column_stack(
        [
            np.ones(len(order4)),
            order4.bus30.to_numpy(float),
            order4[column].to_numpy(float) / 1000.0,
            order4.replaced_inertia_fraction.to_numpy(float),
            order4.min_scr.to_numpy(float),
            order4.gscr.to_numpy(float),
        ]
    )
    model = logistic(design, order4.unstable.to_numpy(float), terms)
    i = model.index("bus30")
    c, e = float(model.coefficients[i]), float(model.standard_errors[i])
    return {
        "adjusted_for": column,
        "n": int(len(order4)),
        "bus30_odds_ratio": float(np.exp(c)),
        "odds_ratio_ci": [
            float(np.exp(np.clip(c - 1.96 * e, -50, 50))),
            float(np.exp(np.clip(c + 1.96 * e, -50, 50))),
        ],
        "bus30_p_value": model.p_of("bus30"),
        "converged": bool(model.converged),
        "pseudo_r2": model.goodness,
        f"coef_{column}_per_1000": float(model.coefficients[model.index(column)]),
        f"p_{column}": model.p_of(column),
    }


def main() -> int:
    pin_blas_threads()
    exp = UCExperiment(
        name="UC02_baseline_reaudit",
        question="Does true removed active power remain a useful baseline?",
        config={
            "frozen_e39": sha256(FROZEN_E39),
            "predictors": sha256(PREDICTORS),
            "census": sha256(CENSUS),
            "permutations": PERMUTATIONS,
            "seeds": {"reproduction": 20260917, "new": 20260918},
        },
    )
    started = time.time()
    data = pd.read_csv(PREDICTORS)
    uc01 = pd.read_csv(RESULTS / "UC01" / "UC01_census_quantities.csv")
    q = uc01[
        [
            "members",
            "replaced_sn_mva",
            "replaced_pg_mw",
            "replaced_pmax_mw",
            "replaced_q_mvar",
        ]
    ]
    data = data.merge(q, on="members", how="left")
    assert data.replaced_pg_mw.notna().all()
    assert np.allclose(data.replaced_sn_mva, data.replaced_mw)

    # 1. reproduce the frozen E39 table exactly
    repro = pd.DataFrame(
        score_rows(data, E39_EXPECTATION, np.random.default_rng(20260917), "repro")
    )
    frozen = pd.read_csv(FROZEN_E39)
    cols = ["roc_auc", "pr_auc", "precision_at_10", "recall_at_10", "permutation_p"]
    joined = repro.merge(frozen, on=["size", "predictor"], suffixes=("", "_frozen"))
    repro_err = float(
        max(np.abs(joined[c] - joined[f"{c}_frozen"]).max() for c in cols)
    )

    # 2. corrected predictors
    new = pd.DataFrame(
        score_rows(data, NEW, np.random.default_rng(20260918), "corrected")
    )
    table = pd.concat([repro, new], ignore_index=True)
    table.to_csv(OUT / "UC02_baseline_scores.csv", index=False)

    # 3. matched analyses
    census = pd.read_csv(CENSUS)[["members", "size", "unstable", "spectral_abscissa"]]
    full = census.merge(
        uc01[["members", "replaced_sn_mva", "replaced_pg_mw"]], on="members"
    )
    full = full.merge(data[["members", "min_scr", "gscr"]], on="members", how="left")
    # the census stores unstable as text with one blank (the size-9 row); an
    # object column would make ~unstable a bitwise NOT on Python bools
    full["unstable"] = full.unstable.map(
        {True: True, False: False, "True": True, "False": False}
    )
    full = full[full.unstable.notna()].copy()
    full["unstable"] = full.unstable.astype(bool)
    e14 = {c: matched_e14(full, c) for c in ("replaced_sn_mva", "replaced_pg_mw")}
    frozen_n4 = pd.read_csv(FROZEN_TABLES / "E14_negative_controls_N4_matched.csv")
    frozen_members = set(
        frozen_n4.get("members", frozen_n4.get("core", pd.Series(dtype=str))).astype(
            str
        )
    )
    e14["overlap_pg_with_frozen"] = len(
        set(e14["replaced_pg_mw"]["members"]) & frozen_members
    )
    pairs = {c: nearest_pairs(full, c) for c in ("replaced_sn_mva", "replaced_pg_mw")}
    for c, frame in pairs.items():
        frame.to_csv(OUT / f"UC02_nearest_pairs_{c}.csv", index=False)
    e18_frame, e18 = e18_rematch(full)
    e18_frame.to_csv(OUT / "UC02_E18_pg_matched_contrast.csv", index=False)
    merged = data.merge(
        pd.read_csv(CENSUS)[["members", "critical_family"]], on="members", how="left"
    )
    e38 = [e38_logistic(merged, "replaced_mw"), e38_logistic(merged, "replaced_pg_mw")]

    def auc(size, predictor, tag):
        r = table[
            (table["size"] == size)
            & (table.predictor == predictor)
            & (table["pass"] == tag)
        ]
        return float(r.roc_auc.iloc[0])

    summary = {
        "frozen_E39_reproduction_max_abs_error": repro_err,
        "auc_old_sn_vs_pg": {
            s: {
                "sn_mva_old_column": auc(s, "replaced_mw", "repro"),
                "pg_mw": auc(s, "replaced_pg_mw", "corrected"),
                "pmax_mw": auc(s, "replaced_pmax_mw", "corrected"),
                "q_mvar": auc(s, "replaced_q_mvar", "corrected"),
                "inertia_fraction": auc(s, "replaced_inertia_fraction", "repro"),
                "min_scr_raw": auc(s, "min_scr", "repro"),
            }
            for s in SIZES
        },
        "spearman_with_inertia_fraction": {
            s: {
                c: float(
                    data[data["size"] == s][c]
                    .rank()
                    .corr(data[data["size"] == s].replaced_inertia_fraction.rank())
                )
                for c in ("replaced_sn_mva", "replaced_pg_mw", "replaced_pmax_mw")
            }
            for s in SIZES
        },
        "e14_matched": {
            k: (
                {kk: vv for kk, vv in v.items() if kk != "members"}
                if isinstance(v, dict)
                else v
            )
            for k, v in e14.items()
        },
        "nearest_pairs_within_5pct": {
            c: int((f.relative_gap <= 0.05).sum()) for c, f in pairs.items()
        },
        "nearest_pairs_median_gap": {
            c: float(f.relative_gap.median()) for c, f in pairs.items()
        },
        "e18_pg_matched": e18,
        "e38": e38,
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "UC02_summary.json", summary)
    exp.finish("COMPUTED", **summary)
    pd.set_option("display.width", 250)
    print(
        table[
            [
                "pass",
                "size",
                "predictor",
                "roc_auc",
                "pr_auc",
                "pr_auc_sign_flipped",
                "precision_at_10",
                "recall_at_10",
                "precision_at_k_unstable",
                "permutation_p",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:7.4f}")
    )
    import json

    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
