"""E38 - overnight 10. Is bus 30 individually ordinary but collectively enabling?

Uses the already-computed 512-portfolio census (E12) and the conventional
predictors (E13), plus one quantity computed here: the participation of the base
inter-area mode in each machine, which is the natural "is this machine part of
the mode" measure and was not in either table.

The hypothesis under test is that bus 30 is unremarkable on every individual
measure and yet turns up disproportionately in the portfolios that fail. The
discipline the protocol imposes is that this may NOT be claimed if bus 30 loses
significance once megawatts and inertia are controlled for.
"""

from __future__ import annotations

import time
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

from _bootstrap import RESULTS
from _overnight import Experiment, pin_blas_threads
from _v2c_common import CORE, base_anchor, nominal_reference
from ibr_cycles.dynamics.modes import eigen_analysis
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case
from ibr_cycles.uncertainty.regression import logistic

NAME = "E38_bus30_collective_role"
FOCUS = 30
PERMUTATIONS = 20000


def parse(members: str) -> tuple[int, ...]:
    return () if members == "BASE" else tuple(int(x) for x in members.split("+"))


def main() -> int:
    pin_blas_threads()
    experiment = Experiment(
        name=NAME,
        question="Is bus 30 individually ordinary but collectively enabling?",
        config={
            "census": "results/tables/E12_compatibility_census_portfolios.csv",
            "predictors": "results/tables/E13_baseline_challenge_predictors.csv",
            "permutations": PERMUTATIONS,
            "rule": "no claim if bus 30 loses significance after controlling MW and inertia",
        },
        workers=1,
    )
    started = time.time()
    rng = np.random.default_rng(20260916)

    census = pd.read_csv(RESULTS / "tables" / "E12_compatibility_census_portfolios.csv")
    predictors = pd.read_csv(
        RESULTS / "tables" / "E13_baseline_challenge_predictors.csv"
    )
    census["parsed"] = census.members.map(parse)
    census["minimal_core"] = census.unstable & (
        census.all_proper_subsets_stable == True  # noqa: E712
    )
    candidates = sorted({b for parsed in census.parsed for b in parsed})
    experiment.note(f"{len(census)} portfolios, {len(candidates)} candidate buses")

    # ---- participation of the base inter-area mode, computed here ----------
    base = solve_case(ReplacementPlan.of({}))
    anchor = base_anchor(eigen_analysis(base.system.A), base, nominal_reference())
    participation = {}
    for bus in candidates:
        indices = [
            i
            for i, name in enumerate(base.system.labels)
            if name in (f"delta_sg{bus}", f"omega_sg{bus}")
        ]
        participation[bus] = float(anchor.participation[indices].sum())

    # ---- per-bus frequencies ------------------------------------------------
    minimal = census[census.minimal_core]
    inter_area = minimal[minimal.critical_family == "sg"]
    pll = census[census.unstable & (census.critical_family == "gfl")]
    singles = census[census["size"] == 1].set_index("members")

    rows = []
    for bus in candidates:
        present = census.parsed.map(lambda p, b=bus: b in p)
        single = singles.loc[str(bus)] if str(bus) in singles.index else None
        block = predictors[predictors.members.map(lambda m, b=bus: b in parse(m))]
        rows.append(
            {
                "bus": bus,
                "in_minimal_cores": int(present[minimal.index].sum()),
                "minimal_cores_total": int(len(minimal)),
                "share_of_minimal_cores": float(present[minimal.index].mean()),
                "in_inter_area_cores": int(present[inter_area.index].sum()),
                "inter_area_cores_total": int(len(inter_area)),
                "share_of_inter_area_cores": float(
                    present[inter_area.index].mean() if len(inter_area) else np.nan
                ),
                "in_pll_cores": int(present[pll.index].sum()),
                "pll_cores_total": int(len(pll)),
                "share_of_pll_cores": float(
                    present[pll.index].mean() if len(pll) else np.nan
                ),
                "individual_abscissa": float(single.spectral_abscissa)
                if single is not None
                else np.nan,
                "individual_unstable": bool(single.unstable)
                if single is not None
                else None,
                "mean_min_scr_when_present": float(block.min_scr.mean()),
                "mean_gscr_when_present": float(block.gscr.mean()),
                "replaced_mw_alone": float(single.replaced_mw)
                if single is not None
                else np.nan,
                "base_mode_participation": participation[bus],
            }
        )
    per_bus = pd.DataFrame(rows).sort_values("share_of_minimal_cores", ascending=False)
    experiment.save_table(per_bus, "E38_bus30_collective_role.csv")

    # ---- order-4 logistic, controlling for MW, inertia and short-circuit ----
    merged = predictors.merge(
        census[["members", "critical_family", "minimal_core"]], on="members", how="left"
    )
    order4 = merged[merged["size"] == 4].copy()
    order4["bus30"] = order4.members.map(lambda m: float(FOCUS in parse(m)))

    def fit(block, label):
        terms = ("const", "bus30", "replaced_mw", "inertia", "min_scr", "gscr")
        design = np.column_stack(
            [
                np.ones(len(block)),
                block.bus30.to_numpy(float),
                block.replaced_mw.to_numpy(float) / 1000.0,
                block.replaced_inertia_fraction.to_numpy(float),
                block.min_scr.to_numpy(float),
                block.gscr.to_numpy(float),
            ]
        )
        model = logistic(design, block.unstable.to_numpy(float), terms)
        i = model.index("bus30")
        coefficient = float(model.coefficients[i])
        error = float(model.standard_errors[i])
        return {
            "population": label,
            "n": int(len(block)),
            "unstable": int(block.unstable.sum()),
            "bus30_coefficient": coefficient,
            "bus30_odds_ratio": float(np.exp(coefficient)),
            "odds_ratio_ci_low": float(np.exp(np.clip(coefficient - 1.96 * error, -50, 50))),
            "odds_ratio_ci_high": float(np.exp(np.clip(coefficient + 1.96 * error, -50, 50))),
            "standard_error": error,
            "bus30_p_value": model.p_of("bus30"),
            "converged": bool(model.converged),
            "pseudo_r2": model.goodness,
            **{
                f"coef_{name}": float(model.coefficients[model.index(name)])
                for name in terms
                if name != "const"
            },
        }

    fits = [fit(order4, "all order-4 portfolios")]
    inter_area_order4 = order4[
        (order4.critical_family == "sg") | (~order4.unstable)
    ]
    fits.append(fit(inter_area_order4, "order-4, inter-area family only"))
    fit_frame = pd.DataFrame(fits)

    # ---- Fisher exact, unadjusted ------------------------------------------
    with30 = order4[order4.bus30 == 1.0]
    without30 = order4[order4.bus30 == 0.0]
    contingency = [
        [int(with30.unstable.sum()), int((~with30.unstable).sum())],
        [int(without30.unstable.sum()), int((~without30.unstable).sum())],
    ]
    odds, fisher_p = fisher_exact(contingency)

    # ---- matched pairs: order-4 portfolios sharing three members -----------
    lookup = {frozenset(parse(m)): bool(u) for m, u in zip(order4.members, order4.unstable, strict=True)}
    wins = losses = ties = 0
    for triple in combinations(sorted(set(candidates) - {FOCUS}), 3):
        with_focus = frozenset(triple) | {FOCUS}
        if with_focus not in lookup:
            continue
        for other in sorted(set(candidates) - {FOCUS} - set(triple)):
            alternative = frozenset(triple) | {other}
            if alternative not in lookup:
                continue
            a, b = lookup[with_focus], lookup[alternative]
            if a and not b:
                wins += 1
            elif b and not a:
                losses += 1
            else:
                ties += 1
    # Matched pairs call for McNemar's exact test, a binomial on the discordant
    # pairs. A 2x2 Fisher table is the wrong test here: the pairs are not
    # independent samples, they are the same three co-members with one swap.
    matched_p = (
        float(binomtest(wins, wins + losses, 0.5).pvalue)
        if (wins + losses)
        else float("nan")
    )

    # ---- permutation test on the bus-30 indicator --------------------------
    labels = order4.unstable.to_numpy(bool)
    indicator = order4.bus30.to_numpy(bool)
    observed = labels[indicator].mean() - labels[~indicator].mean()
    null = np.empty(PERMUTATIONS)
    for i in range(PERMUTATIONS):
        null[i] = (
            labels[rng.permutation(indicator)].mean()
            - labels[~rng.permutation(indicator)].mean()
        )
    # a cleaner null: permute the labels, keep the indicator fixed
    null = np.empty(PERMUTATIONS)
    for i in range(PERMUTATIONS):
        shuffled = rng.permutation(labels)
        null[i] = shuffled[indicator].mean() - shuffled[~indicator].mean()
    permutation_p = float((np.abs(null) >= abs(observed)).mean())

    summary = {
        "minimal_cores": int(len(minimal)),
        "inter_area_cores": int(len(inter_area)),
        "pll_cores": int(len(pll)),
        "bus30_share_of_minimal_cores": float(
            per_bus[per_bus.bus == FOCUS].share_of_minimal_cores.iloc[0]
        ),
        "bus30_rank_by_minimal_share": int(
            per_bus.reset_index(drop=True).index[per_bus.reset_index(drop=True).bus == FOCUS][0]
        )
        + 1,
        "bus30_individual_abscissa": float(
            per_bus[per_bus.bus == FOCUS].individual_abscissa.iloc[0]
        ),
        "bus30_individual_abscissa_rank": int(
            per_bus.sort_values("individual_abscissa", ascending=False)
            .reset_index(drop=True)
            .index[
                per_bus.sort_values("individual_abscissa", ascending=False)
                .reset_index(drop=True)
                .bus
                == FOCUS
            ][0]
        )
        + 1,
        "bus30_participation": participation[FOCUS],
        "bus30_participation_rank": int(
            per_bus.sort_values("base_mode_participation", ascending=False)
            .reset_index(drop=True)
            .index[
                per_bus.sort_values("base_mode_participation", ascending=False)
                .reset_index(drop=True)
                .bus
                == FOCUS
            ][0]
        )
        + 1,
        "fisher_contingency": contingency,
        "fisher_odds_ratio": float(odds),
        "fisher_p": float(fisher_p),
        "matched_pairs_wins": wins,
        "matched_pairs_losses": losses,
        "matched_pairs_ties": ties,
        "matched_pairs_p_mcnemar_exact": matched_p,
        "permutation_difference": float(observed),
        "permutation_p": permutation_p,
        "adjusted_odds_ratio": fits[0]["bus30_odds_ratio"],
        "adjusted_p": fits[0]["bus30_p_value"],
        "survives_adjustment": bool(fits[0]["bus30_p_value"] < 0.05),
    }
    experiment.save_table(fit_frame, "E38_bus30_logistic.csv")

    print()
    print(
        per_bus[
            [
                "bus",
                "share_of_minimal_cores",
                "share_of_inter_area_cores",
                "share_of_pll_cores",
                "individual_abscissa",
                "base_mode_participation",
                "replaced_mw_alone",
                "mean_min_scr_when_present",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:9.4f}")
    )
    print()
    print(fit_frame.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))
    print()
    for key, value in summary.items():
        print(f"  {key:34s} {value}")

    verdict = (
        "SUPPORTED"
        if summary["survives_adjustment"] and summary["fisher_p"] < 0.05
        else "NOT SUPPORTED AFTER ADJUSTMENT"
    )
    experiment.finish(verdict, summary=summary, fits=fits, elapsed_s=time.time() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
