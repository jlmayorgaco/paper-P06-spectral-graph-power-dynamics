"""UC03: planning and mitigation results re-scored on separate, correct axes.

Every frozen planning result that reported "PV MW kept / forgone" used machine
ratings Sn [MVA]. Here each option keeps its frozen stability outcome (spectral
abscissa, closure; nothing re-optimized) and its costs are re-expressed on four
axes that are never merged:

    pg_mw_forgone            converter active dispatch given up [MW], MEASURED
                             from the solved plan (flagship Pg minus plan Pg)
    sn_mva_forgone           converter rating given up [MVA] (the frozen column)
    synchronous_mva_added    machine or condenser rating kept or added [MVA]
    controller_change_norm   unchanged

Inputs (read only): results/tables/E21_minimal_repair_strategies.csv,
results/tables/E23_closure_and_frontier_frontier.csv,
results/tables/E11_replacement_atlas_limits.csv, and the E34 frontier table of
the overnight run. The E34 Pareto front is recomputed per target with the
frozen dominance rule on (pg_mw_forgone, synchronous_mva_added,
controller_change_norm) and compared with the frozen front.
"""

from __future__ import annotations

import re
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _uc import (
    CORE,
    FROZEN_TABLES,
    OVERNIGHT,
    RESULTS,
    WORKERS,
    UCExperiment,
    out_dir,
    sha256,
    write_json,
)

from _overnight import pin_blas_threads  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.units import measure  # noqa: E402

OUT = out_dir("UC03")
E34 = OVERNIGHT / "E34_mitigation_frontier" / "E34_mitigation_frontier.csv"
E34_FRONT = OVERNIGHT / "E34_mitigation_frontier" / "E34_pareto_front.csv"
E21 = FROZEN_TABLES / "E21_minimal_repair_strategies.csv"
E23 = FROZEN_TABLES / "E23_closure_and_frontier_frontier.csv"
E11 = FROZEN_TABLES / "E11_replacement_atlas_limits.csv"


def plan_of(strategy: str, detail: str) -> dict[int, float] | None:
    """The replacement map of an E34 option whose PLAN (not controller) changes P."""

    full = {b: 1.0 for b in CORE}
    m = re.fullmatch(r"machine (\d+) kept", detail)
    if strategy.startswith("M4") and m:
        full.pop(int(m.group(1)))
        return full
    m = re.fullmatch(r"keep ([0-9.]+) of machine (\d+)", detail)
    if strategy.startswith("M5") and m:
        full[int(m.group(2))] = 1.0 - float(m.group(1))
        return full
    return None  # M1/M2/M3: all four converters kept at full dispatch


def pg_of(mapping_items) -> float:
    return measure(solve_case(ReplacementPlan.of(dict(mapping_items)))).replaced_pg_mw


def pareto(block, axes):
    values = block[axes].to_numpy()
    keep = []
    for i in range(len(values)):
        dominated = np.all(values <= values[i], axis=1) & np.any(
            values < values[i], axis=1
        )
        keep.append(not dominated.any())
    return block[np.array(keep)]


def main() -> int:
    exp = UCExperiment(
        name="UC03_planning_rescore",
        question="Which planning statements change once MW and MVA are separated?",
        config={"inputs": {p.name: sha256(p) for p in (E34, E34_FRONT, E21, E23, E11)}},
        workers=WORKERS,
    )
    started = time.time()
    flag_pg = pg_of(tuple((b, 1.0) for b in CORE))
    uc01 = pd.read_csv(RESULTS / "UC01" / "UC01_per_machine.csv")
    m39 = uc01[uc01.benchmark == "IEEE-39"].set_index("bus")

    # ---- E34: measure every plan-changing option ----
    e34 = pd.read_csv(E34)
    plans = {}
    for s, d in zip(e34.strategy, e34.detail, strict=True):
        mapping = plan_of(s, d)
        if mapping is not None:
            plans[(s, d)] = tuple(sorted(mapping.items()))
    with Pool(
        min(WORKERS, max(1, len(set(plans.values())))), initializer=pin_blas_threads
    ) as pool:
        unique = sorted(set(plans.values()))
        measured = dict(zip(unique, pool.map(pg_of, unique), strict=True))
    e34["pg_mw_forgone"] = [
        flag_pg - measured[plans[(s, d)]] if (s, d) in plans else 0.0
        for s, d in zip(e34.strategy, e34.detail, strict=True)
    ]
    e34["sn_mva_forgone"] = e34.pv_mw_forgone  # the frozen column was MVA
    e34["pg_mw_kept"] = flag_pg - e34.pg_mw_forgone
    e34.to_csv(OUT / "UC03_E34_rescored.csv", index=False)
    axes_new = ["pg_mw_forgone", "synchronous_mva_added", "controller_change_norm"]
    axes_old = ["pv_mw_forgone", "synchronous_mva_added", "controller_change_norm"]
    fronts_new, fronts_old = [], []
    for name in e34.target.unique():
        block = e34[(e34.target == name) & e34.meets_target].reset_index(drop=True)
        fronts_new.append(pareto(block, axes_new).assign(target=name))
        fronts_old.append(pareto(block, axes_old).assign(target=name))
    front_new = pd.concat(fronts_new)
    front_old = pd.concat(fronts_old)
    frozen_front = pd.read_csv(E34_FRONT)
    front_new.to_csv(OUT / "UC03_E34_pareto_front_pg.csv", index=False)
    key = ["target", "strategy", "detail"]
    old_keys = set(map(tuple, frozen_front[key].astype(str).to_numpy()))
    new_keys = set(map(tuple, front_new[key].astype(str).to_numpy()))
    recomputed_old = set(map(tuple, front_old[key].astype(str).to_numpy()))

    # ---- E21: re-expressed ----
    e21 = pd.read_csv(E21)
    restore = {f"machine {b} kept": b for b in CORE}
    rows = []
    for r in e21.itertuples():
        if r.detail == "base case, all machines":
            pg_kept = 0.0
        elif r.detail in restore:
            pg_kept = flag_pg - float(m39.loc[restore[r.detail], "Pg_MW"])
        else:
            pg_kept = flag_pg
        rows.append(
            {
                "strategy": r.strategy,
                "detail": r.detail,
                "converter_sn_mva_kept": r.pv_mw_retained,
                "converter_sn_mva_forgone": r.pv_mw_forgone,
                "converter_pg_mw_kept": pg_kept,
                "converter_pg_mw_forgone": flag_pg - pg_kept,
                "synchronous_mva_added": r.synchronous_mva_required,
                "relative_parameter_change": r.relative_parameter_change,
                "spectral_abscissa": r.spectral_abscissa,
            }
        )
    e21_new = pd.DataFrame(rows)
    e21_new.to_csv(OUT / "UC03_E21_reexpressed.csv", index=False)

    # ---- E23: restore choice by min Sn (frozen) vs min Pg forgone ----
    ra = e21_new[e21_new.detail.isin(restore)].copy()
    ra["bus"] = ra.detail.map(restore)
    e23 = pd.read_csv(E23)
    e23_rows = []
    for r in e23.itertuples():
        ok = ra[ra.spectral_abscissa <= r.target_abscissa]
        by_sn = ok.loc[ok.synchronous_mva_added.idxmin()]
        by_pg = ok.loc[ok.converter_pg_mw_forgone.idxmin()]
        e23_rows.append(
            {
                "target_abscissa": r.target_abscissa,
                "frozen_restore_bus": r.restore_bus,
                "restore_bus_min_sn": int(by_sn.bus),
                "restore_bus_min_pg_forgone": int(by_pg.bus),
                "same_choice": int(by_sn.bus) == int(by_pg.bus) == int(r.restore_bus),
                "restore_sn_mva": float(by_sn.synchronous_mva_added),
                "restore_pg_mw_forgone": float(by_pg.converter_pg_mw_forgone),
                "restore_converter_pg_mw_kept": float(by_pg.converter_pg_mw_kept),
                "retune_converter_pg_mw_kept": flag_pg,
                "retune_converter_sn_mva_kept": r.retune_pv_retained,
                "condenser_mva": r.condenser_mva,
            }
        )
    e23_new = pd.DataFrame(e23_rows)
    e23_new.to_csv(OUT / "UC03_E23_restore_choice.csv", index=False)

    # ---- E11: replaceable capacity per bus ----
    e11 = pd.read_csv(E11)
    e11["pg_mw_replaceable"] = [
        float(r.rho_star) * float(m39.loc[int(r.bus), "Pg_MW"])
        for r in e11.itertuples()
    ]
    e11["sn_mva_replaceable"] = e11.mw_replaceable  # frozen column was MVA
    e11.to_csv(OUT / "UC03_E11_limits_reexpressed.csv", index=False)

    summary = {
        "flagship_pg_mw": flag_pg,
        "e34_options": int(len(e34)),
        "e34_measured_plans": int(len(unique)),
        "e34_frozen_front_points": int(len(frozen_front)),
        "e34_front_recomputed_with_frozen_axes_equals_frozen": recomputed_old
        == old_keys,
        "e34_front_pg_points": int(len(front_new)),
        "e34_points_added_by_pg_axis": sorted(new_keys - old_keys),
        "e34_points_removed_by_pg_axis": sorted(old_keys - new_keys),
        "e34_m4_ever_pareto_efficient_pg": bool(
            front_new.strategy.str.startswith("M4").any()
        ),
        "e34_zero_pg_forgone_strategies": sorted(
            set(e34[e34.meets_target & (e34.pg_mw_forgone.abs() < 1e-6)].strategy)
        ),
        "e23_restore_choice_unchanged": bool(e23_new.same_choice.all()),
        "e21": e21_new.to_dict("records"),
        "e23": e23_new.to_dict("records"),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "UC03_summary.json", summary)
    exp.finish(
        "COMPUTED", **{k: v for k, v in summary.items() if k not in ("e21", "e23")}
    )
    pd.set_option("display.width", 250)
    print(
        front_new[
            [
                "target",
                "strategy",
                "detail",
                "pg_mw_forgone",
                "sn_mva_forgone",
                "synchronous_mva_added",
                "controller_change_norm",
                "spectral_abscissa",
            ]
        ].to_string(index=False)
    )
    print(e21_new.round(2).to_string(index=False))
    print(e23_new.to_string(index=False))
    print({k: v for k, v in summary.items() if k not in ("e21", "e23")})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
