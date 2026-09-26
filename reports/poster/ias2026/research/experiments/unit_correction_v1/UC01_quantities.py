"""UC01: canonical MW / MVA / Mvar quantities, measured, for every audited portfolio.

(a) Per machine, three benchmarks: Pg, Qg (base power flow), Sn, documented Pmax.
(b) The 511 E12 census portfolios (rho = 1, matched reactive policy, the frozen
    E12 plans): replaced_sn_mva, replaced_pg_mw (MEASURED from the converter
    injections at the solved equilibrium), replaced_q_mvar, replaced_pmax_mw,
    and two checks: replaced_sn_mva against the legacy replaced_mw column, and
    measured Pg against the base-dispatch sum.
(c) The flagship and every named mitigation / planning case: unrepaired, restore
    SG30/SG33/SG35/SG37 (E21 RA, E34 M4), RB/RC converter retunes, RD condenser
    at 25 percent of each retired Sn (E21), the E34 single-bus condensers on the
    Pareto front (ratings 0.16/0.20/0.24/0.26 at bus 30).
    RB/RC are solved with the retuned controller to show that controller gains
    change no power quantity (P_ref = dispatch).

Nothing frozen is rewritten; the legacy replaced_mw is only read.
"""

from __future__ import annotations

import time
from dataclasses import replace
from multiprocessing import Pool

import pandas as pd
from _uc import (
    CORE,
    FROZEN_TABLES,
    ROOT,
    WORKERS,
    UCExperiment,
    out_dir,
    parse,
    sha256,
    write_json,
)

from _overnight import pin_blas_threads  # noqa: E402
from ibr_cycles.models.ieee39_case import (  # noqa: E402
    InfeasibleReplacement,
    ReplacementPlan,
    solve_case,
)
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402
from ibr_cycles.units import measure, per_machine_table  # noqa: E402

OUT = out_dir("UC01")
NETS = {
    "IEEE-39": ROOT / "configs" / "ias2026" / "ieee39_network.json",
    "Kundur": ROOT / "configs" / "kundur" / "kundur_network.json",
    "IEEE-68": ROOT / "configs" / "ieee68" / "ieee68_network.json",
}
CENSUS = FROZEN_TABLES / "E12_compatibility_census_portfolios.csv"
# E21 RC minimum-norm retune (results/tables/E21_minimal_repair_strategies.csv)
RC = {
    "kp_p": 1.226,
    "ki_p": 0.460,
    "kp_q": 1.173,
    "ki_q": 1.273,
    "kp_i": 1.115,
    "ki_i": 0.960,
}
RB_SCALE = 0.4273


def census_row(members: tuple[int, ...]) -> dict:
    plan = ReplacementPlan.of({b: 1.0 for b in members}, q_policy="matched")
    try:
        q = measure(solve_case(plan)).as_dict()
        status = "OK"
    except InfeasibleReplacement as error:  # pragma: no cover - reported, not hidden
        q, status = {}, f"INFEASIBLE: {error}"
    return {"members": "+".join(map(str, members)) or "BASE", "status": status, **q}


def named_cases() -> list[tuple[str, dict, ConverterParameters | None]]:
    base = ConverterParameters()
    rb = replace(base, kp_p=base.kp_p * RB_SCALE, ki_p=base.ki_p * RB_SCALE)
    rc = replace(base, **{k: getattr(base, k) * v for k, v in RC.items()})
    cases = [
        ("flagship unrepaired (30+33+35+37)", {"rho": {b: 1.0 for b in CORE}}, None),
        ("RB weaken P loop x0.427", {"rho": {b: 1.0 for b in CORE}}, rb),
        ("RC minimum-norm retune", {"rho": {b: 1.0 for b in CORE}}, rc),
        (
            "RD condenser 0.25 x Sn at each retired bus",
            {"rho": {b: 1.0 for b in CORE}, "condenser": {b: 0.25 for b in CORE}},
            None,
        ),
    ]
    for bus in CORE:
        cases.append(
            (f"restore SG{bus}", {"rho": {b: 1.0 for b in CORE if b != bus}}, None)
        )
    for rating in (0.16, 0.20, 0.24, 0.26):
        cases.append(
            (
                f"E34 condenser {rating:.2f} x Sn at bus 30",
                {"rho": {b: 1.0 for b in CORE}, "condenser": {30: rating}},
                None,
            )
        )
    return cases


def named_row(spec) -> dict:
    label, plan_spec, converter = spec
    plan = ReplacementPlan.of(plan_spec["rho"], condenser=plan_spec.get("condenser"))
    case = solve_case(plan, converter=converter) if converter else solve_case(plan)
    q = measure(case).as_dict()
    return {"case": label, "legacy_replaced_mw_column": case.replaced_mw, **q}


def main() -> int:
    exp = UCExperiment(
        name="UC01_quantities",
        question="What are the MW, MVA and Mvar of every audited portfolio, measured?",
        config={
            "census": str(CENSUS),
            "census_sha256": sha256(CENSUS),
            "networks": {k: sha256(v) for k, v in NETS.items()},
        },
        workers=WORKERS,
    )
    started = time.time()
    machines = []
    for name, path in NETS.items():
        for row in per_machine_table(load_network(path)):
            machines.append({"benchmark": name, **row})
    machines = pd.DataFrame(machines)
    machines.to_csv(OUT / "UC01_per_machine.csv", index=False)

    census = pd.read_csv(CENSUS)
    members = [parse(m) for m in census.members]
    with Pool(WORKERS, initializer=pin_blas_threads) as pool:
        rows = pool.map(census_row, members, chunksize=4)
        named = pool.map(named_row, named_cases(), chunksize=1)
    quantities = pd.DataFrame(rows)
    table = census[
        ["members", "size", "unstable", "spectral_abscissa", "replaced_mw"]
    ].merge(quantities, on="members", how="left")
    table = table.rename(columns={"replaced_mw": "legacy_replaced_mw_column_MVA"})
    table["sn_vs_legacy"] = (
        table.replaced_sn_mva - table.legacy_replaced_mw_column_MVA
    ).abs()
    table["pg_vs_base_dispatch"] = (
        table.replaced_pg_mw - table.pg_mw_base_dispatch
    ).abs()
    table["legacy_over_pg"] = table.legacy_replaced_mw_column_MVA / table.replaced_pg_mw
    table.to_csv(OUT / "UC01_census_quantities.csv", index=False)
    named = pd.DataFrame(named)
    named.to_csv(OUT / "UC01_named_cases.csv", index=False)

    m39 = machines[machines.benchmark == "IEEE-39"].set_index("bus")
    flag = m39.loc[list(CORE), ["Pg_MW", "Pmax_MW", "Qg_Mvar", "Sn_MVA"]]
    flag.loc["total"] = flag.sum()
    flag.to_csv(OUT / "UC01_flagship_per_bus.csv")
    nonbase = table[table.members != "BASE"]
    summary = {
        "census_portfolios": int(len(nonbase)),
        "census_status": nonbase.status.value_counts().to_dict(),
        "max_abs_sn_minus_legacy_MVA": float(nonbase.sn_vs_legacy.max()),
        "max_abs_pg_measured_minus_base_dispatch_MW": float(
            nonbase.pg_vs_base_dispatch.max()
        ),
        "legacy_over_pg_range": [
            float(nonbase.legacy_over_pg.min()),
            float(nonbase.legacy_over_pg.max()),
        ],
        "spearman_legacy_vs_pg": float(
            nonbase.legacy_replaced_mw_column_MVA.rank().corr(
                nonbase.replaced_pg_mw.rank()
            )
        ),
        "flagship_per_bus": flag.reset_index().to_dict("records"),
        "named_cases": named.to_dict("records"),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "UC01_summary.json", summary)
    exp.finish("COMPUTED", **{k: v for k, v in summary.items() if k != "named_cases"})
    pd.set_option("display.width", 250)
    print(flag.round(2).to_string())
    print(
        named[
            [
                "case",
                "replaced_sn_mva",
                "replaced_pg_mw",
                "replaced_q_mvar",
                "replaced_pmax_mw",
                "sg_sn_mva_left_at_replaced_buses",
                "condenser_sn_mva",
            ]
        ]
        .round(2)
        .to_string(index=False)
    )
    print(
        {
            k: v
            for k, v in summary.items()
            if k not in ("named_cases", "flagship_per_bus")
        }
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
