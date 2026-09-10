"""BC00-A: MW, MVA and bases, traced column by column from the source data.

Question (theory note, section 1.2): PHASE_E_REPORT.md calls 4270.7 "MW of PV",
and 1040 + 1174.8 + 1085.7 + 970.2 = 4270.7 are the machine ratings Sn in MVA.

Trace
    ReplacementCase.replaced_mw (src/ibr_cycles/models/ieee39_case.py) returns
    sum(slot.weight * SYSTEM_BASE_MVA) over converter slots. A converter's weight
    is rho * Sn / 100, so the property is the CONVERTER RATING IN MVA, not active
    power. E12, E13, E14, E18, E21, E23, E34 and E39 all read it.

Output (new tables; nothing earlier is rewritten)
    BC00_units_per_machine.csv   per benchmark and machine: base, Sn, converter
                                 rating, nameplate MW (PV.pmax of the source case
                                 where it exists), P/Q dispatch, utilization
    BC00_units_per_portfolio.csv the IEEE-39 core portfolios and the E21 repairs:
                                 active MW actually displaced vs the MVA figure
                                 reported as MW
    BC00_E39_ranker_corrected.csv the E39 "replaced MW" ranker recomputed with the
                                 displaced active power (the old column was MVA)
    BC00_units_summary.json
"""

from __future__ import annotations

import json
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from _bc import CONFIGS, ROOT, BCExperiment, out_dir, sha256, write_json

sys.path.insert(0, str(ROOT / "experiments"))
from ibr_cycles.models.ieee39_network import (  # noqa: E402
    load_network,
    solve_power_flow,
)
from ibr_cycles.uncertainty.classification import roc_auc  # noqa: E402

OUT = out_dir("BC00")
NETS = {
    "IEEE-39": ROOT / "configs" / "ias2026" / "ieee39_network.json",
    "Kundur": ROOT / "configs" / "kundur" / "kundur_network.json",
    "IEEE-68": ROOT / "configs" / "ieee68" / "ieee68_network.json",
}
CORE = (30, 33, 35, 37)
BASE_MVA = 100.0


def machine_table() -> pd.DataFrame:
    rows = []
    for name, path in NETS.items():
        payload = json.loads(path.read_text(encoding="utf-8"))
        net = load_network(path)
        pf = solve_power_flow(net)
        pv = {int(r["bus"]): r for r in payload.get("pv", [])}
        slack = {int(r["bus"]): r for r in payload.get("slack", [])}
        loading = float(payload.get("model", {}).get("converter_loading", 0.0))
        for m in payload["machines"]:
            bus = int(m["bus"])
            gen = pf.injection(bus, net.ybus) + net.loads.get(bus, 0j)
            s_mva = abs(gen) * BASE_MVA
            sn = float(m["Sn"])
            src = pv.get(bus) or slack.get(bus) or {}
            pmax = src.get("pmax")
            nameplate = (
                float(pmax) * BASE_MVA
                if pmax is not None and name != "IEEE-68" and float(pmax) < 90
                else np.nan
            )
            inverter = s_mva / loading if loading > 0 else sn
            v = abs(pf.at(bus))
            rows.append(
                {
                    "benchmark": name,
                    "bus": bus,
                    "slack": bus in slack,
                    "S_base_MVA": BASE_MVA,
                    "Sn_machine_MVA": sn,
                    "Sn_meaning": "per-unit base, not a rating"
                    if name == "IEEE-68"
                    else "machine rating",
                    "inverter_rating_MVA": inverter,
                    "inverter_rating_rule": f"|S_gen|/{loading:g}"
                    if loading > 0
                    else "= Sn (weight rho*Sn/100)",
                    "active_nameplate_MW": nameplate,
                    "nameplate_source": "PV/Slack.pmax of the source case"
                    if nameplate == nameplate
                    else "none valid",
                    "P_dispatch_pu": gen.real,
                    "P_dispatch_MW": gen.real * BASE_MVA,
                    "Q_dispatch_Mvar": gen.imag * BASE_MVA,
                    "S_dispatch_MVA": s_mva,
                    "available_PV_MW": np.nan,
                    "available_PV_note": "NOT MODELLED: P_ref = dispatch",
                    "retired_active_MW_if_replaced": gen.real * BASE_MVA,
                    "P_util_of_inverter": gen.real * BASE_MVA / inverter,
                    "Q_util_of_inverter": abs(gen.imag) * BASE_MVA / inverter,
                    "current_util_of_inverter": s_mva / (v * inverter),
                    "S_util_of_Sn": s_mva / sn,
                    "V_pu": v,
                }
            )
    return pd.DataFrame(rows)


def portfolio_table(machines: pd.DataFrame) -> pd.DataFrame:
    m39 = machines[machines.benchmark == "IEEE-39"].set_index("bus")
    rows = []
    candidates = [b for b in m39.index if not m39.loc[b, "slack"]]
    for r in range(1, 5):
        for s in combinations(CORE, r):
            rows.append(("core subset", "+".join(map(str, s)), s, 0.0))
    # E21 repairs of the flagship (PHASE_E_REPORT.md table), condenser at 25 % of Sn
    rows += [
        ("E21 unrepaired", "30+33+35+37", CORE, 0.0),
        ("E21 RA restore 30", "33+35+37", (33, 35, 37), 0.0),
        ("E21 RA restore 37", "30+33+35", (30, 33, 35), 0.0),
        ("E21 RB / RC retune", "30+33+35+37", CORE, 0.0),
        ("E21 RD condenser 25 %", "30+33+35+37", CORE, 0.25),
    ]
    out = []
    for kind, label, members, cond in rows:
        sub = m39.loc[list(members)]
        mva = float(sub.Sn_machine_MVA.sum())
        out.append(
            {
                "benchmark": "IEEE-39",
                "kind": kind,
                "portfolio": label,
                "inverter_rating_MVA": mva,
                "reported_as_MW_by_replaced_mw": mva,
                "displaced_active_MW": float(sub.P_dispatch_MW.sum()),
                "displaced_reactive_Mvar": float(sub.Q_dispatch_Mvar.sum()),
                "active_nameplate_MW": float(sub.active_nameplate_MW.sum()),
                "condenser_MVA": cond * mva,
                "ratio_reported_to_active": mva / float(sub.P_dispatch_MW.sum()),
            }
        )
    frame = pd.DataFrame(out)
    frame.attrs["candidates"] = candidates
    return frame


def e39_corrected(machines: pd.DataFrame) -> pd.DataFrame:
    table = pd.read_csv(
        ROOT / "results" / "tables" / "E13_baseline_challenge_predictors.csv"
    )
    m39 = machines[machines.benchmark == "IEEE-39"].set_index("bus")
    table["displaced_active_MW"] = [
        float(sum(m39.loc[int(b), "P_dispatch_MW"] for b in str(members).split("+")))
        for members in table.members
    ]
    table["rating_check"] = [
        float(sum(m39.loc[int(b), "Sn_machine_MVA"] for b in str(members).split("+")))
        for members in table.members
    ]
    rows = []
    for size in (4, 5, 6):
        g = table[table["size"] == size]
        y = g.unstable.astype(bool).to_numpy()
        if y.all() or not y.any():
            continue
        for col, label in (
            ("replaced_mw", "old column (MVA rating, labelled MW)"),
            ("displaced_active_MW", "displaced active power (MW)"),
        ):
            rows.append(
                {
                    "size": size,
                    "n": int(len(g)),
                    "unstable": int(y.sum()),
                    "predictor": label,
                    "roc_auc_higher_is_more_unstable": float(
                        roc_auc(y, g[col].to_numpy())
                    ),
                }
            )
    frame = pd.DataFrame(rows)
    frame.attrs["rating_matches_old_column"] = bool(
        np.allclose(table.rating_check, table.replaced_mw)
    )
    return frame


def main() -> int:
    exp = BCExperiment(
        name="BC00_units_audit",
        question="Are MW and MVA kept apart in every reported quantity?",
        config={
            "networks": {k: str(v) for k, v in NETS.items()},
            "sha256": {k: sha256(v) for k, v in NETS.items()},
        },
    )
    machines = machine_table()
    machines.to_csv(OUT / "BC00_units_per_machine.csv", index=False)
    portfolios = portfolio_table(machines)
    portfolios.to_csv(OUT / "BC00_units_per_portfolio.csv", index=False)
    auc = e39_corrected(machines)
    auc.to_csv(OUT / "BC00_E39_ranker_corrected.csv", index=False)
    flag = portfolios[portfolios.portfolio == "30+33+35+37"].iloc[0]
    summary = {
        "replaced_mw_is_mva": bool(auc.attrs["rating_matches_old_column"]),
        "flagship_reported_MW": flag.reported_as_MW_by_replaced_mw,
        "flagship_displaced_active_MW": flag.displaced_active_MW,
        "flagship_active_nameplate_MW": flag.active_nameplate_MW,
        "flagship_ratio_reported_to_active": flag.ratio_reported_to_active,
        "E21_rows": portfolios[portfolios.kind.str.startswith("E21")]
        .drop(columns=["benchmark"])
        .to_dict("records"),
        "E39_ranker": auc.to_dict("records"),
        "config_dir": str(CONFIGS),
    }
    write_json(OUT / "BC00_units_summary.json", summary)
    exp.finish("COMPUTED", **{k: v for k, v in summary.items() if k != "E21_rows"})
    pd.set_option("display.width", 220)
    print(
        machines[
            [
                "benchmark",
                "bus",
                "Sn_machine_MVA",
                "inverter_rating_MVA",
                "active_nameplate_MW",
                "P_dispatch_MW",
                "Q_dispatch_Mvar",
                "current_util_of_inverter",
            ]
        ]
        .round(1)
        .to_string(index=False)
    )
    print(portfolios.round(1).to_string(index=False))
    print(auc.to_string(index=False))
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k not in ("E21_rows", "E39_ranker")},
            indent=1,
            default=str,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
