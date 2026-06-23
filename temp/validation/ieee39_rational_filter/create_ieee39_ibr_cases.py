"""Create and test IEEE 39-bus ANDES cases with SG -> IBR replacements.

The packaged ANDES ``ieee39_full.xlsx`` case contains synchronous GENROU
machines. This script creates derived XLSX cases where selected static
generators keep their power-flow injections but their dynamic GENROU, governor,
exciter, and PSS devices are disabled and replaced by renewable models:

* GFL: REGCP1 + REECA1 + REPCA1 + PLL1
* GFM: REGF1

The intent is benchmark construction for the rational graph-filter damping
paper. The script is deliberately conservative: it keeps the slack generator
as a synchronous machine and starts from low-penetration replacement cases.
"""

from __future__ import annotations

import argparse
import contextlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


GFL_TEMPLATE_CASE = "ieee14/ieee14_regcp1.xlsx"
BASE_CASE = "ieee39/ieee39_full.xlsx"


@dataclass(frozen=True)
class CaseSpec:
    name: str
    gfl: tuple[int, ...]
    gfm: tuple[int, ...]


def default_case_specs() -> list[CaseSpec]:
    """Conservative conversion plan.

    Generator 10 is the slack machine in the packaged case and is retained as
    SG for all default cases. The ordering starts at generators electrically
    near the eastern side of the New England benchmark, but the exact order is
    not claimed to be optimal.
    """

    return [
        CaseSpec("ieee39_ibr_gfl20", gfl=(9, 8), gfm=()),
        CaseSpec("ieee39_ibr_gfm20", gfl=(), gfm=(9, 8)),
        CaseSpec("ieee39_ibr_mix20", gfl=(9,), gfm=(8,)),
        CaseSpec("ieee39_ibr_mix40", gfl=(9, 7), gfm=(8, 6)),
        CaseSpec("ieee39_ibr_mix60", gfl=(9, 7, 5), gfm=(8, 6, 4)),
    ]


def read_xlsx(case_path: str | Path) -> dict[str, pd.DataFrame]:
    return pd.read_excel(case_path, sheet_name=None)


def next_uid(df: pd.DataFrame) -> int:
    if df.empty or "uid" not in df:
        return 0
    return int(pd.to_numeric(df["uid"], errors="coerce").max()) + 1


def find_static_gen_row(sheets: dict[str, pd.DataFrame], gen_idx: int) -> tuple[pd.Series, str]:
    for sheet_name in ("PV", "Slack"):
        df = sheets.get(sheet_name)
        if df is None:
            continue
        match = df[df["idx"].astype(int) == int(gen_idx)]
        if not match.empty:
            return match.iloc[0].copy(), sheet_name
    raise KeyError(f"static generator {gen_idx} not found in PV/Slack")


def generator_bus_map(sheets: dict[str, pd.DataFrame]) -> dict[int, int]:
    out: dict[int, int] = {}
    for sheet_name in ("PV", "Slack"):
        df = sheets.get(sheet_name)
        if df is None:
            continue
        for _, row in df.iterrows():
            out[int(row["idx"])] = int(row["bus"])
    return out


def busfreq_map(sheets: dict[str, pd.DataFrame]) -> dict[int, str]:
    df = sheets["BusFreq"]
    return {int(row["bus"]): str(row["idx"]) for _, row in df.iterrows()}


def line_for_bus(sheets: dict[str, pd.DataFrame], bus: int) -> str:
    line = sheets["Line"]
    match = line[(line["bus1"].astype(int) == bus) | (line["bus2"].astype(int) == bus)]
    if match.empty:
        return str(line.iloc[0]["idx"])
    # Prefer the transformer branch attached to the generator terminal bus.
    with contextlib.suppress(Exception):
        trans = match[match["trans"].astype(int) == 1]
        if not trans.empty:
            return str(trans.iloc[0]["idx"])
    return str(match.iloc[0]["idx"])


def disable_synchronous_dynamics(sheets: dict[str, pd.DataFrame], converted_gens: Iterable[int]) -> None:
    """Remove dynamic SG devices superseded by renewable models.

    In the ANDES IEEE39 case, setting ``u=0`` is not sufficient for all
    inherited controllers: disabled TGOV1N rows can still produce initialization
    residuals. For derived benchmark cases we therefore remove the dynamic rows
    attached to converted machines while keeping the static PV/Slack injection.
    """

    converted = {int(g) for g in converted_gens}
    if not converted:
        return

    syn_ids: set[str] = set()
    genrou = sheets.get("GENROU")
    if genrou is not None:
        mask = genrou["gen"].astype(int).isin(converted)
        syn_ids = {str(v) for v in genrou.loc[mask, "idx"]}
        sheets["GENROU"] = genrou.loc[~mask].reset_index(drop=True)

    for sheet_name in ("TGOV1N", "TGOV1"):
        df = sheets.get(sheet_name)
        if df is not None and "syn" in df:
            sheets[sheet_name] = df.loc[~df["syn"].astype(str).isin(syn_ids)].reset_index(drop=True)

    disabled_avrs: set[str] = set()
    ieeex1 = sheets.get("IEEEX1")
    if ieeex1 is not None and "syn" in ieeex1:
        mask = ieeex1["syn"].astype(str).isin(syn_ids)
        disabled_avrs = {str(v) for v in ieeex1.loc[mask, "idx"]}
        sheets["IEEEX1"] = ieeex1.loc[~mask].reset_index(drop=True)

    ieeest = sheets.get("IEEEST")
    if ieeest is not None and "avr" in ieeest:
        sheets["IEEEST"] = ieeest.loc[~ieeest["avr"].astype(str).isin(disabled_avrs)].reset_index(drop=True)

    toggler = sheets.get("Toggler")
    if toggler is not None and {"model", "dev", "u"}.issubset(toggler.columns):
        mask = (toggler["model"].astype(str).str.upper() == "GENROU") & toggler["dev"].astype(str).isin(syn_ids)
        sheets["Toggler"] = toggler.loc[~mask].reset_index(drop=True)


def append_gfl_models(sheets: dict[str, pd.DataFrame], gen_indices: Iterable[int]) -> None:
    import andes  # type: ignore

    template = read_xlsx(andes.get_case(GFL_TEMPLATE_CASE))
    for name in ("REGCP1", "REECA1", "REPCA1", "PLL1"):
        if name not in sheets:
            sheets[name] = pd.DataFrame(columns=template[name].columns)

    regcp_template = template["REGCP1"].iloc[0]
    reeca_template = template["REECA1"].iloc[0]
    repca_template = template["REPCA1"].iloc[0]
    pll_template = template["PLL1"].iloc[0]

    bf = busfreq_map(sheets)
    rows_regcp: list[pd.Series] = []
    rows_reeca: list[pd.Series] = []
    rows_repca: list[pd.Series] = []
    rows_pll: list[pd.Series] = []

    reg_uid = next_uid(sheets["REGCP1"])
    ree_uid = next_uid(sheets["REECA1"])
    rep_uid = next_uid(sheets["REPCA1"])
    pll_uid = next_uid(sheets["PLL1"])

    for local_idx, gen_idx in enumerate(gen_indices, start=1):
        static, _ = find_static_gen_row(sheets, int(gen_idx))
        bus = int(static["bus"])
        rid = int(local_idx)
        # Avoid idx collisions if an input file already has renewable models.
        while rid in set(pd.to_numeric(sheets["REGCP1"].get("idx", pd.Series(dtype=float)), errors="coerce").dropna().astype(int)):
            rid += 1
        pll_idx = f"PLL1_G{int(gen_idx)}"

        pll = pll_template.copy()
        pll["uid"] = pll_uid
        pll["idx"] = pll_idx
        pll["u"] = 1
        pll["name"] = pll_idx
        pll["bus"] = bus
        pll["Kp"] = 1.0
        pll["Ki"] = 0.2
        pll["Tf"] = 0.05
        pll["Tp"] = 0.05
        pll["fn"] = 60
        rows_pll.append(pll)
        pll_uid += 1

        reg = regcp_template.copy()
        reg["uid"] = reg_uid
        reg["idx"] = rid
        reg["u"] = 1
        reg["name"] = f"REGCP1_G{int(gen_idx)}"
        reg["bus"] = bus
        reg["gen"] = int(gen_idx)
        reg["Sn"] = float(static["Sn"])
        reg["Iolim"] = -999.0
        reg["Iqrmax"] = 999.0
        reg["Iqrmin"] = -999.0
        reg["pll"] = pll_idx
        rows_regcp.append(reg)
        reg_uid += 1

        ree = reeca_template.copy()
        ree["uid"] = ree_uid
        ree["idx"] = rid
        ree["u"] = 1
        ree["name"] = f"REECA1_G{int(gen_idx)}"
        ree["reg"] = rid
        ree["Vref0"] = float(static["v0"])
        ree["Vref1"] = float(static["v0"])
        ree["QMax"] = float(static["qmax"])
        ree["QMin"] = float(static["qmin"])
        ree["PMAX"] = float(static["pmax"])
        ree["PMIN"] = float(static["pmin"])
        ree["Imax"] = 10.0
        rows_reeca.append(ree)
        ree_uid += 1

        rep = repca_template.copy()
        rep["uid"] = rep_uid
        rep["idx"] = rid
        rep["u"] = 1
        rep["name"] = f"REPCA1_G{int(gen_idx)}"
        rep["ree"] = rid
        rep["line"] = line_for_bus(sheets, bus)
        rep["busf"] = bf.get(bus, next(iter(bf.values())))
        rep["PLflag"] = 0
        rep["Fflag"] = 0
        rows_repca.append(rep)
        rep_uid += 1

    if rows_regcp:
        sheets["REGCP1"] = pd.concat([sheets["REGCP1"], pd.DataFrame(rows_regcp)], ignore_index=True)
        sheets["REECA1"] = pd.concat([sheets["REECA1"], pd.DataFrame(rows_reeca)], ignore_index=True)
        sheets["REPCA1"] = pd.concat([sheets["REPCA1"], pd.DataFrame(rows_repca)], ignore_index=True)
        sheets["PLL1"] = pd.concat([sheets["PLL1"], pd.DataFrame(rows_pll)], ignore_index=True)


REGF1_COLUMNS = [
    "uid",
    "idx",
    "u",
    "name",
    "bus",
    "gen",
    "Sn",
    "rf",
    "xf",
    "Vdip",
    "Tfrz",
    "PQFLAG",
    "fn",
    "dwmax",
    "dwmin",
    "wdrp",
    "Qdrp",
    "Tr",
    "Te",
    "KPi",
    "KIi",
    "KPv",
    "KIv",
    "Pmax",
    "Pmin",
    "KPplim",
    "KIplim",
    "Qmax",
    "Qmin",
    "KPqlim",
    "KIqlim",
    "Tpm",
    "gammap",
    "gammaq",
]


def append_gfm_models(sheets: dict[str, pd.DataFrame], gen_indices: Iterable[int]) -> None:
    if "REGF1" not in sheets:
        sheets["REGF1"] = pd.DataFrame(columns=REGF1_COLUMNS)

    rows: list[dict[str, object]] = []
    uid = next_uid(sheets["REGF1"])
    for gen_idx in gen_indices:
        static, _ = find_static_gen_row(sheets, int(gen_idx))
        rows.append(
            {
                "uid": uid,
                "idx": f"REGF1_G{int(gen_idx)}",
                "u": 1,
                "name": f"REGF1_G{int(gen_idx)}",
                "bus": int(static["bus"]),
                "gen": int(gen_idx),
                "Sn": float(static["Sn"]),
                "rf": max(float(static.get("ra", 0.0)), 0.0),
                "xf": max(float(static.get("xs", 0.2)), 0.02),
                "Vdip": 0.8,
                "Tfrz": 0.0,
                "PQFLAG": 0,
                "fn": 60,
                "dwmax": 75.0,
                "dwmin": -75.0,
                "wdrp": 0.033,
                "Qdrp": 0.045,
                "Tr": 0.005,
                "Te": 0.005,
                "KPi": 0.5,
                "KIi": 20.0,
                "KPv": 3.0,
                "KIv": 10.0,
                "Pmax": float(static["pmax"]),
                "Pmin": float(static["pmin"]),
                "KPplim": 5.0,
                "KIplim": 30.0,
                "Qmax": float(static["qmax"]),
                "Qmin": float(static["qmin"]),
                "KPqlim": 0.1,
                "KIqlim": 1.5,
                "Tpm": 0.025,
                "gammap": 1.0,
                "gammaq": 1.0,
            }
        )
        uid += 1

    if rows:
        sheets["REGF1"] = pd.concat([sheets["REGF1"], pd.DataFrame(rows, columns=REGF1_COLUMNS)], ignore_index=True)


def create_case(base_case: str | Path, spec: CaseSpec, out_dir: Path) -> Path:
    sheets = read_xlsx(base_case)
    converted = tuple(spec.gfl) + tuple(spec.gfm)
    overlap = set(spec.gfl) & set(spec.gfm)
    if overlap:
        raise ValueError(f"{spec.name}: generators cannot be both GFL and GFM: {sorted(overlap)}")
    if 10 in converted:
        raise ValueError(f"{spec.name}: default workflow keeps slack generator 10 synchronous")

    disable_synchronous_dynamics(sheets, converted)
    append_gfl_models(sheets, spec.gfl)
    append_gfm_models(sheets, spec.gfm)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{spec.name}.xlsx"
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    return out_path


def damping_ratio(s: complex) -> float:
    if abs(s) < 1e-12:
        return float("inf")
    return float(-s.real / abs(s))


def pole_summary(ss, osc_imag_tol: float = 0.1) -> dict[str, object]:
    mu = np.asarray(ss.EIG.mu, dtype=complex)
    stable_osc = [s for s in mu if abs(s.imag) > osc_imag_tol and s.real < 0]
    if stable_osc:
        ratios = np.array([damping_ratio(s) for s in stable_osc])
        critical = stable_osc[int(np.argmin(ratios))]
    else:
        critical = complex(np.nan, np.nan)
    pos = [s for s in mu if s.real > 1e-7]
    pos_osc = [s for s in pos if abs(s.imag) > osc_imag_tol]
    dominant = max(pos, key=lambda s: s.real) if pos else complex(np.nan, np.nan)
    return {
        "n_eigs": int(len(mu)),
        "n_positive_real": int(len(pos)),
        "n_positive_oscillatory": int(len(pos_osc)),
        "n_zero": int(np.sum(np.abs(mu) < 1e-7)),
        "n_oscillatory": int(np.sum(np.abs(mu.imag) > osc_imag_tol)),
        "max_real": float(np.max(mu.real)) if len(mu) else float("nan"),
        "small_signal_stable": bool(len(pos) == 0),
        "dominant_positive": {
            "real": float(dominant.real),
            "imag": float(dominant.imag),
            "zeta": damping_ratio(dominant),
            "freq_hz": float(abs(dominant.imag) / (2 * np.pi)) if np.isfinite(dominant.imag) else float("nan"),
        },
        "critical_stable_osc": {
            "real": float(critical.real),
            "imag": float(critical.imag),
            "zeta": damping_ratio(critical),
            "freq_hz": float(abs(critical.imag) / (2 * np.pi)) if np.isfinite(critical.imag) else float("nan"),
        },
    }


def inventory(ss) -> dict[str, int]:
    names = [
        "Bus",
        "Line",
        "PV",
        "Slack",
        "GENROU",
        "TGOV1N",
        "IEEEX1",
        "IEEEST",
        "REGCP1",
        "REGCA1",
        "REECA1",
        "REPCA1",
        "REGF1",
        "REGF2",
        "PLL1",
    ]
    return {name: int(getattr(getattr(ss, name), "n", 0)) for name in names if hasattr(ss, name)}


def test_case(case_path: Path) -> dict[str, object]:
    import andes  # type: ignore

    result: dict[str, object] = {"case": str(case_path), "ok": False}
    try:
        ss = andes.load(str(case_path), setup=False, no_output=True)
        ss.setup()
        pflow_ok = bool(ss.PFlow.run())
        eig_ok = bool(ss.EIG.run()) if pflow_ok else False
        result.update(
            {
                "ok": bool(pflow_ok and eig_ok),
                "pflow_ok": pflow_ok,
                "eig_ok": eig_ok,
                "inventory": inventory(ss),
            }
        )
        if eig_ok:
            result["poles"] = pole_summary(ss)
    except Exception as exc:  # ANDES errors carry useful model context.
        result.update({"error_type": type(exc).__name__, "error": str(exc)})
    return result


def flatten_result(rec: dict[str, object]) -> dict[str, object]:
    test = rec.get("test", {})
    if not isinstance(test, dict):
        test = {}
    poles = test.get("poles", {})
    if not isinstance(poles, dict):
        poles = {}
    crit = poles.get("critical_stable_osc", {})
    if not isinstance(crit, dict):
        crit = {}
    dom = poles.get("dominant_positive", {})
    if not isinstance(dom, dict):
        dom = {}
    inv = test.get("inventory", {})
    if not isinstance(inv, dict):
        inv = {}
    return {
        "name": rec.get("name"),
        "gfl": " ".join(str(x) for x in rec.get("gfl", [])),
        "gfm": " ".join(str(x) for x in rec.get("gfm", [])),
        "ok": test.get("ok"),
        "small_signal_stable": poles.get("small_signal_stable"),
        "max_real": poles.get("max_real"),
        "n_positive_real": poles.get("n_positive_real"),
        "n_positive_oscillatory": poles.get("n_positive_oscillatory"),
        "dominant_positive_real": dom.get("real"),
        "dominant_positive_imag": dom.get("imag"),
        "critical_stable_osc_real": crit.get("real"),
        "critical_stable_osc_imag": crit.get("imag"),
        "critical_stable_osc_zeta": crit.get("zeta"),
        "critical_stable_osc_freq_hz": crit.get("freq_hz"),
        "GENROU": inv.get("GENROU"),
        "REGCP1": inv.get("REGCP1"),
        "REGF1": inv.get("REGF1"),
        "PLL1": inv.get("PLL1"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-cases", type=Path, default=Path("validation") / "ieee39_rational_filter" / "cases")
    parser.add_argument("--out-results", type=Path, default=Path("outputs") / "ieee39_ibr_cases")
    parser.add_argument("--no-test", action="store_true", help="Only create XLSX cases")
    return parser.parse_args()


def main() -> None:
    import andes  # type: ignore

    args = parse_args()
    base = andes.get_case(BASE_CASE)
    args.out_cases.mkdir(parents=True, exist_ok=True)
    args.out_results.mkdir(parents=True, exist_ok=True)

    results = []
    for spec in default_case_specs():
        case_path = create_case(base, spec, args.out_cases)
        rec: dict[str, object] = {
            "name": spec.name,
            "path": str(case_path),
            "gfl": list(spec.gfl),
            "gfm": list(spec.gfm),
        }
        if not args.no_test:
            rec["test"] = test_case(case_path)
        results.append(rec)

    out_json = args.out_results / "case_generation_summary.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    pd.DataFrame([flatten_result(r) for r in results]).to_csv(args.out_results / "case_generation_summary.csv", index=False)
    print(json.dumps(results, indent=2))
    print(f"Wrote {out_json}")


if __name__ == "__main__":
    main()
