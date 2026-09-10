"""F8 analysis: which services move H_Gamma, which only move margins, and synergy.

Reads results/F8/F8_interventions.csv. Every statement is about pairs of
configurations that differ in exactly one factor, so nothing is summarized into
a single importance score.
"""

from __future__ import annotations

import itertools
import json

import numpy as np
import pandas as pd

from _bootstrap import RESULTS
from F7_report import witnesses
from F8_service_attribution import A_FACTORS, R_FACTORS

F8 = RESULTS / "F8"
TOL_ALPHA = 1e-3


def _code(row, factors):
    return tuple(list(factors[f]).index(row[f]) for f in factors)


def flips(frame: pd.DataFrame, factors: dict) -> pd.DataFrame:
    """Every pair of configurations differing in exactly one factor, per point."""

    rows = []
    for point, group in frame.groupby("point"):
        table = {_code(r, factors): r for _, r in group.iterrows()}
        for code, row in table.items():
            for i, f in enumerate(factors):
                levels = factors[f]
                if code[i] + 1 >= len(levels):
                    continue
                other = code[:i] + (code[i] + 1,) + code[i + 1 :]
                if other not in table:
                    continue
                hi = table[other]
                ok = row["status"] == "OK" and hi["status"] == "OK"
                rows.append(
                    {
                        "point": point,
                        "factor": f,
                        "from_level": levels[code[i]],
                        "to_level": levels[code[i] + 1],
                        "context": "".join(map(str, code)),
                        "status_from": row["status"],
                        "status_to": hi["status"],
                        "H_from": row["label"],
                        "H_to": hi["label"],
                        "kappa_from": row.get("kappa"),
                        "kappa_to": hi.get("kappa"),
                        "H_changes": bool(row["label"] != hi["label"]),
                        "kappa_changes": bool(
                            ok and row.get("kappa") != hi.get("kappa")
                        ),
                        "witness_changes": bool(
                            ok
                            and row["label"] != hi["label"]
                            and row["label"] != "EMPTY"
                            and hi["label"] != "EMPTY"
                            and witnesses(row["label"]) != witnesses(hi["label"])
                        ),
                        "d_alpha_IA": hi.get("alpha_IA", np.nan)
                        - row.get("alpha_IA", np.nan),
                        "d_band_alpha": hi.get("band_alpha", np.nan)
                        - row.get("band_alpha", np.nan),
                        "d_m_cl": hi.get("m_cl", np.nan) - row.get("m_cl", np.nan),
                    }
                )
    return pd.DataFrame(rows)


def flip_summary(fl: pd.DataFrame) -> pd.DataFrame:
    ok = fl[(fl.status_from == "OK") & (fl.status_to == "OK")]
    g = ok.groupby(["point", "factor"])
    out = g.agg(
        pairs=("H_changes", "size"),
        H_changes=("H_changes", "sum"),
        kappa_changes=("kappa_changes", "sum"),
        witness_changes=("witness_changes", "sum"),
        margin_only=("d_band_alpha", lambda s: 0),
        median_d_alpha=("d_band_alpha", "median"),
        median_d_m_cl=("d_m_cl", "median"),
    ).reset_index()
    margin_only = (
        ok[~ok.H_changes & (ok.d_band_alpha.abs() > TOL_ALPHA)]
        .groupby(["point", "factor"])
        .size()
    )
    out["margin_only"] = [
        int(margin_only.get((p, f), 0))
        for p, f in zip(out.point, out.factor, strict=True)
    ]
    status_flip = fl.groupby(["point", "factor"]).apply(
        lambda d: int((d.status_from != d.status_to).sum()), include_groups=False
    )
    out["feasibility_or_base_changes"] = [
        int(status_flip.get((p, f), 0))
        for p, f in zip(out.point, out.factor, strict=True)
    ]
    return out


def synergy(frame: pd.DataFrame, factors: dict) -> pd.DataFrame:
    """Pairs of factors whose joint flip changes H while neither single flip does."""

    names = list(factors)
    rows = []
    for point, group in frame.groupby("point"):
        table = {_code(r, factors): r for _, r in group.iterrows()}
        for code, row in table.items():
            if row["status"] != "OK":
                continue
            for i, j in itertools.combinations(range(len(names)), 2):
                if code[i] + 1 >= len(factors[names[i]]) or code[j] + 1 >= len(
                    factors[names[j]]
                ):
                    continue
                ci = code[:i] + (code[i] + 1,) + code[i + 1 :]
                cj = code[:j] + (code[j] + 1,) + code[j + 1 :]
                cij = list(ci)
                cij[j] += 1
                cij = tuple(cij)
                if not all(
                    c in table and table[c]["status"] == "OK" for c in (ci, cj, cij)
                ):
                    continue
                h0, hi, hj, hij = (
                    row["label"],
                    table[ci]["label"],
                    table[cj]["label"],
                    table[cij]["label"],
                )
                if h0 == hi == hj and hij != h0:
                    rows.append(
                        {
                            "point": point,
                            "pair": f"{names[i]} x {names[j]}",
                            "context": "".join(map(str, code)),
                            "H_base": h0,
                            "H_joint": hij,
                        }
                    )
    return pd.DataFrame(rows)


def factorial_effects(frame: pd.DataFrame, factors: dict, outcome: str) -> pd.DataFrame:
    """Main effects and two-factor interactions of a 2^k design (+-1 coding)."""

    two_level = {f: v for f, v in factors.items() if len(v) == 2}
    rows = []
    for point, group in frame.groupby("point"):
        g = group[group.status == "OK"].dropna(subset=[outcome])
        if len(g) < 2 ** len(two_level) * 0.9:
            continue
        x = np.column_stack(
            [np.where(g[f] == v[1], 1.0, -1.0) for f, v in two_level.items()]
        )
        cols, names = [np.ones(len(g))], ["mean"]
        for i, f in enumerate(two_level):
            cols.append(x[:, i])
            names.append(f)
        for i, j in itertools.combinations(range(len(two_level)), 2):
            cols.append(x[:, i] * x[:, j])
            names.append(f"{list(two_level)[i]} x {list(two_level)[j]}")
        design = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(design, g[outcome].to_numpy(), rcond=None)
        resid = g[outcome].to_numpy() - design @ coef
        for n, c in zip(names, coef, strict=True):
            rows.append(
                {
                    "point": point,
                    "term": n,
                    "effect": 2 * c if n != "mean" else c,
                    "outcome": outcome,
                    "residual_rms": float(np.sqrt(np.mean(resid**2))),
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    frame = pd.read_csv(F8 / "F8_interventions.csv", low_memory=False)
    points = json.loads((F8 / "F8_points.json").read_text(encoding="utf-8"))
    r = frame[(frame["class"] == "R") & frame.config.str.match(r"R_\d{7}$")]
    a = frame[frame["class"] == "A"]
    base = frame[frame.config == "R_none"].set_index("point")

    # intervention table against the plain replacement at the same point
    table = frame[
        [
            "point",
            "config",
            "class",
            "status",
            "label",
            "kappa",
            "alpha_IA",
            "band_alpha",
            "f_IA_hz",
            "m_cl",
            "f_port_hz",
        ]
    ].copy()
    table["H_old"] = table.point.map(base.label)
    table["kappa_old"] = table.point.map(base.kappa)
    table["d_alpha_IA"] = table.alpha_IA - table.point.map(base.alpha_IA)
    table["d_m_cl"] = table.m_cl - table.point.map(base.m_cl)
    table.to_csv(F8 / "F8_intervention_table.csv", index=False)

    fr, fa = flips(r, R_FACTORS), flips(a, A_FACTORS)
    fr.to_csv(F8 / "F8_flips_R.csv", index=False)
    fa.to_csv(F8 / "F8_flips_A.csv", index=False)
    sr, sa = flip_summary(fr), flip_summary(fa)
    sr.to_csv(F8 / "F8_flip_summary_R.csv", index=False)
    sa.to_csv(F8 / "F8_flip_summary_A.csv", index=False)
    syn_r, syn_a = synergy(r, R_FACTORS), synergy(a, A_FACTORS)
    syn_r.to_csv(F8 / "F8_synergy_R.csv", index=False)
    syn_a.to_csv(F8 / "F8_synergy_A.csv", index=False)
    eff = pd.concat(
        [factorial_effects(r, R_FACTORS, o) for o in ("band_alpha", "m_cl")]
        + [factorial_effects(a, A_FACTORS, o) for o in ("band_alpha", "m_cl")]
    )
    eff.to_csv(F8 / "F8_factorial_effects.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    print(json.dumps(points, indent=0)[:1500])
    print(base[["label", "kappa", "alpha_IA", "band_alpha", "m_cl"]])
    for name, fr_ in (("R", r), ("A", a)):
        print(f"--- class {name}: labels per point ---")
        print(fr_.groupby("point").label.value_counts().to_string())
    print(sr.to_string(index=False))
    print(sa.to_string(index=False))
    print(
        "synergy R",
        len(syn_r),
        syn_r.groupby(["point", "pair"]).size().to_string() if len(syn_r) else "",
    )
    print(
        "synergy A",
        len(syn_a),
        syn_a.groupby(["point", "pair"]).size().to_string() if len(syn_a) else "",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
