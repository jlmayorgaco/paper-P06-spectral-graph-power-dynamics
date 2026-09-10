"""F8 key-intervention table: old/new H, old/new kappa, d alpha_IA, d m_cl."""

from __future__ import annotations

import pandas as pd

from _bootstrap import RESULTS

F8 = RESULTS / "F8"
KEY = {
    "R_none": "plain replacement (reference)",
    "R_0000000": "EM only: 25 % condenser, 1 % inertia, frozen EMF, field held, "
    "no PSS/D/Q",
    "R_0100000": "EM + inertia (100 % kinetic energy)",
    "R_0001000": "EM + electromagnetic transients",
    "R_0001100": "EM + transients + AVR",
    "R_0001110": "EM + transients + AVR + PSS",
    "R_0000001": "EM + reactive support",
    "R_0101111": "25 % condenser, all services, D = 0",
    "R_1111111": "100 % condenser, all services, D = 2",
    "R_virtual_inertia": "synthetic inertia on converters (no EM)",
    "A_10111": "surviving fleet: native (reference for class A)",
    "A_00111": "surviving fleet: inertia x0.5",
    "A_20111": "surviving fleet: inertia x2",
    "A_11111": "surviving fleet: D = 2",
    "A_10110": "surviving fleet: PSS off",
    "A_10101": "surviving fleet: AVR manual (field held)",
    "A_10011": "surviving fleet: classical (frozen EMF), D = 0",
    "A_11011": "surviving fleet: classical (frozen EMF), D = 2",
}


def main() -> int:
    t = pd.read_csv(F8 / "F8_intervention_table.csv")
    t = t[t.config.isin(KEY)].copy()
    t["intervention"] = t.config.map(KEY)
    inf = lambda k: "inf" if k == -1 else ("-" if pd.isna(k) else str(int(k)))  # noqa: E731
    lines = []
    for point, g in t.groupby("point", sort=False):
        lines.append(f"\n**{point}**\n")
        lines.append(
            "| intervention | old H | new H | old kappa | new kappa "
            "| d alpha_IA | d m_cl |"
        )
        lines.append("|---|---|---|---|---|---|---|")
        for name in KEY:
            r = g[g.config == name]
            if r.empty:
                continue
            r = r.iloc[0]
            new = r.label if r.status == "OK" else r.status
            lines.append(
                f"| {KEY[name]} | `{r.H_old}` | `{new}` | {inf(r.kappa_old)} "
                f"| {inf(r.kappa)} | "
                f"{r.d_alpha_IA:+.3f} | {r.d_m_cl:+.3f} |"
                if r.status == "OK"
                else f"| {KEY[name]} | `{r.H_old}` | {new} "
                f"| {inf(r.kappa_old)} | - | - | - |"
            )
    (F8 / "F8_key_table.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
