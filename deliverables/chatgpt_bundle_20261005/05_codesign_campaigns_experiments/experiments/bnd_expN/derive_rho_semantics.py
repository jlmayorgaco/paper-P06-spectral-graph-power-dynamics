"""Tabulate the declared ExpN scaling contract from the immutable N01 ledger."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
SOURCE = OUT / "TABLE_N01_original_operating_point.csv"
TARGET = OUT / "TABLE_N02_rho_semantics.csv"


def main() -> None:
    actual = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    expected = (SOURCE.with_suffix(SOURCE.suffix + ".sha256")).read_text().strip()
    if actual != expected:
        raise RuntimeError("N01 changed after freeze")
    with SOURCE.open(newline="", encoding="utf-8") as fh:
        originals = list(csv.DictReader(fh))
    rows = []
    for original in originals:
        bus = int(original["bus"])
        sn = float(original["Sn_original_MVA"])
        h = float(original["H_seconds"])
        p = float(original["P_gen_MW"])
        q = float(original["Q_gen_Mvar"])
        for rho in (0.0, 0.10, 0.50, 0.90, 0.99, 1.0):
            eps = 1.0 - rho
            rows.append(dict(bus=bus, rho=rho, epsilon=eps,
                architecture="SG_ONLY" if rho == 0 else "GFL_ONLY" if rho == 1 else "MIXED",
                SG_Sn_MVA=eps * sn, SG_H_seconds=h if eps > 0 else "ABSENT",
                SG_HSn_MVA_s=eps * sn * h,
                SG_target_P_MW=eps * p, SG_target_Q_Mvar=eps * q,
                GFL_port_scale=rho,
                GFL_internal_filter_current_scale=1.0 if rho > 0 else "ABSENT",
                GFL_per_unit_Rf=0.01 if rho > 0 else "ABSENT",
                GFL_per_unit_Xf=0.03 if rho > 0 else "ABSENT",
                GFL_per_unit_Cdc=1.25 if rho > 0 else "ABSENT",
                GFL_aggregate_physical_capacity_scale=rho,
                GFL_target_P_MW=rho * p, GFL_target_Q_Mvar=rho * q,
                total_target_P_MW=p, total_target_Q_Mvar=q))
    with TARGET.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"N02 rows={len(rows)} output={TARGET}")


if __name__ == "__main__":
    main()
