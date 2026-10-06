"""Export the corrected, fixed-operating-point analytical PD port closure."""
from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"/"bnd_model_audit_m"))
from pd_exact_closure import export_and_validate  # noqa: E402


def main()->None:
    out=ROOT/"reports"/"experiment_M"
    rows=[export_and_validate(out/"matrices"/case) for case in
          ("all_SG","ExpG_candidate","ExpK_nominal")]
    with (out/"tables"/"TABLE_M21_corrected_closure_identity.csv").open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for row in rows:print(row)


if __name__=="__main__":main()
