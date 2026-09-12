# ruff: noqa: E501  -- ledger table kept on one line per claim
"""ParaEMT closure (2026-09-12): add the EMT evidence column to the FINAL validation matrix.

Appends one column, `emt_validation`, to results/20260911_FINAL_VALIDATION_MATRIX.csv. Every
existing column (mathematical, synthetic, IEEE-39 / phasor-TDS, nonlinear, ANDES / independent
tool, robustness, baselines, governors, status, wording, evidence) is asserted unchanged.

Only two EMT categories apply to evaluated claims:
  - EMT network / equilibrium: REPRODUCED (ParaEMT V1 gate G1);
  - everything else: device gates failed (UNRESOLVED) or portfolio experiment never run
    (NEVER EVALUATED). Mathematical and refuted claims are NOT APPLICABLE.
Source of the EMT record: docs/20260912_PAREMT_EMT_CLOSURE.md (V1, V2, V3, SQ1).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

RESEARCH = Path(__file__).resolve().parents[3]
MATRIX = RESEARCH / "results" / "20260911_FINAL_VALIDATION_MATRIX.csv"
COL = "emt_validation"

NA_MATH = "NOT APPLICABLE (mathematical result)"
NEVER = "NEVER EVALUATED - EMT UNRESOLVED (portfolio-level EMT experiment never run; ParaEMT line closed at qualification gates before EMT04)"
EMT = {
    "V01": NA_MATH,
    "V02": NA_MATH,
    "V03": NEVER,
    "V04": NEVER,
    "V05": "NOT APPLICABLE (baseline recount on phasor spectra)",
    "V06": NEVER,
    "V07": NEVER,
    "V08": NA_MATH,
    "V09": "NOT APPLICABLE (connected-cumulant diagnostic; no EMT observable)",
    "V10": NEVER,
    "V11": NEVER,
    "V12": "NOT APPLICABLE (negative scope limit of the linear prediction)",
    "V13": NEVER,
    "V14": NEVER,
    "V15": NEVER,
    "V16": NEVER,
    "V17": NEVER,
    "V18": NEVER,
    "V19": "EMT REPRODUCED - network and operating point (ParaEMT V1 gate G1 PASS); the only EMT-reproduced claim",
    "V20": "EMT UNRESOLVED - SG device gate failed (V1 G3 FAIL; V2 G3a FAIL); not reached in V3 or SQ1",
    "V21": "EMT UNRESOLVED - GFL device gate failed (V1 G4 FAIL); V2 voltage interface never run in a network (V2 stopped at G3a; V3 and SQ1 stopped at instrument gates)",
    "V22": NA_MATH,
    "V23": "NOT APPLICABLE (phasor-port diagnostic)",
    "V24": "NOT APPLICABLE (refuted claim; EMT not used to reinterpret failed claims)",
    "V25": "NOT APPLICABLE (phasor-port closure diagnostic)",
    "V26": NA_MATH,
    "V27": "NOT APPLICABLE (refuted claim; EMT not used to reinterpret failed claims)",
    "V28": "NOT APPLICABLE (refuted claim; transfer across converter models not testable with the same GFL)",
    "V29": NEVER,
    "V30": NEVER,
}


def main() -> int:
    with MATRIX.open(newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], rows[1:]
    if COL in head:
        print("column already present; nothing to do")
        return 0
    ids = [r[0] for r in body]
    assert ids == list(EMT), "claim ids differ from the FINAL matrix"
    out = [head + [COL]] + [r + [EMT[r[0]]] for r in body]
    with MATRIX.open("w", newline="", encoding="utf-8") as f:
        csv.writer(f, lineterminator="\n").writerows(out)
    with MATRIX.open(newline="", encoding="utf-8") as f:
        back = list(csv.reader(f))
    assert [r[:-1] for r in back] == rows, "an existing column changed"
    counts: dict[str, int] = {}
    for r in back[1:]:
        key = r[-1].split(" (")[0].split(" - ")[0]
        counts[key] = counts.get(key, 0) + 1
    print(counts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
