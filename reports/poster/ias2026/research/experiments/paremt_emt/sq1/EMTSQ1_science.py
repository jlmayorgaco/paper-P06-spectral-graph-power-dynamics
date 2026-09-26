# ruff: noqa: E501  -- result tables kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - SQ1 stopped at SQ1-2 (blind synthetic qualification FAIL); the ParaEMT line is stopped permanently (prereg SQ1).
"""SQ1-7/8: scientific core EMT04 -> EMT05 -> EMT06 -> EMT07/08 -> EMT10, then the decision gate
(prereg SQ1 sections 7-8). Requires SQ1-2, SQ1-3, G3a, G3b, G4a, G4b PASS (results/EMTSQ1/EMTSQ1_gates.json).

Uses the block functions of experiments/paremt_emt/v3/EMTV3_science.py unchanged (prepared in V3,
never executed there), with the output and raw directories redirected to SQ1, the frozen V3
estimator (hash asserted), the v1 IEEE-39 realization plus the V2 GFL voltage interface.
Usage: EMTSQ1_science.py core | verify
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "v3"))
sys.path.insert(0, str(HERE.parent))

import _emt  # noqa: E402
import EMTV3_science as S  # noqa: E402

FROZEN_SHA = "849472d8804c694eee69701dd251aaa69d506878007a3a4fc88fde1fb3b0817f"
OUT = _emt.RESULTS / "EMTSQ1"
S.OUT = OUT
S.RAW = _emt.REPO / "external" / "paremt_runs" / "sq1_sci"


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "core"
    assert hashlib.sha256((HERE.parent / "v3" / "emtv3_estimator.py").read_bytes()).hexdigest() == FROZEN_SHA, "frozen V3 estimator changed"
    gates = json.loads((OUT / "EMTSQ1_gates.json").read_text())
    assert all(gates.get(f"GATE_{g}") == "PASS" for g in ("G3a", "G3b", "G4a", "G4b")), "all SQ1 device gates must pass first"
    if mode == "verify":
        d = OUT / "verify"
        d.mkdir(parents=True, exist_ok=True)
        df05, s05, amp = S.emt05(None)
        df05.to_csv(d / "p4_16_subsets.csv", index=False)
        df08, s08 = S.emt08(df05)
        df08[df08.label == "grid"].to_csv(d / "g_grid.csv", index=False)
        same = {f: (d / f).read_bytes() == (OUT / ref / f).read_bytes() for f, ref in (("p4_16_subsets.csv", "EMT05"), ("g_grid.csv", "EMT08"))}
        (OUT / "EMTSQ1_science_rerun_identity.json").write_text(json.dumps(same, indent=1), encoding="utf-8")
        print(same)
        return 0
    summary = {"provenance": _emt.provenance()}
    t0 = time.time()
    df04, g5, hold, noev = S.emt04()
    d = OUT / "EMT04"
    d.mkdir(parents=True, exist_ok=True)
    df04.to_csv(d / "emt04_runs.csv", index=False)
    g5.to_csv(d / "emt04_g5.csv", index=False)
    hold.to_csv(d / "emt04_holdouts.csv", index=False)
    noev.to_csv(d / "emt04_noevent.csv", index=False)
    summary["EMT04"] = {"G5_pass": bool(g5["pass"].all()), "G2_tau_m_pass": bool(hold[hold.holdout == "tau_m_0.5ms"]["pass"].all()),
                        "native_stator_descriptive": hold[hold.holdout == "native_stator"].to_dict("records")}
    decision = {"EMT04_G5": summary["EMT04"]["G5_pass"]}
    if not summary["EMT04"]["G5_pass"]:
        summary["STOPPED"] = "EMT04 G5 failed (prereg SQ1 section 7: STOP, no amendment)"
    else:
        df05, s05, amp = S.emt05(df04)
        d = OUT / "EMT05"
        d.mkdir(parents=True, exist_ok=True)
        df05.to_csv(d / "p4_16_subsets.csv", index=False)
        amp.to_csv(d / "amplitude_checks.csv", index=False)
        s05["amplitude_checks_pass"] = bool(amp["pass"].all())
        summary["EMT05"] = s05
        eff, s06 = S.emt06(df05)
        d = OUT / "EMT06"
        d.mkdir(parents=True, exist_ok=True)
        eff.to_csv(d / "contextual_effects.csv", index=False)
        summary["EMT06"] = s06
        df08, s08 = S.emt08(df05)
        d = OUT / "EMT08"
        d.mkdir(parents=True, exist_ok=True)
        df08.to_csv(d / "g_boundary.csv", index=False)
        df08[df08.label == "grid"].to_csv(d / "g_grid.csv", index=False)
        summary["EMT08"] = s08
        df10, s10 = S.emt10(df05, df08, s08.get("g_star_emt"))
        d = OUT / "EMT10"
        d.mkdir(parents=True, exist_ok=True)
        df10.to_csv(d / "newton_sequence.csv", index=False)
        summary["EMT10"] = s10
        decision.update({"EMT05_primary": s05["primary_pass"], "EMT06": s06["status"] == "PASS", "EMT08_primary": s08["primary_pass"],
                         "EMT08_target": s08["target_pass"], "EMT10": s10["pass"]})
        summary["H_EMT_P4_equals_H4"] = s05["primary_pass"]
        summary["g_star_emt"] = s08.get("g_star_emt")
    decision["corroborated"] = bool(all(decision.values()))
    summary["SQ1-8_decision"] = decision
    summary["NEXT"] = "proceed to EMT09, EMT11-EMT16 (prereg SQ1 section 8)" if decision["corroborated"] else "STOP (prereg SQ1 section 8): the result is the disagreement / non-resolution"
    summary["wall_s"] = round(time.time() - t0, 1)
    (OUT / "EMTSQ1_core.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    with (OUT / "EMTSQ1_core_manifests.jsonl").open("w", encoding="utf-8") as f:
        for m in S.MANIFESTS:
            f.write(json.dumps(m, default=float) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "provenance"}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
