# ruff: noqa: E501  -- criteria tables kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - SQ1 stopped at SQ1-2 (blind synthetic qualification FAIL); the ParaEMT line is stopped permanently (prereg SQ1).
"""SQ1-3: the unused phasor-TDS holdout T1-T9 (prereg SQ1 section 4; run only if SQ1-2 passes).

The frozen V3 tooling is used unchanged; only its output directories are redirected so that the V3
record (results/EMTV3) is not touched:
  EMTSQ1_tds.py gen   (.venv/tx3-analysis)  -> EMTV3_tds_traces.main('holdout'), traces to external/paremt_runs/sq1_tds
  EMTSQ1_tds.py eval  (.venv/xtool-paremt)  -> EMTV3_E0_tds.main('holdout') with the frozen V3 estimator, then the SQ1-3
                                                criteria; writes results/EMTSQ1/SQ1-3/sq1_tds.json
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "v3"))

RESEARCH = HERE.parents[2]
REPO = HERE.parents[6]
OUT = RESEARCH / "results" / "EMTSQ1" / "SQ1-3"
FROZEN_SHA = "849472d8804c694eee69701dd251aaa69d506878007a3a4fc88fde1fb3b0817f"


def main() -> int:
    mode = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    if mode == "gen":
        import EMTV3_tds_traces as T

        T.OUT = OUT
        T.BIG = REPO / "external" / "paremt_runs" / "sq1_tds"
        sys.argv = [sys.argv[0], "holdout"]
        return T.main()
    assert hashlib.sha256((HERE.parent / "v3" / "emtv3_estimator.py").read_bytes()).hexdigest() == FROZEN_SHA
    import EMTV3_E0_tds as E
    import pandas as pd

    E.OUT = OUT
    sys.argv = [sys.argv[0], "holdout"]
    E.main()
    df = pd.read_csv(OUT / "E0_tds_holdout.csv")
    big = df[df.alpha_eig.abs() >= 0.01]
    wrong = ((df.verdict == "STABLE") & (df.alpha_eig > 0)) | ((df.verdict == "UNSTABLE") & (df.alpha_eig < 0))
    res = df[df.resolved == True]  # noqa: E712
    da, dfq = (res.alpha_hat - res.alpha_eig).abs(), (res.f_hat - res.f_eig).abs()
    crit = {"no_crash": bool(not df.crash.any()), "no_wrong_sign": bool(not wrong.any()),
            "all_abs_alpha_ge_0.01_resolved": bool(big.resolved.fillna(False).all()), "unresolved_cases": df[df.resolved != True].case.tolist(),  # noqa: E712
            "max_abs_dalpha_resolved": float(da.max()) if len(res) else None, "max_abs_df_resolved": float(dfq.max()) if len(res) else None,
            "errors_within_tolerance": bool((da <= 0.005).all() and (dfq <= 0.005).all()), "no_unexpected_mode": bool(not (df.verdict == "UNEXPECTED_MODE").any())}
    crit["pass"] = bool(all(crit[k] for k in ("no_crash", "no_wrong_sign", "all_abs_alpha_ge_0.01_resolved", "errors_within_tolerance", "no_unexpected_mode")))
    (OUT / "sq1_tds.json").write_text(json.dumps(crit, indent=1, default=float), encoding="utf-8")
    print(json.dumps(crit, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
