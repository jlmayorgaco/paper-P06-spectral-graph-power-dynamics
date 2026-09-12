# ruff: noqa: E501  -- criteria tables kept on one line
"""Gate E0 phasor-TDS instrument holdout evaluation (prereg V3 section 4.2; criteria B1-B5).

Reads the traces of EMTV3_tds_traces.py (manifest results/EMTV3/E0/tds_manifest_<mode>.json), runs the
V3 targeted tracker + sentinel with the adaptive duration (records truncated at 30/60/90 s or used as
available), and compares with the independently computed linear eigenvalues (band_re, band_hz).
Usage: EMTV3_E0_tds.py dev | holdout. Writes results/EMTV3/E0/E0_tds_<mode>.{csv,json}.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import emtv3_estimator as V  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

RESEARCH = HERE.parents[2]
REPO = HERE.parents[6]
OUT = RESEARCH / "results" / "EMTV3" / "E0"


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "dev"
    man = json.loads((OUT / f"tds_manifest_{mode}.json").read_text())
    pred = pd.read_csv(RESEARCH / "results" / "EMT_PRED" / "phasor_predictions.csv").set_index("case_id")
    rows = []
    for cid, m in man.items():
        d = np.load(REPO / m["file"])
        t, y = d["t"], d["y"]
        a_eig, f_eig = float(pred.loc[m["pred_case_id"], "band_re"]), float(pred.loc[m["pred_case_id"], "band_hz"])
        row = {"case": cid, "pred_case_id": m["pred_case_id"], "alpha_eig": a_eig, "f_eig": f_eig, "run_status": m["run_status"], "record_end_s": m["record_end_s"]}
        try:
            def rec(T, t=t, y=y):
                n = int(round(T / (t[1] - t[0])))
                return t[:n], y[:, :n]

            r = V.classify_adaptive(rec, f_eig)
            row.update(crash=False, verdict=r["verdict"], resolved=bool(r["resolved"]), alpha_hat=r["alpha"], f_hat=r["freq"], ci_lo=r["ci_alpha"][0],
                       ci_hi=r["ci_alpha"][1], duration=r["duration"], reasons=";".join(r["reasons"]), resid_median=r["resid_median"],
                       sentinel_flag=bool(r["sentinel"].get("flag", False)), sentinel_alpha=r["sentinel"].get("alpha"), sentinel_f=r["sentinel"].get("freq"),
                       history=json.dumps(r["history"], default=float))
        except Exception as exc:
            row.update(crash=True, error=f"{type(exc).__name__}: {exc}")
        rows.append(row)
        print(cid, row.get("verdict"), row.get("alpha_hat"), row.get("f_hat"), flush=True)
    df = pd.DataFrame(rows)
    ok_crash = not df.crash.any()
    big = df[df.alpha_eig.abs() >= 0.01]
    correct = ((big.verdict == "STABLE") & (big.alpha_eig < 0)) | ((big.verdict == "UNSTABLE") & (big.alpha_eig > 0))
    near = df[df.alpha_eig.abs() < 0.005]
    b3 = bool(not (((near.verdict == "STABLE") & (near.alpha_eig > 0)) | ((near.verdict == "UNSTABLE") & (near.alpha_eig < 0))).any())
    res = df[df.resolved == True]  # noqa: E712
    da = (res.alpha_hat - res.alpha_eig).abs()
    dfreq = (res.f_hat - res.f_eig).abs()
    crit = {"B1_no_crash": bool(ok_crash), "B2_correct_conclusive_all": bool(correct.all()), "B2_detail": dict(zip(big.case, correct.tolist(), strict=True)),
            "B3_near_boundary_not_wrong_sign": b3, "B4_max_abs_dalpha": float(da.max()) if len(res) else None, "B4_max_abs_df": float(dfreq.max()) if len(res) else None,
            "B4_all_within": bool((da <= 0.005).all() and (dfreq <= 0.005).all()), "B5_no_sentinel_flag": bool(not df.sentinel_flag.fillna(False).any())}
    crit["pass"] = bool(crit["B1_no_crash"] and crit["B2_correct_conclusive_all"] and b3 and crit["B4_all_within"] and crit["B5_no_sentinel_flag"])
    df.to_csv(OUT / f"E0_tds_{mode}.csv", index=False)
    (OUT / f"E0_tds_{mode}.json").write_text(json.dumps({"mode": mode, "criteria": crit}, indent=1, default=float), encoding="utf-8")
    print(json.dumps(crit, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
