# ruff: noqa: E501  -- manifest tables kept on one line
"""Run manifests for the EMT02/EMT03 unit tests and diagnostics (prereg section 9).

Re-executes the three xtool-paremt scripts (unit tests, diagnostics D1-D4, IEEE-39 D5), times each,
and writes results/EMT02/EMT02_03_manifest.json with provenance (git HEAD, ParaEMT upstream SHA,
working-copy diff SHA, environment hash), the per-run specification (dt, duration, disturbance,
seed = none) and the sha256 of every output file. Wall time is per script (the individual runs are
seconds long and are not timed separately). The phasor references (tx3-analysis) are listed by hash.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _emt  # noqa: E402

SCRIPTS = ("EMT02_03_unit_tests.py", "EMT02_03_diagnostics.py", "EMT03_D5_ieee39_chatter.py")
SG_EV = "SG Pm x1.02 on [1.0, 1.2) s"
GFL_EV = {"a": "p_ref x1.01 for t >= 1 s", "b": "|E| x0.99 for t >= 1 s", "c": "E phase +0.02 rad for t >= 1 s"}
UT = "unit-test system: T -- EMT R-L line (R 0, X 0.05) -- INF; source 1.0 behind j1e-4; amendment A1"
QS = "unit-test system with the T--INF line as a synchronous-frame algebraic Norton (diagnostic)"


def runs():
    out = [{"run_id": "EMT02_SG_dt50us", "script": SCRIPTS[0], "system": UT, "dt_s": 50e-6, "duration_s": 10.0, "disturbance": SG_EV, "role": "PREREGISTERED G3"}]
    for g in (0.0, 0.03625, 0.25):
        for c in "abc":
            out.append({"run_id": f"EMT03_GFL_g{g:g}_{c}_dt50us", "script": SCRIPTS[0], "system": UT, "dt_s": 50e-6, "duration_s": 5.0, "disturbance": GFL_EV[c], "role": "PREREGISTERED G4"})
    out += [{"run_id": "D1_SG_quasi_static_line_dt50us", "script": SCRIPTS[1], "system": QS, "dt_s": 50e-6, "duration_s": 10.0, "disturbance": SG_EV, "role": "diagnostic"},
            {"run_id": "D2_SG_emt_line_dt25us", "script": SCRIPTS[1], "system": UT, "dt_s": 25e-6, "duration_s": 10.0, "disturbance": SG_EV, "role": "diagnostic"},
            {"run_id": "D3_SG_emt_line_dt50us_vs_dynphasor", "script": SCRIPTS[1], "system": UT, "dt_s": 50e-6, "duration_s": 10.0, "disturbance": SG_EV, "role": "diagnostic"}]
    for g in (0.0, 0.03625, 0.25):
        for c in "abc":
            out.append({"run_id": f"D1_GFL_quasi_static_line_g{g:g}_{c}_dt50us", "script": SCRIPTS[1], "system": QS, "dt_s": 50e-6, "duration_s": 5.0, "disturbance": GFL_EV[c], "role": "diagnostic"})
    out.append({"run_id": "D4_GFL_g0.25_no_event_every_step", "script": SCRIPTS[1], "system": UT, "dt_s": 50e-6, "duration_s": 0.7, "disturbance": "none", "role": "diagnostic"})
    for name in ("P4_BASE", "P4_H4", "GS_H4"):
        out.append({"run_id": f"D5_IEEE39_{name}_no_event_dt50us", "script": SCRIPTS[2], "system": "IEEE-39 preregistered realization", "dt_s": 50e-6, "duration_s": 30.0, "disturbance": "none", "role": "diagnostic (no estimator, no prediction comparison)"})
    for r in out:
        r["seed"] = None
    return out


def main() -> int:
    walls = {}
    for s in SCRIPTS:
        t0 = time.time()
        subprocess.run([sys.executable, str(HERE / s)], check=True, capture_output=True)
        walls[s] = round(time.time() - t0, 1)
    files = sorted(p for d in ("EMT02", "EMT03") for p in (_emt.RESULTS / d).iterdir() if p.name != "EMT02_03_manifest.json")
    man = {**_emt.provenance(), "script_wall_s": walls, "runs": runs(),
           "outputs": {str(p.relative_to(_emt.RESEARCH)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
           "note": "phasor_ref.npz / dynphasor_ref.npz / op.json are produced by the tx3-analysis scripts EMT02_03_phasor_refs.py and EMT02_03_dynphasor_ref.py (single-thread BLAS, bit-reproducible)"}
    (_emt.RESULTS / "EMT02" / "EMT02_03_manifest.json").write_text(json.dumps(man, indent=1), encoding="utf-8")
    print(json.dumps({k: man[k] for k in ("git_head", "paremt_upstream_sha", "paremt_working_copy_diff_sha", "script_wall_s")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
