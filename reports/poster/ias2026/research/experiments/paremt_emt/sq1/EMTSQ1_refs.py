# ruff: noqa: E501  -- test tables kept on one line
# STATUS: PREPARED BUT NEVER EXECUTED - SQ1 stopped at SQ1-2 (blind synthetic qualification FAIL); the ParaEMT line is stopped permanently (prereg SQ1).
"""SQ1-4/6 references and canonical parameter dumps (.venv/tx3-analysis; run only if SQ1-2 and SQ1-3 pass).

B-SG38 (blind): quasi-static and dynamic-line references through the V2 reference functions unchanged.
Canonical parameters for the exact parameter/base checks: SynchronousMachine parameters of buses 30,
36, 38 with the P4 AVR scaling (ka x1.425, ta x1.5) and weight Sn/100; ConverterParameters defaults,
leak and weight for buses 30, 35, 37.
Writes results/EMTSQ1/refs/{B-SG38_qs.npz, B-SG38_dyn.npz, canonical_params.json, refs_manifest.json}.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import asdict, replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "v2"))

import EMTV2_refs as R2  # noqa: E402
import numpy as np  # noqa: E402

OUT = R2.RESEARCH / "results" / "EMTSQ1" / "refs"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    man = {"files": {}, "meta": {}}
    qs, dyn, meta = R2.sg_refs(38, 7.64783381 + 1.23275808j)
    for kind, d in (("qs", qs), ("dyn", dyn)):
        p = OUT / f"B-SG38_{kind}.npz"
        np.savez_compressed(p, **d)
        man["files"][str(p.relative_to(R2.RESEARCH)).replace("\\", "/")] = hashlib.sha256(p.read_bytes()).hexdigest()
    man["meta"]["B-SG38"] = meta
    net = R2.load_network()
    payload = R2._controller_payload(net.config_path or None)
    canon = {"sg": {}, "gfl": {}}
    for bus in (30, 36, 38):
        p = R2._machine_parameters(net, bus, payload)
        p = replace(p, ka=p.ka * 1.425, ta=p.ta * 1.5)
        canon["sg"][str(bus)] = {**{k: v for k, v in asdict(p).items() if isinstance(v, (int, float)) and not isinstance(v, bool)}, "w": net.machines[bus]["Sn"] / 100.0}
    for bus in (30, 35, 37):
        par = R2.ConverterParameters(voltage_control=True, voltage_gain=0.0, voltage_leak=R2.LEAK)
        canon["gfl"][str(bus)] = {**{k: v for k, v in asdict(par).items() if isinstance(v, (int, float)) and not isinstance(v, bool)}, "w": net.machines[bus]["Sn"] / 100.0}
    (OUT / "canonical_params.json").write_text(json.dumps(canon, indent=1), encoding="utf-8")
    (OUT / "refs_manifest.json").write_text(json.dumps(man, indent=1, default=float), encoding="utf-8")
    print("SQ1 references written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
