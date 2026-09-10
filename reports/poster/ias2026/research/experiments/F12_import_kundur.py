"""F12 - freeze the documented Kundur two-area case into the project network format.

Source: ANDES 2.0.0 bundled case ``cases/kundur/kundur_sexs.xlsx`` (Kundur,
Power System Stability and Control, Example 12.6 network and machines, SEXS
excitation). Run with the ANDES interpreter (.venv/tx3-andes). It writes

    configs/kundur/kundur_network.json          same schema as ieee39_network.json
    configs/kundur/kundur_andes_reference.json  ANDES Ybus and power flow, for the cross-check

No parameter is changed. Declared reductions happen downstream, in the model:
GENROU is represented by its two-axis transient model (subtransient constants
0.03-0.05 s dropped), TGOV1 is not represented (as in the IEEE-39 study), and the
case has no stabilizer, so the PSS gain is zero.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import andes
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CASE = Path(andes.__file__).resolve().parent / "cases" / "kundur" / "kundur_sexs.xlsx"
OUT = ROOT / "configs" / "kundur"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    sheets = pd.read_excel(CASE, sheet_name=None)
    sha = hashlib.sha256(CASE.read_bytes()).hexdigest()
    bus = sheets["Bus"]
    line = sheets["Line"]
    pq = sheets["PQ"]
    pv = sheets["PV"]
    slack = sheets["Slack"]
    gen = sheets["GENROU"]
    sexs = sheets["SEXS"]
    machines, avr, pss = [], [], []
    for _, g in gen.iterrows():
        machines.append({"bus": float(g.bus), "Sn": float(g.Sn), "D": float(g.D), "M": float(g.M),
                         "ra": float(g.ra), "xl": float(g.xl), "xd": float(g.xd), "xq": float(g.xq),
                         "xd1": float(g.xd1), "xq1": float(g.xq1), "xd2": float(g.xd2),
                         "xq2": float(g.xq2), "Td10": float(g.Td10), "Td20": float(g.Td20),
                         "Tq10": float(g.Tq10), "Tq20": float(g.Tq20)})
        e = sexs[sexs.syn == g.idx].iloc[0]
        avr.append({"syn": f"GENROU_{int(g.idx)}", "model": "SEXS", "KA": float(e.K),
                    "TE": float(e.TE), "TATB": float(e.TATB), "TB": float(e.TB),
                    "EMIN": float(e.EMIN), "EMAX": float(e.EMAX)})
        pss.append({"avr": f"SEXS_{int(e.idx)}", "model": "none", "KS": 0.0,
                    "T4": 1.0, "T5": 1.0, "T6": 1.0})
    payload = {
        "source": {"file": str(CASE.name), "andes": andes.__version__, "sha256": sha,
                   "note": "Kundur two-area, SEXS excitation; no PSS in the documented case"},
        "buses": [{"idx": float(r.idx), "Vn": float(r.Vn), "v0": float(r.v0), "a0": float(r.a0),
                   "area": float(r.area)} for r in bus.itertuples()],
        "loads": [{"bus": float(r.bus), "p0": float(r.p0), "q0": float(r.q0)} for r in pq.itertuples()],
        "pv": [{"idx": float(r.idx), "bus": float(r.bus), "Sn": float(r.Sn), "p0": float(r.p0),
                "q0": float(r.q0), "qmax": 99.0, "qmin": -99.0, "v0": float(r.v0),
                "pmax": float(r.pmax), "pmin": float(r.pmin)} for r in pv.itertuples()],
        "slack": [{"idx": float(r.idx), "bus": float(r.bus), "Sn": float(r.Sn), "p0": float(r.p0),
                   "q0": float(r.q0), "v0": float(r.v0), "a0": float(r.a0), "pmax": 99.0,
                   "pmin": 0.0} for r in slack.itertuples()],
        "shunts": [],
        "lines": [{"bus1": float(r.bus1), "bus2": float(r.bus2), "r": float(r.r), "x": float(r.x),
                   "b": float(r.b), "g": float(r.g), "trans": float(r.trans), "tap": float(r.tap),
                   "phi": float(r.phi), "u": float(r.u)} for r in line.itertuples()],
        "machines": machines, "avr": avr, "pss": pss, "gov": [],
        "machine_order": [m["bus"] for m in machines],
        "avr_order": [a["syn"] for a in avr], "pss_order": [p["avr"] for p in pss],
    }
    (OUT / "kundur_network.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")

    # independent reference: ANDES power flow and admittance matrix
    system = andes.load(str(CASE), setup=True, no_output=True, default_config=True)
    system.PFlow.run()
    ybus = None
    if hasattr(system.Line, "build_y"):
        y = system.Line.build_y()
        ybus = np.array(__import__("kvxopt").matrix(y))
    reference = {
        "pflow_converged": bool(system.PFlow.converged),
        "bus_idx": [int(i) for i in system.Bus.idx.v],
        "v": [float(v) for v in system.Bus.v.v], "a": [float(a) for a in system.Bus.a.v],
        "ybus_real": None if ybus is None else np.real(np.asarray(ybus)).tolist(),
        "ybus_imag": None if ybus is None else np.imag(np.asarray(ybus)).tolist(),
    }
    (OUT / "kundur_andes_reference.json").write_text(json.dumps(reference, indent=1), encoding="utf-8")
    print(json.dumps({"sha256": sha, "pflow": reference["pflow_converged"], "buses": len(reference["v"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
