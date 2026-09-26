"""Diagnostic: which ANDES static injection is outside Line + Shunt + PV/Slack + PQ?"""

from __future__ import annotations

import sys
from pathlib import Path

import andes
import numpy as np

andes.config_logger(stream_level=40)
case = Path(andes.__file__).resolve().parent / "cases/ieee39/ieee39_full.xlsx"
ss = andes.load(str(case), setup=True, no_output=True)
ss.PFlow.config.tol = 1e-12
ss.PFlow.run()
bus_ids = [int(b) for b in ss.Bus.idx.v]
V = np.asarray(ss.Bus.v.v) * np.exp(1j * np.asarray(ss.Bus.a.v))
Y = np.array(andes.shared.matrix(ss.build_ybus()), dtype=complex)
S_net = V * np.conj(Y @ V)
S_dev = np.zeros(len(bus_ids), complex)
for b, p, q in zip(ss.PV.bus.v, ss.PV.p.v, ss.PV.q.v, strict=False):
    S_dev[bus_ids.index(int(b))] += complex(p, q)
for b, p, q in zip(ss.Slack.bus.v, ss.Slack.p.v, ss.Slack.q.v, strict=False):
    S_dev[bus_ids.index(int(b))] += complex(p, q)
for b, p, q in zip(ss.PQ.bus.v, ss.PQ.p0.v, ss.PQ.q0.v, strict=False):
    S_dev[bus_ids.index(int(b))] -= complex(p, q)
d = S_net - S_dev
for i in np.argsort(-np.abs(d))[:8]:
    print(f"bus {bus_ids[i]:2d}  S_net-S_dev = {d[i]: .3e}  |V|={abs(V[i]):.6f}")

print("\nmodels with PF-active algebraic equations on Bus.a/Bus.v:")
for name, mdl in ss.models.items():
    if mdl.n == 0:
        continue
    flags = mdl.flags
    if getattr(flags, "pflow", False):
        print(f"  {name:10s} n={mdl.n}")
print("\nPQ config:", dict(ss.PQ.config.as_dict()))
print(
    "PQ p0 vs Ppf:",
    float(np.max(np.abs(np.asarray(ss.PQ.p0.v) - np.asarray(ss.PQ.Ppf.v))))
    if hasattr(ss.PQ, "Ppf")
    else "no Ppf",
)
if len(sys.argv) > 1:
    print(ss.PQ.as_df().to_string())
