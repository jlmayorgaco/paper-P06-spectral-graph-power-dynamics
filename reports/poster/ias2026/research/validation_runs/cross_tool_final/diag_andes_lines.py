# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""Diagnostic: ANDES Line PF equations vs the two-port of Line.build_ybus, per line."""

from __future__ import annotations

from pathlib import Path

import andes
import numpy as np

andes.config_logger(stream_level=40)
case = Path(andes.__file__).resolve().parent / "cases/ieee39/ieee39_full.xlsx"
ss = andes.load(str(case), setup=True, no_output=True)
ss.PFlow.config.tol = 1e-12
ss.PFlow.run()
L = ss.Line
v = np.asarray(ss.Bus.v.v)
a = np.asarray(ss.Bus.a.v)
V = v * np.exp(1j * a)
i1, i2 = np.asarray(L.a1.a), np.asarray(L.a2.a)
ysh = L.u.v * (L.g.v + 1j * L.b.v) / 2
y1 = L.u.v * (L.g1.v + 1j * L.b1.v)
y2 = L.u.v * (L.g2.v + 1j * L.b2.v)
y12 = L.u.v / (L.r.v + 1j * L.x.v)
m = L.tap.v * np.exp(1j * L.phi.v)
If = (y12 + y1 + ysh) / L.tap.v**2 * V[i1] - y12 / np.conj(m) * V[i2]
It = -y12 / m * V[i1] + (y12 + y2 + ysh) * V[i2]
Sf = V[i1] * np.conj(If)
St = V[i2] * np.conj(It)
# ANDES equation values: ExtAlgeb .e holds the device contribution (consumption sign)
Pf_andes, Qf_andes = np.asarray(L.a1.e), np.asarray(L.v1.e)
Pt_andes, Qt_andes = np.asarray(L.a2.e), np.asarray(L.v2.e)
err = np.abs(
    np.c_[
        Sf.real - Pf_andes, Sf.imag - Qf_andes, St.real - Pt_andes, St.imag - Qt_andes
    ]
)
print("max |two-port - ANDES line eq| =", err.max())
for k in np.argsort(-err.max(axis=1))[:8]:
    print(
        f"line {k:2d} {int(L.bus1.v[k])}-{int(L.bus2.v[k])} tap={L.tap.v[k]} r={L.r.v[k]} x={L.x.v[k]} "
        f"b={L.b.v[k]} err={err[k].max():.3e}  Sf={Sf[k]:.6f} andes=({Pf_andes[k]:.6f},{Qf_andes[k]:.6f})"
    )
print("\nLine equation strings:")
for name in ("a1", "a2", "v1", "v2"):
    print(" ", name, getattr(L, name).e_str)
for name in ("gh", "bh", "gk", "bk", "ghk", "bhk", "itap", "itap2"):
    obj = getattr(L, name, None)
    if obj is not None:
        print(" ", name, getattr(obj, "v_str", None) or getattr(obj, "e_str", None))
