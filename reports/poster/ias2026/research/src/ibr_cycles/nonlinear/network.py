"""L1 network: branch incidence, complex taps, pi branches, shunts (v2 eq. N.7-N.12).

    Cf, Ct   (n_branch x n_bus) from- and to-bus selectors
    per branch, with series y = 1/(r + jx), total charging b, tap t = |t| e^{j phi}
    on the FROM side (the convention of the frozen loader and of MATPOWER):

        Yff = (y + j b/2)/|t|^2,  Yft = -y/conj(t),  Ytf = -y/t,  Ytt = y + j b/2
        Yf = diag(Yff) Cf + diag(Yft) Ct,   Yt = diag(Ytf) Cf + diag(Ytt) Ct
        Ybus = Cf^T Yf + Ct^T Yt + diag(y_sh)

    branch currents I_f = Yf V, I_t = Yt V; complex power S_f = V_f conj(I_f).
    With a phase shifter Ybus is NOT symmetric; the test is per-branch current,
    power and loss:  Re(S_f + S_t) = r |y (V_f / t - V_t)|^2.

Realification (v2 eq. N.3): STACKED z = [Re V; Im V] or INTERLEAVED
z = [vx1, vy1, vx2, vy2, ...] (the frozen DAE). The permutation between them is
explicit here and the two are never mixed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

Complex = NDArray[np.complex128]


@dataclass(frozen=True)
class BranchNetwork:
    bus_idx: tuple[int, ...]
    cf: NDArray[np.float64]
    ct: NDArray[np.float64]
    yff: Complex
    yft: Complex
    ytf: Complex
    ytt: Complex
    y_series: Complex
    r: NDArray[np.float64]
    tap: Complex
    y_shunt: Complex
    source_path: str

    @property
    def yf(self) -> Complex:
        return self.yff[:, None] * self.cf + self.yft[:, None] * self.ct

    @property
    def yt(self) -> Complex:
        return self.ytf[:, None] * self.cf + self.ytt[:, None] * self.ct

    @property
    def ybus(self) -> Complex:
        return self.cf.T @ self.yf + self.ct.T @ self.yt + np.diag(self.y_shunt)

    def branch_flows(self, v: Complex) -> dict[str, Complex]:
        i_f, i_t = self.yf @ v, self.yt @ v
        v_f, v_t = self.cf @ v, self.ct @ v
        s_f, s_t = v_f * np.conj(i_f), v_t * np.conj(i_t)
        series = self.y_series * (v_f / self.tap - v_t)
        return {
            "I_f": i_f,
            "I_t": i_t,
            "S_f": s_f,
            "S_t": s_t,
            "loss": (s_f + s_t).real,
            "loss_series": self.r * np.abs(series) ** 2,
        }


def load_branch_network(path: str | Path) -> BranchNetwork:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    order = [int(b["idx"]) for b in payload["buses"]]
    pos = {b: k for k, b in enumerate(order)}
    lines = [ln for ln in payload["lines"] if float(ln["u"]) != 0.0]
    nb, nl = len(order), len(lines)
    cf, ct = np.zeros((nl, nb)), np.zeros((nl, nb))
    y, bch, tap, r = (
        np.zeros(nl, complex),
        np.zeros(nl),
        np.zeros(nl, complex),
        np.zeros(nl),
    )
    for k, ln in enumerate(lines):
        cf[k, pos[int(ln["bus1"])]] = 1.0
        ct[k, pos[int(ln["bus2"])]] = 1.0
        z = complex(ln["r"], ln["x"])
        y[k] = 1.0 / z
        bch[k] = float(ln["b"])  # total charging susceptance (g assumed 0 in the data)
        tap[k] = float(ln["tap"]) * np.exp(1j * float(ln["phi"]))
        r[k] = float(ln["r"])
        if float(ln.get("g", 0.0)) != 0.0:
            raise ValueError(
                "branch conductance g != 0 is not in the frozen data format"
            )
    half = 1j * bch / 2.0
    ysh = np.zeros(nb, complex)
    for sh in payload.get("shunts", []):
        ysh[pos[int(sh["bus"])]] += complex(sh["g"], sh["b"])
    return BranchNetwork(
        bus_idx=tuple(order),
        cf=cf,
        ct=ct,
        yff=(y + half) / np.abs(tap) ** 2,
        yft=-y / np.conj(tap),
        ytf=-y / tap,
        ytt=y + half,
        y_series=y,
        r=r,
        tap=tap,
        y_shunt=ysh,
        source_path=str(path),
    )


def interleave_permutation(n_bus: int) -> NDArray[np.int64]:
    """Index array p with z_interleaved = z_stacked[p]."""

    p = np.empty(2 * n_bus, dtype=np.int64)
    p[0::2] = np.arange(n_bus)
    p[1::2] = n_bus + np.arange(n_bus)
    return p


def realify(m: Complex, order: str = "stacked") -> NDArray[np.float64]:
    """[[Re, -Im], [Im, Re]] in stacked order, or its interleaved permutation."""

    out = np.block([[m.real, -m.imag], [m.imag, m.real]])
    if order == "stacked":
        return out
    p = interleave_permutation(m.shape[0])
    return out[np.ix_(p, p)]
