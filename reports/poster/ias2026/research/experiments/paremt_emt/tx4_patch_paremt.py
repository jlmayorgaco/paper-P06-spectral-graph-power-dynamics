# ruff: noqa: E501  -- long diagnostic strings and code-patch literals kept on one line
"""Apply the documented TX4 network patch to the ParaEMT working copy (external/ParaEMT_tx4).

Idempotent; the pristine clone external/ParaEMT_upstream is never touched. Two changes only,
both with stock behaviour as the default:

1. Off-nominal transformer taps (upstream loads pfd.xfmr_k but never uses it).
   An ideal tap t on the FROM side of the series R-L branch: G stamps
   y/t^2 (from-from), y (to-to), -y/t (off-diagonal); branch current i = (v_F/t - v_T)/Req + Ihis,
   injected as -i/t at F and +i at T. The tap is stored in a new column 9 of Init_net_coe0
   (1.0 for every non-transformer branch, so every other branch is unchanged).
   Lib_BW.InitNet passes the taps only if ini.apply_xfmr_tap is True (default False = stock);
   ini.xfmr_tap_side ('to' for the stock PSS/E-derived pfd, where the power flow closes with the
   ratio on the to-bus) swaps from/to so that the kernel always sees a from-side tap.
2. A switch for ParaEMT's numerical damping elements (the parallel resistor of line/load R-L
   branches, Rp = 20/3 * 2L/ts, and the series resistor of capacitive branches, Rs = 0.15 ts/2C):
   net_damping = 1.0 keeps them (stock), 0.0 removes them (pure trapezoidal rule), so that the
   realized continuous-time elements equal the documented R, L, C exactly.
   Lib_BW.InitNet passes ini.net_damping (default 1.0 = stock).
The time-step history update (numba_updateIhis) and the two history loops of Re_Init are made
tap-aware in the same way.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[6]
TX4 = REPO / "external" / "ParaEMT_tx4"
MARK = "# TX4-PATCH v1"


def patch_lib_numba(src: str) -> str:
    if MARK in src:
        return src
    s = src
    s = s.replace(
        "        # OTHER\n        ts,\n        loadmodel_option,\n):\n    nbus = len(bus_num)",
        "        # OTHER\n        ts,\n        loadmodel_option,\n        xfmr_k,\n        net_damping,\n):\n"
        f"    {MARK}: xfmr_k = from-side tap per transformer; net_damping 1.0 stock, 0.0 pure trapezoidal\n"
        "    damp_on = net_damping > 0.5\n    nbus = len(bus_num)",
    )
    s = s.replace(
        "    Init_net_coe0 = np.zeros((nbranch, 9), dtype=np.complex128)",
        "    Init_net_coe0 = np.zeros((nbranch, 10), dtype=np.complex128)\n    Init_net_coe0[:, 9] = 1.0",
    )
    # every 9-column row assignment writes columns 0..8 only
    s = re.sub(r"Init_net_coe0\[([^\]]+?),\s*:\]\s*=", r"Init_net_coe0[\1, :9] =", s)
    # line R-L damping resistor
    s = s.replace(
        "        if X>0:\n            L =  X / ws\n            Rp = damptrap * (20.0 / 3.0 * 2.0 * L / ts)\n            Rp_inv = 1.0 / Rp\n",
        "        if X>0:\n            L =  X / ws\n            Rp = damptrap * (20.0 / 3.0 * 2.0 * L / ts)\n            Rp_inv = 1.0 / Rp\n"
        "            if not damp_on:\n                Rp_inv = 0.0\n",
    )
    # line charging series resistor
    s = s.replace(
        "        Rs = Rs / damptrap\n\n        if Rs == np.inf:",
        "        Rs = Rs / damptrap\n        if (not damp_on) and Rs != np.inf:\n            Rs = 0.0\n\n        if Rs == np.inf:",
    )
    # load R-L damping resistor and load R||C series resistor
    s = s.replace(
        "                L = X / ws\n                Rp = damptrap * (20.0 / 3.0 * 2.0 * L / ts)\n                Rp_inv = 1.0 / Rp\n\n                C = 0.0",
        "                L = X / ws\n                Rp = damptrap * (20.0 / 3.0 * 2.0 * L / ts)\n                Rp_inv = 1.0 / Rp\n"
        "                if not damp_on:\n                    Rp_inv = 0.0\n\n                C = 0.0",
    )
    s = s.replace(
        "                Rs = 0.15 * ts / 2.0 / C*0.001\n",
        "                Rs = 0.15 * ts / 2.0 / C*0.001\n                if not damp_on:\n                    Rs = 0.0\n",
    )
    # shunt and switched-shunt series resistors
    s = s.replace(
        "        Rs = 0.15 * ts / 2.0 / C / damptrap\n",
        "        Rs = 0.15 * ts / 2.0 / C / damptrap\n        if not damp_on:\n            Rs = 0.0\n",
    )
    # transformer stamps with the from-side tap
    old_x = s[
        s.index(
            "        idx = 12*len(line_from) + 12*i\n        numba_set_coo(G0_rows, G0_cols, G0_data, idx, Fidx, Fidx, 1 / Req)"
        ) : s.index("        coe_idx = 9*len(line_from) + 3*i\n")
    ]
    new_x = (
        "        idx = 12*len(line_from) + 12*i\n"
        "        tap = xfmr_k[i]\n"
        "        for ph in range(3):\n"
        "            off = ph * N1\n"
        "            numba_set_coo(G0_rows, G0_cols, G0_data, idx + 4*ph, Fidx+off, Fidx+off, 1 / Req / (tap*tap))\n"
        "            numba_set_coo(G0_rows, G0_cols, G0_data, idx + 4*ph + 1, Tidx+off, Tidx+off, 1 / Req)\n"
        "            numba_set_coo(G0_rows, G0_cols, G0_data, idx + 4*ph + 2, Fidx+off, Tidx+off, -1 / Req / tap)\n"
        "            numba_set_coo(G0_rows, G0_cols, G0_data, idx + 4*ph + 3, Tidx+off, Fidx+off, -1 / Req / tap)\n\n"
        "        # R-L branch (series side of the ideal tap)\n"
        "        iA_temp = (Init_net_Vt[Fidx] / tap - Init_net_Vt[Tidx]) / complex(R, ws * L)\n"
        "        iB_temp = (Init_net_Vt[Fidx + N1] / tap - Init_net_Vt[Tidx + N1]) / complex(R, ws * L)\n"
        "        iC_temp = (Init_net_Vt[Fidx + N2] / tap - Init_net_Vt[Tidx + N2]) / complex(R, ws * L)\n\n"
    )
    s = s.replace(old_x, new_x)
    s = s.replace(
        "        Init_net_coe0[coe_idx+2, :9] = np.array([Fidx + N2, Tidx + N2, Req, icf,\n                                               Gv1, R, L, 0.0, iC_temp])\n",
        "        Init_net_coe0[coe_idx+2, :9] = np.array([Fidx + N2, Tidx + N2, Req, icf,\n                                               Gv1, R, L, 0.0, iC_temp])\n"
        "        Init_net_coe0[coe_idx:coe_idx+3, 9] = tap\n",
    )
    # initial history with the tap
    s = s.replace(
        "            brch_Ihis_temp = Init_net_coe0[i, 3] * Init_brch_Ipre[i] + Init_net_coe0[i, 4] * (Init_net_V[Fidx] - Init_net_V[Tidx])\n"
        "            Init_node_Ihis[Tidx] += brch_Ihis_temp.real\n\n        Init_brch_Ihis[i] = brch_Ihis_temp.real\n        Init_node_Ihis[Fidx] -= brch_Ihis_temp.real\n",
        "            tap = Init_net_coe0[i, 9].real\n"
        "            brch_Ihis_temp = Init_net_coe0[i, 3] * Init_brch_Ipre[i] + Init_net_coe0[i, 4] * (Init_net_V[Fidx] / tap - Init_net_V[Tidx])\n"
        "            Init_node_Ihis[Tidx] += brch_Ihis_temp.real\n            Init_brch_Ihis[i] = brch_Ihis_temp.real\n"
        "            Init_node_Ihis[Fidx] -= brch_Ihis_temp.real / tap\n            continue\n\n"
        "        Init_brch_Ihis[i] = brch_Ihis_temp.real\n        Init_node_Ihis[Fidx] -= brch_Ihis_temp.real\n",
    )
    # time-step history update with the tap
    s = s.replace(
        "        else:\n            brch_Ipre[i] = (Vsol[Fidx] - Vsol[Tidx])/Init_net_coe0[i,2].real + brch_Ihis[i]\n"
        "            brch_Ihis_temp = Init_net_coe0[i,3] * brch_Ipre[i] + Init_net_coe0[i,4] * (Vsol[Fidx] - Vsol[Tidx])\n"
        "            node_Ihis[Tidx] += brch_Ihis_temp.real\n        brch_Ihis[i] = brch_Ihis_temp.real\n        node_Ihis[Fidx] -= brch_Ihis_temp.real\n",
        "        else:\n            tap = Init_net_coe0[i,9].real\n"
        "            brch_Ipre[i] = (Vsol[Fidx] / tap - Vsol[Tidx])/Init_net_coe0[i,2].real + brch_Ihis[i]\n"
        "            brch_Ihis_temp = Init_net_coe0[i,3] * brch_Ipre[i] + Init_net_coe0[i,4] * (Vsol[Fidx] / tap - Vsol[Tidx])\n"
        "            node_Ihis[Tidx] += brch_Ihis_temp.real\n            brch_Ihis[i] = brch_Ihis_temp.real\n"
        "            node_Ihis[Fidx] -= brch_Ihis_temp.real / tap\n            continue\n"
        "        brch_Ihis[i] = brch_Ihis_temp.real\n        node_Ihis[Fidx] -= brch_Ihis_temp.real\n",
    )
    return s


def patch_lib_bw(src: str) -> str:
    if MARK in src:
        return src
    s = src.replace(
        "    def InitNet(self, pfd, ts, loadmodel_option):\n        (self.Init_net_VbaseA,",
        f"    def InitNet(self, pfd, ts, loadmodel_option):\n        {MARK}: optional taps and damping switch (defaults = stock)\n"
        "        xf, xt = pfd.xfmr_from, pfd.xfmr_to\n"
        "        taps = np.ones(len(pfd.xfmr_from))\n"
        "        if getattr(self, 'apply_xfmr_tap', False):\n"
        "            taps = np.real(np.asarray(pfd.xfmr_k, dtype=np.complex128)).astype(np.float64)\n"
        "            if getattr(self, 'xfmr_tap_side', 'to') == 'to':\n"
        "                xf, xt = pfd.xfmr_to, pfd.xfmr_from\n"
        "        (self.Init_net_VbaseA,",
    )
    s = s.replace(
        "            pfd.xfmr_from,\n            pfd.xfmr_to,\n            pfd.xfmr_RX,\n            pfd.load_bus,",
        "            xf,\n            xt,\n            pfd.xfmr_RX,\n            pfd.load_bus,",
    )
    s = s.replace(
        "            ts,\n            loadmodel_option,\n        )\n        return\n\n    def InitMac",
        "            ts,\n            loadmodel_option,\n            taps,\n            float(getattr(self, 'net_damping', 1.0)),\n        )\n        return\n\n    def InitMac",
    )
    # the two Re_Init history loops (used only after a generator trip)
    for c in ("c1*", ""):
        old = (
            f"            else:\n                brch_Ihis_temp = {c}Init_net_coe0[i, 3] * brch_Ipre[i] + "
            f"{'c2*np.real(' if c else ''}Init_net_coe0[i, 4]{')' if c else ''} * (\n"
        )
        if old in s:
            s = s.replace(
                old,
                old.replace(
                    "            else:\n",
                    "            else:\n                tap = Init_net_coe0[i, 9].real\n",
                ),
            )
    s = s.replace(
        "                            Vsol[Fidx] - Vsol[Tidx])\n                node_Ihis[Tidx] += brch_Ihis_temp.real\n",
        "                            Vsol[Fidx] / tap - Vsol[Tidx])\n                node_Ihis[Tidx] += brch_Ihis_temp.real\n                node_Ihis[Fidx] += brch_Ihis_temp.real * (1.0 - 1.0 / tap)\n",
    )
    s = s.replace(
        "                            Vsol[Fidx] - Vsol[Tidx])\n                node_Ihis_out[Tidx] += brch_Ihis_temp.real\n",
        "                            Vsol[Fidx] / tap - Vsol[Tidx])\n                node_Ihis_out[Tidx] += brch_Ihis_temp.real\n                node_Ihis_out[Fidx] += brch_Ihis_temp.real * (1.0 - 1.0 / tap)\n",
    )
    return s


def main() -> int:
    for name, fn in (("lib_numba.py", patch_lib_numba), ("Lib_BW.py", patch_lib_bw)):
        path = TX4 / name
        raw = path.read_bytes().decode("utf-8")
        crlf = "\r\n" in raw
        src = raw.replace("\r\n", "\n")
        new = fn(src)
        if new != src:
            out = new.replace("\n", "\r\n") if crlf else new
            path.write_bytes(out.encode("utf-8"))
            print("patched", name)
        else:
            print("unchanged (already patched)", name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
