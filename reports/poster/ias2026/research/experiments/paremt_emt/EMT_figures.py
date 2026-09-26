# ruff: noqa: E501  -- plotting calls kept on one line
"""Figures of the ParaEMT campaign as executed (prereg v1; campaign BLOCKED after G3/G4).

F1  network / equilibrium equivalence (EMT01, preregistered, PASS).
FD1 SG unit test (EMT02, preregistered FAIL) with diagnostics D1-D3.
FD2 GFL unit test (EMT03, preregistered FAIL) with diagnostics D1 and D4.
F2-F8 need EMT04-EMT16 and are not produced (BLOCKED).
Writes figures/20260911_EMT_{F1,FD1,FD2}_*.{pdf,png} and figures/20260911_EMT_figure_source.csv.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _emt  # noqa: E402
import EMT02_03_diagnostics as D  # noqa: E402
import EMT02_03_unit_tests as U  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import tx4_case as C  # noqa: E402
import tx4_emt as K  # noqa: E402

FIG = _emt.RESEARCH / "figures"
RES = _emt.RESULTS
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e4e3df"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.titlesize": 8.5, "axes.titleweight": "bold", "legend.fontsize": 7, "legend.frameon": False,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 1.4})


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)


def f1(src):
    s1 = json.loads((RES / "EMT01" / "EMT01_summary.json").read_text())
    with (RES / "EMT01" / "equilibrium_bus_voltages.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    bus = np.array([int(r["bus"]) for r in rows])
    dv = np.abs([float(r["dVm"]) for r in rows])
    da = np.abs([float(r["dAngle_rel39"]) for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.5), sharex=True)
    ref = bus != 39  # bus 39 is the angle reference (0 by definition)
    for a, bb, y, tol, lab in ((ax[0], bus, dv, 1e-4, "|ΔV| (pu)"), (ax[1], bus[ref], da[ref], 1e-3, "|Δθ| rel. bus 39 (rad)")):
        a.vlines(bb, 1e-10, y, color=BLUE, linewidth=1.2)
        a.plot(bb, y, "o", color=BLUE, markersize=3.2)
        a.axhline(tol, color=INK2, linestyle="--", linewidth=1.0)
        a.text(40.3, tol * 1.35, f"preregistered tolerance {tol:.0e}", ha="right", va="bottom", color=INK2, fontsize=7)
        a.set_yscale("log")
        a.set_ylim(1e-10, 3e-3)
        a.set_xlabel("bus")
        a.set_ylabel(lab)
        a.set_xlim(0, 40.5)
    ax[0].set_title("(a) voltage-magnitude residual", loc="left")
    ax[1].set_title("(b) voltage-angle residual", loc="left")
    fig.tight_layout()
    src.append(["F1", "network", "ybus_rel_residual", s1["ybus"]["50us"]["rel_residual_continuous"]])
    save(fig, "20260911_EMT_F1_network_equivalence")
    for b, v, a in zip(bus, dv, da, strict=True):
        src.append(["F1", f"bus {b}", "abs_dVm_pu", v])
        src.append(["F1", f"bus {b}", "abs_dAngle_rad", a])


def fd1(src):
    op = json.loads((RES / "EMT02" / "op.json").read_text())
    ref = np.load(RES / "EMT02" / "phasor_ref.npz")
    dyn = np.load(RES / "EMT02" / "dynphasor_ref.npz")
    runs = {"EMT line, 50 µs vs quasi-static phasor (preregistered)": (D.sg_setup(op, "emt", 50e-6), ref, ORANGE, "-", "o"),
            "EMT line, 25 µs vs quasi-static phasor (D2)": (D.sg_setup(op, "emt", 25e-6), ref, YELLOW, "--", "s"),
            "EMT line, 50 µs vs dynamic-phasor line (D3)": (D.sg_setup(op, "emt", 50e-6), dyn, AQUA, "-", "^"),
            "algebraic line, 50 µs vs quasi-static phasor (D1)": (D.sg_setup(op, "qs", 50e-6), ref, BLUE, "-", "D")}
    m0 = ref["t"] >= 1.0 - 1e-12
    exc = np.abs(ref["x"][m0] - ref["x"][0]).max(axis=0)
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1.35, 1]})
    for lab, ((t, x, _fin), r, col, ls, mk) in runs.items():
        err = x[:, 4] - r["x"][: len(t), 4]
        ax[0].plot(t, err * 1e4, color=col, linestyle=ls, linewidth=1.2)
        m = t >= 1.0 - 1e-12
        rat = np.abs(x[m] - r["x"][: len(t)][m]).max(axis=0) / exc
        ax[1].plot(np.arange(7), rat, marker=mk, linestyle="none", color=col, markersize=4.2)
        for j, nm in enumerate(U.SG_NAMES):
            src.append(["FD1", lab, f"ratio_{nm}", rat[j]])
    tol = 0.02 * exc[4]
    ax[0].axhspan(-tol * 1e4, tol * 1e4, color=GRID, alpha=0.6, linewidth=0)
    ax[0].text(9.9, tol * 1e4 * 1.05, "±2 % of the efd excursion", ha="right", va="bottom", fontsize=6.8, color=INK2)
    ax[0].set_xlim(0.5, 10)
    ax[0].set_xlabel("time (s)")
    ax[0].set_ylabel("efd error (1e-4 pu)")
    ax[0].set_title("(a) AVR-state error after the Pm pulse", loc="left")
    ax[1].axhline(0.02, color=INK2, linestyle="--", linewidth=1.0)
    ax[1].text(6.3, 0.023, "limit 0.02", ha="right", va="bottom", fontsize=6.8, color=INK2)
    ax[1].set_yscale("log")
    ax[1].set_ylim(2e-5, 6e-2)
    ax[1].set_xticks(np.arange(7), ["δ", "ω", "e′q", "e′d", "efd", "pss_w", "pss_l"])
    ax[1].set_ylabel("max error / max excursion")
    ax[1].set_title("(b) all seven states", loc="left")
    handles = [matplotlib.lines.Line2D([], [], color=c, linestyle=ls, marker=mk, markersize=4) for (_r, _ref, c, ls, mk) in runs.values()]
    fig.legend(handles, list(runs), loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02), fontsize=6.8)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save(fig, "20260911_EMT_FD1_sg_unit_test")


def alternation(op3, g, mode):
    vt, s, net, gm, ys = D.quasi_static_g(op3) if mode == "qs" else U.common(op3)
    par = dict(C.GFL_DEFAULTS)
    w = op3["meta_g0.25_a"]["w"]
    x, p_ref, q_ref, v_ref = C.gfl_init(par, w, vt, s)
    row = np.zeros(len(K.GFL_COLS))
    for n, v in dict(par, bus=0, w=w, g=g, leak=C.LEAK, p_ref=p_ref, q_ref=q_ref, v_ref=v_ref).items():
        row[K.GF[n]] = v
    U.DT, U.DS = 50e-6, 1
    t, _sg, _gf, vv, fin = U.run_kernel(net, gm, np.zeros((0, 27)), np.zeros(0, np.complex128), np.zeros((0, 9)),
                                        row[None, :], np.array([x]), ys, np.zeros((0, 5)), 0.8)
    vm = np.abs(vv[:, 0, 0] + 1j * vv[:, 0, 1])
    alt = np.abs(np.diff(vm))
    nb = len(alt) // 20
    env = alt[: nb * 20].reshape(nb, 20).max(axis=1)
    return (np.arange(nb) + 1) * 20 * 50e-6, np.maximum(env, 1e-13), bool(fin), float(t[-1])


def fd2(src):
    op3 = json.loads((RES / "EMT03" / "op.json").read_text())
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1.1, 1.25]})
    for g, mode, col, ls, lab in ((0.0, "emt", ORANGE, "-", "g = 0, EMT line"), (0.25, "emt", YELLOW, "--", "g = 0.25, EMT line"),
                                  (0.25, "qs", BLUE, "-", "g = 0.25, algebraic line (D1)")):
        t, env, fin, tend = alternation(op3, g, mode)
        ax[0].plot(t, env, color=col, linestyle=ls, label=lab + ("" if fin else f" (non-finite at {tend:.2f} s)"))
        src.append(["FD2", lab, "max_alternation_pu", float(env.max())])
    ax[0].set_yscale("log")
    ax[0].set_ylim(5e-14, 1e3)
    ax[0].text(0.4, 1.6e-13, "algebraic line: < 1e-13 (plotted at 1e-13)", ha="center", va="bottom", fontsize=6.5, color=INK2)
    ax[0].set_xlabel("time (s), no disturbance")
    ax[0].set_ylabel("step-to-step |Δ|V_T|| (pu, 1-ms max)")
    ax[0].set_title("(a) Nyquist alternation at the GFL terminal", loc="left")
    ax[0].legend(loc="lower right", bbox_to_anchor=(1.0, 0.1), fontsize=6.3)
    worst = {}
    with (RES / "EMT03" / "EMT03_trajectory_comparison.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            k = r["test"].replace("EMT03_", "")
            worst[k] = max(worst.get(k, 0.0), float(r["ratio"]))
    wq = {}
    with (RES / "EMT03" / "EMT03_D1_quasi_static_line.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            k = r["test"].replace("D1_", "")
            wq[k] = max(wq.get(k, 0.0), float(r["ratio"]))
    keys = list(wq)
    xs = np.arange(len(keys))
    top = 1e4
    pre = np.array([worst[k] for k in keys])
    fin = np.isfinite(pre)
    ax[1].plot(xs[fin], pre[fin], "o", color=ORANGE, markersize=4.5, label="preregistered EMT line")
    ax[1].plot(xs[~fin], np.full((~fin).sum(), top), "x", color=ORANGE, markersize=6, markeredgewidth=1.6, label="non-finite (t = 0.76 s)")
    ax[1].plot(xs, [wq[k] for k in keys], "D", color=BLUE, markersize=4, label="algebraic line (D1)")
    ax[1].axhline(0.02, color=INK2, linestyle="--", linewidth=1.0)
    ax[1].text(len(keys) - 0.6, 0.024, "limit 0.02", ha="right", va="bottom", fontsize=6.8, color=INK2)
    ax[1].set_yscale("log")
    ax[1].set_ylim(1e-3, 3e4)
    ax[1].set_xticks(xs, [k.replace("_", " ") for k in keys], rotation=55, ha="right", fontsize=6.5)
    ax[1].set_ylabel("worst-state error / excursion")
    ax[1].set_title("(b) nine preregistered cases", loc="left")
    ax[1].legend(loc="upper right", fontsize=6.3)
    for k in keys:
        src.append(["FD2", k, "worst_ratio_preregistered", worst[k]])
        src.append(["FD2", k, "worst_ratio_D1_algebraic_line", wq[k]])
    fig.tight_layout()
    save(fig, "20260911_EMT_FD2_gfl_unit_test")


def main() -> int:
    src = []
    f1(src)
    fd1(src)
    fd2(src)
    with (FIG / "20260911_EMT_figure_source.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["figure", "series", "quantity", "value"])
        w.writerows([[a, b, c, f"{d:.6e}"] for a, b, c, d in src])
    print("figures written", len(src))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
