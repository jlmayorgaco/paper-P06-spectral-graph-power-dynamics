from pathlib import Path

import json
import numpy as np
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def lmat(w12, w13, w23):
    """Grounded 3-bus susceptance Laplacian after bus 1 is the reference."""
    return np.array([[w12 + w23, -w23], [-w23, w13 + w23]], dtype=float)


def full_state_matrix(w12, w13, w23, d0, ac, bmat, kmat):
    """First-order realization of T(s)=s^2I+s[D0+K(sI-Ac)^-1B]+L."""
    n = 2
    m = ac.shape[0]
    l = lmat(w12, w13, w23)
    a = np.zeros((2 * n + m, 2 * n + m))
    a[:n, n : 2 * n] = np.eye(n)
    a[n : 2 * n, :n] = -l
    a[n : 2 * n, n : 2 * n] = -d0 * np.eye(n)
    a[n : 2 * n, 2 * n :] = -kmat
    a[2 * n :, n : 2 * n] = bmat
    a[2 * n :, 2 * n :] = ac
    return a


def critical_pole(a):
    eig = np.linalg.eigvals(a)
    osc = [s for s in eig if s.imag > 1e-8 and s.real < 1e-8]
    if not osc:
        raise RuntimeError("No stable oscillatory pole found.")
    return min(osc, key=lambda s: -s.real / abs(s))


def zeta(s):
    return -s.real / abs(s)


def network_screen(w12, w13, w23, d0):
    lam = np.linalg.eigvalsh(lmat(w12, w13, w23))
    nu = lam[-1]
    return {
        "nu_max": float(nu),
        "omega_max": float(np.sqrt(nu)),
        "zeta_prop": float(d0 / (2 * np.sqrt(nu))),
    }


def main():
    # Fixed witness found by seeded search.  Bus 1 is a reference, buses 2 and 3
    # are retained dynamic buses, and the 2-state controller has poles
    # -a +/- j wc.
    w12 = 0.8701154312420417
    w13 = 3.6555746846018327
    w23 = 3.912533790408945
    d0 = 0.7886482292277636
    alpha = 0.2216293514211423
    wc = 1.4782870869433227
    ac = np.array([[-alpha, wc], [-wc, -alpha]])
    bmat = np.array([[-0.12534466, -0.24284206], [0.62850042, 0.11412447]])
    kmat = np.array([[-0.62121562, 0.60119178], [0.20600856, -0.89012032]])
    eps = 0.01

    abase = full_state_matrix(w12, w13, w23, d0, ac, bmat, kmat)
    apost = full_state_matrix(w12, w13, w23 * (1 + eps), d0, ac, bmat, kmat)
    s0 = critical_pole(abase)
    s1 = critical_pole(apost)
    mu = np.linalg.eigvals(ac)[0]
    if abs(s0 - np.linalg.eigvals(ac)[1]) < abs(s0 - mu):
        mu = np.linalg.eigvals(ac)[1]
    net0 = network_screen(w12, w13, w23, d0)
    net1 = network_screen(w12, w13, w23 * (1 + eps), d0)

    report = {
        "w12": w12,
        "w13": w13,
        "w23_base": w23,
        "w23_post": w23 * (1 + eps),
        "d0": d0,
        "control_pole": {"real": float(mu.real), "imag": float(mu.imag)},
        "base_pole": {"real": float(s0.real), "imag": float(s0.imag), "zeta": float(zeta(s0))},
        "post_pole": {"real": float(s1.real), "imag": float(s1.imag), "zeta": float(zeta(s1))},
        "delta_zeta": float(zeta(s1) - zeta(s0)),
        "distance_to_control_base": float(abs(s0 - mu)),
        "distance_to_control_post": float(abs(s1 - mu)),
        "network_screen_base": net0,
        "network_screen_post": net1,
        "delta_nu_max": float(net1["nu_max"] - net0["nu_max"]),
    }
    (ROOT / "synthetic_counterexample_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
        }
    )
    fig = plt.figure(figsize=(7.2, 2.85))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.15, 1.15])

    # Panel 1: circuit sketch.
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    p1 = np.array([0.16, 0.78])
    p2 = np.array([0.16, 0.23])
    p3 = np.array([0.82, 0.50])
    for p, label in [(p1, "1\nref"), (p2, "2"), (p3, "3")]:
        circ = plt.Circle(p, 0.075, fc="#eff6ff", ec="#1d4ed8", lw=1.2)
        ax.add_patch(circ)
        ax.text(*p, label, ha="center", va="center", weight="bold")
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#475569", lw=1.6)
    ax.plot([p1[0], p3[0]], [p1[1], p3[1]], color="#475569", lw=1.6)
    ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color="#dc2626", lw=2.7)
    ax.annotate(
        r"line action: $w_{23}\leftarrow1.01w_{23}$",
        xy=(0.50, 0.36),
        xytext=(0.52, 0.08),
        ha="center",
        color="#dc2626",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#dc2626", lw=0.8),
        arrowprops=dict(arrowstyle="->", lw=1.0, color="#dc2626"),
    )
    ax.set_title("Physical action")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Panel 2: s-plane.
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter([s0.real], [s0.imag], color="#2563eb", label="base pole", zorder=3)
    ax.scatter([s1.real], [s1.imag], color="#dc2626", label="after line", zorder=3)
    ax.scatter([mu.real], [mu.imag], marker="x", s=70, color="#7c2d12", label="control pole")
    ax.annotate("", xy=(s1.real, s1.imag), xytext=(s0.real, s0.imag),
                arrowprops=dict(arrowstyle="->", lw=1.5, color="#dc2626"))
    ax.plot([s1.real, mu.real], [s1.imag, mu.imag], "--", color="#7c2d12", lw=1)
    ax.set_xlabel(r"$\Re(s)$")
    ax.set_ylabel(r"$\Im(s)$")
    ax.set_title("Hidden margin shrinks")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6, loc="lower left")

    # Panel 3: metrics.
    ax = fig.add_subplot(gs[0, 2])
    vals = [
        100 * (net1["nu_max"] - net0["nu_max"]) / net0["nu_max"],
        100 * (abs(s1 - mu) - abs(s0 - mu)) / abs(s0 - mu),
        100 * (zeta(s1) - zeta(s0)) / zeta(s0),
    ]
    labels = [r"$\Delta\nu_{\max}$", r"$\Delta d^c$", r"$\Delta\zeta$"]
    colors = ["#16a34a", "#dc2626", "#dc2626"]
    ax.bar(labels, vals, color=colors)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_title("Naive screen misses risk")
    ax.set_ylabel("percent change")
    ax.set_ylim(-0.50, 0.86)
    for i, v in enumerate(vals):
        ax.text(i, v + (0.055 if v >= 0 else -0.055), f"{v:.2f}%",
                ha="center", va="bottom" if v >= 0 else "top", fontsize=7)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(FIG / "fig_synthetic_counterexample.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fig_synthetic_counterexample.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
