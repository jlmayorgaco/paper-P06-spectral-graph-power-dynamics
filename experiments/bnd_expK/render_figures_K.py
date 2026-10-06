"""Render only evidenced ExpK figures; never imply absent validation passed."""
from pathlib import Path
import csv
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports" / "experiment_K"
TABLES = REPORT / "tables"
FIGS = REPORT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)


def rows(name):
    path = TABLES / name
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def save(name):
    plt.tight_layout()
    plt.savefig(FIGS / name, dpi=180)
    plt.close()


def f(row, key):
    try:
        return float(row[key])
    except (ValueError, KeyError, TypeError):
        return float("nan")


candidate = REPORT / "Z_K_NOMINAL_FINAL.toml"
if candidate.exists():
    import tomllib
    with candidate.open("rb") as file:
        c = tomllib.load(file)
    buses = np.arange(30, 40)
    plt.figure(figsize=(9, 3.5))
    plt.bar(buses, c["rho"], color="#176b8a")
    plt.ylim(0, 1.05)
    plt.xticks(buses)
    plt.ylabel("SG→GFL fraction ρ")
    plt.title("Frozen nominal candidate support and replacement")
    save("FIG_K01_nominal_support_rho.png")

    fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    ax[0].bar(buses, c["Kp"], color="#176b8a")
    ax[0].set_ylabel("Kp")
    ax[1].bar(buses, c["Ki"], color="#db7c26")
    ax[1].set_ylabel("Ki")
    ax[1].set_xticks(buses)
    ax[1].set_xlabel("Generator bus")
    fig.suptitle("Frozen nominal PLL gains")
    save("FIG_K02_Kp_Ki_by_bus.png")

hist = rows("TABLE_K07_support_branch_history.csv")
if hist:
    valid = [r for r in hist if np.isfinite(f(r, "retained_SG_MW"))]
    plt.figure(figsize=(8, 4))
    plt.scatter([f(r, "retained_SG_MW") for r in valid],
                [f(r, "alpha_s_inv") for r in valid], s=4, alpha=.25,
                color="#176b8a")
    plt.axhline(-.05, color="#b73535", linestyle="--", linewidth=1)
    plt.xlabel("Retained SG MW")
    plt.ylabel("Complete-spectrum abscissa (s⁻¹)")
    plt.title("Analytical support and corrector evaluations")
    save("FIG_K03_alpha_vs_retained_SG_MW.png")

    path = [r for r in hist if r["phase"] != "uniform_seed"]
    if path:
        plt.figure(figsize=(8, 4))
        plt.plot(range(len(path)), [f(r, "alpha_s_inv") for r in path],
                 marker=".", linewidth=.6)
        plt.axhline(-.05, color="#b73535", linestyle="--", linewidth=1)
        plt.xlabel("Accepted corrector step")
        plt.ylabel("Rightmost pole real part (s⁻¹)")
        plt.title("Critical pole real-part history (mode switches allowed)")
        save("FIG_K04_critical_pole_trajectory.png")

branches = rows("TABLE_K02_branch_catalog.csv")
if branches:
    labels = sorted({r["family"] for r in branches})
    xs = [int(r["support_mask"]) for r in branches]
    ys = [labels.index(r["family"]) for r in branches]
    plt.figure(figsize=(10, 3.5))
    plt.scatter(xs, ys, s=6, alpha=.7, color="#176b8a")
    plt.yticks(range(len(labels)), labels, fontsize=7)
    plt.xlabel("Support mask (0–1023)")
    plt.title("Active pole families at the interior reference point")
    save("FIG_K05_active_mode_switch_map.png")

front = rows("TABLE_K09_robust_frontier.csv")
if front:
    good = [r for r in front if np.isfinite(f(r, "GFL_MW"))]
    if good:
        plt.figure(figsize=(8, 4))
        plt.semilogx([max(f(r, "beta_req"), 1e-8) for r in good],
                     [f(r, "GFL_MW") for r in good], "o-", color="#176b8a")
        plt.xlabel("Normalized full-block β requirement")
        plt.ylabel("GFL MW")
        plt.title("Best sampled analytical robust candidates (uncertified)")
        save("FIG_K06_GFL_MW_vs_beta.png")
        shutil.copyfile(FIGS / "FIG_K06_GFL_MW_vs_beta.png",
                        FIGS / "FIG_K01_GFL_MW_vs_beta.png")

bnd = rows("TABLE_K10_BND_explanation.csv")
if bnd:
    buses = [int(r["parameter_bus"]) for r in bnd]
    direct = [f(r, "direct_dlambda_real") for r in bnd]
    selfe = [f(r, "self_energy_dlambda_real") for r in bnd]
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1]})
    x = np.arange(len(buses))
    ax[0].bar(x - .18, direct, width=.36, label="Direct")
    ax[0].bar(x + .18, selfe, width=.36, label="Collective self-energy")
    ax[0].set_ylabel("Component (s⁻¹)")
    ax[0].legend()
    ax[1].bar(x, np.array(direct) + np.array(selfe), color="#47715b")
    ax[1].set_ylabel("Net derivative (s⁻¹)")
    ax[1].set_xticks(x, buses)
    ax[1].set_xlabel("SG bus parameter")
    fig.suptitle("Active-mode Schur sensitivity cancellation")
    save("FIG_K07_active_mode_self_energy.png")
    shutil.copyfile(FIGS / "FIG_K07_active_mode_self_energy.png",
                    FIGS / "FIG_K02_active_mode_self_energy.png")

trace = rows("TABLE_K15_analytic_transient_trace.csv")
if trace:
    t = [f(r, "time_s") for r in trace]
    fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    ax[0].plot(t, [f(r, "frequency_Hz_per_MW") for r in trace])
    ax[0].set_ylabel("Frequency (Hz/MW)")
    ax[1].plot(t, [f(r, "RoCoF_Hz_s_per_MW") for r in trace])
    ax[1].set_ylabel("RoCoF (Hz/s/MW)")
    ax[1].set_xlabel("Time (s)")
    fig.suptitle("Exact analytical unit-step response at the worst load bus")
    save("FIG_K08_analytic_frequency_RoCoF.png")

tds = rows("TABLE_K16_tds_trajectories.csv")
if tds:
    groups = sorted({(r["candidate"], r["event_bus"], r["pulse_fraction"]) for r in tds})
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for label, bus, pulse in groups:
        rr = [r for r in tds if (r["candidate"], r["event_bus"], r["pulse_fraction"]) == (label, bus, pulse)]
        tt = np.array([f(r, "time_s") for r in rr])
        ff = np.array([f(r, "COI_frequency_Hz") for r in rr])
        axes[0].plot(tt, ff, label=f"{label}, bus {bus}, pulse {pulse}")
        finite = np.isfinite(tt) & np.isfinite(ff)
        if finite.sum() > 2:
            dt = np.diff(tt[finite])
            df = np.diff(ff[finite])
            good = dt > 0
            axes[1].plot(tt[finite][1:][good], df[good] / dt[good])
    axes[0].set_ylabel("COI frequency deviation (Hz)")
    axes[1].set_ylabel("RoCoF (Hz/s)")
    axes[1].set_xlabel("Time (s)")
    fig.suptitle("Post-freeze nonlinear PowerDynamics pulses")
    axes[0].legend(fontsize=7)
    save("FIG_K09_TDS_frequency_RoCoF_comparison.png")

match = rows("TABLE_K18_linear_nonlinear_traces.csv")
if match:
    plt.figure(figsize=(9, 4))
    for label in sorted({r["candidate"] for r in match}):
        rr = [r for r in match if r["candidate"] == label and
              r["event_bus"] == "8" and f(r, "pulse_fraction") == .001]
        if not rr:
            continue
        tt = [f(r, "time_s") for r in rr]
        plt.plot(tt, [f(r, "linear_frequency_Hz") for r in rr],
                 linestyle="--", label=f"{label} analytic")
        plt.plot(tt, [f(r, "nonlinear_frequency_Hz") for r in rr],
                 label=f"{label} PowerDynamics")
    plt.xlabel("Time (s)")
    plt.ylabel("COI frequency deviation (Hz)")
    plt.title("Frozen-model linear vs nonlinear bus-8 pulse")
    plt.legend(fontsize=8)
    save("FIG_K10_linear_nonlinear_frequency.png")
