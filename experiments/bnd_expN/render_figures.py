"""Render Experiment N evidence figures from frozen tables."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import solve

from compare_validation_cases import load_admittance, transfer

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
CASES = OUT / "validation_cases"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 180})


def matrix(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, delimiter=","))


def read(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)


def port_error_curve(case: str, bus: int) -> tuple[np.ndarray, np.ndarray]:
    path = CASES / case
    anmap = read(path / "AN_state_map.csv")
    ids = np.array([int(r["index"])-1 for r in anmap if int(r["bus"]) == bus])
    yi = np.array([2*bus-2, 2*bus-1])
    a, b, c, d = [matrix(path / f"AN_{n}.csv") for n in
                  ("A_local", "B_port", "C_port", "D_port")]
    mp, ap, bp, cp, dp = [matrix(path / f"PD_bus{bus}_{n}.csv") for n in
                          ("M", "A", "B", "C", "D")]
    original = {int(r["bus"]): r for r in read(
        OUT / "TABLE_N01_original_operating_point.csv")}
    adm = load_admittance(bus, original)
    omega = np.geomspace(1e-4, 1e3, 51)
    errors = []
    for w in omega:
        s = 1j*w
        yan = transfer(np.eye(len(ids)), a[np.ix_(ids, ids)],
                       b[np.ix_(ids, yi)], c[np.ix_(yi, ids)],
                       d[np.ix_(yi, yi)], s) + adm
        ypd = -np.linalg.inv(transfer(mp, ap, bp, cp, dp, s))
        errors.append(np.linalg.norm(ypd-yan)/max(np.linalg.norm(yan), 1e-15))
    return omega, np.array(errors)


def figure_02() -> None:
    fig, ax = plt.subplots(figsize=(6.8, 3.7))
    for case, bus in [("bus31_r099_nom",31),("bus38_r099_nom",38)]:
        w, err = port_error_curve(case, bus)
        ax.loglog(w, err, marker=".", markersize=3, label=f"bus {bus}, ρ=0.99")
    ax.axhline(1e-8, color="crimson", linestyle="--", label="identity gate")
    ax.set(xlabel="Angular frequency ω (rad/s)",
           ylabel="Relative GFL/mixed port error",
           title="Independent PD vs analytical port identity")
    ax.legend(loc="upper right")
    save(fig,"FIG_N02_GFL_port_error_vs_frequency.png")


def figure_03() -> None:
    path=OUT/"postfreeze_pd"
    pdp=pd.read_csv(path/"PD_poles.csv")
    anp=pd.read_csv(path/"AN_poles.csv")
    pdp=pdp.drop(index=np.argmin(np.hypot(pdp.real,pdp.imag)))
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for ax in axs:
        ax.scatter(anp.real,anp.imag,s=17,facecolors="none",edgecolors="#0072B2",
                   label="analytic")
        ax.scatter(pdp.real,pdp.imag,s=7,c="#D55E00",alpha=.6,label="PD")
        ax.axvline(-.05,color="black",linestyle="--",linewidth=.8)
        ax.set(xlabel="Real part (s⁻¹)",ylabel="Imaginary part (s⁻¹)")
    axs[0].set_title("Complete finite spectrum")
    axs[1].set_xlim(-.3,.1);axs[1].set_ylim(-3,3)
    axs[1].set_title("Critical-region detail")
    axs[1].legend()
    save(fig,"FIG_N03_PD_vs_analytic_poles.png")


def figure_04() -> None:
    d=pd.read_csv(OUT/"TABLE_N09_single_support_scan.csv")
    d=d[(d.gainset=="lowKp_highKi") & (d.retained_fraction>0)]
    fig,ax=plt.subplots(figsize=(7.4,4.2))
    for bus in (31,32,34,36,38,39):
        x=d[d.bus==bus].sort_values("retained_SG_MW")
        ax.semilogx(x.retained_SG_MW,x.alpha,marker="o",markersize=2,
                    label=f"SG bus {bus}")
    ax.axhline(-.05,color="black",linestyle="--",linewidth=.9,label="required margin")
    ax.set(xlabel="Retained SG dispatch (MW)",ylabel="Complete spectral abscissa (s⁻¹)",
           ylim=(-.18,.16),title="Support-aware spectral continuation")
    ax.legend(ncol=2,fontsize=7)
    save(fig,"FIG_N04_spectral_abscissa_continuation.png")


def figure_05() -> None:
    d=pd.read_csv(OUT/"TABLE_N10_all_KKT_candidates.csv")
    d=d[d.gainset=="lowKp_highKi"].sort_values("bus")
    fig,ax=plt.subplots(figsize=(7.4,3.8))
    ax.bar(d.bus.astype(str),d.retained_SG_MW,
           color=["#D55E00" if b==38 else "#7AA6C2" for b in d.bus])
    ax.set_yscale("log")
    ax.set(xlabel="Single retained-SG support",ylabel="Minimum corrected retained MW",
           title="Fixed-support active-set boundary candidates")
    save(fig,"FIG_N05_retained_SG_MW_vs_active_branch.png")


def figure_06() -> None:
    d=pd.read_csv(OUT/"TABLE_N09_single_support_scan.csv")
    d=d[(d.gainset=="lowKp_highKi") & (d.retained_fraction>0)]
    fig,axs=plt.subplots(1,2,figsize=(10,3.8))
    for bus in (30,33,36,38):
        x=d[d.bus==bus].sort_values("retained_fraction")
        axs[0].semilogx(x.retained_fraction,x.rightmost_real,marker=".",label=str(bus))
        axs[1].semilogx(x.retained_fraction,abs(x.rightmost_imag),marker=".",label=str(bus))
    axs[0].axhline(-.05,color="black",linestyle="--",linewidth=.8)
    axs[0].set(xlabel="Retained fraction ε",ylabel="Re λ (s⁻¹)",
               title="Rightmost eigenvalue trajectory")
    axs[1].set(xlabel="Retained fraction ε",ylabel="|Im λ| (s⁻¹)",
               title="Mode-family switches")
    axs[1].legend(title="SG bus",fontsize=7)
    save(fig,"FIG_N06_active_eigenvalue_trajectories.png")


def figure_07() -> None:
    d=pd.read_csv(OUT/"TABLE_N15_self_energy_explanation.csv")
    choices=[("rho",38),("Kp",36),("Kp",38),("Ki",36),("Ki",38)]
    q=pd.concat([d[(d.parameter==p)&(d.parameter_bus==b)] for p,b in choices])
    labels=[f"{p}{b}" for p,b in choices]
    direct=np.abs(q.direct_dlambda_real.to_numpy())
    collective=np.abs(q.self_energy_dlambda_real.to_numpy())
    total=np.maximum(direct+collective,1e-30)
    fig,ax=plt.subplots(figsize=(6.5,3.7))
    x=np.arange(len(labels))
    ax.bar(x,direct/total,label="direct/local",color="#0072B2")
    ax.bar(x,collective/total,bottom=direct/total,label="collective/self-energy",
           color="#D55E00")
    ax.set_xticks(x,labels)
    ax.set(ylabel="Share of |direct| + |self-energy|",ylim=(0,1.05),
           title="Active-mode sensitivity at BND pivot bus 36")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5,-0.33), ncol=2)
    save(fig,"FIG_N07_self_energy_direct_vs_collective.png")


def figure_08() -> None:
    d=pd.read_csv(OUT/"TABLE_N14_trajectories.csv")
    fig,axs=plt.subplots(2,1,figsize=(8,5.7),sharex=True)
    colors={8:"#D55E00",16:"#0072B2",29:"#009E73"}
    for bus in (8,16,29):
        for pulse,ls in [(1e-5,"-"),(2e-5,"--")]:
            q=d[(d.event_bus==bus)&np.isclose(d.pulse_fraction,pulse)]
            scale=q.delta_P_MW.iloc[0]
            axs[0].plot(q.time_s,q.COI_frequency_deviation_Hz/scale,
                        ls,color=colors[bus],linewidth=1.1,label=f"bus {bus}, {pulse:g}")
            axs[1].plot(q.time_s,q.RoCoF_Hz_s/scale,
                        ls,color=colors[bus],linewidth=1.1)
    axs[0].set(ylabel="COI Δf / ΔP (Hz/MW)",
               title="Nonlinear PD pulse response and amplitude scaling")
    axs[1].set(xlabel="Time (s)",ylabel="RoCoF / ΔP (Hz/s/MW)")
    axs[0].legend(ncol=3,fontsize=7)
    for ax in axs:
        ax.axvspan(1.0,1.1,color="grey",alpha=.15)
    save(fig,"FIG_N08_nonlinear_frequency_RoCoF_TDS.png")


def main() -> None:
    for fn in (figure_02,figure_03,figure_04,figure_05,
               figure_06,figure_07,figure_08):
        fn()
    print("N_FIGURES_DONE 7; FIG_N01 pre-existing")


if __name__ == "__main__":
    main()
