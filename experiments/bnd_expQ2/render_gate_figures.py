from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_Q2" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 160})

# Q2-01: continuous architecture-invariant grid-frequency measurement as a
# vanishing GFL share is introduced at generator bus 38.
f0 = pd.read_csv(ROOT / "reports/experiment_Q2/F0/TABLE_Q2_F0_frequency_metric_continuity.csv")
e0 = pd.read_csv(ROOT / "reports/experiment_Q2/F0/TABLE_Q2_F0_continuity_errors.csv")
b38 = f0[f0.bus == 38].copy()
e38 = e0[e0.bus == 38].copy()
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.2), constrained_layout=True)
for h, label, color in [(2,"0.04 s window","#0072B2"),(5,"0.10 s window","#D55E00"),
                        (10,"0.20 s window","#009E73"),(20,"0.40 s window","#CC79A7")]:
    d = b38[b38.half_window == h].sort_values("rho_GFL")
    ax[0].plot(d.rho_GFL, 100*d.F_peak_Hz, marker="o", ms=3, label=label, color=color)
ax[0].set_xscale("symlog", linthresh=1e-8)
ax[0].set_xlabel("GFL dispatch share $\\rho_{38}$")
ax[0].set_ylabel("100 MW grid-frequency peak (Hz)")
ax[0].set_title("Same passive bus-voltage measurement")
ax[0].ticklabel_format(axis="y",style="plain",useOffset=False)
ax[0].legend(frameon=False, fontsize=8)
for key, label, color in [("SAVGOL_HALF_WINDOW_5:F_peak_Hz","peak frequency","#0072B2"),
                          ("SAVGOL_HALF_WINDOW_5:R_peak_Hz_s","RoCoF","#D55E00"),
                          ("SAVGOL_HALF_WINDOW_5:F_inf_Hz","steady offset","#009E73")]:
    d=e38[e38.estimator == key].sort_values("rho_GFL")
    ref=float(d.metric_reference.iloc[0])
    ax[1].plot(d.rho_GFL, (d.metric_rho-ref).abs()/max(abs(ref),1e-30),
               marker="o", ms=3, label=label, color=color)
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlabel("GFL dispatch share $\\rho_{38}$")
ax[1].set_ylabel("Relative change from $\\rho=0$")
ax[1].set_title("Vanishing-share continuity (0.10 s estimator)")
ax[1].legend(frameon=False, fontsize=8)
fig.suptitle("FIG_Q2_01 — Architecture-invariant grid-frequency continuity")
fig.savefig(OUT / "FIG_Q2_01_frequency_metric_continuity.png", dpi=220)
plt.close(fig)

# Q2-02: the frozen ExpN candidate is spectrally stable but fails grid-frequency
# excursion and every finite RoCoF bandwidth tested; ideal unfiltered RoCoF has
# a distributional impulse due to its nonzero algebraic phase jump.
base = pd.read_json(ROOT / "reports/experiment_Q2/BASELINE/BASELINE_RESULTS.json", typ="series")
f2 = pd.read_csv(ROOT / "reports/experiment_Q2/F2/TABLE_Q2_F2_ExpN_bandwidth_sensitivity.csv")
fig, ax = plt.subplots(2, 2, figsize=(10.2, 7.0), constrained_layout=True)
alpha=float(base.ExpN_alpha_analytic)
ax[0,0].bar(["ExpN $\\alpha$"],[alpha],color="#009E73")
ax[0,0].axhline(-0.05,color="#222222",ls="--",label="required $-0.05$ s$^{-1}$")
ax[0,0].set_ylabel("s$^{-1}$");ax[0,0].set_title("Complete-spectrum margin: PASS")
ax[0,0].legend(frameon=False,fontsize=8)

labels=["sampled 100 Hz","0.04 s","0.10 s","0.20 s","0.40 s"]
vals=100*f2.R_peak_Hz_s.to_numpy()
colors=["#D55E00" if x>0.5 else "#009E73" for x in vals]
ax[0,1].bar(labels,vals,color=colors)
ax[0,1].axhline(.5,color="#222222",ls="--",label="project limit 0.5 Hz/s")
ax[0,1].set_ylabel("grid RoCoF (Hz/s), 100 MW")
ax[0,1].set_title("RoCoF decision varies with estimator")
ax[0,1].tick_params(axis="x",rotation=25);ax[0,1].legend(frameon=False,fontsize=8)
ax[0,1].text(.02,.98,"Ideal continuous unfiltered RoCoF: impulse",transform=ax[0,1].transAxes,
              va="top",fontsize=8,color="#8B0000")

fpeak=float(f2.F_peak_Hz.iloc[2]);finf=float(f2.F_inf_Hz.iloc[2])
ax[1,0].bar(["peak","steady offset"],[fpeak,finf],color=["#D55E00","#D55E00"])
ax[1,0].axhline(.5,color="#222222",ls="--",label="project limit 0.5 Hz")
ax[1,0].set_ylabel("grid-frequency deviation (Hz)")
ax[1,0].set_title("0.10 s estimator: frequency requirement FAIL")
ax[1,0].legend(frameon=False,fontsize=8)

beta_req=1.6991206999182038e-6; beta_w=float(base.ExpN_beta_pointwise_upper)
ax[1,1].bar(["pointwise $\\beta_*$","required $\\beta$"],[beta_w,beta_req],
             color=["#D55E00","#777777"])
ax[1,1].set_yscale("log");ax[1,1].set_ylabel("normalized uncertainty radius")
ax[1,1].set_title("Frozen ExpN robustness witness: FAIL")
fig.suptitle("FIG_Q2_02 — Stable is not secure for the declared event")
fig.savefig(OUT / "FIG_Q2_02_ExpN_stable_not_frequency_secure.png", dpi=220)
plt.close(fig)
print("WROTE", OUT / "FIG_Q2_01_frequency_metric_continuity.png")
print("WROTE", OUT / "FIG_Q2_02_ExpN_stable_not_frequency_secure.png")
