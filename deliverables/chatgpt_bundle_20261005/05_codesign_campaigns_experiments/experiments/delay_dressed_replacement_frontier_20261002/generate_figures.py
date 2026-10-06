"""Publication-style evidence figures; blocked outputs are shown as not estimated."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = Path(__file__).resolve().parent
BASE = HERE / "baseline_reproduction"
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "axes.titlesize":12,
                     "figure.dpi":150, "savefig.dpi":220, "axes.spines.top":False,
                     "axes.spines.right":False})
BLUE="#174A6E"; TEAL="#28A6A1"; AMBER="#E39C35"; RED="#BD4B4B"; GREY="#657786"

def save(fig, name):
    fig.savefig(HERE/name, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def blocked(name, title, reason, evidence):
    fig, ax = plt.subplots(figsize=(9,4.5)); ax.axis("off")
    ax.text(.5,.83,title,ha="center",va="center",fontsize=17,weight="bold",color=BLUE,transform=ax.transAxes)
    ax.text(.5,.59,"NOT ESTIMATED",ha="center",va="center",fontsize=22,weight="bold",color=RED,transform=ax.transAxes)
    ax.text(.5,.39,reason,ha="center",va="center",fontsize=12,color="#263746",wrap=True,transform=ax.transAxes)
    ax.text(.5,.13,evidence,ha="center",va="center",fontsize=9,color=GREY,wrap=True,transform=ax.transAxes)
    save(fig,name)

# F1: evidence-gated pipeline, not a claimed scientific mechanism.
fig,ax=plt.subplots(figsize=(13,3.6)); ax.set_xlim(0,13); ax.set_ylim(0,3.5); ax.axis("off")
labels=[("IEEE-39 nonlinear\nDAE + physical graph","D0 PASS",TEAL),
        ("Exact PLL-error\nretarded DDE","\u03c4=0 identity PASS",TEAL),
        ("Rightmost DDE\nroot completeness","BLOCKED",RED),
        ("Delay co-design +\nall-event validation","NOT RUN",GREY),
        ("Replacement frontier +\nbounds","NOT ESTABLISHED",GREY)]
xs=[.3,2.85,5.4,7.95,10.5]
for x,(lab,st,c) in zip(xs,labels):
    p=FancyBboxPatch((x,.95),2.15,1.5,boxstyle="round,pad=0.08,rounding_size=.12",fc="white",ec=c,lw=2,transform=ax.transData)
    ax.add_patch(p);ax.text(x+1.075,1.8,lab,ha="center",va="center",fontsize=10,weight="bold",color=BLUE)
    ax.text(x+1.075,1.12,st,ha="center",va="center",fontsize=9,color=c,weight="bold")
for x in [2.48,5.03,7.58,10.13]:
    ax.add_patch(FancyArrowPatch((x,1.7),(x+.32,1.7),arrowstyle="-|>",mutation_scale=14,color=GREY,lw=1.5))
ax.text(6.5,3.06,"Experiment pipeline and present evidence gate",ha="center",fontsize=16,weight="bold",color=BLUE)
ax.text(6.5,.35,"No delay-dependent maximum or feasibility claim crosses the blocked root-completeness gate.",ha="center",fontsize=10,color=GREY)
save(fig,"FIG_D01_THEORY_PIPELINE.png")

# F2: 27 analytic-vs-centered-FD local root sensitivities.
d=pd.read_csv(BASE/"TABLE_D02_DERIVATIVE_VALIDATION.csv")
fig,axs=plt.subplots(1,3,figsize=(12,3.8))
for ax,(param,sub) in zip(axs,d.groupby("parameter",sort=False)):
    lo=min(sub.analytic_real.min(),sub.fd_real.min(),sub.analytic_imag.min(),sub.fd_imag.min())
    hi=max(sub.analytic_real.max(),sub.fd_real.max(),sub.analytic_imag.max(),sub.fd_imag.max())
    ax.scatter(sub.analytic_real,sub.fd_real,s=28,alpha=.8,label="real")
    ax.scatter(sub.analytic_imag,sub.fd_imag,s=28,alpha=.8,label="imag",marker="s")
    span=max(hi-lo,1e-14);pad=.06*span;ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],"--",color=GREY,lw=1)
    ax.set_title(f"{param}: max rel.err {sub.relative_error.max():.2g}")
    ax.set_xlabel("analytic derivative");ax.grid(alpha=.2)
axs[0].set_ylabel("centered finite difference");axs[0].legend(frameon=False)
fig.suptitle("Local simple-root sensitivities: 3 delay patterns × 3 buses",color=BLUE,weight="bold")
fig.tight_layout();save(fig,"FIG_D02_ANALYTIC_VS_FD.png")

# F3-F7 and F9 are explicit not-estimated states, not simulated data.
blocked("FIG_D03_OPERATOR_APPROXIMATION_ERROR.png","Simple PLL graph-operator approximation",
        "No justified Bc/Kp/C reduced map was identified for the full lossy IEEE-39 model.",
        "Exact Schur-slope computations are reported separately; this figure does not substitute them for the requested approximation test.")
blocked("FIG_D04_REPLACEMENT_VS_GRAPH_MODES.png","Replacement versus graph-mode gain dimension",
        "q = 0…5 and full-nodal delayed co-design were not run.",
        "The exact DDE rightmost-root completeness gate is blocked; no replacement values exist to plot.")
blocked("FIG_D05_REPLACEMENT_FRONTIER_VS_DELAY.png","Maximum validated replacement versus delay",
        "P_GFL^max(τ) was not established.",
        "A candidate must pass complete delayed spectrum and all six frozen nonlinear events; neither condition was demonstrated for delayed designs.")
blocked("FIG_D06_SAME_DELAY_DIFFERENT_PLACEMENT.png","Same delay multiset, different bus placement",
        "104 placements were frozen, but no placement-specific maximum replacement was computed.",
        "The graph descriptors are available in TABLE_D05_GSP_DELAY_METRICS.csv; they cannot be correlated with an unavailable frontier.")
blocked("FIG_D07_DELAY_SHADOW_PRICE_MAP.png","Marginal replacement cost of site latency",
        "No validated delayed optimum or stable KKT active set exists.",
        "MW/ms shadow prices are therefore not estimated.")
blocked("FIG_D09_OPTIMALITY_GAP.png","Feasible lower frontier and optimistic upper bound",
        "Neither R_F nor a rigorous R_U is available.",
        "Do not interpret a sampled or locally optimized point as a certified bound.")

# F8: D0 parity residuals normalized by their preregistered tolerances.
p=pd.read_csv(BASE/"TABLE_D00_FULL_SPECTRUM_PARITY.csv")
eq=pd.read_csv(BASE/"TABLE_D00_BASELINE_REPRODUCTION.csv")
metrics=[("equilibrium",eq.equilibrium_rhs_inf.max(),1e-8),
         ("pole match",p.max_spectrum_match_error.max(),1e-5),
         ("rightmost pole",p.alpha_abs_error.max(),1e-6),
         ("frequency peak",p.F_abs_error_Hz.max(),1e-4),
         ("RoCoF peak",p.RoCoF_abs_error_Hz_s.max(),1e-4),
         ("voltage extrema",max(p.Vmin_abs_error_pu.max(),p.Vmax_abs_error_pu.max()),1e-4)]
fig,ax=plt.subplots(figsize=(9,4.5));vals=[m[1]/m[2] for m in metrics]
ax.bar([m[0] for m in metrics],vals,color=[TEAL if v<=1 else RED for v in vals]);ax.axhline(1,color=RED,ls="--",lw=1,label="declared parity tolerance")
ax.set_yscale("log");ax.set_ylabel("maximum error / absolute tolerance");ax.set_title("No-delay ReducedDAE ↔ PowerDynamics parity (D0 PASS)",color=BLUE,weight="bold")
ax.tick_params(axis="x",rotation=20);ax.grid(axis="y",alpha=.2);ax.legend(frameon=False);fig.tight_layout();save(fig,"FIG_D08_FULL_MODEL_VALIDATION.png")

# Additional actual evidence figures.
t=pd.read_csv(HERE/"TABLE_D03_TAYLOR_NETWORK_VALIDATION.csv")
fig,ax=plt.subplots(figsize=(7.5,4.6))
for col,label,color,marker in [("first_order_relative_error","1st order",BLUE,"o"),("second_order_relative_error","2nd order",AMBER,"s"),("third_order_relative_error","3rd order",TEAL,"^")]:
    med=t.groupby("amplitude_rad")[col].median();ax.loglog(med.index,med.values,marker=marker,label=label,color=color,lw=2)
ax.set_xlabel("angle perturbation amplitude [rad]");ax.set_ylabel("median relative branch-power error");ax.set_title("Lossless branch-power Taylor validation (60 tests)",color=BLUE,weight="bold")
ax.grid(which="both",alpha=.2);ax.legend(frameon=False);fig.tight_layout();save(fig,"FIG_D10_TAYLOR_VALIDATION.png")

r=pd.read_csv(BASE/"TABLE_D01_LOCAL_DDE_ROOTS.csv")
fig,ax=plt.subplots(figsize=(7.5,4.4));ax.scatter(r.chi_tau,r.tracked_pair_real_shift_from_tau0*1e9,s=70,c=[AMBER,TEAL,BLUE][:len(r)])
for _,row in r.iterrows(): ax.annotate(row.delay_pattern_id,(row.chi_tau,row.tracked_pair_real_shift_from_tau0*1e9),xytext=(5,5),textcoords="offset points",fontsize=8)
ax.axhline(0,color=GREY,lw=1);ax.set_xlabel("$\\chi_\\tau$ (lossless graph)" );ax.set_ylabel("tracked-pair real-part shift [10$^{-9}$ s$^{-1}$]")
ax.set_title("Local tracked roots only — not a complete DDE spectrum",color=RED,weight="bold");ax.grid(alpha=.2);fig.tight_layout();save(fig,"FIG_D11_LOCAL_ROOT_TRACKING.png")

blocked("FIG_D12_SCHUR_SLOPE_DIAGNOSTIC.png","Gauge-deflated exact Schur slope",
        "The hidden block is singular at s = 0, so D_G and L_tau are undefined for this partition.",
        "Pre-gauge numbers are preserved as a discarded diagnostic and are excluded from interpretation.")

print("generated evidence plots and explicitly blocked frontier panels")
