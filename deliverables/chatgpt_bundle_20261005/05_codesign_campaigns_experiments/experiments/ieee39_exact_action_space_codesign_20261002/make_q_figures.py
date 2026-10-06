from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT = HERE
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.titlesize":12,"axes.labelsize":10,"figure.dpi":150,"savefig.dpi":240})

# Q3: preserve the sampled reduced-determinant locus as an explicitly unverified diagnostic.
d = pd.read_csv(HERE/"TABLE_Q03_REDUCED_NYQUIST_DIAGNOSTIC.csv")
fig, axs = plt.subplots(1,2,figsize=(10.4,4.6),constrained_layout=True)
for ax, tau in zip(axs,[20.0,40.0]):
    q=d[d.tau_ms==tau].sort_values("contour_index")
    zlog=q.D_logabs.to_numpy(float)
    ph=q.D_phase.to_numpy(float)
    mag=np.exp(np.clip(zlog-np.nanmax(zlog),-700,0))
    z=mag*np.exp(1j*ph)
    ax.plot(z.real,z.imag,color="#245b8f",lw=.75,alpha=.6)
    ax.scatter(z.real,z.imag,s=3,c=np.arange(len(z)),cmap="viridis",rasterized=True)
    ax.scatter([z.real[0]],[z.imag[0]],s=36,marker="o",facecolor="#f0a202",edgecolor="black",zorder=5,label="sampled start")
    ax.scatter([0],[0],s=34,marker="x",color="#bb2f3b",zorder=6,label="origin")
    ax.axhline(0,color="0.75",lw=.6);ax.axvline(0,color="0.75",lw=.6)
    ax.set_aspect("equal",adjustable="datalim")
    ax.set_title(f"Uniform delay {tau:.0f} ms")
    ax.set_xlabel("Re{Dτ}/max|Dτ| (rescaled for display)")
    ax.set_ylabel("Im{Dτ}/max|Dτ|")
    ax.legend(frameon=False,fontsize=8,loc="best")
fig.suptitle("Sampled reduced determinant locus — diagnostic only",fontweight="bold")
fig.text(.5,-.02,"Float64 contour samples; winding indices are unverified and are not accepted as DDE root counts.",ha="center",fontsize=9,color="#9c2834")
fig.savefig(OUT/"FIG_Q01_REDUCED_NYQUIST.png",bbox_inches="tight");plt.close(fig)

# Q4: the exact conditional equation places a prescribed target. The branch-local Newton solve starts there.
q=pd.read_csv(HERE/"TABLE_Q04_PI_CLOSED_FORM_VALIDATION.csv")
q["gain_norm"]=np.hypot(q.Kp_closed_form,q.Ki_closed_form)
q["actual_local_frequency_Hz"]=q.target_frequency_Hz+q.root_imag_boundary_error/(2*np.pi)
q["inside_bounds"]=q.gains_within_frozen_bounds.astype(str).str.lower().eq("true")
colors={0.0:"#2b6cb0",20.0:"#dd6b20",40.0:"#319795"}
fig,axs=plt.subplots(1,2,figsize=(10.4,4.5),constrained_layout=True)
ax=axs[0]
for tau,g in q.groupby("tau_ms"):
    ax.scatter(g.target_frequency_Hz,g.actual_local_frequency_Hz,s=25,alpha=.75,color=colors[tau],label=f"τ={tau:g} ms")
lo=min(q.target_frequency_Hz.min(),q.actual_local_frequency_Hz.min());hi=max(q.target_frequency_Hz.max(),q.actual_local_frequency_Hz.max())
ax.plot([lo,hi],[lo,hi],"k--",lw=1,label="identity")
ax.set(xlabel="Prescribed target frequency (Hz)",ylabel="Local Newton endpoint frequency (Hz)",title="Conditional target-root check")
ax.legend(frameon=False,fontsize=8)
ax.text(.03,.97,"Newton initialized at prescribed target; not a nearest-pole search",transform=ax.transAxes,va="top",fontsize=8,color="#9c2834",wrap=True)
ax=axs[1]
for tau,g in q.groupby("tau_ms"):
    for inside,gg in g.groupby("inside_bounds"):
        ax.scatter(gg.sigma_min_G,gg.gain_norm,s=31,alpha=.8,color=colors[tau],marker="o" if inside else "x")
handles=[Line2D([0],[0],marker="o",color="0.3",lw=0,label="inside both gain bounds"),Line2D([0],[0],marker="x",color="0.3",lw=0,label="outside gain box")]
for tau,c in colors.items(): handles.append(Line2D([0],[0],marker="o",color=c,lw=0,label=f"τ={tau:g} ms"))
ax.set_yscale("log")
ax.set(xlabel="PI authority σmin(G)",ylabel="Conditional gain norm √(Kp²+Ki²)",title="Selected-grid conditioning and gain size")
ax.legend(handles=handles,frameon=False,fontsize=7,ncol=2)
fig.suptitle("Closed-form PI boundary law: conditional validation",fontweight="bold")
fig.savefig(OUT/"FIG_Q02_PREDICTED_VS_ACTUAL_BOUNDARY.png",bbox_inches="tight");plt.close(fig)

# Q5: show that the exact quadratic candidates do not yield physical retention in the tested grid.
r=pd.read_csv(HERE/"TABLE_Q05_RHO_CLOSED_FORM_VALIDATION.csv")
fig,axs=plt.subplots(1,2,figsize=(10.2,4.6),constrained_layout=True)
for ax in axs:
    sc=ax.scatter(r.epsilon_root,r.epsilon_imag,c=r.target_frequency_Hz,cmap="viridis",s=19,alpha=.78,edgecolor="none")
    ax.axhline(0,color="black",lw=.8);ax.axvline(0,color="black",lw=.8);ax.axvline(1,color="black",lw=.8)
    ax.axvspan(0,1,color="#4caf50",alpha=.11,label="physical real ε ∈ [0,1]")
    ax.set_xlabel("Re{ε root}");ax.set_ylabel("Im{ε root}")
axs[0].set_title("All conditional roots")
axs[1].set_title("Physical interval detail")
axs[1].set_xlim(-1.5,1.25);axs[1].set_ylim(-1.5,1.5)
cb=fig.colorbar(sc,ax=axs);cb.set_label("Target frequency (Hz)")
axs[1].text(.03,.97,"0 / 120 roots are real and physical",transform=axs[1].transAxes,va="top",fontsize=9,color="#9c2834")
axs[1].legend(frameon=False,loc="lower right",fontsize=8)
fig.suptitle("Conditional rank-two SG-retention roots",fontweight="bold")
fig.savefig(OUT/"FIG_Q03_ANALYTIC_RHO_BOUNDARY.png",bbox_inches="tight");plt.close(fig)
print("wrote Q figures 01-03")
