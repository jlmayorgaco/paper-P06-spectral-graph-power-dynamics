from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/pll_sector_contract"
s = json.loads((OUT/"synthesis.json").read_text(encoding="utf-8"))
v = json.loads((OUT/"verification.json").read_text(encoding="utf-8"))
h = np.load(OUT/"nonlinear_histories.npz")
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10,
                     "axes.spines.top":False, "axes.spines.right":False,
                     "axes.titleweight":"bold", "pdf.fonttype":42})
fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6))
colors = ["#85919e", "#007f89"]
theta = np.linspace(0, np.pi/2, 300)
for c, label, color in zip([v["certificates"][0], v["certificates"][-1]],
                            ["Nominal PLL gains", "Retuned PLL gains"], colors):
    xx = 100*c["epsilon_verified"]*np.cos(theta)
    yy = 1000*c["reference_frequency_envelope_hz"]*np.sin(theta)
    axes[0].plot(xx, yy, color=color, lw=2.1, label=label)
    axes[0].fill_between(xx[::-1], yy[::-1], color=color, alpha=.08)
axes[0].set(xlabel="Terminal voltage perturbation norm (% of 1 pu)",
            ylabel="Reference-frequency magnitude (mHz)",
            title="A  Joint input envelope: all waveforms", xlim=(0,.102), ylim=(0,10.4))
axes[0].legend(frameon=False, loc="upper right", fontsize=9)
gain=100*(v["certificates"][-1]["epsilon_verified"]/v["certificates"][0]["epsilon_verified"]-1)
axes[0].text(.006, 1.1, f"+{gain:.2f}% certified radius", color=colors[1], weight="bold")
for key, label, color, style in [
    ("V0.95_rocof_boundary_state_adversary", "Adversarial input; boundary start", "#007f89", "-"),
    ("V0.95_origin_chirp", "Chirp input; equilibrium start", "#bc6826", "--"),
]:
    a=h[key]
    axes[1].plot(a[:,0], a[:,1], label=label, color=color, ls=style, lw=1.8)
axes[1].axhline(1, color="#303b45", lw=1.2, ls=":")
axes[1].text(2.94,1.02,"Certified boundary",ha="right",fontsize=9)
axes[1].set(xlabel="Time (s)", ylabel=r"Storage $V(z)$", title="B  Nonlinear implementation checks",
            xlim=(0,3), ylim=(0,1.17))
axes[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1,.88), fontsize=8.5)
for ax in axes:
    ax.grid(alpha=.18)
fig.suptitle("A verified nonlinear PLL contract", x=.075, ha="left", fontsize=17, weight="bold")
fig.text(.075,.91,"Conditional component result; current, DC, SG and network interconnection remain uncertified.",
         fontsize=10, color="#52606d")
fig.text(.075,.025,"Guarantee: sine-sector proof + exact rational matrix checks. Time traces audit the implementation.",
         fontsize=9, color="#52606d")
fig.subplots_adjust(left=.075,right=.98,bottom=.18,top=.78,wspace=.28)
for ext in ("png","pdf","svg"):
    fig.savefig(OUT/f"pll_contract.{ext}",dpi=200,facecolor="white")
print(OUT/"pll_contract.png")
