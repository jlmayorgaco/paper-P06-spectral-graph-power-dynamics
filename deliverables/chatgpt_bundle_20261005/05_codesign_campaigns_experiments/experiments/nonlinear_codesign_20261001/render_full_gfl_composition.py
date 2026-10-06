"""Plot verified component results and the explicit algebraic closure witness."""
from pathlib import Path
import csv
import hashlib
import json
import tomllib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/full_gfl_sector_contract"
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
ver = json.loads((OUT / "verification_final.json").read_text(encoding="utf-8"))
wit = json.loads((OUT / "network_closure_witness.json").read_text(encoding="utf-8"))
net = tomllib.loads((OUT / "network_closure.toml").read_text(encoding="utf-8"))
assert wit["verification_sha256"] == sha(OUT / "verification_final.json")
assert wit["network_closure_sha256"] == sha(OUT / "network_closure.toml")
with (OUT / "joint_closure_witness.csv").open(encoding="utf-8", newline="") as f:
    states = list(csv.DictReader(f))
buses = net["buses"]
eps = np.array([c["epsilon_actual_model_verified"] for c in ver["certificates"]])
assert [c["bus"] for c in ver["certificates"]] == buses
Z = np.array(net["current_to_voltage"])
dv = []
for k, row in enumerate(states):
    assert int(row["donor_bus"]) == buses[k]
    t = net["theta0"][k]
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    injected = net["rho"][k] * rot @ np.array([float(row["i_d"]), float(row["i_q"])])
    dv.append((Z[:, 2*k:2*k+2] @ injected).reshape(-1, 2))
alone = np.array([np.max(np.linalg.norm(v, axis=1) / eps) for v in dv])
joint_by_bus = np.linalg.norm(sum(dv), axis=1) / eps
receiver = buses.index(wit["joint_witness"]["receiver_bus"])
joint = joint_by_bus[receiver]
assert abs(alone.max() - wit["joint_witness"]["largest_single_input_usage"]) < 1e-10
assert abs(joint - wit["joint_witness"]["joint_input_usage"]) < 1e-10

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.4), gridspec_kw={"width_ratios": [1, 1.45]})
fig.subplots_adjust(left=.067, right=.98, bottom=.23, top=.75, wspace=.28)
fig.suptitle("Local nonlinear guarantees need a network compatibility condition",
             fontsize=16, fontweight="bold", x=.067, ha="left", y=.98)
fig.text(.067, .88, "Full GFL model: PLL + current loop + DC dynamics; filtered RoCoF instrument.",
         fontsize=11, color="#415466")
x = np.arange(len(buses))
ax[0].bar(x, eps * 1000, color="#197c80", width=.68)
ax[0].set(xticks=x, xticklabels=buses, ylim=(0, 1.04), xlabel="GFL bus (test instance)",
          ylabel=r"Verified joint-input scale $\epsilon_i$ ($10^{-3}$ pu)")
ax[0].set_title("A  Ten conditional invariant regions", loc="left", fontsize=11, pad=17)
ax[0].text(.02, .94, "32 exact-rational vertex checks per GFL", transform=ax[0].transAxes,
           fontsize=9, color="#415466", va="top")
xx = np.r_[x, 11.5]
ax[1].bar(xx, np.r_[alone, joint], color=["#197c80"] * 10 + ["#d75b32"], width=.67)
ax[1].axhline(1, ls="--", lw=1.2, color="#283b4a")
ax[1].text(11.5, joint + .07, f"{joint:.3f}", ha="center", fontsize=11, fontweight="bold")
ax[1].text(-.42, 1.09, "Input assumption boundary", fontsize=9, color="#283b4a")
ax[1].set(xticks=xx, xticklabels=[str(b) for b in buses] + ["Joint"], ylim=(0, 3.2),
          xlabel="Each state variation alone | All ten together",
          ylabel=r"Voltage-input usage $\|\Delta v_i\|/\epsilon_i$")
ax[1].set_title("B  Explicit network closure counterexample", loc="left", fontsize=11, pad=17)
ax[1].text(.03, .96, f"Each alone: largest receiver usage ≈ {alone.max():.3f}\n"
                    f"Joint: receiver {buses[receiver]}; every local storage < 0.020",
           transform=ax[1].transAxes, fontsize=9, color="#415466", va="top")
for a in ax:
    a.set_axisbelow(True)
    a.grid(axis="y", alpha=.13)
fig.text(.067, .10, "An exact algebraic slice of the nonlinear DAE violates the assumed input envelopes of the product region.",
         fontsize=10, color="#283b4a")
fig.text(.067, .054, "Scope: a failure of this certificate composition. Instability and reachability from nominal operation are not established.",
         fontsize=9, color="#415466")
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT / f"full_gfl_composition.{ext}", dpi=200, facecolor="white")
data = dict(script_sha256=sha(Path(__file__)),
            source_hashes={f: sha(OUT / f) for f in ("verification_final.json",
                "network_closure_witness.json", "network_closure.toml", "joint_closure_witness.csv")},
            each_alone_maximum_receiver_usage=alone.tolist(), joint_receiver_usage=float(joint),
            joint_receiver_bus=buses[receiver], joint_all_receiver_usage=joint_by_bus.tolist())
(OUT / "figure_data.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"figure": str(OUT / "full_gfl_composition.png"),
                  "largest_single_usage": float(alone.max()), "joint_usage": float(joint)}))
