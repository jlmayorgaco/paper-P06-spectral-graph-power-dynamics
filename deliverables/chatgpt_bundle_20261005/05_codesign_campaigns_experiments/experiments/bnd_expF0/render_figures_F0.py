from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "reports" / "experiment_D" / "inputs"
OUT = ROOT / "reports" / "experiment_F0"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def read_matrix(path: Path) -> np.ndarray:
    return pd.read_csv(path).iloc[:, 1:].to_numpy(float)


branches = pd.read_csv(DATA / "branch.csv")
buses = pd.read_csv(DATA / "bus.csv")
ports = set(buses.loc[buses.has_gen.astype(bool), "bus"].astype(int))
G = nx.Graph()
G.add_nodes_from(buses.bus.astype(int))
for row in branches.itertuples():
    G.add_edge(int(row.src_bus), int(row.dst_bus), transformer=bool(row.transformer))
pos = nx.spring_layout(G, seed=3906, iterations=300, k=0.6)
fig, ax = plt.subplots(figsize=(10, 7))
line_edges = [(u, v) for u, v, d in G.edges(data=True) if not d["transformer"]]
xfmr_edges = [(u, v) for u, v, d in G.edges(data=True) if d["transformer"]]
nx.draw_networkx_edges(G, pos, edgelist=line_edges, ax=ax, width=1.05, alpha=0.55, edge_color="#667085")
nx.draw_networkx_edges(G, pos, edgelist=xfmr_edges, ax=ax, width=1.5, style="dashed", alpha=0.85, edge_color="#e07a35")
nx.draw_networkx_nodes(G, pos, nodelist=[n for n in G if n not in ports], ax=ax,
                       node_size=260, node_color="#dbe5ef", edgecolors="#344054", linewidths=0.7)
nx.draw_networkx_nodes(G, pos, nodelist=sorted(ports), ax=ax,
                       node_size=420, node_color="#f2b84b", edgecolors="#6b4611", linewidths=1.0)
nx.draw_networkx_labels(G, pos, ax=ax, font_size=7, font_color="#172b4d")
ax.set_title("IEEE-39 passive branch network\nGenerator ports retained for the Kron graph")
ax.text(0.01, 0.01, "Gold: generator buses 30–39   |   Dashed orange: transformer branches",
        transform=ax.transAxes, fontsize=8, color="#475467")
ax.axis("off")
save(fig, "F0_01_ieee39_generator_graph")

spec = pd.read_csv(OUT / "tables" / "TABLE_F02_graph_spectrum.csv")
fig, ax = plt.subplots(figsize=(8.5, 4.5))
colors = ["#d94f45" if z else "#276c91" for z in spec.zero_mode]
ax.stem(spec["mode"], spec.eigenvalue_pu_conductance, linefmt="#8b9aaa", markerfmt="o", basefmt=" ")
ax.scatter(spec["mode"], spec.eigenvalue_pu_conductance, c=colors, zorder=3, s=36)
ax.set(xlabel="Mode index", ylabel="Eigenvalue [p.u. conductance]",
       title="Spectrum of the passive generator-port graph")
ax.grid(axis="y", alpha=0.25)
save(fig, "F0_02_Lc_spectrum")

audit = pd.read_csv(OUT / "tables" / "TABLE_F03_kron_reduction.csv").set_index("metric").value
resid = float(audit["direct_vs_kron_relative_current_residual"])
fig, ax = plt.subplots(figsize=(7.5, 3.8))
ax.barh(["Direct port current vs. Kron map"], [max(resid, 1e-18)], color="#17806d", height=0.48)
ax.axvline(1e-11, color="#b54708", linestyle="--", linewidth=1.4, label="predeclared tolerance: $10^{-11}$")
ax.set_xscale("log")
ax.set_xlabel("Relative current residual")
ax.set_title("Kron reduction identity (7 deterministic voltage trials)")
ax.grid(axis="x", alpha=0.25)
ax.legend(frameon=False, loc="lower right")
save(fig, "F0_03_kron_reduction_residual")

U = read_matrix(OUT / "matrices" / "Lc_eigenvectors.csv")
fig, ax = plt.subplots(figsize=(8.5, 5.2))
im = ax.imshow(U, aspect="auto", cmap="RdBu_r", vmin=-np.max(abs(U)), vmax=np.max(abs(U)))
ax.set(xticks=np.arange(U.shape[1]), xticklabels=np.arange(1, U.shape[1] + 1),
       yticks=np.arange(U.shape[0]), yticklabels=sorted(ports),
       xlabel="Conductance graph mode", ylabel="Generator bus",
       title="Euclidean eigenvectors of $L_c=\\Re(Y_{port})$")
fig.colorbar(im, ax=ax, label="Eigenvector entry")
save(fig, "F0_04_graph_modes")

mode_path = OUT / "tables" / "TABLE_F05_graph_vs_real_modes.csv"
if mode_path.exists():
    modes = pd.read_csv(mode_path)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(modes.pole_id, modes.top_mode_energy, color="#276c91", label="largest single graph-mode energy")
    ax.scatter(modes.pole_id, modes.top3_graph_mode_energy, color="#e07a35", marker="D", s=28,
               label="energy in three strongest graph modes")
    ax.set(ylabel="M-weighted q-subvector energy", xlabel="Frozen ExpC critical pole",
           ylim=(0, 1.04), title="Physical critical-mode content in the preregistered graph basis")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    save(fig, "F0_05_graph_vs_real_mode_overlap")
else:
    fig, ax = plt.subplots(figsize=(7.5, 3.5))
    ax.axis("off")
    ax.text(0.5, 0.6, "Critical-mode overlap unavailable", ha="center", fontsize=15)
    ax.text(0.5, 0.42, "See REPORT_EXP_F0.md for the recorded data availability reason.", ha="center")
    save(fig, "F0_05_graph_vs_real_mode_overlap")

print(f"F0 figures written to {FIG}")
