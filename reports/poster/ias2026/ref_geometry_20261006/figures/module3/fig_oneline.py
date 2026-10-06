# One-line schematic of IEEE-39: real branch list, force-directed layout (positions are schematic only).
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["font.family"] = "Segoe UI"
REPO = "../../../../../../"
br = pd.read_csv(REPO + "reports/experiment_D/inputs/branch.csv")
mc = pd.read_csv(REPO + "reports/experiment_D/inputs/machine.csv")
gens = [int(b) for b in mc["bus"]]
H4 = [30, 33, 35, 37]
edges = [(int(a), int(b)) for a, b in zip(br.src_bus, br.dst_bus)]
buses = sorted({b for e in edges for b in e}); idx = {b: i for i, b in enumerate(buses)}
n = len(buses); print("buses", n, "branches", len(edges), "gens", gens)
rng = np.random.default_rng(7); pos = rng.uniform(-1, 1, (n, 2))
k = 1.0 / np.sqrt(n) * 1.6
for it in range(600):
    t = 0.1 * (1 - it / 600) + 0.002
    d = pos[:, None, :] - pos[None, :, :]; r = np.linalg.norm(d, axis=2) + 1e-9
    f = (d / r[..., None] ** 2 * k * k).sum(1)
    for a, b in edges:
        i, j = idx[a], idx[b]; v = pos[i] - pos[j]; l = np.linalg.norm(v) + 1e-9
        f[i] -= v * l / k; f[j] += v * l / k
    ln = np.linalg.norm(f, axis=1, keepdims=True) + 1e-9
    pos += f / ln * np.minimum(ln, t)
pos -= pos.mean(0)
W, Hh = 133 * 0.787597 / 25.4, 150 * 0.787299 / 25.4
fig = plt.figure(figsize=(W, Hh)); ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
for a, b in edges:
    p, q = pos[idx[a]], pos[idx[b]]
    ax.plot([p[0], q[0]], [p[1], q[1]], color="#B9CCD2", lw=1.0, zorder=1)
OFF = {30: (0.0, 0.26), 33: (-0.4, 0.0), 35: (0.0, -0.36), 37: (0.42, 0.0)}
for b in buses:
    p = pos[idx[b]]
    if b in H4:
        ax.scatter(*p, s=130, c="#C73E3E", ec="#7a1f1f", lw=1.0, zorder=4)
        dv = p / (np.linalg.norm(p) + 1e-9) * 0.0 + np.array(OFF[b])
        ax.text(p[0] + dv[0], p[1] + dv[1], str(b), ha="center", va="center", color="#7a1f1f", fontsize=15.5, fontweight="bold", zorder=5)
    elif b in gens:
        ax.scatter(*p, s=55, c="#2C7A70", ec="white", lw=0.6, zorder=3)
    else:
        ax.scatter(*p, s=14, c="#53625E", zorder=2)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(pos[:,0].min()-0.42, pos[:,0].max()+0.42); ax.set_ylim(pos[:,1].min()-0.36, pos[:,1].max()+0.36)
fig.savefig("oneline.pdf", transparent=True)
