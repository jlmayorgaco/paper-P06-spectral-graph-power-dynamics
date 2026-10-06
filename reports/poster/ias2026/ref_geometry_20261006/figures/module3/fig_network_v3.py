# M3 v3: network lines + load buses only; generator glyphs are drawn in TikZ at the printed coordinates.
# M3 v2: IEEE-39 network illustration (hand-tuned layout of the poster network figure; real branch list).
# source: ../../../generated/figures/network_ieee39.tikz (46 branches, bus coordinates from the frozen poster network figure)
#         ../../../../../experiment_D/inputs is not present in this checkout; topology identical (46 branches).
import re, sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, font_manager
from matplotlib.patches import Circle
import glob, os
FD = os.path.expanduser("~/AppData/Roaming/MiKTeX/fonts/opentype/public/fira/")
for f in ["FiraSans-Regular.otf", "FiraSans-Medium.otf", "FiraSans-SemiBold.otf", "FiraSans-Bold.otf"]:
    font_manager.fontManager.addfont(FD + f)
rcParams["font.family"] = "Fira Sans"
rcParams["pdf.fonttype"] = 42
SRC = "../../../generated/figures/network_ieee39.tikz"
txt = open(SRC, encoding="utf8").read()
pos = {int(a): (float(b), float(c)) for a, b, c in re.findall(r"coordinate \(b(\d+)\) at \(([\d.]+)mm,([\d.]+)mm\)", txt)}
edges = [(int(a), int(b)) for a, b in re.findall(r"\(b(\d+)\)--\(b(\d+)\)", txt)]
assert len(pos) == 39 and len(edges) == 46, (len(pos), len(edges))
W_PX, H_PX = float(sys.argv[1]), float(sys.argv[2])
PXW, PXH = 0.787597 / 25.4, 0.787299 / 25.4
# tikz (x, y) -> picture (px = y, py = -x): same orientation as the PNG network figure
P = {b: (y, x) for b, (x, y) in pos.items()}
xs = [v[0] for v in P.values()]; ys = [v[1] for v in P.values()]
x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
pad = 9.0
sc = min((W_PX) / (x1 - x0 + 2 * pad), (H_PX) / (y1 - y0 + 2 * pad))
fig = plt.figure(figsize=(W_PX * PXW, H_PX * PXH))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
ax.set_xlim(cx - W_PX / sc / 2, cx + W_PX / sc / 2)
ax.set_ylim(cy - H_PX / sc / 2, cy + H_PX / sc / 2)
ax.set_aspect("equal")
H4 = {30, 33, 35, 37}
GENS = set(range(30, 40))
GREEN, TEAL, RED, MIST = "#03534A", "#2C7A70", "#B83A3A", "#A9C4BC"
for a, b in edges:
    (xa, ya), (xb, yb) = P[a], P[b]
    ax.plot([xa, xb], [ya, yb], color="#8FB6AC", lw=3.2, solid_capstyle="round", zorder=1)
for b, (x, y) in P.items():
    if b in GENS:
        continue
    ax.add_patch(Circle((x, y), 2.9 / sc, fc="#DCEBE5", ec=TEAL, lw=1.8, zorder=2))
fig.savefig(sys.argv[3], transparent=True)
left, top = cx - W_PX / sc / 2, cy + H_PX / sc / 2
for b in sorted(GENS):
    x, y = P[b]
    print(b, round((x - left) * sc, 1), round((top - y) * sc, 1))
