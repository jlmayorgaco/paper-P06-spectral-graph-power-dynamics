# Shared style for module 8 plots (DESIGN_V2 section 5). Fira Sans registered from MiKTeX, pdf.fonttype 42, exact slot size.
import sys, matplotlib
matplotlib.use("Agg")
from matplotlib import rcParams, font_manager as fm
FD = "C:/Users/walla/AppData/Roaming/MiKTeX/fonts/opentype/public/fira/"
for f in ["FiraSans-Regular", "FiraSans-Medium", "FiraSans-SemiBold", "FiraSans-Bold"]:
    fm.fontManager.addfont(FD + f + ".otf")
PXW, PXH = 0.787597 / 25.4, 0.787299 / 25.4   # inches per reference px
rcParams.update({"font.family": "Fira Sans", "font.weight": "medium", "font.size": 17.5, "pdf.fonttype": 42,
                 "axes.linewidth": 1.4, "xtick.major.width": 1.2, "ytick.major.width": 1.2,
                 "xtick.major.size": 5, "ytick.major.size": 5, "xtick.direction": "out", "ytick.direction": "out",
                 "mathtext.fontset": "custom", "mathtext.rm": "Fira Sans", "mathtext.it": "Fira Sans:italic", "mathtext.default": "regular"})
GREEN, NAVY, GOLD, RED, TEAL, GREY = "#03534A", "#1E4F8A", "#D9A520", "#B83A3A", "#2C7A70", "#8A9894"
TEXT, MUTED, GRID = "#1C2B27", "#4F5D59", "#E3EBE7"
DER = "../../../../../../research/nhop_ieee39_20261006/derived/"
