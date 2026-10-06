"""Extra figures for poster v5: s-plane small multiples (T16), raw/band traces (T15), PLL loop block diagram."""
import pathlib, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.signal import butter, sosfiltfilt
CH = pathlib.Path(__file__).resolve().parent
OUT = CH.parents[1] / 'reports' / 'poster' / 'ias2026' / 'one_mode_v2_20261006' / 'generated' / 'figures'
BLUE, GREEN, GOLD, RED, GREY, INK = '#176494', '#003B2D', '#C99A20', '#B3363A', '#5b6b73', '#17372D'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 18, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 1.8})
def save(fig, n):
    fig.savefig(OUT / f'{n}.png', dpi=250, bbox_inches='tight'); fig.savefig(OUT / f'{n}.pdf', bbox_inches='tight'); plt.close(fig)

# 1) s-plane small multiples: base (hollow) -> 44 ms (filled), the ~5 Hz PLL cluster
d = pd.read_csv(CH / 'T16_splane_roots.csv'); base = d[d.strategy == 'base40']
fig, axs = plt.subplots(1, 4, figsize=(24, 5.2), sharey=True)
for ax, (k, lab, col) in zip(axs, [('none', 'no retune', GREY), ('lowfreq', 'low-frequency rule', GOLD), ('node491', 'protect 4.91 Hz', BLUE), ('node104', 'protect 1.04 Hz', GREEN)]):
    g = d[(d.strategy == k) & (d.im > 25)]
    ax.axvspan(0, 2, color=RED, alpha=.09, lw=0); ax.axvline(0, color=RED, lw=1.6)
    for _, r in g.iterrows():
        ax.annotate('', xy=(r.re, r.im / 6.2832), xytext=(r.origin_re, r.origin_im / 6.2832), arrowprops=dict(arrowstyle='-|>', color=col, lw=1.8, alpha=.85))
    b = base[base.im > 25]
    ax.scatter(b.re, b.im / 6.2832, s=70, facecolor='white', edgecolor=GREY, lw=1.6, zorder=3)
    ax.scatter(g.re, g.im / 6.2832, s=90, color=col, zorder=4)
    n = int((g.re > 0).sum())
    ax.set_title(f'{lab}\n{n} unstable pairs' if n else f'{lab}\nall stable', fontsize=17, color=col, fontweight='bold')
    ax.set_xlim(-5.6, 1.6); ax.set_xlabel(r'$\Re s$ [s$^{-1}$]'); ax.grid(alpha=.2)
    if k.startswith('node'):
        prot = 4.91 if k == 'node491' else None
        if prot: ax.scatter([-0.8016], [30.8505 / 6.2832], s=420, facecolor='none', edgecolor=BLUE, lw=2.6, zorder=5); ax.text(-0.75, 4.97, 'held', color=BLUE, fontsize=15, fontweight='bold')
axs[0].set_ylabel('frequency [Hz]'); axs[0].set_ylim(4.05, 5.05)
fig.text(0.5, -0.04, 'open circles: 40 ms base     filled: 44 ms     arrows: root motion (continuation of 8 catalogued PLL-cluster roots)', ha='center', fontsize=15, color=GREY)
fig.tight_layout(); save(fig, 'splane4')

# 2) band-limited bus frequency, 4 small multiples on a common linear scale
def band(name):
    f = CH / f'T15_{name}_trajectory.csv'; t = pd.read_csv(f)
    x = sosfiltfilt(butter(4, [4.2, 5.6], btype='band', fs=100, output='sos'), t['f_bus30_Hz'].values) * 1e3
    return t.time_after_event_s.values, x
fig, axs = plt.subplots(1, 4, figsize=(24, 3.4), sharey=True)
for ax, (k, lab, col, end) in zip(axs, [('none44s', 'no retune', GREY, 7.0), ('lowfreq44t', 'low-frequency rule', GOLD, 4.0), ('node491_44', 'protect 4.91 Hz', BLUE, 8), ('node104_44', 'protect 1.04 Hz', GREEN, 8)]):
    t, x = band(k); m = t <= end
    ax.plot(t[m], x[m], color=col, lw=1.4); ax.set_xlim(0, 8); ax.set_ylim(-6, 6); ax.grid(alpha=.2)
    ax.set_title(lab, fontsize=17, color=col, fontweight='bold'); ax.set_xlabel('time [s]')
    if k in ('none44s', 'lowfreq44t'): ax.scatter([t[m][-1]], [0], marker='X', s=180, color=RED, zorder=5); ax.text(t[m][-1], 4.6, 'diverges', color=RED, ha='right', fontsize=14)
axs[0].set_ylabel('5 Hz part of\nbus-30 freq. [mHz]')
fig.tight_layout(); save(fig, 'traces4')

# 3) PLL + grid loop block diagram
fig, ax = plt.subplots(figsize=(10, 3.6)); ax.set_xlim(0, 100); ax.set_ylim(0, 36); ax.axis('off')
def blk(x, y, w, h, t, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.5', fc=fc, ec=ec, lw=2.2)); ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=19, color=INK)
def arr(a, b, c=INK, ls='-'): ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=22, lw=2.4, color=c, ls=ls))
blk(2, 20, 17, 11, r'$e^{-s\tau}$' + '\ndelay', '#FFF4D6', GOLD)
blk(26, 20, 21, 11, r'$k_p+\frac{k_I}{s}$' + '\nPI', '#EAF2ED', BLUE)
blk(54, 20, 18, 11, r'$\frac{1}{s(1+t_f s)}$' + '\nPLL', '#EAF2ED', BLUE)
blk(26, 1, 46, 13, r'grid $G(s;\rho)$' + '\n39 buses, SGs, other PLLs', '#FFF4D6', GOLD)
arr((19.8, 25.5), (25.6, 25.5)); arr((47.8, 25.5), (53.6, 25.5)); arr((72.8, 25.5), (88, 25.5))
ax.text(89, 25.5, r'$\theta$', fontsize=22, va='center'); ax.plot([86, 86], [25.5, 7.5], color=INK, lw=2.4); arr((86, 7.5), (73.2, 7.5))
ax.plot([25.2, 10.5], [7.5, 7.5], color=INK, lw=2.4); arr((10.5, 7.5), (10.5, 19.6))
ax.text(12, 13, 'phase error', fontsize=16, color=GREY)
ax.text(36, 33.5, r'$\kappa(s)=(k_p s+k_I)\,e^{-s\tau}$', fontsize=21, color=BLUE, ha='center')
save(fig, 'pll_loop')
print('ok')
