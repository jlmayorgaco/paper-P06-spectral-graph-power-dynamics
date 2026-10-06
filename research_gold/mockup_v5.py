"""Mock-up v5 (36x48 in, coordinates in mm): diagram- and simulation-led layout inspired by the reference photo."""
import itertools, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle
import matplotlib.image as mpimg

FIG = pathlib.Path(r'C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\reports\poster\ias2026\one_mode_v2_20261006\generated\figures')
OUT = pathlib.Path(r'C:\Users\walla\Documents\Github\paper-P06-spectral-graph-power-dynamics\output\pdf\MOCKUP_v5_36x48.png')
DG, MG, DEEP, GOLD, PALE, PGOLD, BLUE, RED, GREY, INK, NAVY = '#003B2D', '#00664B', '#00291F', '#C99A20', '#EAF2ED', '#FFF4D6', '#176494', '#B3363A', '#5b6b73', '#17372D', '#0B2A4A'
W, H = 914, 1219
fig = plt.figure(figsize=(36, 48), dpi=50); K = 36 / 9.14 / 1.42
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis('off')
def box(x, y, w, h, fc='white', ec='#CFDDD5', lw=1.2, r=3, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw, zorder=z))
def txt(x, y, s, fs=9, c=INK, w='normal', ha='left', va='top', **k):
    ax.text(x, y, s, fontsize=fs * K, color=c, weight=w, ha=ha, va=va, zorder=6, **k)
def panel(x, y, w, h, num, title, fs=11, bar=NAVY):
    box(x, y, w, h); ax.add_patch(Rectangle((x, y), w, 14, fc=bar, ec=bar, zorder=2))
    ax.add_patch(Rectangle((x, y + 13), w, 1.2, fc=GOLD, ec=GOLD, zorder=2))
    if num is not None:
        ax.add_patch(Rectangle((x, y), 13, 14, fc=GOLD, ec=GOLD, zorder=3)); txt(x + 6.5, y + 7.2, str(num), fs, DG, 'bold', ha='center', va='center')
    txt(x + (16 if num is not None else 5), y + 7.2, title, fs, 'white', 'bold', va='center')
def img(name, x, y, w, h):
    im = mpimg.imread(FIG / f'{name}.png'); ih, iw = im.shape[:2]; s = min(w / iw, h / ih); dw, dh = iw * s, ih * s
    ax.imshow(im, extent=(x + (w - dw) / 2, x + (w + dw) / 2, y + (h + dh) / 2, y + (h - dh) / 2), zorder=3)
def words(x, y, w, h, s, fs=8.5):
    box(x, y, w, h, fc=PGOLD, ec=PGOLD, r=1.5, z=2); ax.add_patch(Rectangle((x, y), 1.6, h, fc=GOLD, ec=GOLD, zorder=3)); txt(x + 4, y + 2.5, s, fs)
def arrow(a, b, c=INK, lw=2.5, ms=16): ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=ms, lw=lw, color=c, zorder=5))
L, R0 = 23, 891; GW = R0 - L

# ---------------- header (white logo band + green title band)
txt(30, 14, 'IEEE IAS', 20, NAVY, 'bold'); txt(30, 30, 'Annual Meeting 2026  ·  Vancouver, Canada', 9.5, NAVY)
txt(885, 14, 'Universidad de', 13, INK, ha='right'); txt(885, 28, 'los Andes', 18, INK, 'bold', ha='right')
ax.add_patch(Rectangle((L, 44), GW, 82, fc=DG, ec=DG))
txt(W / 2, 52, 'BEYOND NODAL DAMPING', 11, GOLD, 'bold', ha='center')
txt(W / 2, 64, 'PRESERVING ONE MODE DOES NOT PROTECT THE GRID', 27, 'white', 'bold', ha='center')
txt(W / 2, 96, 'Exact PLL delay compensation, its network-wide price, and how to choose what to protect', 13, GOLD, 'bold', ha='center')
txt(W / 2, 111, 'Jorge Luis Mayorga Taborda  |  Universidad de los Andes, Colombia', 11, 'white', ha='center')

# ---------------- HOOK band: idea in one picture
y = 132; ax.add_patch(Rectangle((L, y), GW, 13, fc=NAVY, ec=NAVY))
txt(L + 4, y + 6.5, 'Hook', 10, GOLD, 'bold', va='center'); txt(L + 22, y + 6.5, 'A LOCAL RETUNE IS EXACT. THE WHOLE GRID PAYS ITS PRICE.', 12, 'white', 'bold', va='center')
box(L, y + 13, GW, 150, r=0)
# (a) block diagram
img('pll_loop', L + 4, y + 18, 270, 100)
txt(L + 8, y + 122, 'Each inverter sees the grid through a PLL whose\nmeasurement arrives τ late. Its only knobs: two gains.', 9.5)
# (b) s-plane cartoon
cx, cy = L + 300, y + 20; ax.add_patch(Rectangle((cx, cy), 270, 125, fc='white', ec='white'))
ax.plot([cx + 200, cx + 200], [cy + 8, cy + 118], color=RED, lw=2.5, zorder=4); ax.add_patch(Rectangle((cx + 200, cy + 8), 60, 110, fc='#F7E6E7', ec='none', zorder=2))
txt(cx + 205, cy + 100, 'unstable', 9, RED, 'bold')
ax.plot([cx + 10, cx + 255], [cy + 108, cy + 108], color=GREY, lw=1.5, zorder=3); txt(cx + 245, cy + 111, 'Re s', 8.5, GREY)
lam = (cx + 165, cy + 48); ax.add_patch(Circle(lam, 6, fc=BLUE, ec='white', lw=1.5, zorder=6)); ax.add_patch(Circle(lam, 11, fc='none', ec=BLUE, lw=2.2, ls='--', zorder=6))
txt(lam[0] - 40, lam[1] + 14, 'protected λ: held exactly', 9.5, BLUE, 'bold')
for (x0, y0, dx, dy) in [(cx + 95, cy + 62, 55, -6), (cx + 120, cy + 82, 95, 4), (cx + 60, cy + 90, 30, 6), (cx + 170, cy + 70, 45, -2)]:
    ax.add_patch(Circle((x0, y0), 4.5, fc='white', ec=GREY, lw=1.6, zorder=6)); arrow((x0 + 4, y0), (x0 + dx, y0 + dy), RED if x0 + dx > cx + 200 else GREY, 2.4)
txt(cx + 12, cy + 12, 'other modes μ move by', 10, INK, 'bold')
txt(cx + 12, cy + 24, r'$\frac{d\mu}{dh}=-k_p(\mu-\lambda)(\mu-\bar\lambda)\,\frac{\partial\mu}{\partial k_I}$', 13, BLUE)
txt(cx + 12, cy + 113, 'PLL: quadratic factor.  NETWORK: residue ∂μ/∂k_I', 8.6, GREY)
# (c) punchline
px = L + 590; box(px, y + 20, 270, 128, fc=PALE, ec=PALE)
txt(px + 135, y + 26, 'Same delay increase, 40 → 44 ms, ten PLLs', 10, DG, 'bold', ha='center')
rows = [('no retune', GREY, '8 unstable', 'aborted'), ('low-frequency rule (usual)', GOLD, '12 unstable', 'aborted'), ('protect 4.91 Hz', BLUE, '0', '2 / 5 events'), ('protect 1.04 Hz', DG, '0', '5 / 5 events')]
for k, (n, c, r, e) in enumerate(rows):
    yy = y + 46 + k * 17; ax.add_patch(Rectangle((px + 6, yy - 1), 4, 12, fc=c, ec=c, zorder=4))
    txt(px + 14, yy, n, 9.5, c, 'bold'); txt(px + 150, yy, r.replace(' unstable',''), 9.5, RED if 'unstable' in r else DG, 'bold'); txt(px + 205, yy, e, 9.5, RED if e == 'aborted' or e.startswith('2') else DG, 'bold')
txt(px + 150, y + 38, 'unstable roots', 8, GREY); txt(px + 205, y + 38, 'events', 8, GREY)
txt(px + 135, y + 118, 'The usual compensation is WORSE than none.\nWhich mode you protect decides everything.', 10.5, RED, 'bold', ha='center')

# ---------------- row 2: law + price
y = 300; w2 = (GW - 8) / 2
panel(L, y, w2, 168, 2, 'THE EXACT LAW: TWO GAINS, ONE MODE')
txt(L + 8, y + 20, 'Delay grows τ → τ + h. Keep the loop factor κ(s) = (k_p s + k_I)e^{-sτ} unchanged at λ = α + jω:', 9.5)
txt(L + 30, y + 36, r'$k_p(h)\,\lambda+k_I(h)=(k_p\lambda+k_I)\,e^{\lambda h}$', 19, BLUE)
txt(L + 8, y + 66, r'$k_p(h)=e^{\alpha h}\left[k_p\cos\omega h+\frac{k_I+\alpha k_p}{\omega}\sin\omega h\right]$', 12.5, BLUE)
txt(L + 8, y + 88, r'$k_I(h)=e^{\alpha h}\left[k_I\cos\omega h-\frac{|\lambda|^2k_p+\alpha k_I}{\omega}\sin\omega h\right]$', 12.5, BLUE)
txt(L + 8, y + 112, 'What is left at every other frequency (exact):', 9.5, INK, 'bold')
txt(L + 8, y + 122, r'$\Delta\kappa(s)=-(s-\lambda)(s-\bar\lambda)\,e^{-s\tau}\int_0^h k_p(t)e^{-st}dt$', 13, BLUE)
words(L + 6, y + 146, w2 - 12, 18, 'No network model needed for the gains. Budget: h < φ_PI(ω)/ω  (8.99 ms at 4.91 Hz).', 9)
x2 = L + w2 + 8
panel(x2, y, w2, 168, 3, 'THE GRID SETS THE PRICE: CHOOSE THE MODE')
txt(x2 + 8, y + 20, r'$\frac{\partial\mu}{\partial k_I}=-\frac{\ell^H Z_{k_I}(\mu)\,r}{\ell^H Z_s(\mu)\,r}$' + '   reduced network operator Z(s), left/right vectors ℓ, r', 11, BLUE)
img('price_prediction', x2 + 4, y + 42, 255, 110)
txt(x2 + 268, y + 44, 'Selection rule', 11, DG, 'bold')
for k, s_ in enumerate(['One eigen-analysis\nof the base design', 'Predict the worst root\nfor every candidate λ', 'Among safe ones,\nkeep the most k_I', 'Verify: root count\nAND nonlinear events']):
    yy = y + 60 + k * 24; ax.add_patch(Circle((x2 + 273, yy + 5), 4.2, fc=GOLD, ec='white', zorder=6)); txt(x2 + 273, yy + 5.3, str(k + 1), 8.5, DG, 'bold', ha='center', va='center')
    txt(x2 + 281, yy, s_, 8.8)
txt(x2 + 8, y + 155, 'Predicted vs exact worst root, 11 candidates: rank corr. 0.99 (same order).', 8.5, GREY)

# ---------------- row 3: s-plane simulations
y = 476; panel(L, y, GW, 172, 4, 'SIMULATION 1 — WHERE THE ROOTS GO (IEEE-39, FULL DELAYED DAE, 204 STATES)')
img('splane4', L + 6, y + 18, GW - 12, 150)
# ---------------- row 4: time simulations
y = 656; panel(L, y, GW, 132, 5, 'SIMULATION 2 — WHAT THE GRID DOES IN TIME (BUS 16 +100 MW, 44 ms)')
img('traces4', L + 6, y + 17, GW - 12, 95)
txt(L + 8, y + 116, '5 Hz part of the bus-30 frequency. Unprotected, the oscillation grows until the delayed nonlinear integration diverges; protected, it dies out within 2 s.\nFull-spectrum count at 44 ms: 8 / 12 / 0 / 0 unstable roots.  Events (5 load steps, limits on Δf, RoCoF, V, SG reserve): aborted / aborted / 2 of 5 / 5 of 5.', 8.6)

# ---------------- row 5: four consequence panels
y = 796; w4 = (GW - 24) / 4
panel(L, y, w4, 196, 6, 'PASS BOTH TESTS', 10.5)
img('central', L + 4, y + 16, w4 - 8, 120)
txt(L + 6, y + 140, 'Roots clean, events fail: 4.91 Hz spends k_I\n(246.7 → 137.8), SG reserve ≈ 0.\nEvents clean, roots unstable: no retune at\n42 ms meets the limits for 30 s with 4 unstable\nroots — 0.5 s windows miss slow 5 Hz growth.', 8.3)
x = L + w4 + 8; panel(x, y, w4, 196, 7, 'REACH: SHARE × DELAY', 10.5)
img('share_delay_map', x + 4, y + 16, w4 - 8, 118)
txt(x + 6, y + 140, 'Rule keeps the margin at 44/48 ms wherever\n40 ms is feasible. No higher maximum share:\n90 % fails events even at 40 ms.\nSecond design (85 %), rule unchanged: 5/5.\nNine-bus bank: +0.886 → −0.131 s⁻¹.', 8.3)
x = L + 2 * (w4 + 8); panel(x, y, w4, 196, 8, 'BEYOND NODAL DAMPING', 10.5)
txt(x + 6, y + 18, r'$D_H=\frac{1}{2}(Z+Z^H)$,  $\delta D_H=\frac{1}{2}(av^H+va^H)$', 10.5, BLUE)
img('rank_two2', x + 4, y + 34, w4 - 8, 98)
txt(x + 6, y + 140, 'One local PLL change is never pure\ndamping: one + and one − direction\n(60/60). Counting law: r bad directions\nneed ≥ r retuned PLLs. Nodal +38.7 vs\ncollective −737.5.', 8.3)
x = L + 3 * (w4 + 8); panel(x, y, w4, 196, 9, 'SAFE ALONE ≠ TOGETHER', 10)
cx = x + w4 / 2; nodes = {}
for size in range(5):
    subs = list(itertools.combinations([30, 33, 35, 37], size))
    for j, sset in enumerate(subs): nodes[sset] = (cx + (j - (len(subs) - 1) / 2) * 24, y + 112 - size * 20)
for a in nodes:
    for b in nodes:
        if len(b) == len(a) + 1 and set(a) <= set(b): ax.plot([nodes[a][0], nodes[b][0]], [nodes[a][1], nodes[b][1]], color='#c9d6cf', lw=0.8, zorder=3)
for s_, (px_, py_) in nodes.items():
    full = len(s_) == 4; ax.add_patch(Circle((px_, py_), 4.5 if full else 3, fc=(RED if full else BLUE), ec='white', lw=1, zorder=4))
txt(cx, y + 18, 'H4 = {30,33,35,37}: 15 stable / 1 unstable', 9, RED, 'bold', ha='center')
txt(x + 6, y + 122, r'$\det(I+Q_H)=\det(I+Q_{RR})\,\det(I-R_{i|R})$', 9.5, BLUE)
txt(x + 6, y + 140, 'Pair retunes too: buses 30 & 37 hold\n−0.0657 alone, −0.0414 together; interaction\n+0.0243 ± 4e-7 (interval-certified), repaired\nto −0.060 in 11 steps. Blind H4 prediction\n16/16 (4 buses), 395/512 (9 buses).', 8.3)

# ---------------- KPI strip
y = 1000; ax.add_patch(Rectangle((L, y), GW, 13, fc=NAVY, ec=NAVY)); txt(W / 2, y + 6.5, 'EVIDENCE AT A GLANCE', 10.5, 'white', 'bold', ha='center', va='center')
box(L, y + 13, GW, 46, r=0)
kp = [('12 → 0', 'unstable roots at 44 ms\nusual rule → our rule'), ('5 / 5', 'delayed nonlinear\nevents passed'), ('0.99', 'rank correlation of the\npredicted price'), ('3.9×10⁻⁵', 'max error of the spillover\nidentity, 219 mode pairs'), ('485', 'interval panels certifying\none root cluster'), ('204', 'states in the full\ndelayed IEEE-39 model')]
for k, (v, l) in enumerate(kp):
    xx = L + 6 + k * (GW / 6); txt(xx + GW / 12 - 3, y + 17, v, 20, RED if k == 0 else DG, 'bold', ha='center'); txt(xx + GW / 12 - 3, y + 41, l, 8.3, GREY, ha='center')
    if k: ax.plot([xx - 3, xx - 3], [y + 17, y + 55], color='#CFDDD5', lw=1.2)
# ---------------- banner
y = 1066; ax.add_patch(Rectangle((L, y), GW, 24, fc=NAVY, ec=NAVY))
txt(W / 2, y + 12, 'CHOOSE THE MODE.  COMPUTE THE GAINS.  CHECK THE WHOLE GRID.', 15, 'white', 'bold', ha='center', va='center')
# ---------------- conclusions / refs / contact
y = 1096; w3 = (GW - 16) / 3
panel(L, y, w3, 104, 10, 'THREE RULES FOR THE ENGINEER', 10)
for k, s_ in enumerate(['Do not compensate PLL delay at low frequency.', 'Choose the protected mode with the spillover law.', 'Check the spectrum AND the events.']):
    ax.add_patch(Circle((L + 9, y + 26 + k * 16), 4, fc=GOLD, ec='white', zorder=6)); txt(L + 9, y + 26.3 + k * 16, str(k + 1), 8.5, DG, 'bold', ha='center', va='center'); txt(L + 17, y + 22 + k * 16, s_, 9.3)
txt(L + 6, y + 74, 'Not claimed: max share, optimal mode, superiority\nover other compensators, novelty of the algebra.', 8, GREY)
x = L + w3 + 8; panel(x, y, w3, 104, 11, 'REFERENCES', 10)
txt(x + 6, y + 20, '[1] U. Markovic et al., IEEE TPWRS 36(5), 2021.\n[2] L. Huang et al., arXiv:1903.05489 (PLL and grid structure).\n[3] W. Michiels, S. Gumussoy, arXiv:2003.05496.\n[4] F. Dorfler et al., IEEE TPWRS 29, 2014.\n[5] D. Bindel, A. Hood, SIAM J. Matrix Anal. Appl., 2013.\n[6] F. Johansson, Arb, IEEE Trans. Comput., 2017.', 8.2)
x = L + 2 * (w3 + 8); panel(x, y, w3, 104, 12, 'CONTACT', 10)
txt(x + 8, y + 24, 'Jorge Luis Mayorga Taborda', 11, DG, 'bold'); txt(x + 8, y + 40, 'jl.mayorga@uniandes.edu.co\nUniversidad de los Andes, Bogota, Colombia\nIEEE IAS Annual Meeting 2026, Vancouver', 9.3)
txt(W / 2, H - 4, 'MOCK-UP v5 (layout preview, 36 x 48 in) — figures are real results; typesetting is not final', 8, GREY, ha='center', va='bottom')
fig.savefig(OUT, dpi=50); print(OUT)
