"""Text-block extents (reference px) in the reference and in the render, for tuning type sizes."""
import numpy as np, pathlib
from PIL import Image
H = pathlib.Path(__file__).resolve().parent
ref = np.array(Image.open(H.parents[3] / 'resources2' / 'Póster científico_ más allá de la amortiguación nodal.png').convert('RGB')).astype(int)
out = np.array(Image.open(H / 'render_ref_scale.png').convert('RGB').resize((1161, 1355))).astype(int)
def box(a, x0, x1, y0, y1, pred):
    s = a[y0:y1, x0:x1]; m = pred(s); ys, xs = np.where(m)
    return (xs.min() + x0, xs.max() + x0, ys.min() + y0, ys.max() + y0) if len(xs) else None
dark = lambda s: (s[..., 0] < 60) & (s[..., 1] < 110) & (s[..., 2] < 90)
textish = lambda s: s.sum(axis=2) < 420
white = lambda s: s.min(axis=2) > 225
yellow = lambda s: (s[..., 0] > 190) & (s[..., 1] > 160) & (s[..., 2] < 130)
checks = {'title': (240, 910, 8, 64, dark), 'subtitle': (190, 925, 66, 86, textish), 'author': (240, 925, 89, 106, textish),
          'foot line1': (240, 950, 1264, 1298, white), 'foot line2': (190, 965, 1300, 1322, yellow), 'tagline': (1070, 1161, 8, 60, textish)}
for k, (x0, x1, y0, y1, p) in checks.items():
    print(f'{k:11s} ref {box(ref, x0, x1, y0, y1, p)}   out {box(out, x0, x1, y0, y1, p)}')
