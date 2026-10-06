"""Compare the LaTeX render with the reference PNG (both at 1161 x 1355): geometry detection + blended overlay + side by side."""
import sys, pathlib
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
H = pathlib.Path(__file__).resolve().parent
REF = H.parents[3] / 'resources2' / 'Póster científico_ más allá de la amortiguación nodal.png'
ref = Image.open(REF).convert('RGB'); out = Image.open(H / 'render_ref_scale.png').convert('RGB').resize(ref.size)
def bars(im):
    a = np.array(im).astype(int); r, g, b = a[..., 0], a[..., 1], a[..., 2]
    m = (r < 40) & (g > 60) & (g < 110) & (b > 40) & (b < 90)
    lab, n = ndi.label(ndi.binary_closing(m, iterations=2)); res = []
    for s in ndi.find_objects(lab):
        h = s[0].stop - s[0].start; w = s[1].stop - s[1].start
        if w > 150 and 15 < h < 60: res.append((s[0].start, s[0].stop, s[1].start, s[1].stop))
    return sorted(res)
rb, ob = bars(ref), bars(out)
print('reference bars:', len(rb), ' render bars:', len(ob))
for x, y in zip(rb, ob): print('ref y%4d-%4d x%4d-%4d | out y%4d-%4d x%4d-%4d' % (x + y))
a = np.array(ref).astype(float); b = np.array(out).astype(float)
Image.fromarray(((a * .5 + b * .5)).astype('uint8')).save(H / 'compare_overlay.png')
w, h = ref.size; sbs = Image.new('RGB', (2 * w, h), 'white'); sbs.paste(ref, (0, 0)); sbs.paste(out, (w, 0)); sbs.save(H / 'compare_side_by_side.png')
# panel border / footer rows
def rows(im):
    a = np.array(im).astype(int); dark = (a.sum(axis=2) < 330); frac = dark.mean(axis=1); return [y for y in range(1150, 1355) if frac[y] > .8]
print('footer dark rows ref', rows(ref)[:1], rows(ref)[-1:], ' out', rows(out)[:1], rows(out)[-1:])
