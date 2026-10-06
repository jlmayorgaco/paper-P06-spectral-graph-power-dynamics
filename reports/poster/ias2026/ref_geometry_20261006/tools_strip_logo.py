"""Make a clean, vector-only copy of the IAS 2026 logo: drop every image XObject (the hidden photos of the Call for Papers)
and crop the page to the logo. usage: python tools_strip_logo.py in.pdf out.pdf  (needs pikepdf)"""
import sys
import pikepdf
from pikepdf import Name, Array

src, dst = sys.argv[1], sys.argv[2]
pdf = pikepdf.open(src)
page = pdf.pages[0]
xo = page.Resources.get('/XObject', {})
images = {k for k, v in xo.items() if v.get('/Subtype') == Name.Image}
ops = []
removed = 0
for operands, op in pikepdf.parse_content_stream(page):
    if str(op) == 'Do' and operands and str(operands[0]) in {str(k) for k in images}:
        removed += 1
        continue
    ops.append((operands, op))
page.Contents = pdf.make_stream(pikepdf.unparse_content_stream(ops))
for k in list(images):
    del xo[k]
H = float(page.MediaBox[3])
# logo box measured on the cairo render (pt, origin top-left): x 238-358, y 351-452
box = Array([238, H - 452, 358, H - 351])
page.MediaBox = box
page.CropBox = box
pdf.remove_unreferenced_resources()
pdf.save(dst)
print('removed image draws:', removed, 'remaining images:', len([1 for v in page.Resources.get('/XObject', {}).values() if v.get('/Subtype') == Name.Image]))
