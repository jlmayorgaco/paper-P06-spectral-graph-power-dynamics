# Print preflight — 36 x 42 in (2026-10-06)

| Check | Result | Status |
|---|---|---|
| Page size | 2592 x 3024 pt = 36.00 x 42.00 in, 1 page; Media/Trim/Bleed boxes equal | OK |
| Fonts | all embedded and subset (Fira Sans, Latin Modern Math, Bahnschrift) | OK |
| Equations / art | vector (TikZ + matplotlib PDF, fonttype 42) | OK |
| Raster images | none: the IEEE IAS logo is now the official vector (cropped from the AM 2026 Call for Papers, images stripped with tools_strip_logo.py); Uniandes logo vector | OK |
| Side margins | panels and (now) header rule + footer inset at 9.5 mm; only the faint header hills/skyline run to the edge (decorative, safe if trimmed) | OK |
| Top / bottom margin | 8.5 mm / 6.1 mm | OK (above the usual 5 mm unprintable edge) |
| Row gutters | were 6.7 / 2.4 mm (row 2 → row 3 too tight); fixed to ≈ 6 mm everywhere (row 3, lower strip +5 px, footer +4 px) | FIXED |
| Column gutters | ≈ 6.3–7 mm; row-3 column edges sit 5–6 mm left of rows 1–2 (inherited from the reference grid) | minor, cosmetic |
| Text size | body 20–24 pt, captions 17–19 pt, pills 15.5 pt; icon-internal labels (PI, θ, Δ, τ) 10.5–15 pt — decorative only | OK |
| Colour | RGB document; deep green #03534A, navy #1E4F8A, red #B83A3A, gold #EABD38 are inside typical CMYK gamut; gold may print slightly duller | ask the printer for an RGB-accepting workflow or a proof |
| Pale tints | panel body #FCFEFC and grid lines #E3EBE7 are near-white: they may vanish on some printers (intended to be subtle) | OK |
| QR | removed on request (no reserved space) | OK |

Send to the printer: `main.pdf` (vector, fonts embedded). Ask for a 100 % scale print with no "fit to page" and a colour proof.
