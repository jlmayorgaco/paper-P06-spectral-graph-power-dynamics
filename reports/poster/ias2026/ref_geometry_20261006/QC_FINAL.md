# QC_FINAL (geometry corrections + Module 1)

## Structural corrections
1. Numbering: scientific panels are 1-9 only. The lower strip is unnumbered (no gold circle): MAIN CONTRIBUTIONS | SCOPE OF THE RESULTS | REFERENCES / CONTACT. PASS.
   Note: the TAKEAWAY banner stays BELOW this strip, as in the accepted geometry and the reference image (the brief says "above"; the accepted geometry was not changed).
2. QR: `assets/repo_qr.png` decodes (OpenCV QRCodeDetector) to `https://github.com/jlmayorgaco/paper-P06-spectral-graph-power-dynamics`, the working repository (= git origin), not a verified final paper/code/data URL. Removed; blank space reserved at x 1058-1146, y 1143-1239 (reference px). PASS.
3. Title: Bahnschrift variable font, wght=700 / wdth=75 (real condensed bold), natural proportions, no anisotropic scaling. Panel titles: natural size, uniform scale-down only when wider than the bar. Header text uses uniform scaling only. PASS.
   (Natural title is slightly shorter than the reference lettering, which is a narrower display face.)

## Module 1 checklist
- compile LuaLaTeX: OK (2 passes, 1 page, 2592 x 3024 pt = 36 x 42 in, fonts embedded).
- render 75 dpi: render_75dpi.png (2700 x 3150). Comparison: compare_side_by_side.png, compare_overlay.png.
- Panel 1 crop: crop_panel1_150dpi.png.
- No clipping; equation hierarchy: objective is the first anchor after the opening statement; flow arrows aligned.
- Smallest text: 17 pt (legend, labels, note). Opening statement 2.3 lines (17.5 pt).
- Prose ~15 % (opening statement + legend + two notes).
- No 40/44 ms, H4, Q, W_i, "optimal" claims; K described as architecture AND gains, not diagonal.
- Remaining cosmetic points: the four flow icons are small relative to their boxes; a blank area right of the legend.

MODULE_1_PASS

## Refinement pass (2026-10-06, "más pixel perfect")
- Header: title Bahnschrift Bold Condensed at natural proportions (140 pt), x 335-820 vs reference 339-820. Subtitle (Bahnschrift SemiBold SemiCondensed 31 pt), author line (28.8 pt) and tagline (23 pt, letter-spaced) placed within 1-4 px of the reference extents (measure.py).
- Footer: both lines in Bahnschrift condensed, extents 284-881 / 211-947 vs reference 284-881 / 207-951; em dash fixed; footer art redrawn in TikZ (pylons, mountains, pines, Vancouver block) without overlapping the text.
- Skyline: thin, low, light buildings + distant range + pines at the far right, behind the tagline only (no overlap with the author line).
- Gold numbers now have a thin deep-gold ring; inner step markers are dark-green discs with a gold ring, as in the reference.
- Module 1 recomposed on the reference layout (statement, SG -> grid -> replace -> GFL schematic, portfolio variable with its meaning, engineering objective, red scope note, four numbered steps, green callout). No optimality claim.
- Known difference: the IEEE IAS logo is the repository asset, whose artwork differs from the logo drawn in the reference image.

MODULE_1_PASS (refined)
