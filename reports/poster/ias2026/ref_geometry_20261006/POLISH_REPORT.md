# Polish report (2026-10-06)

Build: `build_polish/main.pdf` (2592 x 3024 pt, 1 page, 2 LuaLaTeX passes, no errors, no overfull boxes). Crops (150 dpi): `crops_polish/hdr.png, p1..p9.png, strip_a.png, strip_b.png, ftr.png`, detail crops `p5b.png`, `p9b.png`. Helper: `crop.sh`. Note: panels 7-9 and the strip are drawn unshifted at y 801.5+; only their module content has the +5 px shift. No numbers, `% source:` comments, palette, header/footer wording or panel geometry were changed.

## A. Pending items
1. P3 lattice caption (`module3.tex`): now "dot positions: 30 33 (top), 35 37 (bottom)". Crop p3.png.
2. P3 generator glyph sine (`module3.tex`): cause found, the plot used `\x`, which the outer `\foreach` also owns. The plot now uses its own variable `\tt` and stays inside the inner ring. The sine renders (p3.png).
3. P9 failure causes (`module9.tex`): "SG reserve" next to the 52 ms n=0 red dot; "voltage band / integration stop" next to the 52 ms n=3 dots; "|Δf| > 0.5 Hz" was already there. To make room, the red note box now starts at x=1060 and its text is 19.5 pt. Source comment: E6_events.csv. Crop p9b.png.
4. P7 spacing (`module7.tex`): ring illustration moved up 4 px (centre y 918, still 5.5 px below the bar). NoteBox, chip, HeroBox, PLL card and strip all moved up 4-6 px. Gaps are now HeroBox to card 5 px, card to strip 5 px, strip to panel border 3 px (was 2). p7.png.
5. P5 sketch (`module5.tex`): a small network with two A-B paths. The teal upper path has an arrow and the label "damping ↑"; the gold lower path has an arrow and the label "damping ↓". The old +/- discs are gone. p5b.png.
6. Stray fonts (`refstyle.tex`): `\setsansfont` is now Fira Sans (the sans fallback was LM Sans). `pdffonts` now lists only Fira Sans, Bahnschrift and Latin Modern Math.

## B. Print checks
7. Distance to borders: P6 chip and "abort" note were within 1-2 px of the panel border. The evidence cards were shortened to 736-776.5 and the chip moved to y=784.5 (p6.png). P7 fixed as in item 4. Everything else is at least 3 px away.
8. Text size: no text below 15.5 pt apart from icon internals. The only 15.5 pt uses are pills, and module2 lines 42-43 are tag-like.
9. Hairlines: art.tex line 32 and skyline.tex lines 6 and 13 (0.4 pt / 0.45 pt) raised to 0.5 pt. All other strokes are at least 0.7 pt.
10. Colour on colour: the P3 contract pill was white on PGrey; its fill is now PMuted. The P9 "--" placeholder was PMist and is now PMuted. Gold text on mint or cream ("negative", "damping ↓") is now PGoldDeep!78!black. The rest was already compliant.
11. Margins: boxes had side margins of 8-11 px. P2, P3 and P7 boxes and chips were moved to 10 px (P2 415/746, P3 775/1139). P4, P5 and P6 were already 9-10. P1 (11/8), P8 and P9 (11/9-11) were left alone, because their content is built against those edges. NextChip y was aligned to about 9 px above the panel bottom (P4 785.5 to 784.5, P6 787 to 784.5). Chips in P5, P7 and P8 sit inside cards or columns by design.
12. Spelling: "destabilising" to "destabilizing" (P4) and "neighbours" to "neighbors" (P8). No "linearise" found. Units: "0.6223 Hz" in P3 is now a thin-space text unit. "4.91 Hz", "44 ms" etc. already use thin spaces.
13. Plots: all labels are at least 17.5 pt. P9 plot (`figures/module9/make_m9_plot.py`): the "unstable" label sat on the dashed crossing line and is now to its right, regenerated at 231x104. P5 plot: added "symlog axis" to its subtitle (the axis has no unit), regenerated at 331x108. P8 plots were fine (the m=1 / m=2,3 labels do not overlap).

## Not changed / open
- Official high-resolution IAS logo (127 ppi; see PRINT_PREFLIGHT.md).
- Strip texts start about 6 px from their panel edges, which is acceptable.
- A real-size proof is still needed for the print colours.
