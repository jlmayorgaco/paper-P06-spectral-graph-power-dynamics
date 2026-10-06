# Final polish checklist (reviewer, 2026-10-06) — print-oriented
Open items left from earlier reviews + new print checks. Science, numbers and `% source:` comments must not change.

## A. Pending from earlier reviews
1. P3 lattice caption: say which dot position = which bus (e.g. "dot positions: 30 33 / 35 37"), currently "grid 30 33 / 35 37" is cryptic.
2. P3 generator glyphs: the internal sine of the generator symbols does not render — fix or remove the sine cleanly.
3. P9 events matrix: label the failure causes in small text next to the red dots: 52 ms n=0 "SG reserve", 52 ms n=3 "voltage band / integration stop", 44 ms 90 % "|Δf| > 0.5 Hz" (already), source research/nhop_ieee39_20261006/derived/E6_events.csv.
4. P7: the equation HeroBox and the PLL-law card are cramped — add ≥ 4 px vertical breathing room (shrink the ring illustration slightly if needed).
5. P5: the + / − path sketch is plain — make it a small but clear network with a teal arrow "damping ↑" and a gold arrow "damping ↓".
6. Stray fonts: pdffonts shows LMSans12/LMSans17 (Latin Modern Sans text). Find where text falls back to the default text font (likely \textrm/\mathrm/\text or a font switch inside math) and make it Fira Sans.

## B. Print checks (do on 150 dpi crops of every panel + strip + header + footer)
7. Nothing closer than 3 px (2.4 mm) to a panel border or the title bar; no text touching icons, plot frames or other text.
8. Min text 17 pt for any readable word (pills 15.5 pt OK); icon-internal labels may stay smaller.
9. Hairlines: no stroke thinner than 0.5 pt except plot grids (thin lines disappear on large-format printers).
10. Colour on colour: text on mint/cream boxes must be PText/PGreen/PNavy/PRed (no grey lighter than PMuted on tinted fills); white text only on PGreen/PRed/PNavy/PTeal fills.
11. Consistent spacing: same inner margin (≈11 px) on all panels; NextChips aligned to the same bottom-right offset in every panel; HeroBoxes aligned to the panel text edge.
12. Typos/consistency: "linearise" vs "linearize" — use US spelling everywhere (IEEE); "SG→GFL" arrows consistent; units with thin spaces (44 ms, 0.6223 Hz, s⁻¹).
13. Plots: axis labels with units, tick labels ≥ 17 pt, no overlapping annotations.

Deliverable: fixed module files, a list of every change, and before/after crops for anything visible.
