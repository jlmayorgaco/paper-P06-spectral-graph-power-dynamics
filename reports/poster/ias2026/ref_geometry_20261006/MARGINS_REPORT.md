# Margins report (build_margins/main.pdf copied to main.pdf)

Measured on the 150 dpi render (non-white extent), 1 page, 914.4 x 1066.8 mm:
before: left 0.0 (hills/skyline at the page edge), right 0.0, top 8.8, bottom 6.3 mm.
after: left 15.9, right 15.9, top 16.8, bottom 14.1 mm.

## A) Outer margins
- main.tex: everything inside OnPage except the white background fill is in one scope
  `scale around={0.985:(580.5,677.5)}` with `every node/.append style={transform shape}` (text scales too; smallest 17 pt text becomes 16.7 pt).
- Header hills and skyline are clipped to x in [12,1149] (main.tex scope around \HillsLeft/\SkylineRight).

## B) Column alignment (row 3 = row 2 edges)
- Panels 7/8/9 boxes now 13-396, 405-756, 765-1148; panel 8 title start 444 / width 306, panel 9 title start 805.
- Module 7: widened by 6 px (right edge 385): section rule, note, hero box, PLL card, tag, four architecture cards (width 88.75, pitch 90.75).
- Module 8: left text edge 416 (was 409), width 329 (was 336): plots regenerated at 329x58 and 329x102 px (figures/module8), legend, statement, repair card (parbox 314), hero box and Refs re-positioned.
- Module 9: box +1 px on the left; content unchanged (no clipping).

## C) Air above the Refs lines (ink to ink, 150 dpi, reference px)
before: P1 4.4, P2 9.2, P3 5.9, P4 5.2, P5 4.6, P6 5.0, P7 3.3 (0.2 incl. bracket), P8 3.1, P9 4.8.
after: P1 8.7, P2 10.0, P3 8.7, P4 9.8, P5 8.7, P6 8.7, P7 9.0, P8 9.8, P9 8.7.
Changes: Refs lines moved down 1 px in all panels; P1 scene/hero/chips/note lifted 3-6 px, note box 1.5 px shorter, symbol key raised to 17 pt;
P3 hero box 2 px shorter, scope note up 2; P4 hero/audit blocks up 2-3, closing sentence up 3.5; P5 closing card up 3; P6 evidence tiles 39 px tall (was 40.5);
P7 PLL card shorter, architecture cards up 5; P8 spillover plot 102 px high (was 108), authority group up 6; P9 red note box 50 px tall (was 54), dot matrix pitch 9.8 px (was 10.8).
No numbers or science text were removed.

## Notes
- Existing sub-17 pt items not touched: icon-internal glyphs and some `\SubHead` chips at 15.5 pt (tag pills).
- Bottom margin (14.1 mm) is the tightest; scale 0.985 is the minimum satisfying 14 mm on all sides.
