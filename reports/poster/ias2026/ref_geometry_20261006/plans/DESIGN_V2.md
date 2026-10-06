# DESIGN SYSTEM v2 — "designer-made, print-strong, not generic"
Reviewer diagnosis of v1: (1) Latin Modern math is hairline and reads as default LaTeX; (2) too many tinted boxes (gold, mint, pink) with 0.5 pt outlines — the classic AI-template look; (3) icons and diagram strokes too thin for print; (4) plots in default matplotlib style; (5) small grey text too light; (6) uneven internal spacing.
Goal: one coherent typographic family, fewer and calmer containers, heavier strokes, exact alignment — like a layout made in Illustrator/InDesign. The header, footer, geometry, panel bars, palette and science stay as they are.

## 1. Typography (all real fonts, no stretching)
- Text family: **Fira Sans** (MiKTeX, OTF): Regular = body paragraphs; Medium = every text ≤ 19 pt (captions, labels, legends) so it survives print; SemiBold = subheads and emphasis; Bold = callout sentences and tags. Load with fontspec by file name (FiraSans-Regular.otf, -Italic, -Medium, -MediumItalic, -SemiBold, -Bold, -ExtraBold).
- Math: **Fira Math** (`\setmathfont{FiraMath-Regular.otf}`) — monoline sans math that matches Fira Sans. Math colour **navy #1E4F8A** (replace PBlue usage in equations; keep a PBlue alias pointing to the navy).
- Display: keep Bahnschrift condensed (title, panel titles, footer). Replace Arial Narrow in panel titles by Bahnschrift wght 700 / wdth 75, uppercase, letter-spacing +1.5 %.
- Sizes (physical pt): body 21–22.5; captions/labels 17.5–19 (Medium); subheads 23–24 SemiBold; equations 25–32, hero 34–38; tags 15.5 Bold all caps tracked +6 %. Nothing lighter than Regular; no grey lighter than #4F5D59.

## 2. Colour roles (no new colours)
Primary green #03534A, deep green #01463B, gold #EABD38 (accent only), navy #1E4F8A (math, data series 2), teal #2C7A70 (evidence, stable), red #B83A3A (ONLY "unstable / fails / limit"), text #1C2B27, muted #4F5D59, mint #EAF4EF, cream #FBF5E3, grid #E3EBE7.

## 3. Containers — only three kinds, NO hairline outlines
1. `\HeroBox{x0}{y0}{x1}{y1}` — cream fill, no outline, 3 px gold bar on the left, radius 3 px. At most one per panel (the panel's key equation or result).
2. `\CardBox{...}` — mint fill, no outline, radius 3 px. Step cards, evidence blocks.
3. `\NoteBox{...}{colour}` — white fill, no outline, 2.2 px left bar in colour (red = limit/scope, teal = evidence, navy = definition), heading in the same colour (Fira SemiBold).
Delete every pink/red-filled box and every 0.5 pt outline box inside panels (plots keep their axes). Keep the old macro names (\EquationCallout → HeroBox, \GreenCallout → CardBox, \WarningNote → NoteBox red) so existing modules compile, then restyle each module explicitly.
Tags: filled pills — EXACT (teal fill, white text), PROPOSED (gold fill, dark text), NUMERICAL (navy fill, white text), EXPLORATORY (grey #6B7773 fill, white text). Radius 4 px, padding 2 px.
Subhead style: Fira SemiBold 23 pt green, preceded by a 9 px gold square chip with the letter (A, B, C) in deep green — one consistent pattern for every lettered section.

## 4. Diagrams (TikZ)
Strokes: primary 1.8 pt, secondary 1.1 pt, never below 0.8 pt (except plot grids). Round caps/joins. Arrowheads Stealth 7 pt. Node discs ≥ 12 px, labels Fira Medium ≥ 17 pt. Fills only mint / cream / white; accents gold; semantic red only for failing elements. Flow strips = filled mint chips (no outline) joined by gold chevrons, labels Bahnschrift SemiBold caps.

## 5. Plots (matplotlib, generated at the exact slot size, 1:1)
rcParams: font Fira Sans (register the OTF files with matplotlib.font_manager), size 19, axes linewidth 1.4, spines top/right off, ticks outward 5 pt width 1.2, y-grid only colour #E3EBE7 width 0.9, lines 2.6, markers 8 with white edge 1.0, legend frameless; colours green #03534A, navy #1E4F8A, gold #D9A520, red #B83A3A, teal #2C7A70, grey #8A9894; `pdf.fonttype = 42`. Direct labelling preferred over legends. Units on every axis.

## 6. Spacing and alignment
Inner panel margin 11 px (reference px) left and right, 7 px below the bar. Block gap 6 px. All left edges of a panel align to one x. Captions 3 px below their figure, Medium 17.5 pt. Equations centred only inside HeroBox; elsewhere left-aligned on the panel's text edge.

## 7. Copy
Sentence case, active voice, no filler adjectives ("novel", "robust", "leverage"), no exclamation marks, at most two lines per caption, hyphenation off for headings.

## 8. Module-specific fixes
P1: icon strokes to 1.8 pt; δ meaning line in Medium; symbol key Medium 17.5 pt.
P2: flow strip → mint chips + gold chevrons; delay note → NoteBox navy.
P3: "15/15 …" → NoteBox red (white fill); Hasse colour map mint→deep teal for stable (darker = more negative α⊥), H4 node red; numbers Fira SemiBold; caption adds "∅ = no replacement"; one-line strokes 1.6 pt, generator discs 9 px, H4 discs red with white label.
P4: arrow diagram must connect ΔY_i×T₀⁻¹ to both "self return M_ii" and the Q_S box; audit box → NoteBox teal; LOCAL/COLLECTIVE boxes → CardBox (mint) and NoteBox red heading, no outlines.
P5: plot restyled per §5 with y label "eig δD_H [MW]" (check the unit in the CSV header/origin script; if unknown, no unit claim — write "eigenvalues of δD_H").
P6: step-2 derivation line gets 4 px more space; delay strip → NoteBox red; evidence table with Fira Medium, no shaded header bar (thin green rule instead).
