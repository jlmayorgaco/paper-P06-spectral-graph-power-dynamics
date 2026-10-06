# DESIGN v3 — from "robotic template" to "illustrated editorial poster"
Author verdict on v2: "looks robotic, generic, lacks detail". Reviewer diagnosis (from the full 75 dpi render):
1. Uniformity: every panel uses the same mint card + cream box pattern at the same weight → no rhythm, no focal point; the eye has nowhere to land.
2. Flat, schematic drawings: circles, boxes and dots; no recognisable power-system objects (generators, inverters, lines, PLL loops).
3. Too much small text of equal weight; captions and legends look like a data dump.
4. Plots read as default charts dropped into boxes (thin bars, cramped labels, legends in corners).
5. No visual story across panels: nothing links 1 → 9.
Keep: page geometry, header, footer, palette, fonts (Fira Sans / Fira Math / Bahnschrift), all science and numbers, all `% source:` comments.

## 1. Icon library (new file icons.tex — one consistent hand-drawn-quality set, 1.6–2.2 pt strokes, round joins, two-tone: deep green outline + mint/cream fill + one gold accent)
\IconSG (turbine–generator: circle with rotor + shaft + sine), \IconGFL (inverter cabinet with IGBT switch + DC capacitor + output sine), \IconPLL (small loop block: detector ⊗ → PI → ∫ → θ with feedback), \IconBus (bar), \IconLine (transmission line with pylon), \IconNetwork (mesh of 6 buses), \IconClock (delay τ), \IconSpectrum (complex plane with roots and margin line), \IconCheck / \IconCross / \IconWarn (filled discs), \IconWave (frequency event), \IconGraphHop (node + rings), \IconMatrix, \IconFunnel (Schur elimination), \IconLock. Each macro takes (x, y, scale). Use them everywhere instead of plain circles/boxes.

## 2. Visual hierarchy per panel (one focal element, then support)
Every panel gets exactly ONE "hero visual" occupying 35–50 % of the body (illustration or plot), ONE HeroBox (key equation/result), and at most THREE supporting elements. Remove or merge anything beyond that (move fine print to `% notes` comments, not the poster).
- P1 hero: large illustrated SG → GFL replacement scene (real icons, network between them, PLL on the GFL), objective in the HeroBox.
- P2 hero: a vertical "zoom" illustration: full grid → linearised Jacobian blocks → eliminated network (funnel) → transverse quotient (plane with gauge arrow); equations sit beside each stage.
- P3 hero: IEEE-39 drawn with generator ICONS at buses 30–39 (H4 red, others green) + the lattice "cliff" (glyphs rising to the red top); one number in HeroBox.
- P4 hero: two GFL icons with their own PLL loops (teal arcs) and the network between them carrying the red closing loop — make it a real scene, bigger, with labels on the arrows.
- P5 hero: the δD_H plot redrawn as a "diverging lollipop" per site (stems from zero, teal up, gold down) with site numbers on the axis; tiny network sketch showing + and − paths.
- P6 hero: a pipeline of four illustrated stations (open port icon → row equation → W_i dial → gains Kp, KI knobs) joined by a thick gold flow arrow.
- P7 hero: the hop rings drawn ON the IEEE-39-like communication graph with GFL icons at nodes; ring labels n = 0, 1, 2, 3 on a curved path.
- P8 hero: spillover plot as large bars with direct annotations; authority plot smaller; one-line repair result.
- P9 hero: the delay root locus with the crossing highlighted (gold marker + label), nonlinear events as a clean dot matrix with row/column headers.

## 3. Rhythm and depth
- Alternate body backgrounds subtly: odd panels white #FFFFFF, even panels #FAFCFB (barely visible) — gives column rhythm.
- Section labels inside panels: small caps gold-chip letter + SemiBold heading, with a 0.8 pt hairline rule extending to the right edge (editorial look).
- HeroBox: keep cream + left gold bar, add 6 px inner padding and a soft drop shadow (TikZ `drop shadow={opacity=0.12, shadow xshift=0.6pt, shadow yshift=-0.6pt}`) — ONLY on HeroBox.
- Big-number callouts: when a panel's key result is a number (P3 +0.127, P8 11/12, P9 41.57 ms), typeset it as a large Bahnschrift numeral (60–72 pt) with a small Fira label beside it — magazine style.
- Story thread: a small gold "→ next" connector chip at the bottom-right of panels 1–8 stating the bridge in ≤ 5 words (e.g. P3 "Why? → closure (4)", P4 "Who repairs it? → gains (6)").

## 4. Plots
Thicker everything (lines 3 pt, bars width 0.7, markers 10), Fira Sans 19 pt, direct labels on data, no boxed legends, light horizontal grid only, highlight colour only on the data point that matters (others grey/teal).

## 5. Text discipline
Body ≥ 21 pt; captions 18 pt Medium; no paragraph longer than 2 lines; no captions that repeat what the illustration says.
