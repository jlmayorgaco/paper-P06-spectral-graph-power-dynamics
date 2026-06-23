# IEEE IAS Poster Session Package

This folder contains a one-page A0 portrait IAS poster using Option B:
the main story is the line-ranking damping-margin witness, with the modal
substitutability result inserted as a focused "full lever space" block.

## Main Message

Robustness planning is choosing in a structured space of levers. A frequency
ranking can pick a line that raises critical modal stiffness while reducing
damping margin; the structural block then shows where line actions sit in the
larger lever space of existing lines, new lines, and shunt/virtual-inertia
actions.
The structural statement used in the block is: virtual-inertia perturbations
live in existing-line directions plus a shunt/grounding residual, not in the
new-topology subspace.  The shunt residual is what makes inertia or a
synchronous condenser hard to replace at low-substitutability nodes.

## Outputs

- `poster_A0_ieee_session.pdf` - compiled A0 portrait poster (`841 x 1189 mm`).
- `poster.tex` - editable LaTeX source.
- `build/poster_preview.png` - rendered preview used for visual checks.
- `figures/fig_witness_bar_ranking.*` - dominant line-ranking witness figure.
- `figures/fig_modal_substitutability.*` - modal substitutability map used in
  the inserted structural block.
- `figures/fig_opposite_levers.*` - secondary modal-substitutability figure,
  generated for reuse but not currently dominant in the Option-B poster.
- `logos/ias_annual_2026_logo.png` - official IEEE IAS Annual Meeting 2026
  logo used in the poster header.

## Rebuild

From this folder:

```powershell
.\compile.ps1
```

Manual equivalent:

```powershell
python .\scripts\generate_demo_and_figures.py
python .\scripts\generate_enhanced_poster_assets.py
python .\scripts\generate_modal_substitutability_assets.py
python .\scripts\verify_inertia_lever_space.py
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build poster.tex
Copy-Item -LiteralPath .\build\poster.pdf -Destination .\poster_A0_ieee_session.pdf -Force
```

## Line-Ranking Witness Data

- Master seed: `20260608`.
- Case seed: `727986944`.
- Network: 6 buses, 10 tested line reinforcements.
- Frequency-top line: `1-3`.
- Critical stiffness movement: `Delta nu_c = +0.0931127997`.
- Full-QEP damping margin movement:
  `zeta_min: 0.0420405828 -> 0.0406764855`,
  `Delta zeta_min = -0.0013640973`.
- Damping-aware top line: `2-5`, with
  `Delta zeta_min = +0.0004510725`.
- Reversal count: `5/10` tested lines raise stiffness while lowering damping
  margin.

## Modal-Substitutability Block Data

- Test system: 6 buses.
- Critical modal stiffness: `nu_c = 7.51`.
- Modal substitutability:
  `s = [0.54, 0.14, 0.58, 0.36, 0.90, 0.50]`.
- Classifier spread: `8.9x`.
- Witness contrast:
  - Node 4: `s = 0.90`, topology-first candidate.
  - Node 1: `s = 0.14`, inertia/SynCon candidate.
- Structural placement check:
  `scripts/verify_inertia_lever_space.py` tests 40 random connected graphs
  (`n = 4..11`) and verifies that inertia-induced perturbations are representable
  by existing-edge Laplacians plus diagonal shunts, with no new-edge component.
  Current maximum relative residual: `5.670191096847199e-15`.

## Claim Discipline

The poster does not claim the algebraic edge/diagonal decomposition as new. It
is presented as a taxonomy. The useful structural result is that virtual-inertia
perturbations live in existing-line directions plus a shunt/grounding residual,
not in the new-topology subspace.

The poster also does not claim final calibrated IEEE-39 / ANDES IBR validation.
The network-inertia structure is scoped to the reduced subsystem before control
self-energy effects dominate.
Multi-graph robustness of the substitutability classifier remains a next gate.
The current poster reports the six-bus classifier witness and the multi-graph
structural placement check.

## Design Notes

The visual hierarchy is:

1. Title and hook: one planning question.
2. Left column: model and why frequency-only line ranking fails.
3. Center column: dominant witness figure and interpreted numbers.
4. Right column: full lever space, modal substitutability map, and workflow.
5. Bottom: what is shown, what is not claimed, and the 30-second story.

Color semantics:

- Green: acceptable topology/line lever.
- Red: rejected/hard-to-replace lever, inertia/SynCon where relevant.
- Black: hook and 30-second story bands.
- IAS green: institutional headers, footer, rules, and logo alignment. The
  dominant green is sampled from the official IAS Annual Meeting 2026 logo
  (`RGB 0,128,64`).

## Branding and Size Checks

- Poster size verified with `pdfinfo`:
  `2383.94 x 3370.39 pts (A0)`, 1 page.
- IAS 2026 official site confirms the event dates, Vancouver location, and
  student poster competition context.
- IAS Branding & Communications confirms that official IAS logos are available
  for print/event use and should be used consistently with IAS branding.

Sources checked:

- 2026 IEEE IAS Annual Meeting:
  `https://ias-am.ieee.org/2026/`
- Student Poster Competition Program:
  `https://ias-am.ieee.org/2026/student-poster-competition-program/`
- IEEE IAS Branding & Communications:
  `https://ias.ieee.org/about/branding-communications/`
