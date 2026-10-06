# Jury review + author requests (reviewer, 2026-10-06)
Author: "no QR"; "remove SCOPE OF THE RESULTS and use the space for better things"; "remove the final phrase — too catchy, too ChatGPT"; earlier: "floating big numbers with no meaning look AI-generic" (fixed in P3 only).
Rule for the whole poster: factual, technical tone. No slogans, no marketing taglines, no rhetorical questions. Science, numbers and `% source:` comments unchanged.

## 1. Footer (main.tex + art.tex)
Remove both footer sentences ("PRESERVING ONE MODE…" and "More information gives…"). Replace the 84 px green band by a slim 30–34 px deep-green band at the page bottom carrying only factual identity text in Fira Medium 19 pt white/pale: left "IEEE IAS Annual Meeting 2026 · Vancouver, Canada · October 4–8, 2026"; right "Jorge Luis Mayorga Taborda · Universidad de los Andes · jl.mayorga236@uniandes.edu.co". Keep the thin gold top rule. Remove the footer mountains/pylons/trees art (FooterArtLeft/Right) — they read as stock decoration. Give the freed height to the lower strip.

## 2. Lower strip (module10 + the three \RefSection bars in main.tex)
Remove SCOPE OF THE RESULTS. New strip (taller, using the freed footer height):
- MAIN CONTRIBUTIONS (left, ≈ 45 % width): the five items, one line each where possible, icons kept.
- KEY RESULTS (middle, ≈ 30 % width): a compact 4-row table, Fira Medium 18 pt, columns "model | finding | evidence": TX4 IEEE-39 | H4 blocker, 15/15 proper subsets stable, α⊥(H4) = +0.127 s⁻¹ | full-spectrum, fixed policy; full IEEE-39, exact delay | first margin crossing 41.57 ms; uniform K_p −14 % restores margin, 5/5 events | argument-principle counts + Julia events; n-hop, full IEEE-39 | n_min = 0,1,1,2 for m = 1–4; 11/12 exact assignments leave roots beyond the margin | modal equations + full counts; δD_H | 60/60 cases with one positive and one negative eigenvalue | Julia validation. Take every number from the panel that already prints it (copy its `% source:` comment).
- REFERENCES / CONTACT (right): as now, full width of its section (no QR space).

## 3. Remove catchy story devices
- Delete every \NextChip (all panels).
- Replace slogan-like sentences by plain technical statements (same meaning), e.g. P4 "Nothing failed locally. The network closed the destabilizing loop." → "All local factors stay regular; the instability enters through det(I+Q_S)." P5 hero line → "One local change raises damping along one direction and lowers it along another; damping at a single frequency is not a stability test." P8 "More modal authority ≠ a safer grid." → "Exact assignment of selected modes does not guarantee full-spectrum stability." P6 "Preserving one mode does not preserve the others." keep as "Preserving one mode exactly does not preserve the others." P1 bottom callout keep (it is the problem statement) but phrase as "Safe SG→GFL replacement is a joint portfolio and controller design problem."
## 4. Floating big numbers (same issue the author flagged in P3)
- P8 "11/12": remove the BigNumber; put the fact inside the plot caption or a short statement line ("11 of 12 exact assignments leave roots beyond the −0.05 s⁻¹ margin (no retune: 10).").
- P9 "41.57 ms": remove the BigNumber; annotate the crossing directly on the plot (gold marker + label "41.57 ms") and keep one caption line.
## 5. Jury detail fixes
- P2: the STABLE / UNSTABLE pills collide with the HeroBox edge — give them padding/space or write "α⊥ < 0 stable, α⊥ > 0 unstable" in the box as text.
- P4: the definition lines at the top are dense; increase line spacing by ≥ 3 px and keep the audit note text ≥ 18 pt.
- P5: lollipop y axis: make the symlog ticks readable (e.g. −1000, −10, 0, 10) and label the axis "eigenvalues of δD_H (symlog)".
- P7: left-side caption texts ≥ 18 pt; ensure ring labels do not touch nodes.
- Everywhere: badge/pill padding ≥ 2 pt (already set globally); no text closer than 3 px to box edges.
