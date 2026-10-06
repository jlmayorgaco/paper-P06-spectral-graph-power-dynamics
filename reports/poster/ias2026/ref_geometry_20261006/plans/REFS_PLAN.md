# References plan (reviewer, 2026-10-06) — all entries verified online
Author request: drop the extra/weak citation, add 5 more real references in an additional two-column block, and print under every module which references it uses.

## Final reference list (numbering is final)
[1] P. Kundur, *Power System Stability and Control*, McGraw-Hill, 1994.
[2] F. Dörfler, M. R. Jovanović, M. Chertkov, F. Bullo, "Sparsity-promoting optimal wide-area control of power networks," IEEE Trans. Power Syst., vol. 29, no. 5, pp. 2281–2291, 2014.
[3] L. Huang, H. Xin, W. Dong, F. Dörfler, "Impacts of grid structure on PLL-synchronization stability of converter-integrated power systems," arXiv:1903.05489.
[4] W. Michiels, S. Gumussoy, "Eigenvalue based algorithms and software for the design of fixed-order stabilizing controllers for interconnected systems with time-delays," arXiv:2003.05496.
[5] U. Marković, O. Stanojev, P. Aristidou, E. Vrettos, D. S. Callaway, G. Hug, "Understanding small-signal stability of low-inertia systems," IEEE Trans. Power Syst., vol. 36, no. 5, pp. 3997–4017, 2021.
[6] X. Wang, M. G. Taul, H. Wu, Y. Liao, F. Blaabjerg, L. Harnefors, "Grid-synchronization stability of converter-based resources — an overview," IEEE Open J. Ind. Appl., vol. 1, pp. 115–134, 2020.
[7] F. Dörfler, F. Bullo, "Kron reduction of graphs with applications to electrical networks," IEEE Trans. Circuits Syst. I, vol. 60, no. 1, pp. 150–163, 2013.
[8] S. Skogestad, I. Postlethwaite, *Multivariable Feedback Control: Analysis and Design*, 2nd ed., Wiley, 2005.
[9] A. Ortega, P. Frossard, J. Kovačević, J. M. F. Moura, P. Vandergheynst, "Graph signal processing: overview, challenges, and applications," Proc. IEEE, vol. 106, no. 5, pp. 808–828, 2018.
Removed: "J. L. Mayorga T., Beyond nodal damping, in prep." (not a citable source).

## Layout
- REFERENCES section: widen it and set the list in TWO COLUMNS (≈ 5 entries per column), Fira Medium 17 pt, hanging indent for the [n] label, each entry ≤ 3 lines; contact line (name · Universidad de los Andes · jl.mayorga236@uniandes.edu.co) stays at the bottom of the section. Rebalance the three section widths of the lower strip as needed (e.g. contributions ≈ 38 %, key results ≈ 28 %, references ≈ 34 %) and, if necessary, increase the strip height by taking up to 14 px from the empty space of the footer/bottom margin (keep ≥ 6 mm bottom margin).
- Under EVERY panel 1–9 add one muted line at the bottom of the panel body, Fira Medium 17 pt, colour PMuted, left-aligned on the panel text edge: `Background: [..], [..]  ·  results: this work` — use exactly:
  P1 Background: [1] [5] [6] · P2 [1] [7] · P3 [5] · P4 [7] [8] · P5 [8] · P6 [3] [6] · P7 [2] [9] · P8 [2] · P9 [4] [9].
  Format "Refs: [1] [5] [6]" (no "results: this work" text needed). These lines need ~14 px of vertical room per panel: take it from the panel's existing slack, never by shrinking text below 17 pt (move the panel's bottom-most element up or tighten spacing). Panels whose body is full (P1, P6, P8, P9) may need small recomposition.
- The lines must not collide with any NextChip (all NextChips were removed in JURY_FIXES) or with the panel border (≥ 3 px).

Source/consistency rule: keep every number and `% source:` comment unchanged; do not change science.
