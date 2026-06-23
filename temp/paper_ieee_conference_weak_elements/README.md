# IEEE Conference Paper: Weak Nodes and Weak Links

This folder contains a focused 6-page IEEE conference draft:

**Weak Nodes and Weak Links Are Not Graph-Centrality Objects in IBR Grids**

Scope:
- Uses the validated Phase-0D Schur-NEP bridge.
- Uses the registered Phase Braess-NEP line audit.
- Claims only the controller-aware weak-element certificate and the modal-local control-resonant Braess mechanism.
- Does not claim full estimator/planning validation.
- Includes a native TikZ mechanism diagram in `main.tex` explaining how line reinforcement can move a network-family pole toward a PLL/control pole and reduce the hidden damping margin.

Build:

```powershell
python make_application_figures.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Note: On this Windows/MiKTeX install, `pdflatex` writes `main.pdf` correctly but returns a nonzero process code because MiKTeX cannot write its global log under `AppData\Local\MiKTeX\miktex\log`. The local `main.log` has no undefined references or overfull boxes after the final build.

Figures are copied from:
- `outputs/phase0d_nep_contour_solver/fig_phase0d_s_plane.pdf`
- `outputs/phase_braess_nep/fig_braess_nep_s_plane.pdf`

Application figures are generated from:
- `outputs/phase_braess_nep/phase_braess_nep_line_decomposition.csv`
- `outputs/phase_braess_nep/phase_braess_nep_status.json`

Generated application artifacts:
- `figures/fig_application_workflow.pdf`
- `figures/fig_audit_mechanism_summary.pdf`
- `figures/fig_application_weaklink_report.pdf`
- `application_report_weaklink_line20_34.md`
