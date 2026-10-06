# IAS26-020 — final F1 audit

Status: **PASS**

RUN_ID: `20260926T161116Z_d0fecb32_ias26_020_f1_audit_v1`. The exact source NPZ was reused byte-for-byte: 16 portfolios, 1,216 transverse eigenvalues and eigenvectors. Decision threshold `tau_dec=1e-08 s^-1` was read from the renderer at the clean source commit, not selected after inspecting the spectrum.

The full-spectrum gate gives 15 stable proper subsets, 1 unstable H4, and 0 indeterminate portfolios. `alpha_perp` and `alpha_Omega` agree for all 16. H4 has alpha `0.127006467828 s^-1`, frequency `0.62227967 Hz`, 86 DAE states and 84 transverse coordinates.

The eigenpair residuals were checked against reconstructed transverse operators without re-running an eigensolver. The figure is explicitly Python-only for the V4 lattice; stored Python–Julia parity covers H4 GFL11 only. V9 is not included in the P1 claim. No F2, F3, F4, Monte Carlo, or TDS work was performed.
