# REPORT: T1, E1, E2, E5 (2026-10-06). Not committed. Preregistration untouched; one dated amendment (AMENDMENT_01.md, written before E1/E2/E5).
Code: `src/nhop.py` (T1 library), `src/t1_tests.py`, `src/e1_e2.py`, `src/e5_spillover.py`, `src/figs.py`. Logs `logs/`. Tables `derived/`. Figures `figures/`.
(Bug caught and fixed during the run: `fm.py` defines its own `CAMP`; first outputs landed in the feedback_cycle folder as untracked files and were moved here; nothing tracked there was changed. Results unaffected.)

## T1 tests (`derived/T1_tests.json`, `logs/t1_tests.log`)
- Parity diag(Kp,KI)=baseline vs `Model.delta`: max rel. err 1.5e-18 (10 random s, tau 40 and 44 ms). PASS (<=1e-12). Refinement identical to `Model.refine` on the 44 ms critical roots.
- Counter (same bands as t3): 40 ms -> N_unstable 0, N_margin 0; 44 ms -> 8 unstable / 10 beyond -0.05 (gauge root excluded). PASS, non-integer winding <= 7e-12. ~43 s per (unstable+margin) count.
- Open-PLL operator: omega and xi do NOT enter the rest dynamics (A0 rows of rest, cols omega/xi: max 0.0; C cols omega/xi: 0.0; PLL rows from rest: 0.0). G(s)=C_theta + C_r (sI-A_rr)^-1 A_r,theta (174 rest states).
- Determinant identity, 24 random s, random dense Kp, KI, random tau: det Delta = (600 pi)^10 * det(sI-A_rr) * det(s^2(1+t_f s)I - (s Kp+KI) E G); constant c = t_f^-10 derived analytically (Schur on the rest block, then on omega/xi); max rel. err 1.1e-11. PASS (<=1e-8).
- Y.csv block convention: 2x2 blocks have form [[a,-b],[b,a]] (violation 1.4e-14) -> a+jb; Yc symmetric to 7e-15; Z = Yc^-1 (rcond 0.11). r_ij uses a modulus so it is invariant to the sign choice of b.

## G_c (`derived/GRAPH_Gc.csv`, `figures/Gc.pdf`) -- PROBLEM
The frozen rule gives a DISCONNECTED graph (13 edges): components {30,31,32,37,38,39} and {33,34,35,36}; diam undefined (infinite). See AMENDMENT_01: primary = frozen graph with n = 0..3 (largest finite hop distance, within-component); sensitivity = "bridge" (adds the min-r cross edge 30-33, `GRAPH_Gc_bridge.csv`), diam 4. All tables carry a `graph` column (frozen | bridge).

## E1 (`derived/E1_targets.csv`, `E1_residuals.csv`, `E1_nmin.csv`, `figures/E1_residual_vs_hops.pdf`)
Targets (44 ms; PLL-family = PLL participation >= 0.5 in the Delta null vector; rightmost 5): Re/f = +0.969/4.649 Hz, +0.929/4.681, +0.677/4.713, +0.254/4.566, -0.018/4.637; all participation >= 0.975.
n_min(m) (eps<=1e-8), frozen: m=1:0, m=2:1, m=3:1, m=4:2, m=5: NOT reached (eps 4.8e-4 at n=3). Bridge: 0,1,1,2,2.
eps(n) non-increasing in n (nesting) in all cases; below threshold it is ~1e-15 (machine level) or >=1e-3, no intermediate values.
Interpretation / doubt: the sharp transition is dimension counting (2m real equations vs 2|N_i| unknowns per site; N_i smallest sizes 1,3,4,4 for n=0..3 on the frozen graph), not a subtle spectral effect. "n_min increases with m for >=2 m": SUPPORTED (0,1,1,2 / frozen), but it is largely generic counting; the frozen graph's cut leaves m=5 unreachable.
## E2 (`derived/E2_architectures.csv`, `figures/E2_architectures.pdf`, shown for m=4)
Free n-hop reaches ~1e-15 at n>=2 (m=4). Shared polynomial (b) stays at 1e-2..7e-2 for all n and m (never <=1e-8; bridged and frozen). Node-varying polynomial (c) equals free n-hop only in low cases (m=1 all n; reaches 1e-15 for m=2 n>=1, bridge m=3 n>=2, m=4 n>=3, m=5 n=4 bridge), otherwise 1e-3..2e-2 (e.g. frozen m=3 n>=2 3.6e-3, frozen m=5 n=3 5e-3). Polynomials are subsets of free n-hop, consistent.

## E5 (`derived/E5_spillover.csv`, `figures/E5_spillover.pdf`; `_bridge` variants too)
62 exact assignments (eps<=1e-8; 2 variants each). Variant "absolute" = literal min-norm row solution (baseline diagonal NOT retained; gains up to 3.9e3 for KI); variant "delta" = baseline + min-norm change (my reading of "smallest" intervention; the preregistration is ambiguous). Targets verified exact in every case (sigma_min of reduced char. at lambda* <= 5.5e-13, refine dev 0).
Roots beyond -0.05 left (baseline 10, unstable 8): 56/62 assignments leave roots beyond the margin ("safer grid" REFUTED in these cases); N_margin = 0 only for: m=1,n=0 (both graphs, both variants) and bridge m=5,n=3,4 (delta variant). Frozen graph delta: best N_margin 2 (m=4). Exact assignment of more modes with absolute gains often worsens spectrum (up to 18 beyond margin, up to 11 unstable). Counts floating-point, not certified.

## Caveats
Counter not certified; roots refined with exact exponentials only (no Pade). E1/E2 conclusions depend on the amended graph treatment and on q taken from the baseline null vector at fixed shape.
