# PREREGISTRATION — information-radius (n-hop) PLL co-design on the full IEEE-39 model
Written 2026-10-06 before any n-hop computation. Frozen by `PREREGISTRATION.sha256`. Deviations only as dated `AMENDMENT_xx.md` files written BEFORE the affected run.

## Model (frozen, identical to research/feedback_cycle_fullmodel_20261005)
- Exported PowerDynamics/ReducedDAE linearisation, `experiments/graph_gsp_codesign_20261003/model/*.csv`, reader `experiments/interaction_decision_20261004/model.py`; 204 states; GFL sites = buses 30..39 (index i = bus-30); rho_i = 0.875; nominal Kp0 = 2*pi*5 = 31.4159, KI0 = (2*pi*5)^2/4 = 246.7401; baseline Kp_ii = 28.2743, KI_ii = 246.7401; t_f = 1/(600*pi) s; exact exponential delays (never Pade in reported roots).
- Unit injections: `M['Bp'][:,i]` (= 1/t_f at omega_i), `M['Bi'][:,i]` (= 1 at xi_i); detector rows `Model.C` (10 x 204).
- **Distributed PLL (new):** gain MATRICES Kp, KI (10x10, row = receiving PLL i, column = measured detector j):
  Delta(s) = sI - A0 - (Bp Kp + Bi KI) E(s) C,  E(s) = diag(exp(-s*(tau_j + d_ij))) -> with d_ij = 0 (baseline) E = diag(exp(-s tau_j)).
  Diagonal Kp, KI must reproduce `Model` exactly (parity test, rel. error <= 1e-12).
- Margin sigma_req = 0.05 s^-1. Counts: banded argument principle (`research_gold/checks/t3_full_spectrum_count.py` generalised), floating point, NOT certified.
- Stress cases: uniform tau = 44 ms (baseline gains: 8 unstable / 10 beyond margin), 48 ms (16/16), 52 ms.

## Communication graph G_c (frozen rule)
- Complex 10x10 port admittance Yc from the real 20x20 Kron matrix `model/Y.csv` (2x2 blocks [[a,-b],[b,a]] -> a+jb; agent must verify the block convention and record it). Z = Yc^{-1} (if singular: use the pseudo-inverse and record it).
- Electrical distance r_ij = |Z_ii + Z_jj - 2 Z_ij|. Edge i~j iff j is one of the 2 nearest sites of i OR i is one of the 2 nearest of j. Unweighted Laplacian L of G_c. Hop distance dist(i,j) on G_c; N_i^(n) = {j: dist(i,j) <= n}, n = 0..diam(G_c).
- Report the graph (edges, diameter) before using it.

## Experiments
E1 Modal authority (free n-hop). Targets: the m = 1..5 rightmost PLL-family roots at tau = 44 ms (baseline gains). For each target r: lambda_r* = -0.30 + j*Im(lambda_r); q_r = PLL-angle components (theta_i) of the right null vector of Delta at lambda_r (unit 2-norm). Open-PLL operator G(s) (theta -> e with every PLL loop open, grid closed) derived from the partition of the 204 states into PLL states (theta, omega, xi) x 10 and the rest; verify det Delta(s) = det(sI - A_rr) * det(s^2(1+t_f s) I - (s Kp + KI) E(s) G(s)) / normalisation at >= 20 random s (rel. error <= 1e-8). y_r = G(lambda_r*) q_r. Row conditions per site i:
  sum_{j in N_i^(n)} ( lambda z_ijr kp_ij + z_ijr kI_ij ) = b_ir,  z_ijr = exp(-lambda_r* tau_j) (y_r)_j,  b_ir = lambda_r*^2 (1 + t_f lambda_r*) (q_r)_i.
  Stack real/imag rows -> C_i^(n) k_i = b_i; residual eps_i(n) = ||(I - C C^+) b_i||_2 / ||b_i||_2; report eps(n) = max_i eps_i(n) and n_min(m) = smallest n with eps(n) <= 1e-8.
E2 Architecture comparison, same targets: (a) free n-hop (E1); (b) shared graph polynomial Kp = sum_{l<=n} a_l L^l, KI = sum b_l L^l (global LS); (c) node-varying polynomial Kp = sum diag(a_l) L^l, KI likewise. Residual vs n. Graph polynomials are NOT claimed to equal free n-hop.
E3 Closed-loop repair (predictor-corrector engine v2 of research/feedback_cycle_fullmodel_20261005/src/repair.py, generalised to matrix gains). Variables: diagonal entries in log scale within [0.25, 4] x nominal; off-diagonal entries kp_ij = Kp0*u_ij, kI_ij = KI0*v_ij with |u|,|v| <= 1, only for j in N_i^(n). Step cap 0.1, <= 60 iterations, step acceptance + best iterate. Success = full count N_margin = 0. Grid: tau in {44, 48, 52} ms, n = 0..diam. Report success, iterations, ||Delta K||_F (normalised), worst root.
E4 Replacement share: tau = 44 ms, rho uniform in {0.875, 0.90, 0.925, 0.95}; for each n the largest rho whose E3 repair succeeds (N_margin = 0). Label every point "largest verified point in the declared sweep, not a certified optimum".
E5 Spillover: apply the exact E1 assignments (n with eps <= 1e-8) directly (no corrector) and count roots beyond the margin. Expected possibility: exact assignment of m modes does not stabilise the spectrum.
E6 Nonlinear: for each successful E3/E4 design, the five frozen events (`DelayedEvents.CASES`, 60 s, guards 0.5 Hz, 0.5 Hz/s, V in [0.9,1.1], SG slack >= 0.002) with a distributed-gain Julia runner; diagonal-only parity with the existing runner first.

## Falsification statements
- "More information radius gives more modal authority" is SUPPORTED if eps(n) is non-increasing in n for free n-hop (guaranteed by nesting) AND n_min(m) increases with m for at least two m. 
- "More authority gives a safer grid" is REFUTED for a case if an exact assignment (E5) leaves roots beyond the margin.
- "Information radius extends the repairable region" is SUPPORTED only if some (tau or rho) case succeeds for n >= 1 and fails for n = 0 under identical engine settings.
Every outcome is reported, including negative ones.
