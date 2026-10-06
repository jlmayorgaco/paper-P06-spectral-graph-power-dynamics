# REPORT: E3 and E4 (2026-10-06). Not committed. Preregistration untouched; no AMENDMENT_02 needed (no setting changed after any result).

Code: `src/repair_nhop.py` (engine v2 generalised to matrix gains; `repair.py` of the feedback-cycle campaign unmodified), `src/run_e3_e4.py` (runner), `src/agg_e3_e4.py` (tables, figures).
Outputs: `derived/E3_repair.csv`, `derived/E4_replacement.csv`, `derived/E4_largest_rho.csv`, `figures/E3_repair.pdf`, `figures/E4_replacement.pdf`, `raw/runs/*.json` (per-iteration history), `raw/designs/<label>.toml` (34 successful designs; keys rho, tau, Kp, KI).
Counts are floating point, not certified. "worst Re" = worst catalogued root (Re > -1.5, Im > 0.3) of the best iterate. Success is decided only by the full-spectrum counter (N_margin = 0, gauge root excluded).

## Setup as run
Variables: log Kp_ii, log KI_ii in [0.25,4] x nominal; u_ij, v_ij in [-1,1] (kp_ij = Kp0 u_ij, kI_ij = KI0 v_ij), j in N_i^(n), j != i. Start = baseline (28.2743 / 246.7401, off-diagonal 0). Step cap 0.1 per variable (log units for diagonals, u/v units for off-diagonals), <= 60 iterations, trust-region halving, acceptance iff the worst catalogued root improves by > 1e-6, best iterate, stop at Re <= -0.0501 or radius < 1e-3. Effort = ||dKp||_F/Kp0 + ||dKI||_F/KI0 (full matrices, diagonal change included).
Sensitivities: exact simple-root formula, dlam/dK_ij = l^H B_i e^{-lam tau_j} (C r)_j / (l^H Delta'(lam) r). Finite-difference check (central, h = 1e-5, refine on perturbed model, at the rightmost root of the baseline, one diagonal Kp, one diagonal KI, one off-diagonal u, one v): max relative error over all 54 runs 9.9e-8 (limit 1e-5) -> PASS (`fd_max_rel_err` column).
Engine note: the n = 0 run reproduces diagonal-only repair (same code path with 0 off-diagonal variables). Cases (44 ms, 0.875) are shared by E3 and E4 (run once). 54 unique runs (27 E3 + 27 additional E4 rho != 0.875); 35-92 s each, about 68 min summed, 331 s wall time with 12 parallel workers.

## E3 (rho = 0.875): tau in {44, 48, 52} ms
| graph | n | tau ms | rho | success | N_unst | N_margin | iters | worst Re | effort | runtime s |
|---|---|---|---|---|---|---|---|---|---|---|
| bridge | 0 | 44 | 0.875 | yes | 0 | 0 | 2 | -0.0529 | 0.249 | 65 |
| bridge | 1 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0529 | 0.208 | 56 |
| bridge | 2 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0512 | 0.204 | 73 |
| bridge | 3 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0519 | 0.203 | 74 |
| bridge | 4 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0509 | 0.203 | 55 |
| bridge | 0 | 48 | 0.875 | yes | 0 | 0 | 4 | -0.0530 | 0.908 | 64 |
| bridge | 1 | 48 | 0.875 | yes | 0 | 0 | 7 | -0.0530 | 0.869 | 59 |
| bridge | 2 | 48 | 0.875 | yes | 0 | 0 | 5 | -0.0530 | 0.932 | 57 |
| bridge | 3 | 48 | 0.875 | yes | 0 | 0 | 5 | -0.0524 | 0.903 | 84 |
| bridge | 4 | 48 | 0.875 | yes | 0 | 0 | 5 | -0.0518 | 0.900 | 72 |
| bridge | 0 | 52 | 0.875 | yes | 0 | 0 | 5 | -0.0528 | 1.512 | 47 |
| bridge | 1 | 52 | 0.875 | yes | 0 | 0 | 13 | -0.0530 | 1.455 | 65 |
| bridge | 2 | 52 | 0.875 | yes | 0 | 0 | 7 | -0.0536 | 1.798 | 84 |
| bridge | 3 | 52 | 0.875 | yes | 0 | 0 | 8 | -0.0590 | 2.015 | 83 |
| bridge | 4 | 52 | 0.875 | **FAIL** | 10 | 10 | 13 | +0.8020 | 2.179 | 72 |
| frozen | 0 | 44 | 0.875 | yes | 0 | 0 | 2 | -0.0529 | 0.249 | 55 |
| frozen | 1 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0530 | 0.209 | 35 |
| frozen | 2 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0515 | 0.208 | 57 |
| frozen | 3 | 44 | 0.875 | yes | 0 | 0 | 3 | -0.0514 | 0.209 | 54 |
| frozen | 0 | 48 | 0.875 | yes | 0 | 0 | 4 | -0.0530 | 0.908 | 57 |
| frozen | 1 | 48 | 0.875 | yes | 0 | 0 | 6 | -0.0553 | 0.904 | 57 |
| frozen | 2 | 48 | 0.875 | yes | 0 | 0 | 6 | -0.0551 | 0.998 | 86 |
| frozen | 3 | 48 | 0.875 | yes | 0 | 0 | 10 | -0.0529 | 1.025 | 60 |
| frozen | 0 | 52 | 0.875 | yes | 0 | 0 | 5 | -0.0528 | 1.512 | 72 |
| frozen | 1 | 52 | 0.875 | yes | 0 | 0 | 12 | -0.0517 | 1.431 | 62 |
| frozen | 2 | 52 | 0.875 | **FAIL** | 18 | 18 | 21 | +0.6956 | 1.086 | 76 |
| frozen | 3 | 52 | 0.875 | yes | 0 | 0 | 18 | -0.0529 | 1.668 | 84 |


Summary: n = 0 (diagonal-only) already succeeds at 44, 48 and 52 ms. Failures: **frozen n = 2, tau = 52 ms** (N_unstable 18, N_margin 18, stopped at radius < 1e-3 after 21 iterations, worst catalogued Re +0.696) and **bridge n = 4, tau = 52 ms** (N_unstable 10, N_margin 10, 13 iterations, worst +0.802). Both fail while the smaller-radius sets (nested subspaces, n = 0 and n = 1 and n = 3 frozen / n = 0..3 bridge) succeed: the failure is a local-search / trust-region stall of the engine, not infeasibility. Not rerun, not retuned.
Effort at 52 ms: n = 0 1.512; n = 1 1.431 (frozen), 1.455 (bridge); n = 3 frozen 1.668, bridge 2.015: extra freedom does not reduce effort, the engine returns the first feasible point reached. At 44 ms effort falls slightly with n (0.249 -> 0.203-0.209).

## E4 (tau = 44 ms): rho in {0.875, 0.90, 0.925, 0.95}
rho = 0.875 rows are the E3 44 ms rows above. Further rows:

| graph | n | tau ms | rho | success | N_unst | N_margin | iters | worst Re | effort | runtime s |
|---|---|---|---|---|---|---|---|---|---|---|
| bridge | 0 | 44 | 0.9 | yes | 0 | 0 | 2 | -0.0529 | 0.244 | 45 |
| bridge | 1 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0539 | 0.203 | 72 |
| bridge | 2 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0507 | 0.199 | 80 |
| bridge | 3 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0517 | 0.198 | 74 |
| bridge | 4 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0505 | 0.197 | 69 |
| bridge | 0 | 44 | 0.925 | **FAIL** | 0 | 1 | 2 | -0.0529 | 0.240 | 44 |
| bridge | 1 | 44 | 0.925 | **FAIL** | 0 | 1 | 2 | -0.0592 | 0.198 | 71 |
| bridge | 2 | 44 | 0.925 | **FAIL** | 0 | 1 | 3 | -0.0509 | 0.193 | 71 |
| bridge | 3 | 44 | 0.925 | **FAIL** | 0 | 1 | 3 | -0.0519 | 0.192 | 82 |
| bridge | 4 | 44 | 0.925 | **FAIL** | 0 | 1 | 3 | -0.0507 | 0.192 | 53 |
| bridge | 0 | 44 | 0.95 | **FAIL** | 0 | 3 | 2 | -0.0530 | 0.235 | 43 |
| bridge | 1 | 44 | 0.95 | **FAIL** | 0 | 3 | 2 | -0.0570 | 0.192 | 56 |
| bridge | 2 | 44 | 0.95 | **FAIL** | 0 | 3 | 3 | -0.0519 | 0.188 | 72 |
| bridge | 3 | 44 | 0.95 | **FAIL** | 0 | 3 | 3 | -0.0524 | 0.186 | 59 |
| bridge | 4 | 44 | 0.95 | **FAIL** | 0 | 3 | 3 | -0.0518 | 0.186 | 69 |
| frozen | 0 | 44 | 0.9 | yes | 0 | 0 | 2 | -0.0529 | 0.244 | 55 |
| frozen | 1 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0546 | 0.204 | 53 |
| frozen | 2 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0526 | 0.202 | 92 |
| frozen | 3 | 44 | 0.9 | yes | 0 | 0 | 3 | -0.0527 | 0.202 | 79 |
| frozen | 0 | 44 | 0.925 | **FAIL** | 0 | 1 | 2 | -0.0529 | 0.240 | 76 |
| frozen | 1 | 44 | 0.925 | **FAIL** | 0 | 1 | 2 | -0.0527 | 0.198 | 71 |
| frozen | 2 | 44 | 0.925 | **FAIL** | 0 | 1 | 3 | -0.0530 | 0.196 | 58 |
| frozen | 3 | 44 | 0.925 | **FAIL** | 0 | 1 | 3 | -0.0530 | 0.196 | 68 |
| frozen | 0 | 44 | 0.95 | **FAIL** | 0 | 3 | 2 | -0.0530 | 0.235 | 75 |
| frozen | 1 | 44 | 0.95 | **FAIL** | 0 | 3 | 2 | -0.0555 | 0.193 | 78 |
| frozen | 2 | 44 | 0.95 | **FAIL** | 0 | 3 | 3 | -0.0530 | 0.190 | 59 |
| frozen | 3 | 44 | 0.95 | **FAIL** | 0 | 3 | 3 | -0.0530 | 0.190 | 53 |

Largest verified rho per n (every n, both graphs): **0.90**. Each point is the largest verified point in the declared sweep, not a certified optimum. rho = 0.925 and 0.95 fail for every n, with N_margin 1 and 3 respectively, although the catalogued worst root is already about -0.053 to -0.059: the remaining roots beyond the margin are not in the Pade-seeded catalogue (or sit at Re > -1.5 with Im < 0.3), which only the full counter exposes. The engine therefore cannot see the blocking mode; this is a limitation of the preregistered engine, not evidence that no gain matrix exists.

## Falsification statements
- "Information radius extends the repairable region" (SUPPORTED only if some tau or rho case succeeds for n >= 1 and fails for n = 0 under identical engine settings): **NOT SUPPORTED / refuted in this grid**. No case succeeds for n >= 1 and fails for n = 0: n = 0 succeeds at every tau (44, 48, 52) and at rho 0.875 and 0.90, and fails at rho 0.925 and 0.95 exactly as every n >= 1 does. The only discordant cases go the other way (n = 2 frozen and n = 4 bridge fail at 52 ms where n = 0 succeeds).
- Largest verified rho does not increase with n (0.90 for all n, both graphs). Replacement share unchanged by information radius in this sweep.
- The other two statements (modal authority; spillover) are from REPORT_E1_E2_E5.md and are not changed by E3/E4.
- Interpretation: closed-loop repair of this fault is achievable by diagonal retuning (n = 0) up to the sweep limits; off-diagonal PLL cross-coupling is not needed for the repair.

## All failures (complete list)
frozen n2 tau52 rho0.875; bridge n4 tau52 rho0.875; rho 0.925 and 0.95 at n = 0..3 (frozen) and n = 0..4 (bridge), i.e. 18 runs. Total: 54 unique runs, 34 succeed, 20 fail (2 at rho 0.875, 18 at rho 0.925/0.95).

## Caveats
Counter and catalogue are floating-point, uncertified. Each case is one deterministic engine run from the baseline; no restarts. The 52 ms failures show the engine result is path dependent, so "fails" means "this engine did not reach N_margin = 0 within its limits". Bridge variant is the AMENDMENT_01 deviation (extra edge 30-33). Julia nonlinear E6 not started.
