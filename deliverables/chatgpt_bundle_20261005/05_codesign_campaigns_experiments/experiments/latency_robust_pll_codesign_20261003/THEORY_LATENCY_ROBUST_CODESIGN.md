# Exact gain action, delay-margin gradient, and multimode retuning

## Model and scope

At the frozen interior replacement vector \(\rho_H\), each of ten GFL PLLs receives the *exogenous* measurement delay \(\tau_i\) on its phase-detector error. Its gains are \(K_{p,i},K_{i,i}\); no other converter path is delayed. The reduced physical DDE characteristic is

\[
\Delta(s;K,\tau)=sI-A_0-\sum_{i=1}^{10}b_i(K)c_i^T e^{-s\tau_i},
\qquad b_i(K)=K_{p,i}b_{p,i}+K_{i,i}b_{I,i}.
\]

The descriptor implementation retains 282 variables (203 physical differential modes after gauge deflation). For a *fixed* replacement vector, operating point, topology, and PLL channel support, \(A_0\), \(c_i\), \(b_{p,i}\), and \(b_{I,i}\) do not change with K. The gain difference from any reference K₀ is therefore **exactly**

\[
\boxed{T(s;K,\tau)-T(s;K_0,\tau)=U_K(K,K_0)H(s,\tau),\quad \operatorname{rank}\le10,}
\]

where column i of \(U_K\) inserts \(\Delta K_{p,i}b_{p,i}+\Delta K_{i,i}b_{I,i}\) in the PLL state rows, and row i of \(H\) is \(-e^{-s\tau_i}c_i^T\) (including descriptor voltage coordinates before algebraic elimination). Thus **20 free gain variables act through ten measured PLL error channels**. The joint \((\rho,K_p,K_i)\) descriptor update from the earlier campaign has rank at most 30; the fixed-ρ gain-only bound is sharper. At nonsingular reference points,

\[
\det T(K)=\det T(K_0)\det[I_{10}+H T(K_0)^{-1}U_K].
\]

The algebraic rank statement is an exact identity for the declared model. `L1_ACTION_SPACE_VALIDATION.csv` checks 48 descriptor cases (both stored gain sets and 20 prespecified gain-box random vectors), including delays near 37–40 ms. Random vectors were **not** screened by nonlinear events. An exact factorization does not by itself make its floating-point root count interval-certified.

## Local latency-margin derivative

For a simple root \(\Delta(\lambda)v=0\), \(w^H\Delta(\lambda)=0\), with nonzero \(w^H\Delta_s v\),

\[
\boxed{\lambda_p=-\frac{w^H\Delta_p v}{w^H\Delta_s v}.}
\]

At fixed ρ,

\[
\Delta_{K_{p,i}}=-e^{-s\tau_i}b_{p,i}c_i^T,\qquad
\Delta_{K_{i,i}}=-e^{-s\tau_i}b_{I,i}c_i^T,
\]

\[
\Delta_{\tau_i}=s e^{-s\tau_i}b_i c_i^T,\qquad
\Delta_s=I+\sum_i\tau_i e^{-s\tau_i}b_i c_i^T.
\]

For a common delay \(\tau_i=\tau\), set \(\Delta_\tau=\sum_i\Delta_{\tau_i}\). If one root alone limits the margin and \(\alpha_\tau=\operatorname{Re}\lambda_\tau>0\), the implicit-function theorem gives

\[
\boxed{\frac{\partial\tau_{c}}{\partial K_j}
=-\frac{\operatorname{Re}\lambda_{K_j}}{\operatorname{Re}\lambda_\tau}.}
\]

The derivative is local. It fails as a single smooth derivative of the *minimum* delay margin when two families become co-critical. The `L2` tables compare these formulas with independently relinearized finite perturbations and explicitly label mode switches.

## Analytical step and exact correction

In log coordinates \(u_j=\log K_j\), a single-mode Euclidean trust ball \(\|\Delta u\|_2\le\kappa\) has first-order support direction \(\Delta u^*=\kappa g/\|g\|_2\), where \(g=\nabla_u\tau_c\). For co-critical roots \(m\), the actual local objective is the lower envelope \(\min_m\tau_{c,m}(u)\). The implemented small LP maximizes \(t\) subject to

\[
t\le\tau_{c,m}(u)+g_m^T\Delta u\quad(\text{all discovered near-active modes}),
\qquad |\Delta u_j|\le0.01,\quad\|\Delta u\|_1\le0.05,
\]

and the source gain bounds. The trust sizes were frozen in `L3_OPTIMIZATION_CONTRACT.md`. This is an analytical **predictor**, not an optimizer of the full DDE/events. Every proposal is relinearized; the exact exponential characteristic is solved; a complete numerical DDE contour checks SAFE below and UNSAFE above the candidate boundary. Zero-delay equilibrium, full physical spectrum, and all five frozen nonlinear events determine acceptance. If the active SG-actuator margin fails, the step is rejected and backtracked. Event sensitivities inferred from a secant along a tested direction are only local approximations; the executed event remains the feasibility oracle.

Neither an interval certificate of every DDE root nor a global optimization bound is available. Report \(K^*,\tau_c^*\) only as *best fully validated found* when all required checks pass. No positive-delay nonlinear event claim follows from this spectral experiment.
