# Q5–Q6 — replacement and joint action-space structure

## Terminal-row update for one replacement fraction

For a fixed interior SG/GFL architecture, keep the local state ordering fixed. The repository closure is a descriptor system with device differential rows and network algebraic KCL rows. At candidate bus `i`, the KCL contribution is affine in `epsilon_i=1-rho_i`: the retained SG Norton/output map has weight `epsilon_i`, and the GFL port-current map has weight `1-epsilon_i`. All terms changed by this scalar enter the two real KCL equations at that terminal. Therefore the augmented descriptor pencil has an affine update of rank at most two:

`T(s,epsilon_i) = T_-i(s) + epsilon_i U_i D_i(s) V_i^H`, with `rank(U_i D_i V_i^H) <= 2`.

When `T_-i(s)` is nonsingular, define `N_i(s)=D_i(s)V_i^H T_-i(s)^(-1)U_i`. The determinant lemma gives

`det T = det T_-i det(I_2 + epsilon_i N_i)`

`= det T_-i [1 + epsilon_i tr(N_i) + epsilon_i^2 det(N_i)]`.

Thus conditional candidates solve `1 + epsilon_i tr(N_i) + epsilon_i^2 det(N_i)=0`; when `det(N_i) != 0`,

`epsilon_i = [-tr(N_i) ± sqrt(tr(N_i)^2 - 4 det(N_i))]/[2 det(N_i)]`, and `rho_i=1-epsilon_i`.

A candidate is physically admissible only when `epsilon_i` is real and lies in `[0,1]`. It is a conditional root at the selected `s`, with other device actions frozen—not a closed-loop eigenvalue after equilibrium re-solution. The coefficients and physical roots must be evaluated from the actual repository descriptor at each operating point. The Kron-reduced state matrix is generally rational in `epsilon_i`, because the network closure is inverted. The theorem also excludes architecture endpoints where state dimension changes, singular eliminated blocks, and a changed equilibrium/scaling contract.

The structural audit found a rank-two row update at each candidate terminal, zero row-support residual, maximum affine second-difference residual `1.27e-21`, and maximum augmented-Schur discrepancy `2.16e-18`. Across 120 candidates (10 buses, 3 target frequencies, 2 delays, 2 quadratic roots), zero roots were real and inside the physical interval `[0,1]`. The exact quadratic law therefore did not yield a usable SG-retention boundary in this tested grid.

## Joint gain and replacement updates

For ten candidate terminals, each replacement action updates at most two descriptor rows, giving rank at most 20. The PLL gain update at each GFL changes the proportional and integral state rows through one shared phase-detector row; jointly changing both gains at one device is one rank-one state-row action, giving rank at most 10 over ten devices. In this fixed-coordinate augmented descriptor representation, the joint action-space correction has dimension bounded by

`r <= 20 + 10 = 30`.

Fixed delays change the scalar entries of the PLL update through `exp(-s tau_i)`; they do not add action directions. This is a rank bound, not an optimization result or root-count certificate. The numerical audit reconstructed the full descriptor difference with maximum residual `9.07e-15` over six complex frequency/delay cases; maximum numerical difference rank was 30. This pass is restricted to the fixed interior architecture, fixed equilibrium/port scaling, and the tested descriptor coordinates.

**Claim class:** Q5 is a reduced-descriptor theorem under the stated fixed-coordinate port-scaling contract plus numerical validation. Q6 is an exact finite-dimensional determinant factorization when the reference pencil is nonsingular, with numerical reconstruction validation. Neither establishes a replacement frontier.
