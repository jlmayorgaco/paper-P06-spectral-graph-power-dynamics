# Q6 — unified finite-dimensional action-space descriptor

For the fixed interior SG/GFL architecture and fixed coordinate ordering, write the augmented descriptor characteristic pencil as `T_ref(s,tau)` plus the actions in the ten PLL channels and ten SG/GFL terminal ports. Each device's joint `(Kp, Ki)` change is one state-row direction because proportional and integral terms share the same scalar PLL phase-detector row; across ten devices this contributes rank at most 10. A change in retained synchronous share at one terminal changes only its two real KCL rows, so all ten replacement actions contribute rank at most 20. Delays enter those same PLL update directions through `exp(-s tau_i)` and add no new action coordinates.

Consequently there are matrices `U`, `V`, and an action coefficient matrix `Theta(y,s,tau)` such that

`T(s;y,tau) = T_ref(s,tau) + U Theta(y,s,tau) V^H`, with `r = number_of_columns(U) <= 30`.

Whenever `T_ref` is nonsingular, the determinant lemma gives

`det T = det T_ref det(I_r + Theta V^H T_ref^{-1} U)`.

This gives the exact compact action-space determinant conditional on the fixed descriptor and nonsingularity of the reference pencil. The tested maximum numerical action rank was 30, and the factorization reconstruction residual was at most `9.07e-15` over six complex frequency/delay samples (`TABLE_Q06_MASTER_ACTION_SPACE.csv`). It is not an optimization result, a global root oracle, or an equivalence at singular eliminated/reference blocks. State-dimension-changing endpoints and re-solved equilibria require rebuilding the descriptor and are outside this rank statement.
