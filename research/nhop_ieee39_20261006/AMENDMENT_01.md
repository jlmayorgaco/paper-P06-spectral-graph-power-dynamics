# AMENDMENT_01 (2026-10-06, written before any E1/E2/E5 run)
Finding (T1): the frozen G_c rule (2 nearest electrical neighbours, symmetric union) yields a DISCONNECTED graph with 2 components
{30,31,32,37,38,39} and {33,34,35,36} (sites index 0..9 = bus-30). Hop distance across components is infinite; "n = 0..diam(G_c)" is
therefore undefined as written.
Resolution (no change to the frozen rule or primary results):
1. PRIMARY analysis uses the frozen G_c as is. N_i^(n) = {j : finite dist(i,j) <= n}; n runs 0..D_f where D_f = largest FINITE hop distance (= 3).
   Cross-component information is never available in the primary analysis; eps(n) may therefore never reach 1e-8 and n_min may be undefined. That is a reported result.
2. SENSITIVITY variant "Gc+bridge": add the single minimum-electrical-distance cross-component edge (r_ij minimal over i,j in different components) so G_c is connected; n = 0..diam of the bridged graph. Labelled as a deviation in every table (column graph = frozen | bridge).
Block convention of Y.csv: [[a,-b],[b,a]] -> a+jb verified (violation 1.4e-14); r_ij = |Z_ii+Z_jj-2Z_ij| is invariant under complex conjugation, so the graph does not depend on the sign convention.
