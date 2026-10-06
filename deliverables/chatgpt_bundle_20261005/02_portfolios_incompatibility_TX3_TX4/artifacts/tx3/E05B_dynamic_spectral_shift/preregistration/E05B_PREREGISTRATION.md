# E05B — Dynamic spectral-shift decomposition preregistration

**Freeze:** TX3-E05B-HOLDOUT/NUMERICS/GATES-1.0. **Git SHA before freeze:** `24dd4019678c896085500459c62928284dbd20a3`.

E05B tests a new claim and does not rescue E05. `C3a=REJECTED`, its damping threshold
`|Delta_S zeta| >= 0.0025`, and the accepted C1/C2 decisions are immutable. E05C, E06,
TDS surgery, and new ParaEMT science are outside scope.

For every E03 pair and triple at all 24 prior operating points, and subsequently at
16 newly generated independent Sobol points, the descriptor determinant is separated as

`log|det T| = log|det(-gy)| + log|det M| + log|det(sI-Ad)|`,

where `Ad=M^-1(fx-fy gy^-1 gx)`. The three finite Mobius components are evaluated at
identical coalition vertices. Finite poles are computed independently by dense eigensolution;
the frozen E03 full-descriptor vertices and a pointwise sparse-LU audit test the factorization.

Only the 24 prior points are development data. They select two pair and two triple dynamic
candidates by the frozen reproducibility/rank rule; A7-A8, A3-A6, and A2-A7-A8 remain mandatory
diagnostics. Candidate identities and OP00 contour families are hashed and committed before
the new holdout may be evaluated.

The connected dynamic determinant is represented by
`Xi_S(s)=sum_U (-1)^(|S|-|U|) tr[(sI-Ad(U))^-1]`. Circular contours are determined solely by
the local baseline pole gap. Exact residue moments and a 128-node numerical Xi integral must
agree. D1 is exact-factor closure, D2 is reproducible finite-pole externality, and D3 is
localized connected pole motion. These are not an engineering damping-materiality gate.
