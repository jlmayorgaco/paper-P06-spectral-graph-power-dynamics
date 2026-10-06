# TX3 theory freeze for E05B dynamic spectral shift

**Freeze ID:** `TX3-TF-E05B-1.0`

For the index-one finite DAE operator
`T(s)=[[sM-fx,-fy],[-gx,-gy]]`, nonsingular `gy` and strictly positive diagonal
`M` give the exact Schur identity

`det T(s)=det(-gy) det(M) det(sI-Ad)`,

with `Ad=M^-1(fx-fy gy^-1 gx)`. Therefore every finite Mobius externality separates
linearly into algebraic, mass-scaling, and finite-pole terms. Algebraic and mass factors
are independent of `s`; they vanish from the logarithmic derivative.

For coalition S, define the meromorphic connected logarithmic derivative directly as
`Xi_S(s)=sum_{U subset S} (-1)^(|S|-|U|) tr[(sI-Ad(U))^-1]`. No branch of a complex
logarithm is needed. In this finite-dimensional, generally nonnormal setting the identity
follows from Jacobi's determinant formula. Perturbation-determinant/spectral-shift literature
is conceptual background, not a claim that self-adjoint operator theorems apply unchanged.

For a contour avoiding all poles,
`mu_{S,k}=(2 pi i)^-1 integral s^k Xi_S(s) ds` is the Mobius sum of enclosed pole moments.
`mu_0` is connected pole count and `mu_1` is connected enclosed pole sum. A single isolated
positive-imaginary pole is used; its conjugate carries the conjugate moment. Stable one-pole
occupancy is required at every coalition vertex. Contour results establish localized pole
motion, not material damping-margin consequence.
