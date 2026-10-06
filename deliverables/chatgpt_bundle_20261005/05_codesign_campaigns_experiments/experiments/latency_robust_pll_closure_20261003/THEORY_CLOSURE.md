# Analytical closure of the fixed-replacement latency problem

The replacement vector `rho` is frozen to the parent experiment's ten values.
Only the ten PLL proportional and ten integral gains vary. The delay is an
exogenous common latency on the PLL measurement/error channel. For a fixed
operating point, the quotient linearization has the exact retarded form

\[
\Delta(s;K,\tau)=sI-A_0(K)-\sum_{i=1}^{10}A_i(K)e^{-s\tau},
\qquad A_i=B_i(K)C_i^T.
\]

The full characteristic is evaluated with the exponential factors. No Padé
approximation is used for the spectral gate. At fixed `rho`, the model's gain
dependence satisfies

\[
A(K)-A(K_0)=[B(K)-B(K_0)]C^T,
\qquad \operatorname{rank}(A(K)-A(K_0))\le10.
\]

`TABLE_F11_GAIN_ACTION_RANK.csv` checks this identity along the evaluated
trajectory; the largest relative reconstruction error is of order `1e-14`.
The observed rank of an individual *difference* can be smaller than ten.

For a simple delayed root, let `v` and `w` be right and left null vectors of
`Delta`. Differentiating `Delta(lambda,p)v=0` gives

\[
\lambda_p=-\frac{w^H\Delta_p(\lambda,p)v}
                    {w^H\Delta_s(\lambda,p)v}.
\]

At a boundary `Re(lambda_m(K,tau_m))=-sigma_req`, provided the crossing is
transverse (`alpha_{m,tau}>0`), the implicit-function theorem gives

\[
\nabla_{\log K}\tau_m
=-\frac{\nabla_{\log K}\alpha_m}{\alpha_{m,\tau}}.
\]

The security latency margin is a lower envelope of all modal crossings:
`tau_crit(K)=min_m tau_m(K)`. Therefore it is nonsmooth at family switches. A
single-mode gradient cannot describe its ascent at a multimode balance.

The implemented local predictor maximizes `t` over log-gain change `du` with
one inequality per discovered mode,

\[
t\le \tau_m(K)+\nabla\tau_m(K)^Tdu,
\]

and a measured bus-16 `+100 MW` actuator-surrogate inequality,

\[
s_{\rm act}(K)+E_{\rm act}^Tdu\ge s_{\rm required}.
\]

The trust box and `L1` step budget restrict extrapolation. Positive
multipliers on several modal constraints and the actuator inequality are KKT
multipliers of this **local linear program**. They do not certify the KKT
conditions of the complete nonlinear multi-event problem. Every proposed
point requires a new equilibrium, full exponential-characteristic contour,
all relevant root refinements, and the five frozen zero-delay events.

At the five-family equalization point, the relevant physical explanation is
minimax competition among delayed PLL families under a nearly active SG
governor constraint. This is a numerical interpretation; global optimality,
robust safety, and nonlinear positive-delay safety require separate evidence.

For the largest fully event-validated design in this closure,
`full20_step15_medium`, the five corrected crossings span
`43.797191610–43.798385186 ms` (spread `0.001193576 ms`). The local LP
predictor assigned positive multipliers to all five modal constraints and
to the bus-16 actuator surrogate. Its own KKT residual was `4.44e−15`.
The complete nonlinear/full-characteristic optimization is still improving,
so the term *equalization* describes the observed mode competition and does
not identify a proved optimum. The earliest crossing is for the required
security line `Re(s)=−0.05 s⁻¹`, not necessarily the first pole with
`Re(s)=0`.
