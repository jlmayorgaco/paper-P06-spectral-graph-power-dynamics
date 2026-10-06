# ExpQ2B derivation and numerical method

## Frozen model and contract

The ExpN analytical source freeze is reused unchanged (SHA-256 model id
`e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a`). The
generator objective weights are the initialized component active powers
`Pgen0_i`; ZIP loads at buses 31 and 39 remain separate. For each bus,
`epsilon_i` is the retained SG fraction and `rho_i=1-epsilon_i`; the existing
component P/Q split, SG rating convention, and gain box are inherited from
ExpN. No PowerDynamics call is made by the design evaluator.

## Causal phase measurement

For the passive bus-voltage measurement, the reduced model gives the algebraic
voltage jump `Delta theta_k(0+)` and the post-event bus-frequency trajectory
`f_k(t)=C_k exp(A t) B_l d + D_k d` (with the corresponding step integral for
the state). The trajectory is integrated once to obtain unwrapped phase
`theta_k(t)`. For a frozen window `T`,

\[
f_{T,k}(t)=\frac{\theta_k(t)-\theta_k(t-T)}{2\pi T},\qquad
R_{T,k}(t)=\frac{\theta_k(t)-2\theta_k(t-T)+\theta_k(t-2T)}{2\pi T^2},
\]

and `theta(t)=0` for pre-event time in the phase-deviation coordinates. The
implementation uses the same phase channels in every support, includes the
algebraic jump, and applies no feedback. The sweep uses `T={0.2,0.5,1,2}s`;
`T=0.5s` is the nominal event-study window. A uniform time lattice is used for
candidate extrema with `dt=0.05s` during SQP and `dt=0.01s` for the final
re-evaluation. This is a numerical peak estimate; grid convergence and
one-sided window-edge checks remain part of the validation record.

## DAE reduction and steady output

For a fixed physical architecture, ExpN's exact descriptor is reduced through
the regular algebraic block, and the validated rotational vector is quotiented.
The step output is evaluated from the quotient model. The DC output is

\[
H_0=-C A^{-1}B+D,
\]

so `F_inf` is the largest absolute generator-bus DC entry times the event
magnitude. PLL gains have no direct DC contribution under the inherited
SimpleGFLDC trim; the optimizer sets this derivative to zero and finite-difference
checks the frequency-domain result near the candidate.

## Spectrum and robustness

The complete finite spectrum is recalculated from the fixed-architecture
quotient for every accepted design step. The required nominal condition is
`alpha<=-0.05 s^-1`. The ExpG robustness definition is the full-complex additive
state-matrix block at the shifted boundary. During local SQP, the corrector
tracks the zero-frequency singular value and singular values at the three
least-damped pole frequencies to identify the active robust branch. This is a
search estimate, not a full-band certificate. Any final beta claim requires
the ExpG full resolvent refinement and a separate lower-bound audit.

The robust derivative uses the envelope identity for a simple minimum
singular value at its active frequency,
`d beta/dz = Re(uᴴ (-A_z) v)`, with `A_z` formed by centered descriptor
perturbations in a fixed gauge basis. Its singular gap and conditioning are
stored with the derivative audit. The active modal derivative uses the ExpN
left/right eigenvector sensitivity inside a fixed support. The DC derivative
uses descriptor finite differences only in retained-SG coordinates; gain
derivatives are exactly zero. If the peak frequency is the DC plateau, its
gradient equals the DC gradient; otherwise a complete finite-window
re-evaluation is used.

## Active-set SQP and certification scope

Variables are dimensionless: `epsilon` in `[0,1]` and Kp/Ki scaled to the
frozen ExpK boxes. The primary objective is linear in retained SG MW. The
feasible-start active-set QP uses a diagonal positive Hessian, bound and trust
constraints, a deterministic exact-model line search, and complete constraint
re-evaluation for each candidate. The first executed branch is the ten-bus
mixed architecture seeded adjacent to all-SG. Disturbance continuation uses
step halving when the current design is not feasible at the proposed event size.

This branch-local search is not an exhaustive support search. A small
stationarity residual is not assumed from a feasible point; LICQ and SOSC are
reported only when their explicit audits pass. Full robust globality and
support completeness are not inferred from the sampled corrector or from the
frequency-only values.
