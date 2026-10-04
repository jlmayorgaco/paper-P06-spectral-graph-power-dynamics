# Mathematical review: network moments and control authority

Independent read-only review of `moments.py`, `export_inputs.jl`, the protocol,
generated moment tables, and the multimode corrector. Only this review file was
added. No broad simulations were run.

**Finding.** The transfer equations, Laurent recurrence, difference formula and
moment signs are correct for the stated symmetry-preserving model. There are
useful exact restrictions on spatial moment shaping, and a larger space of
moment-preserving actions than sitewise compensation. Modal feasibility remains
a separate requirement. The present computations validate a symmetry-restored
local transfer model; they are not an interval certificate for independent exact
binary64 matrix entries.

## 1. Input/output equations and assumptions

With all PLL triples removed, let the remaining state be `h`, PLL angles `q`,
and the physical disturbance `u`. The exported equations are

\[
 (sI-A_{hh})h=A_{h\theta}q+B_{u,h}u,\quad
 e=C_hh+C_\theta q+e_u u,\quad
 \phi=O_hh+O_\theta q+\phi_u u.
\]

Thus

\[
\begin{aligned}
 G&=C_\theta+C_h(sI-A_{hh})^{-1}A_{h\theta},&
 g&=e_u+C_h(sI-A_{hh})^{-1}B_{u,h},\\
 T&=O_\theta+O_h(sI-A_{hh})^{-1}A_{h\theta},&
 Z&=\phi_u+O_h(sI-A_{hh})^{-1}B_{u,h}.
\end{aligned}
\]

These match `moments.py`. In the exporter, `A0=Fx-B*C` and
`Binput=F_u-B*detector_input` remove the instantaneous PLL contribution before
reintroducing its delayed version. The algebraic phase feedthrough is retained
in `Z`. The removed PLL rows of exported `Binput` are at most 3.56e-15 in
absolute value, consistent with their structural zero. The recorded physical
input finite-difference errors are 1.35e-9--2.16e-9. These are numerical checks,
not interval bounds.

For diagonal controllers, elimination gives exactly

\[
 K(s)q=g(s)u,\qquad
 K(s)=\operatorname{diag}\!\left[
 \frac{s^2(1+t_fs)e^{s\tau}}{sK_p+K_i}\right]-G(s).
\]

The input cancellation requires disturbances in the plant/detector path and
zero perturbation history. It need not hold for inputs injected into PLL states
or after the delayed measurement. Assume positive integral gains; fixed plant,
replacement, equilibrium, input/output maps and filters; analyticity of the
opened plant near zero; and one simple rotational singularity:

\[
 K_0r=0,\quad \ell^*K_0=0,\quad
 \operatorname{rank}K_0=n-1,\quad d=\ell^*K_1r\ne0.
\]

Then `P=r ell*/d` is the inverse residue: `K(s)^-1=P/s+O(1)`.
For calibrated bus phases, rotational equivariance gives `r=1` and `T(0)r=1`.
The implemented bordered recurrence additionally uses `ell*r != 0`; that
condition is satisfied numerically here but is not required by the general
Laurent theorem. A different border is needed if that pairing vanishes.

## 2. Exact finite moment laws

Write the windowed frequency transfer as

\[
 H(s)=\frac{1-e^{-sW}}{2\pi W}\{T(s)K(s)^{-1}g(s)+Z(s)\}
      =H_0+sH_1+s^2H_2+\cdots.
\]

The matrices `G,T,g,Z` are independent of PLL gains and delays. Since the
controller contribution starts at degree two, `K0` and `K1` are unchanged by
any such tuning, including finite changes of integral gain. Consequently

\[
 \boxed{H_0=T(0)Pg(0)/(2\pi)}
\]

is invariant under all these tunings. This is an algebraic DC statement;
interpreting it as a final response requires stability after cancellation of
the rotational pole.

Define real network weights `w_i=conj(ell_i) r_i/d`. For finite integral-gain
changes, allowing simultaneous proportional-gain/delay changes,

\[
 \kappa_2=\sum_i w_i\left(K_{i,i}^{\prime-1}-K_{i,i}^{-1}\right),\qquad
 \boxed{\Delta H_1=-\kappa_2H_0.}
\]

For fixed integral gains,

\[
 \kappa_3=\sum_iw_i\left(\frac{\Delta\tau_i}{K_{i,i}}
                      -\frac{\Delta K_{p,i}}{K_{i,i}^2}\right),\qquad
 \boxed{\Delta H_0=\Delta H_1=0,\quad\Delta H_2=-\kappa_3H_0.}
\]

**Proof.** The controller coefficients are
`s^2/Ki+s^3[(tf+tau)/Ki-Kp/Ki^2]+O(s^4)`.
Use the exact inverse identity
`Knew^-1-Kold^-1=-Knew^-1(Knew-Kold)Kold^-1` and
`PDP=(ell*D r/d)P`. Multiplication by the window factor, whose leading term is
`s/(2pi)`, gives the displayed orders and signs. The recurrence in
`coefficients()` obtains the same coefficients through successive solvability
conditions. Its window coefficients `1,-W/2,W^2/6` are correct. The stable
evaluation in `difference()` also has the correct sign and finite gain/delay
difference formula.

For a step `u0`, let `f_inf=H0 u0`. Under sufficient decay,

\[
 \int_0^\infty\Delta f\,dt=-\kappa_2 f_\infty.
\]

At fixed integral gains,

\[
 \int_0^\infty\Delta f\,dt=0,\qquad
 \int_0^\infty t\,\Delta f\,dt=\kappa_3 f_\infty.
\]

Here `H2` is the coefficient, not the second derivative: `H''(0)=2H2`.
No conclusion about nadir, settling, RoCoF, or stability follows from these
integral identities alone.

## 3. Spatial and bounded-gain consequences

For each nonzero final response, define the normalized signed area

\[
 A_{j,u}=\int_0^\infty(1-f_j(t)/f_{\infty,j})\,dt
        =-(H_1u)_j/(H_0u)_j.
\]

Then `Delta A_{j,u}=kappa2`, independently of bus and input. Thus every
difference `A_{j,u}-A_{k,v}` is invariant under all PLL gain/delay changes at
fixed plant, whenever the normalizations exist. Equivalently, for bus outputs
and any spatial contrast `Pi*1=0`, `Pi*H1` is controller-invariant. At fixed
integral gains, `Pi*H2` is also invariant under proportional-gain/delay tuning.
These are restrictions on low-frequency spatial shaping, not on complete
response waveforms. If a DC response is zero, use the unnormalized matrix laws.

There is an exact bounded-integral-gain test for matching specified area
targets. Put `x_i=1/Ki_i` and
`x_i in [1/Ki_max,i,1/Ki_min,i]=[xlo_i,xhi_i]`.
Relative to reference `x0`, every requested area change must be the same scalar
`a`. That scalar must lie in

\[
 \left[\sum_i\min\{w_i(xlo_i-x0_i),w_i(xhi_i-x0_i)\},\quad
       \sum_i\max\{w_i(xlo_i-x0_i),w_i(xhi_i-x0_i)\}\right].
\]

These conditions are necessary and sufficient for matching the algebraic
moment targets by bounded integral gains. They are only necessary for a stable,
otherwise feasible controller. The observed signed weights preclude a general
monotonic interpretation of increasing integral gain.

At fixed integral gains, finite delay compensation needs only

\[
 c^T\Delta K_p=b,\quad c_i=w_i/K_{i,i}^2,\quad
 b=\sum_iw_i\Delta\tau_i/K_{i,i}.
\]

Sitewise `Delta Kp_i=Ki_i Delta tau_i` is one solution. Cross-site solutions are
equally exact for these moments. For gain-change bounds `[l,u]`, existence is
equivalent to
`b in [sum min(c_i l_i,c_i u_i),sum max(c_i l_i,c_i u_i)]`.

## 4. Moment and modal compatibility

Let `R>0` specify gain effort, and let `a=grad_Kp Re(lambda)` for a simple root.
For the exact moment equality `c^T k=b`, the minimum-effort point and remaining
inverse metric are

\[
 k_M=\frac{R^{-1}cb}{c^TR^{-1}c},\qquad
 Z_R=R^{-1}-\frac{R^{-1}cc^TR^{-1}}{c^TR^{-1}c}.
\]

The remaining first-order modal authority is
`chi=a^T Z_R a`. For a linearized modal requirement `a^T k<=h`, including the
known delay contribution in `h`, let `nu=a^T k_M-h`.
If `nu<=0`, `k_M` already passes that linearized condition. If `nu>0` and
`chi>0`, the minimum additional correction is

\[
 z_*=-\nu Z_Ra/\chi,\qquad
 \tfrac12z_*^TRz_*=\nu^2/(2\chi).
\]

If `chi=0` and `nu>0`, no first-order correction preserving the moment can
meet that modal requirement. Bounds can invalidate the closed-form solution;
use the equality-constrained bounded QP. With several modal gradient rows `J`,
the authority matrix is `J Z_R J^T`; all active branches must be included.
This is standard constrained sensitivity geometry applied to a physically
derived exact moment constraint. Zero authority at one point is not a finite
stabilization impossibility theorem.

A stronger conditional finite exclusion is available. On a convex region of
the moment hyperplane containing `k_M`, where a root branch is regular, suppose
`grad(alpha)^T Z_R grad(alpha)<=L^2` uniformly and
`||k-k_M||_R<=r`. Integration along the segment gives
`|alpha(k)-alpha(k_M)|<=L*r`.
If `alpha(k_M)-L*r` exceeds the required boundary, every design in that region
fails it. This requires a verified uniform derivative bound. A fixed-contour
cluster mean can replace the simple root across internal collisions, with the
usual necessary-only interpretation of a mean margin.

## 5. A valid graph communication corollary

For a new controller architecture with matrix gains and source measurement
delays before the gains,

\[
 s^2(I+sT_f)\theta=(sK_p+K_i)E_\tau e,
\]

the correct order is

\[
 K=E_\tau^{-1}(sK_p+K_i)^{-1}s^2(I+sT_f)-G.
\]

Its controller coefficients are
`K2=Ki^-1` and
`K3=Ttau Ki^-1+Ki^-1 Tf-Ki^-1 Kp Ki^-1`.
Therefore, at fixed invertible `Ki` and fixed delays,

\[
 \boxed{\Delta K_pK_i^{-1}r=0}
\]

is sufficient to preserve `H0,H1,H2` exactly. If `Ki=ki I`, any zero-row-sum
proportional coupling works. In particular, `Delta Kp=p(L)`, `p(0)=0`, `Lr=0`
works. For heterogeneous integral gains, `Delta Kp=p(L)Ki` works; unscaled
`p(L)` generally does not. A graph polynomial is only one way to satisfy this
annihilation condition.

For a simple finite-frequency root of
`F=s^2(I+sTf)-(sKp+Ki)E_tau G`, a direction `Delta Kp` has sensitivity

\[
 \lambda_\eta=
 \frac{\lambda\,v^*\Delta K_pE_\tau(\lambda)G(\lambda)q}
      {v^*F_s(\lambda)q}.
\]

The numerator can vanish or have either sign. The corollary supplies a
moment-preserving design space, not a guarantee of modal improvement. It also
requires new measurement links, gain bounds and communication-delay semantics;
it is not an available action in the existing diagonal-PLL experiment.
Row-sum cancellation is familiar disagreement-control structure. The useful
specialization is its exact full-network moment consequence with the stated
PLL dynamics, not a new graph-polynomial identity.

## 6. Evidence and claim limits

The gauge corrections are explicit: about 6.68e-12 for `G0` and 1.29e-11 for
`T0`. This is appropriately disclosed. A tiny independent exact-binary gauge
residual would change the literal limit at zero, so the restoration cannot be
silently relabeled an exact result for the unmodified binary export.
The moment identity tables and small-positive-frequency checks are consistent
with the theory, but the source assumptions and rounding are not yet enclosed.

Changing SG/GFL replacement changes the plant, DC gain and weights. Recompute
them at the new replacement before applying these laws. A new design may match
all moment equalities and still destabilize another mode: the preserved failed
single-mode attempt is directly relevant evidence. The multimode QP and its
restoration slack are predictors; endpoint full-spectrum checks remain the
independent feasibility gate. No support maximum, arbitrary-disturbance safety,
or controller superiority follows from this review.

## Final targeted audit of the theorem section and interval assumptions

`THEORY_MOMENT_SECTION.tex` correctly distinguishes finite parameter changes
from finite-frequency accuracy and finite-amplitude nonlinear behavior. Its
Taylor coefficients, signed-area signs, rearrangement bounds, projected modal
correction and matrix ordering are correct. The time-domain statements should
be read as zero-state responses from the same locked pre-event history. A
graph direction can have zero modal sensitivity as well as either sensitivity
sign; its existence alone does not establish useful finite-frequency authority.

The section's sufficient graph direction is sound:

\[
 \Delta K_P=K_I S K_I,\quad Sr=0
 \quad\Longrightarrow\quad
 \Delta K_3^{\rm ctrl}=-K_I^{-1}\Delta K_PK_I^{-1}=-S,
 \quad \Delta K_3^{\rm ctrl}r=0.
\]

This is equivalent to the earlier sufficient right-annihilation condition,
with an invertible change of parameterization. No commutation of `Ki`, the
graph matrix or heterogeneous delays is needed. The delay operator must remain
on the side specified in the controller equation.

The interval script's certificate logic is valid for its stated restored
point model. Exact right multiplication by `I-11^T/10` supplies the structural
nullvector; a nonzero 9-by-9 minor then proves rank exactly nine. The first nine
columns are independent, and their negative sum is the last column. Therefore
solving orthogonality to those columns plus `ell^T 1=1` constructs the full
left nullvector. Enclosing a nonzero `ell^T K1 1` proves a simple analytic zero.
Successful finite ball inversion of the hidden matrix proves its regularity.
The script hash matches the certificate. An independent, read-only Arb192
recalculation reproduced the positive minor, approximately 0.2498901337062252,
and `d`, approximately 0.09286351224410158.

This certifies hidden regularity, the restored gauge rank, the residue
denominator and weight signs. It does not certify physical model uncertainty,
the unmodified export's exact gauge, stability, nonlinear trajectories, or a
finite-frequency remainder. Output equivariance is imposed/documented
separately by the transfer-map restoration. The JSON field
`max_float_weight_error` is a floating comparison, not itself an outward error
bound; independent interval subtraction bounded the maximum discrepancy by
1.0584197516831e-8.

Using the same independently enclosed weights and the actual saved binary64
gains/delays in `globalized/designs`, the moment residuals were enclosed as
follows; the listed intervals are conservatively rounded outward:

| Additional uniform delay | Certified `kappa3` interval, s^2 |
|---|---:|
| 0.1 ms | [-1.245952383e-17, -1.245952382e-17] |
| 0.5 ms | [-2.773374275e-16, -2.773374273e-16] |
| 1.0 ms | [-5.204016787e-16, -5.204016785e-16] |

Thus the computed designs are extremely close to moment-neutral in the
restored model; literal equality to zero should not be claimed for rounded
gains. The coefficient identity bounds their remaining `H2` difference by
`|kappa3|` times `|H0|`, entrywise, without suggesting a frequency-band bound.

The recorded 41 ms candidate has catalog critical real part
-0.0585728980256243/s and status `MAX_ITERATIONS`. It passes the original
-0.05/s catalog guard but misses the optimization target -0.0600001/s.
It is neither a converged solution nor an optimum. Full-spectrum and nonlinear
tests are separate gates; their running status was not converted into a result
in this review.
