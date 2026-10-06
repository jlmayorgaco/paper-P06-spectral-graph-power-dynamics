"""Assemble the open manuscript from preserved derivations and result files."""
from pathlib import Path
import json,re,hashlib,math
import pandas as pd
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
from paths import resolve
TARGET=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
old=(OUT/'THEORY_before_closure.tex').read_text(encoding='utf-8')
def cut(a,b): return old[old.index(a):old.index(b)]
pre=old[:old.index(r'\title{')]
pre+=r"""
\title{Beyond Nodal Damping:\\Collective Interaction Bounds for Delay-Aware SG--GFL Co-Design}
\author{Research manuscript}
\date{3 October 2026}
\begin{document}
\maketitle
\begin{abstract}
Grid-following penetration limits depend on collective device--network
interactions that bounded local PLL retuning may not remove.
We formulate a regional bound problem for synchronous-to-grid-following
replacement with fixed measurement delays and freely moving modal patterns.
Retaining network algebraic variables yields an affine descriptor
perturbation for a declared fixed-equilibrium sharing contract.
A contour identity converts this perturbation into finite changes of
spectral-cluster means. We derive their full parameter Hessian,
an explicit closed-walk remainder, and an error budget for omitting
inter-site interactions. Necessary cluster conditions give linear and
semidefinite outer relaxations of the secure replacement problem;
their bounds require verified contour and numerical errors.
Three frozen IEEE--39 algebraic tests give walk remainders
18--46 times smaller than a scalar-norm bound, without establishing
replacement superiority. Complex ball arithmetic certifies one root
cluster throughout a declared IEEE--39 parameter box.
A separate preserved nonlinear experiment
validates one analytical PLL retuning step at 40 ms uniform delay:
87.6\% GFL, with five declared events passing. Unchanged gains also
pass its spectral test. The theory therefore supplies a conditional
finite-bound method and a validated feasible construction, not a
demonstrated physical replacement maximum. We identify the missing
informative-bound and comparison tests required for that stronger claim.
\end{abstract}
\noindent\textbf{Index terms:} grid-following converters, PLL delay,
descriptor systems, nonlinear eigenvalues, collective interactions,
regional replacement bounds.

\section{Introduction and research question}
Replacing synchronous generation changes inertia, electrical coupling
and converter synchronization simultaneously. Retuning each PLL against
one operating point can preserve an assigned oscillatory mode while
moving other modes toward a security boundary. The design question is
therefore quantitative: within a declared replacement and gain region,
how much replacement can any admissible local PLL tuning support,
at fixed physical delays and under a stated security contract?

Detailed low-inertia models already connect converter control settings
to collective modal behavior \cite{markovic}. Grid structure and PLL
synchronization are also established concerns \cite{huang}.
Stability-constrained planning incorporates frequency and grid-strength
constraints \cite{chuteng}; recent expansion-planning work explicitly
addresses converter-driven modes \cite{jia}. Thus neither PLL tuning,
delay-induced instability, nor the use of network information is by
itself a new contribution.

The methodological distinction examined here is a \emph{necessary
finite regional bound across all allowed local gains}, rather than
only a successful tuning or a fixed-pattern assignment obstruction.
Descriptor lifting \cite{peaucelle}, nonlinear-eigenvalue contour
localization \cite{bindel}, contour moments \cite{beyn}, structured
DDE robustness \cite{borgioli}, and fixed-structure delayed-controller
optimization \cite{michielsdesign} are antecedents, not discoveries
of this work. The proposed synthesis specializes them to the actual
SG/GFL sharing model and makes inter-site curvature and remainder
terms explicit.

The manuscript establishes three conditional mathematical results:
an affine full-model intervention representation, a cluster
curvature/closed-walk bound, and a controller-bounded necessary
replacement relaxation. Its numerical evidence has two distinct
scopes: identity and enclosure tests of the bound machinery, and one
previously frozen full-model feasible retuning step. An informative
IEEE--39 replacement ceiling below the parameter-box ceiling has not
been demonstrated. This limitation prevents a claim of superior
replacement capacity or certified near-optimality.

\paragraph{Organization.}
The main text states the physical contract, derives the finite
regional bounds, and reports their evidence and limitations.
Appendices preserve the physical damping derivation, frequency-moment
law, nonlinear RoCoF argument, exact gain map, and local design
constructions. These auxiliary results have their own assumptions;
they are not combined into an unstated nonlinear stability theorem.
"""
model=cut(r'\section{Scientific contract and actual model}',r'\section{First correction')
model=model[model.index('With fixed device support'):]
model=r'\section{Physical model and optimization contract}\label{sec:model}'+'\n'+model
model+=r"""
For a compact support-preserving region $\mathcal D$, define
\[
 P^\star_{\mathcal D}(\tau)=
 \sup_{(\rho,K_p,K_i)\in\mathcal D}
 (P^0)^T\rho
\]
over designs satisfying the complete required nonrotational spectral
margin and every declared nonlinear-event constraint.
The delays $\tau_i$ are exogenous. A finite event set is a specified
security test, not robustness to arbitrary disturbances.
The initialized dispatch $(P^0)^T\mathbf1$ is the denominator for
the replacement fraction; net load and installed capacity are not
substituted for it.

All identities involving fixed exported matrices apply to the
repository's fixed-support, fixed-equilibrium current-sharing chart.
The internal per-unit device equations and trims are held fixed
while their terminal-current contributions are scaled. Physical
rating/inertia interpretation must use that same normalization.
Redispatch, topology changes, removal of device states at endpoints,
or different nonlinear converter models require a new derivation
or explicitly bounded model terms.
"""
regional=cut(r'\section{A finite regional replacement',r'\section{Claim ledger')
regional=regional[:regional.index(r'\subsection{What must be verified')]
regional=regional.replace(
    regional[regional.index(r'\paragraph{What is being closed.}'):regional.index(r'\subsection{A finite identity')],
    r"""A regional necessary bound must permit the modal patterns to move
with the design. The following construction uses complete analytic
characteristic matrices and bounds root clusters on fixed contours.
It does not impose the eigenvector compatibility conditions used
by the separate gain-construction method.

""")
regional=regional.replace('let $\\Delta_0(s)',r'define $\mathcal D_d=\{d:p_0+d\in\mathcal D\}$, and let $\Delta_0(s)')
regional=regional.replace(r'\sup_{d\in\mathcal D,',r'\sup_{d\in\mathcal D_d,')
regional=regional.replace(r'\sup_{\Gamma,\mathcal D}',r'\sup_{\Gamma,\mathcal D_d}')
regional=regional.replace('If every physical root satisfies', 'If every required nonrotational physical root satisfies')
regional=regional.replace(
    'Let $R_d$ contain the corresponding nonnegative half-widths, repeated\nfor the two occurrences of each replacement coordinate.',
    r"""For increments $l_a\le d_a\le u_a$, put
$r_a=\max(|l_a|,|u_a|)$. Let $R_d$ contain these nonnegative
displacement radii, repeated for the two occurrences of each
replacement coordinate. They equal box half-widths only when
$p_0$ is its center.""")
regional=regional.replace(
    r'\subsection{Stronger route: retain algebraic variables to remove that remainder}',
    r'\subsection{Affine descriptor realization with no parameter remainder}')
regional=regional.replace(
    'A finite covering closes the contour.',
    'A finite covering establishes this contour condition. The displayed constants\napply in these coordinates; any nonunit state/equation scaling requires\nrederiving the perturbation bounds in the scaled coordinates.')

algorithm=r"""
\section{A reproducible analytical design and bound procedure}
\label{sec:algorithm}
The optimization and certification roles are distinct. The exact
all-PLL formula in Eq.~\eqref{eq:allgain} constructs candidates.
The regional outer problem permits arbitrary admissible modal
patterns and limits what any such candidate can achieve.

\begin{enumerate}
\item Freeze the sharing convention, fixed delay vector, physical
gain limits, replacement region, event set, and output definitions.
Solve and audit the equilibrium; verify the algebraic chart.
\item At a candidate center, assemble the exact descriptor with
exponential delay factors. Use numerical roots to propose contours;
verify their counts, boundary regularity and parameter-box envelopes.
Numerical root locations alone do not close this step.
\item Integrate the reference moments and $J_\Gamma$ (and, if used,
$Q_\Gamma$) with error enclosures. Evaluate the walk remainder.
If a panel or algebraic chart cannot be verified, subdivide it or
return an unresolved box; do not fill it by interpolation.
\item Solve the necessary LP or lifted outer problem and verify an
upper objective bound or a strict controller-exclusion witness.
Compare this upper bound against the box-only replacement ceiling.
\item Use the analytical gain map or a trust-region local step to
propose a design, then evaluate the full delayed spectrum and all
events. Failed designs do not supply a feasible lower bound.
\item Report the same-domain feasible and upper values separately,
with their uncertainty. Subdivide parameter regions only under a
declared branch rule if a wider bound is sought.
\end{enumerate}
No global convergence theorem for the predictor/corrector is claimed.
A global upper bound would require a verified cover of the full
design domain, including all support charts and unresolved boxes.
The full event-feasible problem generally has more active constraints
than the two-row reduced LP; its SG-retention structure cannot be
inferred from a two-anchor theorem.
"""
validation=cut(r'\section{Targeted validation',r'% BEGIN REGIONAL CERTIFICATE')
validation=validation.replace('The independent delayed full-model calculation',
    'An independently assembled Jacobian of the same model, followed by the delayed calculation,')
validation=validation.replace(
    'infinity residual $2.807\\times10^{-9}$.',
    'infinity residual $2.807\\times10^{-9}$. This cross-check uses the same\nnonlinear device equations; it is not an independent-simulator or\nphysical-data validation.')
data=json.loads(resolve('experiments/all_pll_gain_map_validation_20261003/SUMMARY.json').read_text())
validation+=r"""
\paragraph{Reproducible gain vector.}
For buses 30--39, every listed replacement fraction is $\rho_i=0.876$.
These are constructed feasible-candidate gains, not optimal gains.
\begin{center}\small
\begin{tabular}{rrr}\toprule
Bus & $K_{p,i}$ & $K_{i,i}$\\\midrule
"""
for i,(kp,ki) in enumerate(zip(data['Kp'],data['Ki']),30):
    validation+=f'{i} & {kp:.8f} & {ki:.8f}'+r'\\'+'\n'
validation+=r'\bottomrule\end{tabular}\end{center}'+'\n'

checks=r"""
\section{Bound checks, falsification and regional evidence}
\label{sec:boundchecks}
\status{NUMERICALLY\_VALIDATED}
The descriptor determinant and affine-factor identities were checked
at the three frozen frequencies 0.5, 5 and 10 Hz, with
$s=-0.05+\ii 2\pi f$, heterogeneous delays uniformly spread from
0 to 40 ms, replacement displacement radii 0.001, and ordinary gain
radii 1\%. This delay pattern is different from the uniform-delay
nonlinear trial; their lower and upper evidence cannot be mixed.
The largest previously recorded normalized identity discrepancy was
$1.011\times10^{-12}$.

The new checks retain those same points and box. They compare the
scalar logarithm-tail bound with the entrywise closed-walk majorant.
The table contains pointwise bounds on the trace logarithm, not
eigenvalue shifts or a sampled replacement certificate.
\begin{center}\small
\begin{tabular}{rrrr}\toprule
$f$ [Hz] & scalar tail & walk tail & reduction factor\\\midrule
"""
df=pd.read_csv(OUT/'TABLE_01_WALK_REMAINDER.csv')
for _,r in df.iterrows():
    checks+=f"{r.frequency_Hz:g} & {r.scalar_norm_remainder:.6g} & {r.walk_remainder_order1:.6g} & {r.scalar_norm_remainder/r.walk_remainder_order1:.2f}"+r'\\'+'\n'
checks+=r"""\bottomrule\end{tabular}\end{center}
The first- and second-order trace predictions and the cross-site
omission budgets passed all three pointwise numerical checks.
This supports improved bound sharpness at these points; it does
not demonstrate that graph terms increase attainable replacement.

For the separate two-state pure-DDE fixture
\[
 \Delta_{\rm toy}(s,p)=
 \diag(s+1,s+2)-pe^{-s/10}
 \begin{bmatrix}1&1\\1&-1\end{bmatrix},\qquad |p|\le 0.01 ,
\]
the circle centered at $-1$ with radius $1/4$ admits the analytic
whole-circle bound $\kappa\le12/175<1$ and one enclosed root.
Implicit differentiation gives
$\lambda'(0)=e^{1/10}$ and $\lambda''(0)=1.8e^{1/5}$.
The contour Hessian reproduces both values to numerical error.
For the five frozen values of $p$ at 256 and 512 quadrature nodes,
the largest quadratic error is $3.140\times10^{-6}$, below the
analytic cubic-tail bound $5.770\times10^{-5}$.
This is a theorem/check on the declared fixture, not on the power grid.
"""
result=OUT/'CONTOUR_RESULT.json'
if result.exists():
    cc=json.loads(result.read_text())
    checks+=r"""
\subsection{Continuous IEEE--39 contour test}
Complex ball arithmetic at 128-bit precision encloses every point
of a square centered at the preserved slow-mode estimate
$-0.12229688588491999+\ii\,0.32408378548166134$, with half-side
0.04. The region uses uniform fixed delay 40 ms,
$\rho_i\in[0.874,0.876]$ and ordinary gains within 1\% of the
baseline. The matrix inputs are the exact binary64 values parsed
from the frozen exports. The fixed delay is the exact decimal
$0.04=1/25$ s enclosed in ball arithmetic, not its binary64 rounding.
Physical and equilibrium uncertainty
are not enclosed by this numerical certificate.
"""
    q_display=math.ceil(cc['q_max_upper']*1e8)/1e8
    checks+=f"An interval inverse verifies $G(\\rho)$ throughout that box.\nThe contour cover uses {cc['panel_count']} accepted panels and gives\n$\\kappa\\le {q_display:.8f}<1$.\n"
    if cc['count_one_proved']:
        checks+=r"""The enclosed argument-principle integral contains exactly one
integer, 1. Consequently one root remains inside this square for
every permitted parameter combination. Its real part is below
$-0.08229\ \mathrm{s}^{-1}$, so this one cluster satisfies the
$-0.05\ \mathrm{s}^{-1}$ margin throughout the box.
This verifies one cluster; it is not a complete-spectrum or nonlinear
certificate.
"""
    else:
        checks+=r"""The continuous nonsingularity and box envelope pass, but the
enclosed argument-principle integral does not uniquely identify the
root count. The count/moment gate therefore remains unresolved.
No count or regional stability certificate is inferred from rounding.
"""
    checks+=r"""All ball enclosures, parameters and code are saved in the
reproduction package. The coarse first-order scalar moment remainder
is bounded by $4.458179\ \mathrm{s}^{-1}$, much larger than the
required margin. This bound is conservative even though the contour
geometry proves this cluster safe. The sharper walk majorant has
not yet been integrated with verified errors around this contour.
"""
else:
    checks+=r"""\subsection{Continuous IEEE--39 contour test}
The preregistered complex-ball contour verification is pending.
No result is claimed until its count and enclosure gates complete.
"""
checks+=r"""
\subsection{A decisive negative result about the chosen regional box}
\status{NEGATIVE\_RESULT}
To the reported numerical precision, the feasible 87.6\% design has
$\rho_i=0.876$ at every candidate and its gains lie inside the
uniform-delay verification box. Its reported objective coincides with
the box-only upper bound
\[
 P_{\rm box}=0.876(P^0)^T\mathbf1=4732.818715\ {\rm MW}.
\]
Conditional on true feasibility of the exact upper-corner design,
this is the optimum \emph{of that box}, by monotonicity of a positive
dispatch objective. No spectral theorem is needed for this fact.
It is not an irreducible SG-support bound and is not a meaningful
near-global-optimality result. A purported upper bound below this
value would contradict the numerical witness and require resolving
that discrepancy. For strict arithmetic bookkeeping, the binary64
literal $0.876$ exceeds the exact decimal endpoint by
$1/1125899906842624000$; the numerical witness is not treated as
an interval-certified corner point. No rigorous lower--upper gap
is claimed from these rounded values.
Thus even perfect verification of this regional box cannot establish
the main application claim of a physical interior replacement limit.
The manuscript does not expand the box after observing this result.
"""
discussion=r"""
\section{Discussion, contribution boundary and limitations}
\begin{center}\small
\begin{tabular}{p{0.29\linewidth}p{0.61\linewidth}}\toprule
Result & Evidentiary scope\\\midrule
Affine descriptor and gain map & Exact identities in the fixed chart.\\
Cluster Hessian and walk tails & Proved under the stated contour/envelope hypotheses.\\
Regional LP / lifted outer bound & Necessary relaxation; verified data and dual bound required.\\
87.6\% candidate, five events & Full-model numerical result for uniform 40 ms delay.\\
18--46-fold tail reduction & Three frozen algebraic points, not a capacity increase.\\
Physical maximum or global gap & Not established.\\
Local PLLs fundamentally insufficient & Not established with freely varying modal patterns.\\
GSP compression or physical cycle law & Not established.\\\bottomrule
\end{tabular}
\end{center}

\paragraph{What is technically added.}
The model-specific lift links replacement and PLL gain increments
to a finite characteristic perturbation without a network-inverse
parameter remainder. Cluster curvature resolves pairwise collective
actions, and the walk budget quantifies the error of dropping them.
The outer relaxation allows all admissible gains and modal patterns;
a strict finite dual witness would therefore be stronger than
failure to preserve a prescribed eigenvector. This is an applied
bound construction assembled from established tools. Bibliographic
priority for the combination has not been proved.

\paragraph{What still prevents a strong application paper.}
An informative verified regional upper bound or controller-exclusion
witness below the box-only ceiling is missing. The preserved trial
also lacks a paired nonlinear fixed-gain comparison and a fair
fixed-structure spectral-abscissa/SQP comparator with identical
constraints. A controller tuned with $I+L+L^2$ terms has not been
implemented or shown necessary. The graph Fourier basis may express
interactions, but no low-frequency truncation theorem follows.

\paragraph{Physical limits.}
The certificate concerns a local linearization of a specific
nonlinear model. Nonlinear event tests remain separate; no global
nonlinear stability, arbitrary-disturbance safety, or unmodeled
converter-current-limit guarantee follows. Delays may be
heterogeneous in the derivations, but the preserved nonlinear
feasible result has uniform delays. The graph moment in the
appendices concerns a declared swing-network class and has not
been transferred to the full IEEE--39 converter model.

\section{Conclusion}
Collective effects can be made explicit at the level of finite
parameter changes: the full delayed network produces a cluster
Hessian and a closed-walk remainder, which lead to necessary
regional replacement bounds. The derivations permit moving modes
and fixed heterogeneous delays without a low-frequency delay
approximation. The available IEEE--39 results validate algebraic
sharpness and one feasible analytical retuning construction.
One continuous complex-ball calculation also certifies retention and
margin stability of a slow cluster over its complete local box.
These results do not yet establish a physically limiting replacement ceiling
or superior co-design performance. The next decisive result is an
informative bound on a preregistered region extending beyond the
current feasible corner, followed by a comparison under the same
nonlinear security contract.

\section*{Reproducibility and research integrity}
The main source is edited in place; pre-edit source bytes and prior
experiment artifacts are preserved. New checks, protocol, exact
input copies, package versions, tables, claim ledger and manifests
are under \path{experiments/regional_paper_closure_20261003/}.
Numerical table entries in the new main text are generated from
the saved result files. The closure performs linear-algebra and
contour checks only; it does not rerun nonlinear trajectories.
The historical event evidence remains explicitly attributed to its
frozen experiment. Three independent reviews audited the equations,
evidence provenance and literature boundary; their corrections and
unresolved objections are recorded in the package.

\clearpage
\appendix
"""
appendices=cut(r'\section{First correction',r'\section{Targeted validation')
appendices=appendices.replace(r'\section{Closing the design problem before simulation}',
                              r'\section{All-PLL construction and its exact scope}')
bib=old[old.index(r'\begin{thebibliography}'):old.index(r'\end{thebibliography}')]
bib+=r"""
\bibitem{markovic} U. Markovi\'c, O. Stanojev, P. Aristidou,
E. Vrettos, D. Callaway, and G. Hug,
\emph{Understanding Small-Signal Stability of Low-Inertia Systems},
IEEE Transactions on Power Systems, 36(5):3997--4017, 2021.
\url{https://doi.org/10.1109/TPWRS.2021.3061434}.
\bibitem{chuteng} Z. Chu and F. Teng,
\emph{Stability Constrained Optimization in High IBR-Penetrated
Power Systems---Part I: Constraint Development and Unification},
author version, arXiv:2307.12151.
\url{https://arxiv.org/abs/2307.12151}.
\bibitem{jia} H. Jia et al.,
\emph{Converter-Driven Stability Constrained Generation Expansion
for High IBG-Penetrated Power Systems},
IEEE Transactions on Power Systems, 2025.
\url{https://doi.org/10.1109/TPWRS.2025.3558879}.
\bibitem{peaucelle} D. Peaucelle and Y. Ebihara,
\emph{Robust stability analysis of discrete-time systems with
parametric and switching uncertainties}, IFAC World Congress, 2014.
\url{https://skoge.folk.ntnu.no/prost/proceedings/ifac2014/media/files/0282.pdf}.
\bibitem{bindel} D. Bindel and A. Hood,
\emph{Localization Theorems for Nonlinear Eigenvalue Problems},
SIAM Journal on Matrix Analysis and Applications, 2013.
\url{https://www.cs.cornell.edu/~bindel/papers/2013-simax.pdf}.
\bibitem{borgioli} F. Borgioli, W. Michiels, D. Lu, and B. Vandereycken,
\emph{A globally convergent method to compute the real stability
radius for time-delay systems}, author manuscript, 2019.
\url{https://www.unige.ch/math/vandereycken/papers/published_Borgioli_MLV.pdf}.
\bibitem{arb} F. Johansson,
\emph{Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval
Arithmetic}, IEEE Transactions on Computers, 66(8):1281--1292, 2017.
\url{https://doi.org/10.1109/TC.2017.2690633}.
\end{thebibliography}
\end{document}
"""
checks=checks.replace('Complex ball arithmetic at',r'Complex ball arithmetic \cite{arb} at')
text=pre+model+regional+(OUT/'COLLECTIVE_BOUNDS.tex').read_text()+algorithm+validation+checks+discussion+appendices+bib
# Self-contained native-editor source: no external input or graphics dependency.
assert r'\input{' not in text and r'\includegraphics' not in text
labels=re.findall(r'\\label\{([^}]+)\}',text)
assert len(labels)==len(set(labels)), 'duplicate labels'
refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)
assert set(refs)<=set(labels), sorted(set(refs)-set(labels))
keys=set(re.findall(r'\\bibitem\{([^}]+)\}',text))
for cite in re.findall(r'\\cite\{([^}]+)\}',text):
    assert set(cite.split(','))<=keys,cite
TARGET.write_text(text,encoding='utf-8')
(OUT/'MANUSCRIPT_BUILD.json').write_text(json.dumps({
    'source_snapshot_sha256':hashlib.sha256((OUT/'THEORY_before_closure.tex').read_bytes()).hexdigest(),
    'written_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),
    'labels':len(labels),'references':len(keys),'unresolved_internal_references':0,
    'continuous_result_integrated':result.exists()},indent=2)+'\n')
print('Updated open source:',TARGET)
