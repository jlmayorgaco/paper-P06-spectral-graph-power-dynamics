"""Append the independently reviewed poster closure evidence to the open article."""
from pathlib import Path
import hashlib,json,shutil
import pandas as pd
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
THEORY=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
OLD=ROOT/'experiments/all_pll_gain_map_validation_20261003'
snapshot=OUT/'THEORY_before_poster_closure.tex'
if not snapshot.exists():shutil.copyfile(THEORY,snapshot)
original=THEORY.read_text(encoding="utf-8-sig")
analytic=pd.read_csv(OLD/'TABLE_08_NONLINEAR_EVENTS.csv')
fixed=pd.read_csv(OUT/'TABLE_01_FIXED_GAIN_EVENTS.csv')
p=analytic.merge(fixed,on=['bus','delta'],suffixes=('_analytic','_fixed'),validate='one_to_one')
assert len(p)==5 and p.pass_analytic.all() and p.pass_fixed.all()
for field in ['F','R','slack','Vmin','Vmax']:
    p[field+'_difference']=p[field+'_analytic']-p[field+'_fixed']
p.to_csv(OUT/'TABLE_03_PAIRED_COMPARISON.csv',index=False)
summary=json.loads((OLD/'SUMMARY.json').read_text(encoding="utf-8-sig"))
ports=pd.read_csv(ROOT/'experiments/graph_gsp_codesign_20261003/model/ports.csv')
gains=pd.DataFrame({'bus':range(30,40),'rho':summary['rho'],'Kp':summary['Kp'],'Ki':summary['Ki'],
    'tau_ms':summary['tau_ms'],'P0_MW':ports.P0,
    'SG_MW':ports.P0*(1-pd.Series(summary['rho'])),
    'GFL_MW':ports.P0*pd.Series(summary['rho'])})
gains.to_csv(OUT/'TABLE_04_DESIGN_VECTOR.csv',index=False)
rows=[]
for r in p.itertuples():
    rows.append(f'{r.bus}, {r.delta:+.0f} & {r.F_fixed:.6f} & {r.F_analytic:.6f} & {r.R_fixed:.6f} & {r.R_analytic:.6f} & Pass / Pass'+r'\\')
section=r"""
\section{Targeted closure: physical damping and matched nonlinear comparison}
\label{sec:posterclosure}
This follow-up adds two declared tests without changing the historical
design, events, gains, or frozen result files. Its protocol and outputs
are stored in \texttt{experiments/poster\_julia\_validation\_20261003/}.
The harmonic tests are algebraic; the event comparator uses Julia 1.11.9.

\subsection{Finite-delay physical damping law in IEEE--39}
\status{NUMERICALLY\_VALIDATED}
For the historical baseline and analytical designs, at each of
$f=0.5,5,10$ Hz, increase the delay at one of buses 30--39 by exactly
the numerical increment 1 ms, starting from uniform 40 ms.
All $2\times3\times10=60$ cases are retained. No case was chosen by
its sign or stability outcome. At each point, independently eliminate
the 194 hidden states in the old and modified 204-state operators and
compare the difference with the finite rank-one formula
in Appendix~\ref{app:physical}.
The physical mass scaling is
$M_i=2H_iS_i^0(1-\rho_i)$, in MW\,s for per-unit speed.

The largest relative error in $\delta Z=\gamma_i u_i v_i^*$ is
$4.85\times10^{-8}$. Every computed Hermitian increment has one
positive and one negative eigenvalue above the declared relative
threshold $10^{-7}\|\delta D_H\|_2$; the other eight eigenvalues
have relative magnitude at most $1.40\times10^{-9}$.
The minimum squared sine between the two propagation paths is
$0.00327$, so none of these sampled paths is numerically collinear.
Hidden-block condition numbers and all signed eigenvalues are
reported in \texttt{TABLE\_02\_FINITE\_DELAY\_RANK\_LAW.csv}.
This supports the finite physical identity at the declared points.
It does not certify an interval of frequencies, imply stability from
the sign of the increment, or improve a replacement-capacity bound.

\begin{figure}[htbp]\centering
\includegraphics[width=.87\linewidth]{reports/poster/ias2026/collective_interaction_20261003/generated/figures/rank_two.pdf}
\caption{Physical damping redistribution at the preregistered illustrative
frequency 5 Hz, analytical design. One PLL delay is increased from
40 to 41 ms at a time. The symmetric-log scale exposes both signs.
This is a harmonic linear-model check, not a nonlinear stability test.}
\label{fig:rankclosure}
\end{figure}

\subsection{Same replacement, unchanged versus analytical gains}
\status{NUMERICALLY\_VALIDATED}
The missing matched comparator is now executed at the same
$\rho_i=0.876$ and $\tau_i=40$ ms as the analytical candidate.
Its gains remain at the common baseline values.
The canonical method-of-steps solver uses Rodas5P, tolerance $10^{-9}$,
maximum step 0.01 s, and a 60 s event horizon. Both phase-derived
frequency and RoCoF metrics use the canonical 0.5 s window over all
39 bus voltages. The five nominal event magnitudes correspond to
the existing constant-impedance perturbation contract.
Equilibrium residual and gain bounds are checked before execution.

\begin{table}[htbp]\centering\small
\caption{All five matched events. $F$ is peak absolute frequency deviation
in Hz; $R$ is peak absolute RoCoF in Hz/s. A denotes the analytical map;
fixed retains the original gains. All original guards also pass.}
\begin{tabular}{rrrrrl}\toprule
Bus, MW & $F_{\rm fixed}$ & $F_{\rm A}$ & $R_{\rm fixed}$ & $R_{\rm A}$ & Guards\\\midrule
""" + '\n'.join(rows) + r"""
\\[-2.5mm]\bottomrule
\end{tabular}\label{tab:pairedclosure}
\end{table}

Both designs pass all five events. Analytical retuning increases
frequency excursion by $0.000177$--$0.000207$ Hz and reduces RoCoF
by $0.0000744$--$0.000435$ Hz/s across these cases. Its SG actuator
slack is also slightly smaller. These are small numerical tradeoffs,
not a demonstrated engineering advantage; no refinement study was
performed to turn these small differences into accuracy-certified
performance claims. The purpose of the exact gain map is achieved:
the assigned PLL eigenpair is preserved after the prescribed
replacement step. This matched trial does not show that retuning
was necessary for feasibility, that it maximized replacement, or
that it dominates fixed gains.

\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{reports/poster/ias2026/collective_interaction_20261003/generated/figures/events.pdf}
\caption{Preserved analytical-design trajectories for every frozen
event, plotted across the full 60 s horizon. Frequency and RoCoF use
the declared 0.5 s phase window. The new comparator's complete
trajectories and paired metrics are included separately in the archive.}
\label{fig:allfrozen}
\end{figure}

\paragraph{Updated contribution boundary.}
The physical finite-delay theorem, the collective contour Hessian,
and the walk remainder provide explicit equations for network-mediated
effects. The numerical evidence now includes all 60 physical-law
checks and a matched five-event comparator. No informative
verified physical replacement maximum follows. The posterior
comparison therefore supports a theory-and-validation poster;
the stronger application claim remains open.
"""
if r'\section{Targeted closure:' in original:
    a=original.index(r'\section{Targeted closure:')
    b=original.index(r'\section{Discussion, contribution boundary and limitations}',a)
    original=original[:a]+original[b:]
text=original
if r'\usepackage{graphicx,microtype}' not in text:
    text=text.replace(r'\hypersetup',r'\usepackage{graphicx,microtype}'+'\n'+
                      r'\setlength{\emergencystretch}{2em}'+'\n'+r'\hypersetup',1)
text=text.replace(r'\author{Research manuscript}',r'\author{Jorge Luis Mayorga Taborda\\Universidad de los Andes, Colombia}')
text=text.replace('Unchanged gains also\npass its spectral test.','Unchanged gains also\npass its spectral test and, in a matched follow-up, all five nonlinear events.')
text=text.replace('the spectral margin. Its nonlinear events were not run, so nonlinear\nsuperiority or equivalence is undetermined.',
                  'the spectral margin. The subsequent matched nonlinear comparison\n'+r'in Section~\ref{sec:posterclosure} finds that both designs pass all five'+'\nevents, with small opposing changes in frequency and RoCoF metrics.\nNeither superiority nor exact performance equivalence is established.')
text=text.replace('Its numerical evidence has two distinct\nscopes: identity and enclosure tests of the bound machinery, and one\npreviously frozen full-model feasible retuning step.',
                  'Its numerical evidence includes identity and enclosure tests of the\nbound machinery, a frozen full-model feasible retuning step, and a\nmatched five-event fixed-gain follow-up.')
text=text.replace('The preserved trial\nalso lacks a paired nonlinear fixed-gain comparison and a fair\nfixed-structure spectral-abscissa/SQP comparator with identical\nconstraints.',
                  'The matched nonlinear fixed-gain comparison now shows both designs\nfeasible, with small opposing metric changes. A fair fixed-structure\nspectral-abscissa/SQP comparator with identical constraints remains missing.')
text=text.replace('IEEE--39 path independence has not been tested in this theory-only work.',
                  r'Section~\ref{sec:posterclosure} checks path independence numerically at 60 frozen IEEE--39 harmonic cases; no frequency interval is certified.')
text=text.replace(r'\section{Physical collective damping and its finite latency law}',
                  r'\section{Physical collective damping and its finite latency law}\label{app:physical}')
text=text.replace(r'\label{app:physical}\label{app:physical}',r'\label{app:physical}')
text=text.replace('Sources inspected are\n'+
                  r'\path{experiments/nonlinear_codesign_20261001/ReducedDAE.jl},'+'\n'+
                  r'\path{src/bnd_model_expN/PDExactDesignN.jl}, and'+'\n'+
                  r'\path{experiments/graph_gsp_codesign_20261003/DelayedEvents.jl}.',
                  r'''Sources inspected include:
\begin{itemize}
\item \path{experiments/nonlinear_codesign_20261001/ReducedDAE.jl};
\item \path{src/bnd_model_expN/PDExactDesignN.jl};
\item \path{experiments/graph_gsp_codesign_20261003/DelayedEvents.jl}.
\end{itemize}''')
text=text.replace(r'''Because $G,N$ are affine in $\rho$, computable bounds are
$a=\sum_j r_{\rho,j}\|G_0^{-1}G_{\rho_j}\|_2$ and
$v=\sum_j r_{\rho,j}\|G_0^{-1}(N_{\rho_j}+G_{\rho_j}V_0)\|_2$.''',
                  r'''Because $G,N$ are affine in $\rho$, computable bounds are
\begin{align*}
a&=\sum_j r_{\rho,j}\|G_0^{-1}G_{\rho_j}\|_2,\\
v&=\sum_j r_{\rho,j}\|G_0^{-1}(N_{\rho_j}+G_{\rho_j}V_0)\|_2.
\end{align*}''')
text=text.replace(r'\section{Discussion, contribution boundary and limitations}',
                  section+'\n'+r'\section{Discussion, contribution boundary and limitations}',1)
text=text.replace('the saved result files. The closure performs linear-algebra and\ncontour checks only; it does not rerun nonlinear trajectories.',
                  'the saved result files. The initial closure performed linear-algebra\nand contour checks. The subsequent poster closure adds the 60 physical-law\nchecks and five new fixed-gain nonlinear comparator trajectories under\n'+r'\path{experiments/poster_julia_validation_20261003/}.')
text=text.replace('The manuscript establishes three conditional mathematical results:\nan affine full-model intervention representation, a cluster\ncurvature/closed-walk bound, and a controller-bounded necessary\nreplacement relaxation.',
                  'The manuscript establishes four conditional mathematical results:\na finite physical damping-redistribution identity, an affine full-model\nintervention representation, a cluster curvature/closed-walk bound,\nand a controller-bounded necessary replacement relaxation.')
aa=text.index(r'\begin{abstract}')+len(r'\begin{abstract}')
bb=text.index(r'\end{abstract}',aa)
text=text[:aa]+r"""
Local PLL tuning and measurement delay act through the whole electrical
network. Starting from the exact delayed linearization, we show that
changing one PLL delay produces a rank-one update of the physical
torque--speed impedance. With independent propagation paths and a
nonzero update, its Hermitian damping increment has one positive
and one negative eigenvalue: latency redistributes harmonic damping
across collective directions.
For a fixed-equilibrium SG/GFL sharing contract, retaining network
algebraic variables also gives an affine descriptor perturbation.
A contour identity yields finite spectral-cluster changes, their full
parameter Hessian, and an explicit closed-walk remainder.
Necessary cluster conditions lead to linear and semidefinite outer
relaxations of secure replacement, conditional on verified contour
and remainder bounds.
In IEEE--39, all 60 declared harmonic tests confirm the finite-delay
identity numerically. Three frozen algebraic points give walk remainder
bounds 18--46 times smaller than a scalar-norm bound. Complex ball
arithmetic certifies one root cluster throughout a declared parameter box.
A preserved analytical gain construction at 40 ms delay yields 87.6\%
GFL and passes five nonlinear events. A new matched Julia comparison
finds that unchanged gains also pass all five, with small opposing
frequency and RoCoF changes. The results establish conditional
interaction laws and a feasible construction; they do not establish
replacement superiority or a physically limiting maximum.
"""+text[bb:]
THEORY.write_text(text,encoding='utf-8',newline='\r\n')
(OUT/'GENERATED_CLOSURE_SECTION.tex').write_text(section,encoding='utf-8')
(OUT/'ARTICLE_UPDATE_PROVENANCE.json').write_text(json.dumps({
    'before_sha256':hashlib.sha256(snapshot.read_bytes()).hexdigest(),
    'after_sha256':hashlib.sha256(THEORY.read_bytes()).hexdigest(),
    'no_historical_result_mutation':True,'comparison_events':5,'rank_law_cases':60},indent=2)+'\n')
print('Article evidence updated; inspect any remaining historical comparator wording.')
