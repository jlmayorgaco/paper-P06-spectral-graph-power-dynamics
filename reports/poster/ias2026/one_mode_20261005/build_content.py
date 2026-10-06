"""Build the new poster from immutable evidence and a preserved modular template."""
from pathlib import Path
import hashlib,json,re,shutil,sys,unicodedata
sys.path.insert(0,str(Path(__file__).resolve().parent/'build/vendor'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import qrcode

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
EXP=ROOT/'experiments/poster_julia_validation_20261003'
OLD=ROOT/'experiments/all_pll_gain_map_validation_20261003'
REG=ROOT/'experiments/regional_paper_closure_20261003'
FIG=HERE/'generated/figures';DATA=HERE/'generated/data'
DATA.mkdir(exist_ok=True)
GREEN='#00664B';BLUE='#176494';GOLD='#C99A20';RED='#BC3030';INK='#17372D'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.labelsize':17,
    'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,
    'text.color':INK,'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,
    'pdf.fonttype':42,'savefig.facecolor':'white'})
sources=[]
def source(p):
    sources.append({'path':str(p.relative_to(ROOT)).replace('\\','/'),
                    'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return p
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.08,
                metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight',pad_inches=.08)
    plt.close(fig)
def write(path,s): path.write_text(s,encoding='utf-8',newline='\n')
S=json.loads(source(OLD/'SUMMARY.json').read_text(encoding="utf-8-sig"))
E=pd.read_csv(source(OLD/'TABLE_08_NONLINEAR_EVENTS.csv'))
W=pd.read_csv(source(REG/'TABLE_01_WALK_REMAINDER.csv'))
R=pd.read_csv(source(EXP/'TABLE_02_FINITE_DELAY_RANK_LAW.csv'))
RS=json.loads(source(EXP/'RANK_LAW_SUMMARY.json').read_text(encoding="utf-8-sig"))
CR=json.loads(source(REG/'CONTOUR_RESULT.json').read_text(encoding="utf-8-sig"))
panels=json.loads(source(REG/'PANELS.json').read_text(encoding="utf-8-sig"))
for p in [OLD/'SUMMARY.json',OLD/'TABLE_08_NONLINEAR_EVENTS.csv',
          REG/'TABLE_01_WALK_REMAINDER.csv',REG/'CONTOUR_RESULT.json',
          EXP/'TABLE_02_FINITE_DELAY_RANK_LAW.csv',EXP/'RANK_LAW_SUMMARY.json']:
    shutil.copyfile(p,DATA/p.name)

# Figure 1: all physical branches, all ten candidate buses, all event buses.
coordinates=source(HERE/'reference_template/network_ieee39.tikz').read_text(encoding="utf-8-sig")
pos={int(b):(float(y),float(x)) for b,x,y in re.findall(
    r'\\coordinate \(b(\d+)\) at \(([\d.]+)mm,([\d.]+)mm\)',coordinates)}
branches=pd.read_csv(source(ROOT/'reports/experiment_D/inputs/branch.csv'))
fig,ax=plt.subplots(figsize=(10.5,6.1))
for _,e in branches.iterrows():
    a,b=pos[int(e.src_bus)],pos[int(e.dst_bus)]
    ax.plot([a[0],b[0]],[a[1],b[1]],color='#a6b9ad',lw=2,zorder=1)
for bus,(x,y) in pos.items():
    color=BLUE if bus>=30 else (RED if bus in [8,16,29] else '#e4ede7')
    ax.scatter(x,y,s=490 if bus>=30 else 365,c=color,edgecolor=INK,lw=.7,
               marker='s' if bus in [8,16,29] else 'o',zorder=3)
    ax.text(x,y,str(bus),ha='center',va='center',fontsize=14,
            color='white' if bus>=30 or bus in [8,16,29] else INK,zorder=4)
ax.set_aspect('equal');ax.axis('off')
ax.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=BLUE,markersize=13,label='SG / GFL sites 30–39'),
                   Line2D([],[],marker='s',color='none',markerfacecolor=RED,markersize=12,label='Frozen event sites')],
          loc='lower center',bbox_to_anchor=(.5,-.12),ncol=2,frameon=False,fontsize=16)
fig.tight_layout();save(fig,'network')

# Figure 2: the finite rank-two damping law at the preregistered 5 Hz.
r=R[(R.design=='analytic')&(R.frequency_Hz==5.)]
fig,ax=plt.subplots(figsize=(9.3,4.1))
ax.axhline(0,c=INK,lw=1)
ax.bar(r.bus-.18,r.lambda_max,width=.36,color=BLUE,label=r'$\lambda_+(\delta D_H)>0$')
ax.bar(r.bus+.18,r.lambda_min,width=.36,color=RED,label=r'$\lambda_-(\delta D_H)<0$')
ax.set_yscale('symlog',linthresh=1)
ax.set_xticks(r.bus);ax.set_xlabel('PLL site receiving +1 ms')
ax.set_ylabel(r'$\delta D_H$ eigenvalue [MW]')
ax.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.24),fontsize=15)
fig.tight_layout();save(fig,'rank_two')

# Figure 3: exact pointwise logdet remainder majorants, not capacity gains.
fig,ax=plt.subplots(figsize=(8.8,4.0))
x=np.arange(len(W));width=.25
ax.bar(x-width,W.scalar_norm_remainder,width,color='#9ea9a2',label='Scalar norm')
ax.bar(x,W.walk_remainder_order1,width,color=BLUE,label='Closed walks')
ax.bar(x+width,W.actual_error_order1,width,color=GOLD,label='Actual error')
ax.set_yscale('log');ax.set_xticks(x,[str(v).rstrip('0').rstrip('.') for v in W.frequency_Hz])
ax.set_xlabel('Frequency [Hz]');ax.set_ylabel('Trace-log remainder')
ax.legend(frameon=False,ncol=3,fontsize=13,loc='upper center',bbox_to_anchor=(.5,1.22))
fig.tight_layout();save(fig,'walk_bounds')

# Figure 4: every frozen event, complete 60 s traces, no trajectory selection.
fig,axes=plt.subplots(1,2,figsize=(13.6,4.0),sharex=True)
colors=[BLUE,GREEN,GOLD,'#9358a6',RED]
for color,row in zip(colors,E.itertuples()):
    t=pd.read_csv(source(OLD/'nonlinear'/f'bus{row.bus}_{int(row.delta)}_trajectory.csv'))
    label=f'{row.bus}: {int(row.delta):+d} MW'
    axes[0].plot(t.time_after_event_s,t.Fmax_Hz,color=color,lw=1.9,label=label)
    axes[1].plot(t.time_after_event_s,t.Rmax_Hz_s,color=color,lw=1.9)
for ax in axes:
    ax.axhline(.5,c=RED,lw=1.4,ls='--');ax.set_xlim(0,60);ax.set_ylim(0,.54)
    ax.set_xlabel('Time after event [s]');ax.grid(axis='y',color='#e4ede7')
axes[0].set_ylabel(r'Max. windowed $|\Delta f|$ [Hz]')
axes[1].set_ylabel('Max. windowed RoCoF [Hz/s]')
fig.legend(*axes[0].get_legend_handles_labels(),ncol=5,frameon=False,fontsize=13,
           loc='upper center',bbox_to_anchor=(.5,1.04))
fig.tight_layout(rect=(0,0,1,.95));save(fig,'events')

# Figure 5: whole-contour interval cover; marker is a nominal estimate only.
center=np.array([-.12229688588491999,.32408378548166134])
segs=[np.array([p['a'],p['b']])+center for p in panels]
fig,ax=plt.subplots(figsize=(6.6,4.0))
lc=LineCollection(segs,array=np.array([p['q_upper'] for p in panels]),cmap='viridis',linewidth=3.2)
lc.set_clim(0,1);ax.add_collection(lc)
ax.plot(*center,'x',ms=10,mew=2,color=INK,label='Nominal root estimate')
ax.axvline(-.05,color=RED,ls='--',lw=1.5)
ax.text(-.048,.284,'Security margin',rotation=90,va='bottom',fontsize=12,color=RED)
ax.set_xlim(-.175,-.025);ax.set_ylim(.273,.377)
ax.set_xlabel(r'Re $s$ [s$^{-1}$]');ax.set_ylabel(r'Im $s$ [rad/s]')
ax.legend(frameon=False,fontsize=11,loc='upper left')
cb=fig.colorbar(lc,ax=ax,pad=.02);cb.set_label('Interval envelope q',fontsize=13)
fig.tight_layout();save(fig,'contour')

url='https://github.com/jlmayorgaco/paper-P06-spectral-graph-power-dynamics/tree/codex/collective-interaction-bounds-20261003'
qrcode.make(url).save(HERE/'assets/repo_qr.png')
ratios=W.scalar_norm_remainder/W.walk_remainder_order1
macros={
'GFLPercent':f'{S["GFL_percent"]:.1f}',
'GFLMW':f'{S["GFL_MW"]:,.2f}',
'SGMW':f'{S["retained_SG_MW"]:,.2f}',
'AddedMW':f'{S["added_GFL_MW"]:.3f}',
'FreqMax':f'{S["max_frequency_Hz"]:.4f}',
'RocofMax':f'{S["max_RoCoF_Hz_s"]:.4f}',
'VoltageMin':f'{S["min_voltage_pu"]:.3f}',
'VoltageMax':f'{S["max_voltage_pu"]:.3f}',
'SlackMin':f'{S["min_SG_actuator_slack"]:.5f}',
'PanelCount':str(CR['panel_count']),
'RankChecks':str(RS['cases']),
'RankError':f'{RS["max_relative_update_error"]:.2e}',
'WalkMin':str(int(np.floor(ratios.min()))),
'WalkMax':str(int(np.ceil(ratios.max()))),
'TargetFreq':f'{S["target_frequency_Hz"]:.3f}'}
write(HERE/'generated/results.tex','% Values generated from hashed evidence. Do not edit by hand.\n'+
      '\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items())+'\n')
write(HERE/'generated/claims.tex',r"""% Scientific scope shared by the modules.
\newcommand{\FeasibleClaim}{A constructed feasible design; replacement was not optimized.}
\newcommand{\RegionalClaim}{One exported-model root cluster, not the complete spectrum.}
\newcommand{\WalkClaim}{Three frozen algebraic points; not a replacement-capacity bound.}
""")

def replace_content(name,body):
    p=HERE/'sections'/name
    old=p.read_text(encoding="utf-8-sig")
    suffix=old[old.index(r'\ifdefined\PosterFull'):]
    write(p,'% Collective-interaction version. Original poster preserved.\n'+
          '\\long\\def\\PosterSectionContent{%\n'+body.lstrip('\n')+'\n}%\n'+suffix)
def panel(name,width,height,num,title,body,ticket=False):
    env='TicketTwoPanel' if ticket else 'PosterPanel'
    replace_content(name,rf"""\begin{{minipage}}[t][\{height}][t]{{\{width}}}\vspace{{0pt}}
  \begin{{{env}}}{{\{height}}}{{{num}}}{{{title}}}
{body}
  \end{{{env}}}
\end{{minipage}}%""")

head=(HERE/'sections/header_green.tex').read_text(encoding="utf-8-sig")
head=head.replace('WHEN STABLE REPLACEMENTS FAIL TOGETHER','LOCAL TUNING. COLLECTIVE EFFECTS.')
head=head.replace(r'Minimal Incompatibility Hypergraphs for Safe SG $\to$ IBR Transition Planning',
                  r'Delay-aware PLL co-design: interaction laws and IEEE-39 validation')
write(HERE/'sections/header_green.tex',head)

panel('01_motivation_question.tex','PosterOneWidth','PosterTopRowHeight',1,'THE DESIGN QUESTION',r"""
    {\PosterParagraph How much SG can a grid replace with delayed, locally controlled GFL?\par}
    \vspace{2mm}
    {\HeroEquation\color{IASBlue}\[
       P_{\mathrm{GFL}}=\sum_i P_i^0\rho_i
    \]}
    {\FigureSize Initialized dispatch $P_i^0$; sharing $\rho_i$ at each generator site.\par}
    \vspace{2mm}
    {\Subhead\color{IASDarkGreen}The physical PLL\par}
    {\EquationSize\[
    \begin{aligned}
    \dot\vartheta_i&=\omega_{p,i}\\
    t_{f,i}\dot\omega_{p,i}&=\xi_i+K_{p,i}e_i^\tau-\omega_{p,i}\\
    \dot\xi_i&=K_{i,i}e_i^\tau
    \end{aligned}
    \]}
    {\FigureSize $e_i^\tau=e_i(t-\tau_i)$ is the delayed voltage phase error. Delays are fixed; only $\rho,K_p,K_i$ are design variables.\par}
""",True)

panel('02_flagship_ieee39.tex','PosterTwoWidth','PosterTopRowHeight',2,'IEEE-39: A VALIDATED CONSTRUCTION',r"""
    \noindent\begin{minipage}[c]{282mm}\centering
      \includegraphics[width=277mm,height=155mm,keepaspectratio]{generated/figures/network.pdf}
    \end{minipage}\hfill%
    \begin{minipage}[c]{115mm}\raggedright
      {\PosterCondensed\bfseries\fontsize{68pt}{72pt}\selectfont\color{IASBlue}\GFLPercent\%\par}
      {\PosterParagraph GFL dispatch share\par}\vspace{5mm}
      {\MetricSize\bfseries\color{IASDarkGreen}\GFLMW\par}
      {\FigureSize MW GFL\par}\vspace{4mm}
      {\MetricSize\bfseries\color{IASGold!70!black}\SGMW\par}
      {\FigureSize MW retained SG\par}\vspace{4mm}
      {\Subhead\color{IASDarkGreen}40 ms\par}
      {\FigureSize at every PLL\par}
    \end{minipage}\par
    \vspace{3mm}
    {\EquationSize\color{IASBlue}
      $\rho_i:0.875\longrightarrow0.876\qquad\Delta P_{\mathrm{GFL}}=\AddedMW\ \mathrm{MW}$\par}
    \vspace{4mm}
    {\FigureSize\FeasibleClaim\ Five frozen external events pass the declared nonlinear security checks.\par}
""",True)

panel('03_target_vs_path.tex','PosterThreeWidth','PosterTopRowHeight',3,'GAINS FROM AN EQUATION',r"""
    {\FigureSize Fix one nonreal pole $\lambda$ and PLL pattern $q$. Eliminate hidden states:\par}
    {\EquationSize\[
       y=G_{\mathrm{PLL}}(\lambda,\rho)q
    \]}
    {\EquationSize\[
       W_i=\frac{\lambda^2(1+t_{f,i}\lambda)e^{\lambda\tau_i}q_i}{y_i}
    \]}
    {\HeroEquation\color{IASBlue}\[
       K_{p,i}=\frac{\Im W_i}{\Im\lambda}
    \]}
    {\EquationSize\color{IASBlue}\[
       K_{i,i}=\Re W_i-\Re\lambda\,K_{p,i}
    \]}
    {\FigureSize Requires nonsingular eliminated blocks, $y_i\ne0$, $\Im\lambda\ne0$ and admissible gains.\par}
    \vspace{2mm}
    {\FigureSize \textbf{Target: \TargetFreq\ Hz.} One mode assigned; full spectrum and events need separate checks.\par}
""",True)

replace_content('result_strip.tex',r"""
\leavevmode{\setlength{\fboxsep}{0pt}%
\colorbox{IASDeepGreen}{\begin{minipage}[c][\PosterStripHeight][c]{\PosterGridWidth}
\begin{minipage}[c]{215.6mm}\MetricCell{\RankChecks/ \RankChecks}{finite-delay identities}{two designs, three frequencies}\end{minipage}%
\begin{minipage}[c]{2mm}\textcolor{IASPaleGreen}{\rule{.7pt}{30mm}}\end{minipage}%
\begin{minipage}[c]{215.6mm}\MetricCell{$\WalkMin$--$\WalkMax\times$}{tighter remainder}{three pointwise comparisons}\end{minipage}%
\begin{minipage}[c]{2mm}\textcolor{IASPaleGreen}{\rule{.7pt}{30mm}}\end{minipage}%
\begin{minipage}[c]{215.6mm}\MetricCell{5/5}{nonlinear events pass}{Julia method of steps}\end{minipage}%
\begin{minipage}[c]{2mm}\textcolor{IASPaleGreen}{\rule{.7pt}{30mm}}\end{minipage}%
\begin{minipage}[c]{215.6mm}\MetricCell{\PanelCount}{verified contour panels}{one cluster over a full box}\end{minipage}
\end{minipage}}}\par%""")

panel('04_network_closure.tex','PosterGridWidth','PosterFourHeight',4,
      'ONE LOCAL DELAY REDISTRIBUTES COLLECTIVE DAMPING',r"""
    \noindent\begin{minipage}[t]{270mm}\vspace{0pt}
      {\Subhead\color{IASDarkGreen}Physical torque--speed ports\par}
      {\EquationSize\[
      \Delta=sI-A_0-\sum_i b_ic_i^T e^{-s\tau_i}
      \]}
      {\EquationSize\[
      Z=M(\Delta_{\nu\nu}-\Delta_{\nu z}\Delta_{zz}^{-1}\Delta_{z\nu})
      \]}
      {\EquationSize\[
      D_H(\Omega)=\operatorname{Herm}Z(\mathrm{j}\Omega)
      \]}
      {\FigureSize $\nu$: per-unit SG speed deviation; $z$: hidden states.
      $M=\operatorname{diag}[2H_iS_i^0(1-\rho_i)]$.
      Harmonic mechanical supply: $\tfrac12\widehat\nu^*D_H\widehat\nu$.\par}
    \end{minipage}\hfill%
    \begin{minipage}[t]{282mm}\vspace{0pt}
      {\Subhead\color{IASDarkGreen}Finite-change identity\par}
      {\HeroEquation\color{IASBlue}\[
      \delta Z=\gamma_i u_i v_i^*
      \]}
      {\HeroEquation\color{IASBlue}\[
      \delta D_H=\tfrac12(a_i v_i^*+v_i a_i^*)
      \]}
      {\FigureSize $a_i=\gamma_i u_i$; $u_i,v_i$ are the two network propagation paths from a single PLL channel.\par}\vspace{4mm}
      {\PosterParagraph If the paths are independent and $\gamma_i\ne0$, the change has one positive and one negative eigenvalue.\par}
    \end{minipage}\hfill%
    \begin{minipage}[t]{288mm}\vspace{0pt}\centering
      {\Subhead\color{IASDarkGreen}IEEE-39: +1 ms at one site\par}
      \includegraphics[width=281mm,height=87mm,keepaspectratio]{generated/figures/rank_two.pdf}\par
      {\FigureSize Figure: 5 Hz, analytic design; symmetric-log axis. All 60 checks (2 designs $\times$ 3 frequencies $\times$ 10 sites) pass; error $<5\times10^{-8}$.\par}
      \vspace{2mm}
      {\FigureSize Requires invertible hidden blocks. These signs alone do not certify stability.\par}
    \end{minipage}\par
""")

panel('05_policy_atlas.tex','PosterHalfWidth','PosterMiddleRowHeight',5,
      'BEYOND ADDING NODAL EFFECTS',r"""
    {\FigureSize Affine descriptor lift $F-F_0=U\Theta(d)V$;
    $T=VF_0^{-1}U$ and $H=\Theta(d)T$ on a root contour $\Gamma$.\par}
    {\EquationSize\[
       z_\Gamma(d)-z_\Gamma(0)
       =-\frac{1}{2\pi\mathrm{j}m}\oint_\Gamma
                \operatorname{tr}\log(I+H)\,ds
    \]}
    {\FigureSize $z_\Gamma$: mean of $m$ enclosed poles. A uniform $\|H\|<1$ preserves the count; the gauge pole is excluded.\par}
    {\HeroEquation\color{IASBlue}\[
       \mathcal Q_{ab}=
       \Re\frac{1}{2\pi\mathrm{j}m}\oint_\Gamma
                    \operatorname{tr}(E_aTE_bT)\,ds
    \]}
    {\FigureSize $E_a$ selects a physical parameter. Mixed curvature captures two interventions coupled through return paths. These are signed interactions on the intervention-port graph, not a Laplacian.}
""")

panel('06_synthetic_envelope.tex','PosterHalfWidth','PosterMiddleRowHeight',6,
      'KEEP THE WALKS, TIGHTEN THE BOUND',r"""
    {\EquationSize\[
      |H|\le\mathcal M,\quad\varrho(\mathcal M)<1,
      \qquad
      \Psi_2(\mathcal M)=-\log\det(I-\mathcal M)-\operatorname{tr}\mathcal M
    \]}
    \noindent\begin{minipage}[c]{274mm}\centering
      \includegraphics[width=271mm,height=97mm,keepaspectratio]{generated/figures/walk_bounds.pdf}
    \end{minipage}\hfill%
    \begin{minipage}[c]{143mm}\raggedright
      {\PosterCondensed\bfseries\fontsize{62pt}{65pt}\selectfont\color{IASBlue}\WalkMin--\WalkMax$\times$\par}
      {\PosterParagraph smaller remainder bound than the scalar norm bound\par}\vspace{4mm}
      {\FigureSize Same parameter box; heterogeneous delays spanning 0--40 ms. All three points reported.\par}
    \end{minipage}\par
    \vspace{2mm}
    {\FigureSize\WalkClaim\ A replacement upper bound additionally needs a verified integral over the complete contour.}
""")

panel('07_retuning.tex','PosterHalfWidth','PosterLowerRowHeight',7,
      'NONLINEAR DELAYED VALIDATION IN JULIA',r"""
    {\FigureSize Full retained-device DAE reduction; true delay by method of steps, Rodas5P; 60 s, tolerance $10^{-9}$. All 39 bus phase frequencies.\par}
    \vspace{2mm}
    \includegraphics[width=419mm,height=106mm,keepaspectratio]{generated/figures/events.pdf}\par
    \vspace{2mm}
    {\MetricSize\bfseries\color{IASBlue}$|\Delta f|_{\max}=\FreqMax$ Hz
       \qquad $|\mathrm{RoCoF}|_{\max}=\RocofMax$ Hz/s\par}
    \vspace{2mm}
    {\FigureSize Dashed limits: 0.5 Hz and 0.5 Hz/s. Both metrics use a 0.5 s phase window. Frozen impedance-load events at the indicated nominal MW.}
""")

panel('08_safe_paths.tex','PosterHalfWidth','PosterLowerRowHeight',8,
      'A CONTINUOUS REGIONAL CERTIFICATE',r"""
    \noindent\begin{minipage}[c]{251mm}\centering
      \includegraphics[width=248mm,height=127mm,keepaspectratio]{generated/figures/contour.pdf}
    \end{minipage}\hfill%
    \begin{minipage}[c]{163mm}\raggedright
      {\Subhead\color{IASDarkGreen}Every point in the box\par}\vspace{4mm}
      {\EquationSize $|\delta\rho_i|\le0.001$\par}
      {\EquationSize $|\delta K_i|/K_i^0\le1\%$\par}
      {\FigureSize for both $K_p,K_i$; center $\rho_i^0=0.875$.\par}\vspace{4mm}
      {\Subhead\color{IASBlue}Exactly one root\par}
      {\EquationSize $\Re\lambda<-0.08229$\par}
      {\FigureSize inside the contour for every parameter vector in the box.\par}\vspace{3mm}
      {\FigureSize \PanelCount\ interval panels; 128-bit complex ball arithmetic; envelope $q<0.98794235<1$.\par}
    \end{minipage}\par
    \vspace{4mm}
    {\FigureSize\RegionalClaim\ The certificate encloses the whole contour and parameter box at exact $\tau=1/25$ s; exported binary64 coefficients are treated as exact data.}
""")

panel('09_evidence_scope.tex','PosterGridWidth','PosterBottomRowHeight',9,
      'WHAT IS ESTABLISHED -- AND WHAT REMAINS OPEN',r"""
    \noindent\begin{minipage}[c]{210mm}\centering
      {\MetricSize\bfseries\color{IASBlue}\VoltageMin--\VoltageMax\ pu\par}
      {\FigureSize voltage range; limits [0.9, 1.1]\par}
    \end{minipage}\hfill%
    \begin{minipage}[c]{193mm}\centering
      {\MetricSize\bfseries\color{IASBlue}\SlackMin\par}
      {\FigureSize minimum SG slack; limit 0.002\par}
    \end{minipage}\hfill%
    \begin{minipage}[c]{213mm}\centering
      {\MetricSize\bfseries\color{IASDarkGreen}Feasibility shown\par}
      {\FigureSize all-root margin checked numerically\par}
    \end{minipage}\hfill%
    \begin{minipage}[c]{227mm}\centering
      {\MetricSize\bfseries\color{IASGold!70!black}Maximum still open\par}
      {\FigureSize no certified replacement ceiling\par}
    \end{minipage}\par
    \vspace{3mm}
    {\FigureSize \textbf{Same-$\rho$ comparison:} fixed gains also pass all five events. Retuning slightly lowers RoCoF but increases frequency excursion; superiority is not established. No converter current-limit claim.}
""")

footer=(HERE/'sections/footer.tex').read_text(encoding="utf-8-sig")
footer=footer.replace('SAFE TRANSITIONS REQUIRE SAFE PREFIXES.','LOCAL DELAY. TWO OPPOSITE DAMPING EFFECTS.')
footer=footer.replace('Check every portfolio prefix against the minimal blockers of its control policy.',
 'Network paths determine the sign; collective bounds quantify the interaction.')
start=footer.index(r'      {[1]')
end=footer.index('\n\n',start)
footer=footer[:start]+r"""      {[1] Markovi\'c et al., IEEE TPWRS, 2021.\\{}
       [2] Bindel \& Hood, SIAM J. Matrix Anal., 2013.\\{}
       [3] Beyn, Linear Algebra Appl., 2012.};
    \node[anchor=north west,inner sep=0pt,outer sep=0pt,
      font={\PosterFooterReferenceType},text=IASWhite,align=left] at (410,24)
      {[4] Johansson, IEEE Trans. Comput., 2017.\\{}
       Theory, proofs, full gain vectors and data: QR.};"""+footer[end:]
write(HERE/'sections/footer.tex',footer)
compilepath=HERE/'compile_section.ps1'
cs=(HERE/'reference_template/compile_section_original.ps1').read_text(encoding="utf-8-sig").replace("'../../..'","'../../../..'").replace('IAS2026_Poster_Pulido','IAS2026_Collective_Damping_20261003')
cs=unicodedata.normalize('NFKD',cs).encode('ascii','ignore').decode('ascii')
write(compilepath,cs)
write(HERE/'SOURCE_MANIFEST.json',json.dumps({'inputs':sources,'qr_target':url,
      'geometry':'36 x 48 inches; inherited from existing user poster, not a claim about conference requirements',
      'copied_template':'reports/poster/ias2026; original modules and main preserved',
      'scientific_scope':'exact conditional identities + pointwise checks + one regional cluster + nonlinear feasible construction'},
      indent=2)+'\n')
print('Poster content and figures generated:',HERE)
