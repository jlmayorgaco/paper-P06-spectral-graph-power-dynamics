"""Generate reports, figures and manuscript values from executed tables."""
from moments import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shutil

def main():
    read=lambda f:pd.read_csv(OUT/f)
    w=read('TABLE_02_NETWORK_WEIGHTS.csv');ident=read('TABLE_03_MOMENT_IDENTITIES.csv')
    finite=read('TABLE_04_FINITE_FREQUENCY_CHECK.csv');orig=read('TABLE_06_CODESIGN.csv')
    final=read('globalized/TABLE_13_GLOBALIZED.csv');vec=read('globalized/TABLE_14_DESIGNS.csv')
    spec=read('TABLE_11_FULL_SPECTRUM.csv');events=read('TABLE_12_NONLINEAR_EVENTS.csv')
    parity=read('TABLE_10_JULIA_PARITY.csv');cert=json.loads((OUT/'MOMENT_ASSUMPTION_CERTIFICATE.json').read_text())
    enclosures=read('TABLE_17_DESIGN_MOMENT_ENCLOSURES.csv')
    assert enclosures[enclosures.design.str.endswith('_collective')].neutral_within_1e_12.all()
    assert len(events)==5 and len(spec)==4
    coll=final[(final.pattern=='uniform_1ms')&(final.method=='collective')].iloc[0]
    modal=final[(final.pattern=='uniform_1ms')&(final.method=='modal_only')].iloc[0]
    moment=final[(final.pattern=='uniform_1ms')&(final.method=='moment_only')].iloc[0]
    unchanged=orig[(orig.pattern=='uniform_1ms')&(orig.method=='unchanged')].iloc[0]
    sitewise=orig[(orig.pattern=='uniform_1ms')&(orig.method=='sitewise')].iloc[0]
    chosen=vec[(vec.pattern=='uniform_1ms')&(vec.method=='collective')].sort_values('bus')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':220,'axes.titleweight':'bold'})
    def save(fig,name):
        fig.tight_layout();fig.savefig(OUT/(name+'.png'),bbox_inches='tight');fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,3.5));ax.bar(w.bus,w.weight,color=['#b45309' if z<0 else '#007c91' for z in w.weight])
    ax.axhline(0,color='black',lw=.8);ax.set(xlabel='PLL bus',ylabel='Network moment weight $w_i$',xticks=w.bus,
        title='One network, signed control weights');save(fig,'FIG_01_NETWORK_WEIGHTS')
    fig,ax=plt.subplots(figsize=(7,3.6))
    for method,color,marker in [('unchanged','#7c7c7c','o'),('sitewise','#b45309','s'),('collective','#007c91','D')]:
        df=final[final.method==method] if method=='collective' else orig[orig.method==method]
        values=[-.06]+[float(df[df.pattern==f'uniform_{t:g}ms'].critical_real.iloc[0]) for t in [.1,.5,1]]
        ax.plot([0,.1,.5,1],values,marker=marker,color=color,label=method.replace('_',' '))
    ax.axhline(-.05,color='black',ls='--',label='Required margin');ax.axhline(0,color='gray',lw=.7)
    ax.set(xlabel='Additional uniform measurement delay (ms)',ylabel='Rightmost catalog real part (1/s)',
           title='Moment preservation alone does not preserve damping');ax.legend(ncol=2,fontsize=8)
    save(fig,'FIG_02_DELAY_AND_DAMPING')
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    dk=chosen.Kp.to_numpy()-P[10:20]
    axes[0].bar(w.bus-.17,P[20:]*.001,width=.34,label='Sitewise',color='#b45309')
    axes[0].bar(w.bus+.17,dk,width=.34,label='Collective',color='#007c91');axes[0].axhline(0,color='black',lw=.7)
    axes[0].set(xlabel='PLL bus',ylabel=r'$\Delta K_p$',xticks=w.bus,title='Same +1ms delay, different tuning');axes[0].legend(fontsize=8)
    methods=['moment only','sitewise','modal only','collective'];objects=[moment,sitewise,modal,coll]
    axes[1].bar(methods,[r.critical_real for r in objects],color=['#a0a0a0','#b45309','#627d98','#007c91'])
    axes[1].axhline(-.05,color='black',ls='--');axes[1].set(ylabel='Catalog real part (1/s)',title='Ablation at 41ms')
    axes[1].tick_params(axis='x',rotation=18);save(fig,'FIG_03_TUNING_AND_ABLATION')
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    for method,color in [('unchanged','#7c7c7c'),('sitewise','#b45309'),('collective','#007c91')]:
        df=final[final.method==method] if method=='collective' else orig[orig.method==method]
        vals=[df[df.pattern==f'site{i}_1ms'].critical_real.iloc[0] for i in range(30,40)]
        axes[0].plot(range(30,40),vals,'o-',color=color,label=method)
    axes[0].axhline(-.05,color='black',ls='--');axes[0].set(xlabel='Site receiving +1ms',ylabel='Catalog real part (1/s)',title='Identical delay multiset');axes[0].legend(fontsize=8)
    for n,d in finite.groupby('case'):
        axes[1].loglog(d.s,d.relative_asymptotic_error,lw=.8,alpha=.6)
    axes[1].set(xlabel='Real Laplace argument $s$ (1/s)',ylabel='Leading-coefficient relative error',title='Asymptotic law, not a broadband surrogate')
    save(fig,'FIG_04_PLACEMENT_AND_VALIDITY')
    fig,axes=plt.subplots(1,2,figsize=(9,3.5));xx=np.arange(5);labels=[f'{r.bus}: {r.delta:+g}' for r in events.itertuples()]
    old=pd.read_csv(OLD/'TABLE_12_NONLINEAR_EVENTS.csv')
    for ax,col,title in [(axes[0],'F','Frequency excursion (Hz)'),(axes[1],'R','Windowed RoCoF (Hz/s)')]:
        ax.bar(xx-.18,old[col],width=.36,label='40ms previous design',color='#627d98')
        ax.bar(xx+.18,events[col],width=.36,label='41ms collective',color='#007c91')
        ax.axhline(.5,color='black',ls='--');ax.set(xticks=xx,xticklabels=labels,ylabel=title,xlabel='Bus: load parameter (MW)',title=title)
    axes[0].legend(fontsize=8);save(fig,'FIG_05_NONLINEAR_EVENTS')
    extra=100*(coll.relative_Kp_norm**2/modal.relative_Kp_norm**2-1)
    status={'experiment_status':'THEORY_AND_TARGETED_VALIDATION_COMPLETE','claim_status':'NUMERICALLY_VALIDATED',
        'GFL_MW':float(PORTS.P0@P[:10]),'GFL_percent':float(100*PORTS.P0@P[:10]/PORTS.P0.sum()),
        'retained_SG_MW':float(PORTS.P0@(1-P[:10])),'replacement_optimized':False,
        'globalized_status':coll.status,'catalog_alpha':float(coll.critical_real),
        'full_spectrum_status':spec[spec.design=='uniform_1ms_collective'].status.iloc[0],
        'nonlinear_five_events_pass':bool(events['pass'].all()),'F_max_Hz':float(events.F.max()),'R_max_Hz_s':float(events.R.max()),
        'worst_frequency_bus':int(events.iloc[events.F.argmax()].bus),'worst_frequency_delta':float(events.iloc[events.F.argmax()].delta),
        'moment_constraint_cost_percent_vs_local_modal_only':float(extra),
        'max_partial_catalog_alpha_collective':float(final[final.method=='collective'].critical_real.max()),
        'journal_ready':False,'global_optimality':'NOT_ESTABLISHED','graph_communication_implemented':False,
        'certificate_scope':cert['scope'],'interval_weight_signs':cert['all_weight_signs_enclosed']}
    dump('STATUS.json',status)
    report=f'''# Leyes de red para el ajuste de PLL — resultado del experimento

## Resultado central

Para este modelo SG/GFL, manteniendo planta, rho y Ki, los cambios finitos de
Kp y retardo cumplen exactamente una ley entre coeficientes:

    H_b(s)-H_a(s) = -s² kappa H0 + O(s³)
    kappa = sum_i w_i (delta_tau_i/Ki_i - delta_Kp_i/Ki_i²)

Los pesos w salen de la singularidad rotacional del retorno de red completo;
no se eligieron como centralidades. Tienen signos distintos. Con estabilidad,
el área firmada de la diferencia de frecuencias es cero y su primer momento
temporal es kappa*f_infinity. No es una cota de nadir, RoCoF ni seguridad.
Cambiar Ki mueve todas las áreas normalizadas en una misma dirección:
sus diferencias entre buses son invariantes dentro de esta arquitectura.

La derivación para control matricial muestra además que
delta_KP=KI*S*KI, S*1=0, preserva H0,H1,H2. Un polinomio de Laplaciano sin
término constante cumple esa condición. Es un resultado teórico condicionado;
ese controlador con comunicación no fue implementado ni comparado.

## Qué se ejecutó en IEEE-39

204 estados, diez PLL filtrados, retardos exponenciales puros en el detector,
red AC con pérdidas,39 salidas, tres entradas linealizadas. Las ecuaciones
se comprobaron con recurrencia de Laurent, transferencia completa204D y
derivadas por diferencias finitas. Los signos de los pesos y las hipótesis de
regularidad se encerraron con aritmética de bolas192 bits para la exportación
binaria con la restauración rotacional explícita. No cubre incertidumbre física.
Los14 diseños colectivos redondeados satisfacen |kappa|<1e-12 mediante encierro
intervalar; no se atribuye un cero literalmente exacto a los gains numéricos.

Conservamos **{status['GFL_MW']:.6f} MW GFL ({status['GFL_percent']:.6f}%)** y
**{status['retained_SG_MW']:.6f} MW SG**. No se optimizó ni aumentó ese porcentaje.
Solo cambió Kp; Ki y rho están congelados en los valores previos.

Al pasar de40 a41ms:

| Método | Parte real crítica del catálogo (1/s) | Norma de cambio relativo Kp |
|---|---:|---:|
| Sin ajuste | {unchanged.critical_real:.9f} | {unchanged.relative_Kp_norm:.7f} |
| Cancelación PLL por PLL | {sitewise.critical_real:.9f} | {sitewise.relative_Kp_norm:.7f} |
| Solo momento, mínimo cuadrático | {moment.critical_real:.9f} | {moment.relative_Kp_norm:.7f} |
| Solo margen modal | {modal.critical_real:.9f} | {modal.relative_Kp_norm:.7f} |
| Momento colectivo + margen modal | {coll.critical_real:.9f} | {coll.relative_Kp_norm:.7f} |

El último estado del algoritmo es **{coll.status}**: no equivale a optimalidad
global ni, si agotó iteraciones, a convergencia local. El coste cuadrático frente
al otro candidato local modal cambia {extra:.3f}%; es una comparación de dos
resultados del algoritmo, no el precio óptimo probado de la restricción.
Además, el candidato colectivo no alcanza el objetivo auxiliar de-0.0600001/s,
mientras que el modal directo sí. Ambos se evalúan con el límite original-0.05/s;
sus costes no permiten inferir superioridad a igual margen modal alcanzado.

Espectro completo: **{status['full_spectrum_status']}**, conteo numérico por
contorno con cota de región y exponenciales exactas; no certificado intervalar.
Cinco eventos no lineales de60s: **{int(events['pass'].sum())}/5 pasan**.
Máximo |delta f|={events.F.max():.9f}Hz; RoCoF={events.R.max():.9f}Hz/s;
tensión mínima={events.Vmin.min():.9f}, máxima={events.Vmax.max():.9f};
holgura SG mínima={events.slack.min():.9f}. Ventana heredada0.5s.
Los eventos son cambios de carga Z parametrizados en MW, no potencia constante
exacta durante toda la trayectoria. No hay evidencia de seguridad de limitador
de corriente de convertidor: no está modelado.

## Resultados negativos que se conservan

- Igualar momentos de baja frecuencia no conserva el amortiguamiento modal:
  la regla por PLL hace inestable el caso41ms aunque kappa sea cero.
- Un solo modo protegido falló; el catálogo completo es necesario aquí.
- El SQP sin globalización falló a41ms; sus resultados permanecen sin cambios.
- La ley es asintótica. A s=.0025/s, el error uniforme de primer término sigue
  siendo {100*finite[(finite.case=='uniform_1ms')&(finite.s==.0025)].relative_asymptotic_error.iloc[0]:.2f}%.
- Las pantallas SDP y la cota causal del experimento hermano no proporcionaron
  una barrera continua de todos los gains ni un piso SG útil.

## Novedad y publicación

La identidad específica y su restricción de autoridad son candidatas a una
aportación útil. Residuos, momentos, cancelación de consenso y QP son clásicos.
No se encontró una coincidencia exacta en la búsqueda acotada, pero eso no
prueba prioridad. Las revisiones independientes están en REVIEW_*.md.

**No se declara listo para journal ni ganador de póster.** Falta demostrar un
beneficio operativo de preservar este momento frente al ajuste modal directo,
transferir el resultado a otro punto/red, y desarrollar/validar el controlador
de grafo si se quiere vender esa arquitectura. La máxima sustitución SG→GFL,
su cota superior y la ventaja global frente a otros métodos siguen abiertas.

## Reproducción y entregables

THEORY_MOMENT_SECTION.tex contiene pruebas; el documento THEORY.tex abierto
integra la sección y los datos generados. TABLE_09 y globalized/TABLE_14 guardan
rho,Kp,Ki,tau completos. TABLE_11 contiene el espectro independiente y TABLE_12
los cinco eventos. Los protocolos, fallos, fuentes, versiones, hashes y revisión
se preservan junto al código. REPRODUCE.md da los comandos.
'''
    for a,b in {'pérdidas,39':'pérdidas, 39','bolas192':'bolas de 192','Los14':'Los 14',
                'de40 a41ms':'de 40 a 41 ms','de60s':'de 60 s','heredada0.5s':'heredada de 0.5 s',
                'objetivo auxiliar de-':'objetivo auxiliar de -','límite original-':'límite original -',
                'a41ms':'a 41 ms','de primer término sigue':'del primer término sigue'}.items():report=report.replace(a,b)
    (OUT/'REPORT_ES.md').write_text(report,encoding='utf8')
    claims=f'''# Claim--evidence ledger

- EXACT_IDENTITY: Theorem network moment hierarchy and collective scalar equality, under its stated full-linearized-model assumptions. Evidence: THEORY_MOMENT_SECTION.tex, independently reviewed proof.
- EXACT_IDENTITY: Spatial normalized signed-area differences are invariant under PLL-only retuning at fixed plant; stability and nonzero DC response required. Not a nonlinear peak bound.
- EXACT_IDENTITY: KI*S*KI proportional communication increments with S*1=0 preserve H0,H1,H2. Conditional architecture, unimplemented; no graph-controller superiority claim.
- NUMERICALLY_VALIDATED: Coefficient checks over39 outputs and3 inputs; maximum cubic relative error {ident[ident.power==2].coefficient_relative_error.max():.3e}; quadratic {ident[ident.power==1].coefficient_relative_error.max():.3e}. These are not independent network replicates.
- NUMERICALLY_VALIDATED: Exported-model hidden regularity/gauge derivative and all signed weights enclosed with Arb; scope is explicitly symmetry-restored binary matrices, not physical uncertainty.
- NUMERICALLY_VALIDATED: All14 rounded collective designs have interval-enclosed |kappa|<1e-12; literal zero is not claimed for numerical gains.
- NEGATIVE_RESULT: Sitewise moment cancellation has critical real part {sitewise.critical_real:.9f}/s at41ms; moment matching does not guarantee stability.
- NUMERICALLY_VALIDATED: At unchanged {status['GFL_percent']:.6f}% replacement, collective41ms candidate has catalog alpha {coll.critical_real:.9f}/s, spectrum status {status['full_spectrum_status']}, and {int(events['pass'].sum())}/5 frozen nonlinear events passing.
- SUPPORTED_LOCAL: Globalized numerical candidate, status {coll.status}; no global or certified near-optimal gain claim.
- BLOCKED: Maximum replacement, global SG floor, robust operating-region safety, publication priority and graph-controller advantage are not established.
'''
    (OUT/'POSTER_CLAIMS.md').write_text(claims,encoding='utf8')
    sci=lambda x:(lambda pair:pair[0]+r'\times10^{'+str(int(pair[1]))+'}')(f'{x:.2e}'.split('e'))
    generated=rf'''
\status{{NUMERICALLY\_VALIDATED}}
The Julia input derivative export agrees with centered differences to
${sci(read('TABLE_01_INPUT_PARITY.csv').input_relative_error.max())}$ relative.
The maximum coefficient-check errors are
${sci(ident[ident.power==2].coefficient_relative_error.max())}$ (cubic PLL change)
and ${sci(ident[ident.power==1].coefficient_relative_error.max())}$ (integral-gain change).
Arb arithmetic encloses the hidden inverse, a nonzero nine-dimensional gauge
minor, a nonzero residue denominator, and all ten weight signs for the explicitly
gauge-restored binary export. This is not a physical uncertainty certificate.
The fourteen rounded collective designs have interval-enclosed
$|\kappa_3|<10^{{-12}}$; literal exact neutrality is not attributed to rounded gains.

At the fixed replacement of {status['GFL_percent']:.6f}\%, or
{status['GFL_MW']:.6f} MW, changing uniform PLL delay from 40 to 41 ms gives:
\begin{{center}}
\begin{{tabular}}{{lrr}}
\toprule Method & Catalog $\alpha$ (s$^{{-1}}$) & $\|\Delta K_p/K_p\|_2$\\
\midrule
Unchanged & {unchanged.critical_real:.6f} & {unchanged.relative_Kp_norm:.6f}\\
Sitewise compensation & {sitewise.critical_real:.6f} & {sitewise.relative_Kp_norm:.6f}\\
Moment only & {moment.critical_real:.6f} & {moment.relative_Kp_norm:.6f}\\
Modal only & {modal.critical_real:.6f} & {modal.relative_Kp_norm:.6f}\\
Collective moment and modes & {coll.critical_real:.6f} & {coll.relative_Kp_norm:.6f}\\
\bottomrule
\end{{tabular}}
\end{{center}}
The collective candidate's solver status is
\texttt{{{coll.status.replace('_',r'\_')}}}; no global optimum is claimed.
It passes the original $-0.05$ guard but misses the auxiliary $-0.0600001$
target. The modal-only candidate reaches that tighter target but its line
search stalls. Their costs do not establish superiority at an equal achieved
modal margin; neither run proves a local optimum.
The complete-region numerical DDE contour gives
\texttt{{{status['full_spectrum_status'].replace('_',r'\_')}}}.
All {int(events['pass'].sum())} of 5 frozen 60 s nonlinear events pass; the
worst excursion is {events.F.max():.6f}Hz and windowed RoCoF is
{events.R.max():.6f}Hz/s. Current-limiter safety is not modeled.
The delay-perturbed candidate retains {status['retained_SG_MW']:.6f} MW of SG;
this experiment does not optimize replacement.

\begin{{figure}}[htbp]\centering
\includegraphics[width=.84\linewidth]{{experiments/network_moment_codesign_20261004/FIG_01_NETWORK_WEIGHTS.pdf}}
\caption{{Network moment weights for the fixed IEEE--39 design. Their signs
are interval-verified for the symmetry-restored exported model.}}\end{{figure}}
\begin{{figure}}[htbp]\centering
\includegraphics[width=.88\linewidth]{{experiments/network_moment_codesign_20261004/FIG_02_DELAY_AND_DAMPING.pdf}}
\caption{{Exact-exponential pole calculations: equal moments do not protect
modal damping. Lines connect the declared samples, not a certified frontier.}}\end{{figure}}
\begin{{figure}}[htbp]\centering
\includegraphics[width=\linewidth]{{experiments/network_moment_codesign_20261004/FIG_03_TUNING_AND_ABLATION.pdf}}
\caption{{Finite tuning increments and algorithm ablations at 41 ms. Gain effort
is compared between local numerical candidates, not global optima.}}\end{{figure}}
\begin{{figure}}[htbp]\centering
\includegraphics[width=\linewidth]{{experiments/network_moment_codesign_20261004/FIG_05_NONLINEAR_EVENTS.pdf}}
\caption{{Five frozen nonlinear events. Frequency and RoCoF use the inherited
0.5 s window over all 39 bus voltage phases.}}\end{{figure}}
'''
    (OUT/'GENERATED_RESULTS.tex').write_text(generated,encoding='utf8')
    theory=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
    backup=OUT/'THEORY_before_moment.tex'
    if not backup.exists():shutil.copyfile(theory,backup)
    source=backup.read_text(encoding='utf8')
    section=(OUT/'THEORY_MOMENT_SECTION.tex').read_text().replace('% GENERATED_MOMENT_RESULTS_PLACEHOLDER',generated)
    source=source.replace('\\begin{thebibliography}',section+'\n\\begin{thebibliography}',1)
    bib=r'''
\bibitem{malladamoments} F. Paganini and E. Mallada,
\emph{Global Analysis of Synchronization Performance for Power Systems:
Bridging the Theory-Practice Gap}, IEEE Transactions on Automatic Control,
65(7):3007--3022,2020. \url{https://doi.org/10.1109/TAC.2019.2942536}.
\bibitem{aalipourmoments} F. Aalipour and T. Das,
\emph{Shaping the Transient Response of Nonlinear Systems to Satisfy a Class
of Integral Constraints}, Advances in Control Applications,4(3):e110,2022.
\url{https://arxiv.org/abs/2012.12493}.
'''
    source=source.replace('\\end{thebibliography}',bib+'\n\\end{thebibliography}')
    source=source.replace('Collective Interaction Bounds for Delay-Aware SG--GFL Co-Design','Collective PLL Interactions and Network Frequency-Response Constraints')
    source=source.replace('therefore quantitative: within a declared replacement and gain region,\nhow much replacement can any admissible local PLL tuning support,\nat fixed physical delays and under a stated security contract?',
        'which features of network frequency response finite PLL retuning can preserve\nwhile controlling collective modes. Section~\\ref{sec:networkmoments} derives\nexact constraints on this authority and tests a targeted delay-compensation\ndesign. The ultimate maximum-replacement problem remains open.')
    source=source.replace('The methodological distinction examined here is a \\emph{necessary\nfinite regional bound across all allowed local gains}, rather than\nonly a successful tuning or a fixed-pattern assignment obstruction.',
        'A separate unresolved route, preserved below, examines a \\emph{necessary\nfinite regional bound across all allowed local gains}. The present executed\nclosure concerns finite interactions and network response moments; it does\nnot claim that the regional-bound route has delivered a replacement ceiling.')
    added=(' A full-network Laurent analysis additionally separates the authority of '
        'integral gains, proportional gains and delays over frequency-response moments. '
        'It yields an exact network-weighted compensation equality and a conditional '
        'moment-neutral graph-control direction. IEEE--39 tests distinguish moment '
        'matching from modal stability, with independent delayed nonlinear validation. '
        'This remains a research manuscript; engineering superiority and global '
        'replacement optimality are not established.\n')
    source=source.replace('\\end{abstract}',added+'\\end{abstract}',1)
    theory.write_bytes(source.replace('\r\n','\n').replace('\n','\r\n').encode('utf8'))
    print(json.dumps(status,indent=2))

if __name__=='__main__':main()
