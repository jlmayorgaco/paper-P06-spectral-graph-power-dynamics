"""Create reports, figures, input snapshot, manifest and a portable research zip."""
from analyze_ports import *
from scipy.optimize import least_squares
import platform, shutil, tomllib, zipfile, datetime

def markdown_table(df,cols=None):
    if cols: df=df[cols]
    header='| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'
    return header+'\n'.join('| '+' | '.join(f'{v:.8g}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'

def fit_growth():
    summary=pd.read_csv(OUT/'TABLE_07_NONLINEAR_DDE.csv');rows=[]
    for row in summary.itertuples():
        filename=f'{row.case}_a{row.amplitude}_{row.tag}.csv'
        # Julia writes 1.0e-5 instead of Python's 1e-05.
        candidates=list((OUT/'nonlinear').glob(f'{row.case}_a*_{row.tag}.csv'))
        path=next(p for p in candidates if np.isclose(float(p.stem.split('_a')[-1].split('_')[0]),row.amplitude,rtol=1e-10,atol=0))
        a=pd.read_csv(path);t=a.time.to_numpy();y=a.modal_real.to_numpy()+1j*a.modal_imag.to_numpy()
        omega=2*np.pi*row.predicted_frequency_hz
        def residual(theta):
            alpha,om=theta
            basis=np.column_stack([np.exp((alpha+1j*om)*t),np.exp((alpha-1j*om)*t),np.ones(len(t))])
            coef=la.lstsq(basis,y)[0];r=(basis@coef-y)/max(la.norm(y),1e-30)
            return np.r_[r.real,r.imag]
        if row.tau_ms>0:
            fit=least_squares(residual,[row.predicted_alpha,omega],xtol=1e-12,ftol=1e-12,gtol=1e-12,
                bounds=([row.predicted_alpha-2,omega-2],[row.predicted_alpha+2,omega+2]))
            fitted=fit.x[0];freq=fit.x[1]/(2*np.pi);err=la.norm(fit.fun)
        else: fitted=freq=err=np.nan # <0.05 cycle in the short zero-delay record.
        rows.append(dict(case=row.case,amplitude=row.amplitude,tag=row.tag,predicted_alpha=row.predicted_alpha,
            fitted_alpha=fitted,fitted_frequency_hz=freq,fit_relative_residual=err,
            relative_trajectory_error=row.relative_trajectory_error))
    df=pd.DataFrame(rows);df.to_csv(OUT/'TABLE_09_NONLINEAR_FITS.csv',index=False)
    fig,axes=plt.subplots(2,2,figsize=(10,6),layout='constrained')
    for ax,(key,tau) in zip(axes.flat,[('N',40),('T',40),('N',45),('T',45)]):
        case=f'{key}_{tau}ms';p=next(p for p in (OUT/'nonlinear').glob(f'{case}_a*_standard.csv') if np.isclose(float(p.stem.split('_a')[-1].split('_')[0]),1e-5,atol=0))
        d=pd.read_csv(p)
        ax.plot(d.time,d.reference_real,'k--',lw=1,label='Exact linear mode')
        ax.plot(d.time,d.modal_real,color='#007f88',lw=1,label='Nonlinear pure-delay model')
        row=df[(df['case']==case)&np.isclose(df.amplitude,1e-5)&(df.tag=='standard')].iloc[0]
        ax.set_title(f'{key}: {tau} ms; growth {row.predicted_alpha:+.4f} / {row.fitted_alpha:+.4f} s$^{{-1}}$')
        ax.set_xlabel('Time (s)');ax.set_ylabel('Modal projection');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Local nonlinear DDE validation — PLL-angle history amplitude 10$^{-5}$ rad')
    fig.savefig(OUT/'FIG_03_NONLINEAR_DDE.png',dpi=180);fig.savefig(OUT/'FIG_03_NONLINEAR_DDE.pdf');plt.close(fig)
    return df

def extra_figures(roots,action):
    r=roots[(roots.design=='T')&(roots.tau_ms==45)].sort_values('real').iloc[-1]
    m=Model('T');s=complex(r.real,r.imag);Z=m.impedance(1j*s.imag,.045);H=(Z+Z.conj().T)/2
    pd.DataFrame(H.real,index=range(30,40),columns=range(30,40)).to_csv(OUT/'MATRIX_T45_DH_REAL.csv')
    pd.DataFrame(H.imag,index=range(30,40),columns=range(30,40)).to_csv(OUT/'MATRIX_T45_DH_IMAG.csv')
    fig,ax=plt.subplots(figsize=(6,5),layout='constrained');v=np.max(abs(H.real))
    im=ax.imshow(H.real,cmap='RdBu_r',vmin=-v,vmax=v);ax.set_xticks(range(10),range(30,40));ax.set_yticks(range(10),range(30,40))
    ax.set_xlabel('Retained SG port bus');ax.set_ylabel('Retained SG port bus')
    ax.set_title('T, 45 ms: real part of physical Hermitian damping\nFrequency-specific operator; not a Laplacian')
    fig.colorbar(im,ax=ax,label='Coefficient in physical torque/speed normalization')
    fig.savefig(OUT/'FIG_04_PHYSICAL_OPERATOR.png',dpi=180);fig.savefig(OUT/'FIG_04_PHYSICAL_OPERATOR.pdf');plt.close(fig)
    a=action[action.parameter=='tau'];fig,ax=plt.subplots(figsize=(7,4),layout='constrained');x=np.arange(2)
    ax.bar(x-.2,a.nodal_real/1000,.2,label='Diagonal action')
    ax.bar(x,a.collective_real/1000,.2,label='Collective action')
    ax.bar(x+.2,a.full_real/1000,.2,label='Full derivative')
    ax.set_xticks(x,a.design);ax.axhline(0,color='k',lw=.8);ax.set_ylabel('d Re(pole) / d delay [(1/s)/ms]')
    ax.set_title('Actual pole sensitivity at 40 ms — exploratory attribution');ax.legend(fontsize=8)
    fig.savefig(OUT/'FIG_05_POLE_ACTION.png',dpi=180);fig.savefig(OUT/'FIG_05_POLE_ACTION.pdf');plt.close(fig)

def provenance():
    files=['Project.toml','Manifest.toml','src/bnd_design/AnalyticSG.jl','src/bnd_design/AnalyticGFLPLL.jl',
        'src/bnd_design_e/CollectiveModel.jl','src/bnd_model_expN/PDExactDesignN.jl',
        'experiments/nonlinear_codesign_20261001/ReducedDAE.jl',
        'experiments/delay_dressed_replacement_frontier_20261002/DelayCharacteristic.jl',
        'experiments/latency_robust_pll_codesign_20261003/designs/N_nominal.toml',
        'experiments/latency_robust_pll_codesign_20261003/designs/Z_zero_delay_tuned.toml',
        'experiments/latency_robust_pll_codesign_20261003/L0_REPRODUCTION.csv',
        'experiments/latency_robust_pll_closure_20261003/designs/full20_step15_medium.toml',
        'experiments/latency_robust_pll_closure_20261003/evaluations/full20_step15_medium/ROOTS.csv',
        'experiments/optimal_latency_margin_20261002/T08_FAST_MODE_PROVENANCE.csv',
        'reports/experiment_A/matrices/bus33_baseline_equilibrium.csv',
        'reports/experiment_N/TABLE_N01_original_operating_point.csv',
        'reports/experiment_N/TABLE_N01_original_operating_point.csv.sha256']
    files += [str(p.relative_to(ROOT)).replace('\\','/') for p in (ROOT/'reports/experiment_D/inputs').glob('*.csv')]
    snapshot=OUT/'source_snapshot';hashes={}
    for rel in files:
        src=ROOT/rel;dst=snapshot/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
        hashes[rel]=hashlib.sha256(src.read_bytes()).hexdigest()
    git=lambda *a:subprocess.check_output(['git',*a],cwd=ROOT,text=True,encoding='utf-8',errors='replace').strip()
    pkg=tomllib.loads((ROOT/'Manifest.toml').read_text())
    versions={k:v[0].get('version','stdlib') for k,v in pkg.get('deps',{}).items() if k in ['PowerDynamics','NetworkDynamics','ForwardDiff','SciMLBase','OrdinaryDiffEqRosenbrock','CSV','DataFrames']}
    import scipy
    manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),git_HEAD=git('rev-parse','HEAD'),
        git_dirty_status=git('status','--short'),source_hashes=hashes,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        julia_manifest_version=pkg.get('julia_version'),julia_packages=versions,
        protocol_sha256=hashlib.sha256((OUT/'PROTOCOL.md').read_bytes()).hexdigest(),
        experiment_scope='physical mechanism and local nonlinear DDE validation; fixed rho and frozen gains',
        frozen_external_results_modified=False,commits_or_pushes=False)
    (OUT/'EXPERIMENT_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    (OUT/'MODEL_PROVENANCE.md').write_text(f'''# Model provenance

Git HEAD: `{manifest['git_HEAD']}`. Repository already dirty; full status and
exact source hashes are in EXPERIMENT_MANIFEST.json. This experiment writes only
inside its own directory. Prior frozen designs/results are copied, not modified.

The model uses the repository's exact algebraic network elimination and full
SG/AVR/governor and GFL/DC/current/PLL differential equations. Physical-supply
DC convention. Pure delay is inserted only at the PLL voltage-error detector.
The original PLL frequency low-pass state remains, independently of pure delay.
No physical line dynamics or current limiter absent from the source model is
added or claimed. The ten retained SG ratings define physical torque-speed ports.

Python: {platform.python_version()}; NumPy {np.__version__}; SciPy {scipy.__version__}.
Julia manifest version: {pkg.get('julia_version')}. Package versions: {versions}.

`source_snapshot/` includes the exact minimal repository dependencies and inputs.
The downloadable zip lays those out at the repository-relative paths required
by the code. Project/Manifest files and original source paths are retained.
''',encoding='utf-8')
    return files

def reports(fit):
    parity=pd.read_csv(OUT/'TABLE_01_PARITY.csv');checks=pd.read_csv(OUT/'TABLE_02_PORT_CHECKS.csv')
    roots=pd.read_csv(OUT/'TABLE_03_TRACKED_MODES.csv');hopf=pd.read_csv(OUT/'TABLE_04_HOPF.csv')
    action=pd.read_csv(OUT/'TABLE_08_ACTION_ATTRIBUTION.csv');nl=pd.read_csv(OUT/'TABLE_07_NONLINEAR_DDE.csv')
    extra_figures(roots,action)
    critical=roots.sort_values('real').groupby(['design','tau_ms']).tail(1).sort_values(['design','tau_ms'])
    critical.to_csv(OUT/'TABLE_10_TRACKED_ENVELOPE.csv',index=False)
    first=hopf.dropna(subset=['tau_ms']).sort_values('tau_ms').groupby('design').head(1).sort_values('design')
    t45=critical[(critical.design=='T')&(critical.tau_ms==45)].iloc[0]
    act=action[action.parameter=='tau'].copy();act['slope_per_ms']=act.full_real/1000
    report=f'''# Experimento ejecutado: amortiguamiento colectivo con puertos físicos

**Resultado: sí existe una contribución colectiva física y material en este
modelo; su signo aislado no proporciona una ley suficiente de estabilidad PLL
ni una nueva solución de reemplazo óptimo.**

Se mantuvo el reemplazo en **88.4551408%**, {parity.GFL_MW.iloc[0]:.6f} MW GFL y
{parity.SG_MW.iloc[0]:.6f} MW SG. No se optimizó rho, Kp o Ki en esta campaña.
N = ganancias nominales; Z = ajuste anterior sin retraso; T = ajuste anterior
de margen de latencia. Sus vectores completos están en models/*/design.toml.

## Resultado físico — EXACT_IDENTITY + NUMERICALLY_VALIDATED

Se derivó la impedancia exacta de par/velocidad Z(s) al eliminar estados ocultos
del modelo completo a frecuencia finita. Su parte hermítica D_H describe trabajo
mecánico incremental para un movimiento armónico prescrito. Se separó la diagonal
nodal del término entre nodos sin llamarlos laplacianos ni matrices constantes.
Las hipótesis y unidades están en THEORY.md.

En T, 45 ms, sobre el movimiento del modo PLL seguido:
- contribución nodal: **{t45.mode_damping_nodal:.6f}**;
- contribución colectiva: **{t45.mode_damping_cross:.6f}**;
- suma: **{t45.mode_damping_total:.6f}**.

Son coeficientes de amortiguamiento en la normalización de velocidad fijada
(norma del fasor de velocidad = 1), no MW de generación reemplazable. Aquí
omitir el término colectivo cambia el signo. Su fracción absoluta es
{100*t45.cross_absolute_fraction:.3f}%. No todos los diseños/frecuencias muestran
el mismo efecto; se incluyen todos los casos del protocolo.

## Polos del modelo con exponenciales exactas — SUPPORTED_LOCAL

Máximo real entre las ramas PLL rastreadas (no todo el espectro infinito):

{markdown_table(critical[critical.tau_ms.isin([40,45])],['design','tau_ms','real','frequency_hz'])}

Primer cruce Re(s)=0 entre esas ramas:

{markdown_table(first,['design','tau_ms','frequency_hz','mode_damping_nodal','mode_damping_cross'])}

Estos cruces de cero son distintos de los límites anteriores de seguridad
Re(s)=-0.05. El barrido no certifica ausencia de todas las demás raíces DDE.
En particular, los modos lentos pueden ser más cercanos al eje que la familia
PLL cuando esa familia está bien amortiguada.

## Validación no lineal con retraso puro — SUPPORTED_LOCAL

Se ejecutó el método de pasos con Rodas5P, usando historial de un modo exacto,
dos amplitudes (1e-5 y 1e-4 rad PLL), seis casos N/T a 0, 40 y 45 ms, y 3 s.
La referencia lineal es la solución analítica Re(a v exp(lambda t)) del DDE.

{markdown_table(fit[(fit.amplitude==1e-5)&(fit.tag=='standard')],['case','predicted_alpha','fitted_alpha','relative_trajectory_error'])}

El ajuste de crecimiento se omite a 0 ms: el modo de ~0.015 Hz no completa un
ciclo en 3 s; la trayectoria sí se compara. La tabla completa mantiene ambas
amplitudes y cualquier refinamiento requerido. Máximo error relativo de
trayectoria entre ejecuciones estándar: {nl[nl.tag=='standard'].relative_trajectory_error.max():.6g}.
Si el error crece con amplitud y no con refinamiento, es efecto no lineal.
**No son los cinco eventos de +/-100 MW, ni una demostración de seguridad robusta.**

## Acción colectiva sobre el polo — seguimiento exploratorio

Se derivó la sensibilidad del polo con Z_p y Z_s, conservando las partes reactiva
y disipativa y los modos izquierdo/derecho. Se separó la acción diagonal y
colectiva, validada contra la característica completa y diferencias centradas.

{markdown_table(act,['design','slope_per_ms','collective_fraction'])}

Para 42 derivadas (retraso uniforme y 20 ganancias por diseño), máximo error
relativo frente a diferencias finitas: {action.fd_relative.max():.6g}; máximo
desacuerdo entre formulación completa y de puertos: {action.port_full_relative.max():.6g}.
La fracción colectiva anterior es la fracción absoluta de las dos contribuciones
a la sensibilidad, no MW/ms de reemplazo. Todos los componentes están en TABLE_08.
Esta parte fue añadida tras el hallazgo del límite del criterio de signo y está
marcada exploratoria. La identidad de sensibilidad es conocida; no se vende
como un nuevo teorema de optimización.

## Resultado negativo principal — NEGATIVE_RESULT

El signo del amortiguamiento armónico proyectado **no clasifica universalmente
la estabilidad**: N a 40 ms tiene un polo seguido con parte real positiva
y amortiguamiento proyectado positivo. La dinámica reactiva, los estados ocultos
y la geometría de los modos importan. Por tanto queda rechazada la sustitución
del problema PLL completo por una condición escalar de signo de D_H.

La reducción anterior en s=0 sigue bloqueada para su partición original; este
experimento no la repara retroactivamente. La presente reducción es dependiente
de frecuencia y solo válida donde el bloque eliminado es invertible. El máximo
condicionamiento observado en los bloques ocultos de las ramas es
{roots.hidden_condition.max():.4g}: se reporta y no se presenta como certificación
en aritmética validada. Máximo residuo de raíz completo {roots.residual.max():.3g};
máximo residuo de Schur {roots.schur_root_residual.max():.3g}. Máximo error de
identidad transferencia/impedancia {checks.impedance_transfer_error.max():.3g}.

## Qué falta para la contribución del póster

1. Un diseño predictivo que use la interacción colectiva y supere un comparador
   con las mismas variables, restricciones, tolerancias y presupuesto de cálculo.
2. Volver a optimizar rho junto con Kp/Ki y ejecutar los cinco eventos congelados
   con retraso no lineal; el mecanismo por sí mismo no levanta el límite de reserva.
3. Verificar heterogeneidad espacial y cualquier compresión GSP en datos reservados.
4. Una cota superior válida antes de decir máximo global o casi óptimo.

**Decisión: resultado de mecanismo útil; aún no justifica un póster que anuncie
máximo reemplazo óptimo mediante Beyond Nodal Damping.** El parentesco con
literatura de impedancias/par complejo y control basado en grafos se documenta
en LITERATURE_POSITION.md; no se afirma novedad por renombrar esas herramientas.
'''
    (OUT/'REPORT_ES.md').write_text(report,encoding='utf-8')
    (OUT/'STATUS.md').write_text(f'''# Status

EXPERIMENT: COMPLETED_MECHANISM_PILOT
BASELINE_PARITY: NUMERICALLY_VALIDATED
FINITE_FREQUENCY_PHYSICAL_PORT_IDENTITY: EXACT_IDENTITY + NUMERICALLY_VALIDATED
COLLECTIVE_DAMPING_MATERIALITY: SUPPORTED_LOCAL
SCALAR_HARMONIC_DAMPING_STABILITY_RULE: NEGATIVE_RESULT
PURE_DELAY_NONLINEAR_MODE_VALIDATION: SUPPORTED_LOCAL
POLE_SENSITIVITY_ATTRIBUTION: NUMERICALLY_VALIDATED, exploratory
FRESH_COMPLETE_DDE_ROOT_COUNT: NOT_RUN
FULL_FROZEN_LARGE_EVENT_SET_WITH_DELAY: NOT_RUN
NEW_REPLACEMENT_OPTIMIZATION: NOT_RUN
GLOBAL_OPTIMALITY_BOUND: BLOCKED
GSP_COMPRESSION: NOT_TESTED
POSTER_READY_FOR_MAXIMUM_REPLACEMENT_CLAIM: NO

Missing gates: constrained replacement co-design, all five nonlinear delayed
events, appropriate full spectral coverage, fair method comparison and a valid
upper bound for any near-optimality claim. Local nonlinear mode tests do not
substitute for those gates. See REPORT_ES.md and the complete tables.
''',encoding='utf-8')
    (OUT/'POSTER_CLAIMS.md').write_text(f'''# Claim ledger

- EXACT_IDENTITY: In the stated physical torque-speed coordinates, an exact
  finite-frequency Schur impedance gives harmonic incremental mechanical work
  through its Hermitian part where the eliminated block is invertible.
- NUMERICALLY_VALIDATED: This identity agrees with the full transfer matrix on
  all 60 prescribed checks (maximum relative error {checks.impedance_transfer_error.max():.3g}).
- SUPPORTED_LOCAL: In frozen design T at 45 ms, the nodal projection is
  {t45.mode_damping_nodal:.4f} and the collective term {t45.mode_damping_cross:.4f};
  omitting the latter changes the projected damping sign.
- NEGATIVE_RESULT: Positive projected harmonic damping is compatible with an
  unstable PLL pole in the nominal design at 40 ms; this scalar is not a complete
  stability criterion.
- SUPPORTED_LOCAL: Actual nonlinear method-of-steps simulations test the six
  specified local histories, at two amplitudes, against exact linear DDE modes.
  Full results and amplitude dependence are in TABLE_07 and TABLE_09.
- NUMERICALLY_VALIDATED (exploratory): Physical-port pole sensitivities agree
  with full-characteristic sensitivities and centered finite differences; worst
  relative finite-difference error {action.fd_relative.max():.3g} for 42 derivatives.

No new penetration optimum, large-event delayed security, novel general theorem,
global guarantee, heterogeneous-delay result, or GSP compression is claimed.
''',encoding='utf-8')

def package(files):
    bundle=OUT/'physical_collective_damping_20261003_bundle.zip'
    relative_exp=OUT.relative_to(ROOT)
    with zipfile.ZipFile(bundle,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for rel in files: z.write(OUT/'source_snapshot'/rel,rel)
        for p in sorted(OUT.rglob('*')):
            if not p.is_file() or 'source_snapshot' in p.parts or '__pycache__' in p.parts or p.suffix=='.zip' or p.name.endswith('.zip.sha256'): continue
            z.write(p,str(relative_exp/p.relative_to(OUT)))
    with zipfile.ZipFile(bundle) as z:
        assert z.testzip() is None
        count=len(z.namelist())
    digest=hashlib.sha256(bundle.read_bytes()).hexdigest()
    bundle.with_suffix('.zip.sha256').write_text(digest+'  '+bundle.name+'\n')
    print('BUNDLE',bundle,'entries',count,'bytes',bundle.stat().st_size,'sha256',digest)

if __name__=='__main__':
    fit=fit_growth();reports(fit)
    if '--no-git' not in sys.argv:
        files=provenance();package(files)
