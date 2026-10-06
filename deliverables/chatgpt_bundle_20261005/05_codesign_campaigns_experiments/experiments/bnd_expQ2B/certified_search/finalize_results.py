"""Assemble measured results; never optimize or edit the frozen candidate."""
from pathlib import Path
import csv, hashlib, json, subprocess, tomllib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'reports/experiment_Q2B/CERTIFIED_SEARCH'
def toml(name):
    return tomllib.loads((OUT/name).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

c=toml('Z_LOCAL_SECURE_FINAL.toml'); k=toml('LOCAL_CERTIFICATE.toml')
ks=toml('LOCAL_CERTIFICATE_FEASIBLE.toml'); br=toml('FULL_BAND_CERTIFICATE.toml')
bt=toml('ALL_TIME_CERTIFICATE.toml'); pdbr=toml('PD/RICCATI_FULL_BAND_CERTIFICATE.toml')
pdsp=toml('PD/SPECTRUM.toml'); events=pd.read_csv(OUT/'PD/EVENTS.csv')
event=events.loc[events.event=='bus16_100MW'].iloc[0]
oos_fail=events.loc[(events.event!='bus16_100MW')&(~events.frequency_pass|~events.rocof_pass)]
oos_text='Ninguno de los eventos adicionales registrados excedió estos dos límites.' if oos_fail.empty else '; '.join(f'{r.event}: F={r.Fpeak:.6f} Hz, R={r.Rpeak:.6f} Hz/s' for r in oos_fail.itertuples())
assert sha(OUT/'Z_LOCAL_SECURE_FINAL.toml')==(OUT/'Z_LOCAL_SECURE_FINAL.toml.sha256').read_text().strip()
assert k['local_KKT_certificate'] and ks['local_KKT_certificate'] and br['strict_negative_LMI'] and bt['pass']
assert pdbr['strict_negative_LMI'] and pdsp['margin_pass']
assert event.frequency_pass and event.rocof_pass

ap=pd.read_csv(OUT/'FROZEN_ANALYTICAL_POLES.csv')
pp=pd.read_csv(OUT/'PD/ALL_PD_POLES.csv')
av=ap.real.to_numpy()+1j*ap.imag.to_numpy();pv=pp.real.to_numpy()+1j*pp.imag.to_numpy()
pv=np.delete(pv,np.argmin(abs(pv)))
ii,jj=linear_sum_assignment(abs(pv[:,None]-av[None,:]))
errs=abs(pv[ii]-av[jj]);assert max(errs)<1e-5
pd.DataFrame(dict(PD_real=pv[ii].real,PD_imag=pv[ii].imag,AN_real=av[jj].real,
                 AN_imag=av[jj].imag,error=errs)).to_csv(OUT/'PD/BIJECTIVE_POLE_MATCHING.csv',index=False)
orig=pd.read_csv(ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv').sort_values('bus')
bus=orig.bus.to_numpy();weights=orig.P_gen_MW.to_numpy();e=np.array(c['epsilon'])
design=pd.DataFrame(dict(bus=bus,epsilon=e,rho=c['rho'],retained_MW=weights*e,Kp=c['Kp'],Ki=c['Ki']))
design.to_csv(OUT/'FINAL_DESIGN_BY_BUS.csv',index=False)
old=tomllib.loads((ROOT/'reports/experiment_Q2B/CORE/Z_Q2B_SECURE_T05_FINAL.toml').read_text())
gain=old['retained_SG_MW']-c['J_MW'];gap=c['J_MW'];gap_pp=100*gap/c['total_original_SG_MW']

ledger=[]
incmask=sum(1<<(b-30) for b in c['support'])
for mask in range(1024):
    ledger.append(dict(mask=mask,support=';'.join(str(b) for b in bus if mask&(1<<int(b-30))),
        rigorous_lower_MW=0.,status='LOCAL_KKT_POINT_FOUND_BRANCH_COMPLETENESS_OPEN' if mask==incmask else 'OPEN',
        globally_resolved=False))
pd.DataFrame(ledger).to_csv(OUT/'GLOBAL_SUPPORT_LEDGER.csv',index=False)
globality=dict(lower_MW=0.,upper_MW=gap,gap_MW=gap,
    replacement_fraction_gap_percentage_points=gap_pp,global_certified=False,
    lower_scope='Nonnegative retained dispatch; exact trivial bound',
    upper_scope='Numerically certified analytical feasibility plus independent primary nonlinear event',
    local_scope='Fixed-support numerical KKT/SOSC; strict one-sided insertion witnesses',
    unresolved_supports=1024,formal_interval_proof=False,
    stronger_Bernstein_bound='Attempted; inconclusive and NOT admitted as rigorous lower bound')
(OUT/'GLOBAL_GAP.json').write_text(json.dumps(globality,indent=2),encoding='utf-8')

# Check immutable model, dependencies, and prior candidate without touching them.
mf=json.loads((ROOT/'reports/experiment_N/MODEL_FREEZE.json').read_text())
checks=[dict(path=p,expected=h,actual=sha(ROOT/p),passed=sha(ROOT/p)==h) for p,h in mf['inputs'].items()]
assert all(r['passed'] for r in checks)
pd.DataFrame(checks).to_csv(OUT/'IMMUTABLE_INPUT_HASH_AUDIT.csv',index=False)
sw=tomllib.loads((ROOT/'reports/experiment_N/SOFTWARE_PROVENANCE.toml').read_text())
assert sha(ROOT/'Project.toml')==sw['project_sha256'] and sha(ROOT/'Manifest.toml')==sw['manifest_sha256']
assert sha(ROOT/'reports/experiment_Q2B/CORE/Z_Q2B_SECURE_T05_FINAL.toml')=='f53250383cdb906a2d0cb3130ca9c1e72c63b472717c5c82aea1d5838ac785ce'
provenance=dict(branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
    parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    model_sha=mf['MODEL_SHA'],candidate_sha=sha(OUT/'Z_LOCAL_SECURE_FINAL.toml'),
    frozen_model_hashes_passed=len(checks),project_unchanged=True,manifest_unchanged=True,
    old_candidate_unchanged=True,commit=False,push=False)
(OUT/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140})
figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
hist=pd.read_csv(OUT/'MULTIPEAK_SOC_HISTORY.csv')
fig,axs=plt.subplots(1,3,figsize=(12,3.6))
axs[0].plot(hist.iteration,hist.retained_MW,color='#126c8c');axs[0].set(ylabel='Retained SG [MW]',xlabel='SQP iteration',title='Joint co-design')
axs[1].semilogy(hist.iteration,np.maximum(hist.stationarity,1e-14),label='Stationarity')
axs[1].semilogy(hist.iteration,np.maximum(hist.primal,1e-14),label='Primal');axs[1].axhline(1e-6,color='grey',ls='--');axs[1].legend();axs[1].set(xlabel='SQP iteration',title='Normalized KKT residuals')
axs[2].bar(['Reduced direction 1','Reduced direction 2'],k['projected_Hessian_eigenvalues'],color='#18977a')
axs[2].set(ylabel='Projected Hessian eigenvalue',title='Second-order local check');fig.tight_layout();fig.savefig(figdir/'FIG_C01_LOCAL_CERTIFICATE.png');plt.close(fig)
fig,axs=plt.subplots(3,1,figsize=(9,7),sharex=True)
for ax,name,label in zip(axs,['retained_MW','Kp','Ki'],['Retained SG [MW]','PLL Kp','PLL Ki']):
    ax.bar(bus,design[name],color=['#126c8c' if z>0 else '#c8cfd2' for z in e]);ax.set(ylabel=label,xticks=bus)
axs[-1].set_xlabel('Generator bus');fig.suptitle('Frozen joint design: 312.523288 MW retained');fig.tight_layout();fig.savefig(figdir/'FIG_C02_DESIGN.png');plt.close(fig)
lin=pd.read_csv(OUT/'FROZEN_LINEAR_TRAJECTORIES.csv');nonlin=pd.read_csv(OUT/'PD/bus16_100MW_SENSORS.csv')
regular=[]
for t in (7.29445,14.62293):
    sl=lin.loc[abs(lin.event_time_s-t)<.035]
    co=np.polyfit(sl.event_time_s-t,np.abs(sl.F_bus34),2)
    regular.append(dict(bus=34,time_s=t,absolute_frequency_second_derivative=2*co[0],
        strict_local_time_maximum=2*co[0]<0,method='local quadratic cross-check of frozen exact trajectory'))
assert all(r['strict_local_time_maximum'] for r in regular)
pd.DataFrame(regular).to_csv(OUT/'ACTIVE_PEAK_REGULARITY.csv',index=False)
scaling=[]
for mw in (1,25,50,100):
    path=OUT/f'PD/bus16_{mw}MW_SENSORS.csv'
    if not path.exists():continue
    trace=pd.read_csv(path);sel=(trace.event_time_s>=.01)&(trace.event_time_s<=60)
    tt=trace.event_time_s[sel].to_numpy();errors=[];denom=[]
    for b in bus:
        pred=np.interp(tt,lin.event_time_s,lin[f'F_bus{b}'])*mw/100
        obs=trace.loc[sel,f'grid_Hz_bus{b}'].to_numpy();errors.append(np.max(abs(pred-obs)));denom.append(np.max(abs(pred)))
    scaling.append(dict(MW=mw,max_frequency_mismatch_Hz=max(errors),relative_frequency_mismatch=max(errors)/max(denom)))
pd.DataFrame(scaling).to_csv(OUT/'PD/LINEAR_NONLINEAR_SCALING.csv',index=False)
fig,axs=plt.subplots(2,1,figsize=(9,6),sharex=True)
for ax,lc,nc,limit,label in [(axs[0],'F_bus','grid_Hz_bus',.5,'Max local |frequency deviation| [Hz]'),(axs[1],'R_bus','rocof_Hz_s_bus',.5,'Max local |RoCoF| [Hz/s]')]:
    ax.plot(lin.event_time_s,np.abs(lin[[lc+str(b) for b in bus]]).max(axis=1),label='Exact linear model')
    ax.plot(nonlin.event_time_s,np.abs(nonlin[[nc+str(b) for b in bus]]).max(axis=1),label='Independent nonlinear PD',ls='--')
    ax.axhline(limit,color='#b44336',ls=':');ax.set(ylabel=label,xlim=(0,60));ax.legend()
axs[-1].set_xlabel('Time after +100 MW bus-16 step [s]');fig.tight_layout();fig.savefig(figdir/'FIG_C03_LINEAR_VS_PD.png');plt.close(fig)
authority=pd.read_csv(OUT/'MARGINAL_SECURITY_VALUES.csv')
fig,ax=plt.subplots(figsize=(8,3.8))
ax.bar(authority.bus.astype(str),authority.peak_frequency,label='Frequency peaks',color='#126c8c')
ax.bar(authority.bus.astype(str),authority.robust,bottom=authority.peak_frequency,label='Robustness',color='#d49a26')
ax.axhline(1,color='black',ls='--');ax.set(xlabel='Retained SG bus',ylabel='Multiplier-weighted marginal value',title='Interior retention law: weighted authority = 1');ax.legend();fig.tight_layout();fig.savefig(figdir/'FIG_C04_SECURITY_VALUE.png');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,2.9));ax.plot([0,gap],[0,0],lw=12,color='#d7e4e8');ax.scatter([0,gap],[0,0],s=80,c=['#777','#126c8c']);ax.set(xlim=(-10,gap+15),yticks=[],xlabel='Minimum retained SG [MW]',title='Global uncertainty remains open (not an estimated distance)');ax.text(0,.025,'Valid trivial lower bound\n0 MW',ha='left');ax.text(gap,.025,f'Feasible upper bound\n{gap:.6f} MW',ha='right');ax.set_ylim(-.03,.085);fig.tight_layout();fig.savefig(figdir/'FIG_C05_GLOBAL_GAP.png');plt.close(fig)

table='\n'.join(f'| {int(r.bus)} | {r.epsilon:.9f} | {r.rho:.9f} | {r.retained_MW:.6f} | {r.Kp:.9f} | {r.Ki:.9f} |' for r in design.itertuples())
report=f'''# Óptimo local numéricamente certificado — cierre de Q2B

## Resultado y alcance

**PASS_LOCAL_SECURE_NUMERICAL.** Se ejecutó co-diseño conjunto de las cinco retenciones presentes y las 20 ganancias PLL, con SQP explícito, múltiples picos activos y corrección de segundo orden. El modelo y el evento permanecen congelados. El óptimo es local para el **problema analítico de seguridad especificado**; la factibilidad no lineal se valida independientemente. Esto no demuestra optimalidad del problema con restricciones no lineales de trayectoria ni globalidad.

- Retención: **{gap:.9f} MW**; GFL: **{c['converted_GFL_MW']:.9f} MW**, **{100*c['GFL_fraction']:.6f}%** del despacho original.
- Soporte: **30, 33, 35, 37, 39**. Mejora frente al candidato congelado Q2B anterior: **{gain:.6f} MW** menos de SG.
- El punto KKT antes del resguardo de redondeo cuesta {k['retained_MW']:.12f} MW. El representante factible añade solamente {gap-k['retained_MW']:.9g} MW.
- No se modificó el candidato anterior ni el modelo ExpN; {len(checks)}/{len(checks)} hashes y ambos archivos de dependencias coinciden. Sin commit, push o cambio de rama.

## Certificado local

| Prueba | Punto KKT | Representante congelado |
|---|---:|---:|
| Residuo primal | {k['primal']:.3g} | {ks['primal']:.3g} |
| Estacionariedad | {k['stationarity']:.6g} | {ks['stationarity']:.6g} |
| Complementariedad | {k['complementarity']:.6g} | {ks['complementarity']:.6g} |
| Rango activo / restricciones | 23 / 23 | 23 / 23 |
| Condición del Jacobiano normalizado | {k['active_Jacobian_condition']:.6g} | {ks['active_Jacobian_condition']:.6g} |

Todos los multiplicadores activos son positivos. Hay dos direcciones tangentes libres. Los autovalores de la Hessiana reducida son {k['projected_Hessian_eigenvalues'][0]:.6f} y {k['projected_Hessian_eigenvalues'][1]:.6f}; la variación entre tres pasos de diferenciación es {k['Hessian_step_variation']:.6f}. La menor curvatura positiva supera ampliamente esa variación. Objetivo dividido por 1000; epsilon y ganancias normalizadas a sus intervalos. Estos residuos no equivalen a una distancia en MW al óptimo global.

Es un **certificado numérico a tolerancias declaradas**, no una prueba con intervalos del cero exacto de las ecuaciones KKT. El pequeño resguardo factible deja residuos no nulos documentados. Se verifican LICQ, signos, complementariedad y SOSC; KKT por sí solo no habría bastado.

Las restricciones activas son beta, dos picos de frecuencia del bus 34 (aproximadamente 7.29445 y 14.62293 s) y las 20 cotas de ganancias. La restricción modal está holgada; su multiplicador es cero. No se presenta un modo modal activo ficticio ni una explicación BND de una restricción que no limita este óptimo.

## Factibilidad analítica completa

| Magnitud | Resultado |
|---|---:|
| Alpha, 143 polos físicos | {c['alpha']:.12f} 1/s |
| Beta observada | {c['beta_observed']:.12g} |
| Beta mínima cubierta por certificado | {br['beta_lower_tested']:.12g} |
| Pico de frecuencia | {c['Fpeak']:.12f} Hz |
| Frecuencia estacionaria | {ks['Fsteady']:.12f} Hz |
| Pico de RoCoF | {c['Rpeak']:.12f} Hz/s |

La robustez se cubre con una desigualdad estricta bounded-real/Riccati a 256 bits: X positivo y residual negativo definido. Residuo de Riccati: {br['Riccati_residual_big']:.5g}; margen mínimo de la desigualdad: {br['negative_LMI_min_eigenvalue_Float64']:.5g}. Cubre toda frecuencia y la cola infinita. No depende de que una malla haya encontrado el máximo. Redondeo al más cercano: **NUMERICALLY_CERTIFIED**, no FORMALLY_CERTIFIED. El intento previo de cubrir la banda por subdivisiones se interrumpió por costo; no se usa como evidencia de certificación.

Para tiempo se combinaron raíces de las derivadas, cotas de Taylor sobre {bt['nodes']} intervalos y una cota modal de la cola desde {bt['tailtime']:.0f} s. La comprobación es numérica; el margen de redondeo se declara en el código. Las mallas densas son contraste y visualización.

## Validación independiente PowerDynamics

- Alpha PD: **{pdsp['alpha_PD']:.12f} 1/s**; error con diseño: {pdsp['alpha_error']:.4g} 1/s.
- Emparejamiento biyectivo de todos los polos: error máximo {max(errs):.4g} 1/s; error relativo de A en las mismas coordenadas físicas: {pdsp['matrix_relative_error']:.4g}.
- Trim: residuo {pdsp['trim_residual']:.4g}; error P {pdsp['max_P_error_pu']:.4g} pu; error Q {pdsp['max_Q_error_pu']:.4g} pu. Cargas 31/39 separadas e invariantes.
- Se obtuvo además el certificado Riccati directamente sobre la matriz reconstruida por PD.
- Evento no lineal bus16/100 MW: **frecuencia {event.Fpeak:.8f} Hz; RoCoF {event.Rpeak:.8f} Hz/s**. Ambos límites de 0.5 pasan. Frecuencia media final: {event.Fsteady:.8f} Hz.
- Error máximo de balance: {event.power_balance_error_MW:.4g} MW; cambio Pdc: {event.Pdc_change:.4g}. Valores y trazas por evento en `PD/EVENTS.csv` y `PD/*SENSORS.csv`.
- No se declara protección por límite de corriente GFL ni seguridad para perturbaciones arbitrarias. El evento cambia Pset del ZIP exactamente como en el proyecto; el consumo efectivo y pérdidas se registran durante el evento. Ventana de medida T=0.5 s, sin escogerla después de ver el resultado.

El contraste local de 1 MW tiene aproximadamente 0.2705% de discrepancia máxima relativa de frecuencia. A 100 MW la discrepancia alcanza 13.7542%: la respuesta no lineal pasa los límites, pero no se afirma identidad lineal para una perturbación grande. Véase `PD/LINEAR_NONLINEAR_SCALING.csv`.

**Fuera del caso de diseño:** {oos_text}. Los eventos de buses 8/29 son validación externa; sus fallos se conservan y limitan el alcance. No se incorporaron como restricciones ocultas ni se retocó el candidato. En el instante exacto del escalón se reconstruye el límite algebraico derecho manteniendo todos los estados diferenciales fijos: la interpolación de PD había combinado inicialmente parámetros posteriores con voltajes anteriores. La corrección afecta sólo la exportación de esa muestra, no la trayectoria integrada.

## Parámetros congelados

| Bus | epsilon | rho | SG MW | Kp | Ki |
|---|---:|---:|---:|---:|---:|
{table}

SHA-256 del candidato: `{provenance['candidate_sha']}`.

## Cambios de arquitectura y globalidad

Se construyeron límites exactos unilaterales para insertar SG en 31, 32, 34, 36 y 38. Los cinco violan estrictamente beta en omega=0; la evaluación del inverso fue a 256 bits. Las matrices con estados adicionales se proyectan mediante una coisometría sobre el sistema anterior, de modo que su norma resolvente no puede ser menor. Esto excluye también inserciones simultáneas infinitesimales por continuidad, bajo el modelo de incertidumbre congelado. Los residuos se conservan en `ADJACENT_EXACT_LIMITS.csv`. La afirmación es local y numérica; no excluye inserciones finitas ni ramas desconectadas.

Para el problema global, la única cota inferior admitida aquí es **L=0 MW**. La cota superior factible es **U={gap:.9f} MW**. Por tanto:

    0 <= J_global* <= {gap:.9f} MW
    0 <= J_candidato - J_global* <= {gap:.9f} MW

La mejora global todavía posible está acotada por **{gap_pp:.6f} puntos porcentuales del despacho original**. Esto es una cota máxima, no una estimación de cuánto falta. No se puede calcular un error porcentual respecto de J_global* usando una cota inferior cero.

El intento de cota mediante eliminación DC exacta y polinomios de Bernstein no dio una cota útil y carece de envolventes con redondeo exterior. No se usa el heurístico afín de 57.2 MW. La identidad estática se comprobó en 100 casos (error máximo 1.64e-11 Hz), pero eso no convierte sus coeficientes Float64 en una prueba global. Los 1024 soportes siguen abiertos en cuanto a completitud global; uno contiene este punto local certificado. Véase `GLOBAL_SUPPORT_LEDGER.csv`.

## Relevancia para el paper

La contribución que estos datos sostienen es un co-diseño físico de retención SG y PLL locales, con restricciones espectrales, robustas y temporales, evidencia de optimalidad local y validación independiente. SQP, la corrección de segundo orden y el lema bounded-real son herramientas conocidas; no se reclaman como teoremas nuevos. El resultado nuevo del caso es la asignación de seguridad y su comprobación reproducible. La novedad frente a toda la literatura requiere una comparación bibliográfica específica.

Los valores marginales ponderados son aproximadamente uno en los cinco buses interiores. Buses 30 y 39: domina la frecuencia; buses 33,35,37: domina beta. La explicación está en `MARGINAL_SECURITY_VALUES.csv`, con signo coherente con minimizar SG retenido.

## Tiempo, obstáculos y reproducción

La ejecución SQP registrada tardó {hist.elapsed_s.iloc[-1]:.3f} s sin compilación, en {len(hist)} iteraciones externas. La certificación robusta tardó {br['elapsed_s']:.3f} s y la temporal {bt['elapsed_s']:.3f} s. Los contadores completos de evaluaciones del SQP no se guardaron en esa ejecución; el código ahora los registra al reproducirla. No se inventa ese dato. Las primeras pruebas fallidas de gradientes, corrección Newton y cobertura por subdivisión se conservan para trazabilidad.

Se conserva la secuencia previa fallida en `LOCAL_CERTIFICATE`; este directorio la supera con evidencia nueva, sin reescribirla. Reproducción: véase `REPRODUCE.md`. Falta una cota global fuerte, completar ramas y una prueba formal con intervalos. Esas limitaciones no anulan el resultado local.

Referencias metodológicas: [Boyd et al., Linear Matrix Inequalities in System and Control Theory](https://web.stanford.edu/~boyd/lmibook/), [Boyd, Balakrishnan y Kabamba, norma H-infinity](https://stanford.edu/~boyd/papers/bisection_hinfty.html), [Gould y Robinson, SQP y corrección de segundo orden](https://optimization-online.org/2008/12/2192/).
'''
report=report.replace('**PASS_LOCAL_SECURE_NUMERICAL.**','**PASS_LOCAL_SECURE — certificación numérica.**')
(OUT/'RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md').write_text(report,encoding='utf-8')
summary=dict(status='PASS_LOCAL_SECURE',certificate_type='NUMERICAL_NOT_FORMAL',retained_MW=gap,GFL_MW=c['converted_GFL_MW'],
    GFL_fraction=c['GFL_fraction'],support=c['support'],primal=ks['primal'],stationarity=ks['stationarity'],
    complementarity=ks['complementarity'],LICQ='PASS',SOSC='PASS_NUMERICAL',beta_full_band='PASS_NUMERICAL_256BIT',
    PD_alpha=pdsp['alpha_PD'],PD_Fpeak=float(event.Fpeak),PD_Rpeak=float(event.Rpeak),
    PD_primary_event='PASS',global_lower_MW=0.,global_upper_MW=gap,global_gap_MW=gap,
    out_of_sample_failures=oos_fail.event.tolist(),
    global_certified=False,formal_interval_certificate=False,push=False,commit=False)
(OUT/'TERMINAL_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
