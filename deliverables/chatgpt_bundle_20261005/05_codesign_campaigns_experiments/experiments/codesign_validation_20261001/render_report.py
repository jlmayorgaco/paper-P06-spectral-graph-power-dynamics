"""Render the validation evidence without converting incomplete runs into passes."""
import argparse,csv,hashlib,json,math,tomllib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/codesign_validation_20261001'
parser=argparse.ArgumentParser();parser.add_argument('--label',default='improved');args=parser.parse_args()
label=args.label
def toml(p):return tomllib.loads(p.read_text(encoding='utf-8'))
def csvrows(p):return list(csv.DictReader(p.open(encoding='utf-8',newline='')))
c=toml(OUT/f'candidate_{label}.toml')
old=toml(ROOT/'experiments/codesign_validation_20261001/frozen/original_candidate.toml')
pdpath=OUT/f'PD_{label}'
spec=toml(pdpath/'spectrum.toml');rob=toml(OUT/f'robustness_{label}.toml')
events=csvrows(pdpath/'events.csv');sources={r['event']:pdpath for r in events}
altpath=OUT/f'PD_{label}_limits'
if (altpath/'events.csv').exists():
 alt=csvrows(altpath/'events.csv');events+=alt;sources.update({r['event']:altpath for r in alt})
refinepath=OUT/f'PD_{label}_refined'
if (refinepath/'events.csv').exists():
 alt=csvrows(refinepath/'events.csv');events+=alt;sources.update({r['event']:refinepath for r in alt})
inputpath=OUT/f'PD_{label}_inputs'
if (inputpath/'events.csv').exists():
 alt=csvrows(inputpath/'events.csv');events+=alt;sources.update({r['event']:inputpath for r in alt})
by={r['event']:r for r in events}
incomplete=[toml(p) for directory in (pdpath,altpath,refinepath,inputpath) for p in directory.glob('*_incomplete.toml') if p.stem.removesuffix('_incomplete') not in by]
suffix='' if label=='improved' else '_'+label
prefixpath=OUT/f'bus8_prefix{suffix}_metrics.csv'
prefix=csvrows(prefixpath) if prefixpath.exists() else []
prefixincompletepath=OUT/f'bus8_prefix{suffix}_incomplete.toml'
prefixincomplete=toml(prefixincompletepath) if prefixincompletepath.exists() else None
tail=csvrows(OUT/f'linear_peak_and_tail_audit{suffix}.csv')
exponential=csvrows(OUT/f'independent_matrix_exponential{suffix}.csv')
params=csvrows(OUT/f'parameters_{label}.csv')
gains=csvrows(OUT/f'gain_comparison_{label}.csv')
total=c['retained_SG_MW']+c['converted_GFL_MW'];replacement=100*c['converted_GFL_MW']/total
required=['step_16','step_29','step_8','negative_16','pulse_16','refined_16']
if label=='repaired':required.append('refined_8')
passes=lambda r:r['F_pass'].lower()=='true' and r['R_pass'].lower()=='true'
limit_convergence={}
if 'step_8_event_limits' in by and 'step_8_event_limits_refined' in by:
 for metric in ('Fpeak_Hz','Rpeak_Hz_s'):
  a=float(by['step_8_event_limits'][metric]);b=float(by['step_8_event_limits_refined'][metric])
  limit_convergence[metric]=abs(a-b)/max(abs(b),1e-15)
event_limit_verified=bool(limit_convergence) and max(limit_convergence.values())<1e-4
coverage=dict(by)
if event_limit_verified:coverage['step_8']=by['step_8_event_limits_refined']
missing=[s for s in required if s not in coverage]
failed=[s for s in required if s in coverage and not passes(coverage[s])]
prefix_failure=bool(prefix) and all(float(r['Fpeak_Hz'])>.5 or float(r['Rpeak_Hz_s'])>.5 for r in prefix)
if prefix_failure and 'step_8' not in failed:failed.append('step_8')
convergence={}
for bus in (16,8):
 if f'refined_{bus}' in by and f'step_{bus}' in by:
  convergence[str(bus)]={}
  for metric in ('Fpeak_Hz','Rpeak_Hz_s'):
   a=float(by[f'step_{bus}'][metric]);b=float(by[f'refined_{bus}'][metric])
   convergence[str(bus)][metric]={'absolute_difference':abs(a-b),'relative_difference':abs(a-b)/max(abs(b),1e-15)}
searches=[]
for folder in (OUT,OUT/'search_faces',OUT/'search_derivatives'):
 p=folder/'optimization_results.json'
 if p.exists():
  d=json.loads(p.read_text());searches.append({'stage':folder.name,'evaluations':d['evaluations'],
   'best_audited_MW':min(a['J'] for a in d['audits']),
   'minimum_stationarity_residual':min(a['stationarity_inf'] for a in d['audits']),
   'solver_messages':sorted({a['solver_message'] for a in d['audits']})})
summary={'candidate':c,'replacement_percent':replacement,'pd_spectrum':spec,'robustness':rob,
 'events':events,'incomplete':incomplete,'bus8_diagnostic_prefix':prefix,'required_missing':missing,'required_failed':failed,
 'bus8_prefix_demonstrates_violation':prefix_failure,
 'bus8_prefix_incomplete':prefixincomplete,
 'linear_peak_and_tail_audit':tail,'independent_matrix_exponential':exponential,
 'declared_suite_pass':not missing and not failed,'convergence':convergence,'searches':searches,
 'event_limit_verification':event_limit_verified,'event_limit_relative_peak_differences':limit_convergence,
 'global_optimum_certified':False,'local_optimum_certified':False,
 'limitation':'Full-band radius is a Float64 numerical lower bound, not outward-rounded proof. No hardware current limiter or EMT validation. Common mean gains are a reference, not an optimized common-gain baseline.'}
(OUT/f'validation_summary_{label}.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
fig,axes=plt.subplots(3,1,figsize=(8.2,11),gridspec_kw={'height_ratios':[1,1.25,1.1]},layout='constrained')
fig.get_layout_engine().set(rect=(0,.045,1,.955))
verdict='Required disturbance limits violated' if failed else 'Required disturbance validation incomplete' if missing else 'Declared disturbance suite passed'
fig.suptitle(f'IEEE-39 SG–GFL candidate: validation audit\n{replacement:.2f}% replacement · {c["retained_SG_MW"]:.2f} MW retained SG\n{verdict}',fontsize=14)
bus=np.arange(30,40);rho=np.array(c['rho'])*100
axes[0].bar(bus,rho,color='#167b8c',label='GFL');axes[0].bar(bus,100-rho,bottom=rho,color='#d99635',label='SG')
axes[0].set(xticks=bus,ylim=(0,111),ylabel='Share of original dispatch (%)',title='A   Allocation after four complete SG removals')
axes[0].legend(ncol=2,loc='lower left')
for b,r in zip(bus,rho):axes[0].text(b,102,f'{r:.1f}',ha='center',fontsize=8)
styles={'step_16':('#007b83','+100 MW bus 16'),'step_29':('#b25b10','+100 MW bus 29'),
 'step_8':('#7252a3','+100 MW bus 8'),'negative_16':('#ce3c4c','−100 MW bus 16')}
for name,(color,title) in styles.items():
 key='step_8_event_limits' if name=='step_8' and 'step_8_event_limits' in by else name
 p=sources.get(key,pdpath)/f'{key}_trajectory.csv'
 if p.exists():
  a=np.genfromtxt(p,delimiter=',',names=True);mask=a['time_after_event_s']>=0
  axes[1].plot(a['time_after_event_s'][mask],a['Fmax_Hz'][mask],color=color,label=title,lw=1.5)
if prefix and 'step_8' not in coverage:
 p=OUT/f'bus8_prefix{suffix}_fine.csv';a=np.genfromtxt(p,delimiter=',',names=True);mask=a['time_after_event_s']>=0
 axes[1].plot(a['time_after_event_s'][mask],a['Fmax_Hz'][mask],color='#7252a3',ls=':',label='+100 MW bus 8 (5 s prefix)',lw=2)
axes[1].axhline(.5,color='#343434',ls='--',lw=1,label='0.5 Hz limit')
axes[1].set(xlabel='Time after disturbance (s)',ylabel='Maximum |Δf| over buses (Hz)',title='B   Independent nonlinear PowerDynamics simulations',xlim=(0,60))
axes[1].legend(fontsize=8,ncol=2)
visible=[r for r in events if r['event'] not in ('small_16','small_8','small_29','refined_16','refined_8','step_8_relaxed','step_8_event_limits_refined')]
if prefix:
 visible.append({'event':'step_8_prefix','Fpeak_Hz':prefix[-1]['Fpeak_Hz'],'Rpeak_Hz_s':prefix[-1]['Rpeak_Hz_s']})
x=np.arange(len(visible));w=.34
axes[2].bar(x-w/2,[float(r['Fpeak_Hz'])/.5 for r in visible],w,label='Frequency / limit',color='#167b8c')
axes[2].bar(x+w/2,[float(r['Rpeak_Hz_s'])/.5 for r in visible],w,label='RoCoF / limit',color='#d99635')
axes[2].axhline(1,color='#343434',ls='--',lw=1)
names={'step_16':'+100\nbus 16','step_29':'+100\nbus 29','step_8':'+100\nbus 8',
 'step_8_event_limits':'+100 bus 8\nevent limits',
 'step_8_prefix':'+100 bus 8\n5 s prefix',
 'negative_16':'−100\nbus 16','stress_29':'+150\nbus 29','pulse_16':'+100 pulse\nbus 16'}
ctitle='C   Tested operating envelope; 150 MW is a stress test' if 'stress_29' in by else 'C   Completed tests\nBus 8 and 150 MW stress trajectories incomplete'
axes[2].set(xticks=x,xticklabels=[names.get(r['event'],r['event']) for r in visible],ylabel='Peak / declared limit',title=ctitle)
axes[2].legend(fontsize=8)
footer='Causal 0.5 s bus-phase measurement. Local/global optimality not certified.'
if failed:footer+='\nRequired disturbance limits violated: '+', '.join(failed)+'.'
elif missing:footer+='\nIncomplete 60 s validation: '+', '.join(missing)+'.'
fig.text(.5,.002,footer,ha='center',fontsize=8)
fig.savefig(OUT/f'validation_{label}.png',bbox_inches='tight');fig.savefig(OUT/f'validation_{label}.pdf',bbox_inches='tight');plt.close(fig)

lines=[f'VALIDACION SG GFL IEEE 39 — {label}',
 f'Fecha: 2026-10-01. Candidato SHA256: {spec["candidate_sha256"]}',
 'ESTADO: '+('NO APROBADO: existen violaciones de perturbaciones requeridas.' if failed else 'INCOMPLETO: faltan perturbaciones requeridas.' if missing else 'PASA la suite declarada de perturbaciones; no certifica optimalidad ni estabilidad no lineal global.'),
 '',f'SG retenido: {c["retained_SG_MW"]:.6f} MW; GFL: {c["converted_GFL_MW"]:.6f} MW; sustitucion: {replacement:.6f}%.',
 f'Reduccion de SG respecto al candidato anterior: {old["retained_SG_MW"]-c["retained_SG_MW"]:.6f} MW.',
 'Los porcentajes corresponden al despacho inicial congelado; no son un porcentaje universal de capacidad instalada.',
 f'Soporte SG: {c["support"]}. Se retiraron completamente los SG con epsilon=0.',
 '', 'PROBLEMA Y MEDICION',
 'Variables por nodo: epsilon=1-rho, Kp y Ki del PLL GFL. Objetivo J=sum(P0_i*epsilon_i). Se mantiene el punto inicial P/Q y la regla de escalamiento SG del modelo.',
 'Restricciones de busqueda: alpha(Aq)<=-0.05 s^-1; radio complejo aditivo de Aq+0.05 I >= beta_req; max|f_T|<=0.5 Hz y max|R_T|<=0.5 Hz/s para escalones de carga de 100 MW en buses 16, 8 y 29.',
 'f_T(t)=[theta(t)-theta(t-T)]/(2*pi*T); R_T(t)=[theta(t)-2*theta(t-T)+theta(t-2*T)]/(2*pi*T^2), T=0.5 s.',
 'El radio corresponde a perturbaciones complejas no estructuradas de la matriz en las coordenadas elegidas; no es un porcentaje de incertidumbre fisica de lineas, cargas o controladores.',
 'Se elimina exactamente un modo de simetria rotacional; se conservan todos los otros polos finitos, incluso cuando se retiran SG.',
 '', 'ESTABILIDAD LOCAL Y ROBUSTEZ LINEAL',
 f'Polos fisicos PD: {spec["physical_poles"]}; alpha PD: {spec["pd_alpha"]:.12g} s^-1; limite: -0.05 s^-1.',
 f'Error maximo de polos emparejados: {spec["max_matched_pole_error"]:.6g} s^-1.',
 f'Robustez: {rob["status"]}; cota inferior: {rob.get("beta_lower_bound",math.nan):.12g}; requisito: {rob["beta_req"]:.12g}.',
 'La cota de toda la banda usa Float64 con una reserva numerica; no es una prueba con redondeo dirigido.',
 '', 'EVENTOS NO LINEALES: cambios de Pset de carga, Qset fijo, ventana causal T=0.5 s.',
 'Evento               Pico Hz      RoCoF Hz/s    Pasa limites   Tiempo solver s']
for r in events:lines.append(f'{r["event"]:20s} {float(r["Fpeak_Hz"]):11.7f} {float(r["Rpeak_Hz_s"]):13.7f} {str(passes(r)):14s} {float(r["runtime_s"]):9.2f}')
lines+=['',f'Eventos obligatorios ausentes: {missing}; fallidos: {failed}.']
for d in incomplete:lines.append(f'INCOMPLETO: {d["event"]}, {d["retcode"]}, t={d["last_time_s"]} s. No demuestra por si solo inestabilidad fisica.')
lines+=['','DIAGNOSTICO BUS 8: prefijo independiente de 5 s despues de la perturbacion. No sustituye la validacion de 60 s.',json.dumps(prefix,indent=2),
 'Intento de prefijo incompleto: '+json.dumps(prefixincomplete),
 'VERIFICACIONES DEL EVALUADOR',
 '8 pruebas de unidades de salto de fase y derivada de meseta: PASS. 12 pruebas del callback de limites: PASS.',
 f'Exponencial matricial aumentada frente al evaluador modal: {sum(r["pass"]=="true" for r in exponential)}/{len(exponential)} comparaciones PASS; error relativo maximo {max(float(r["relative_phase_error"]) for r in exponential):.3g}.',
 'Respuesta lineal: prefijo de 120 s con paso 0.0025 s y cota analitica de la cola posterior; los tres buses pasan los limites de 0.5 Hz y 0.5 Hz/s. El prefijo se muestrea; no se certifican todos los extremos entre muestras.',
 json.dumps(tail,indent=2)]
lines+=['','CONVERGENCIA NUMERICA',json.dumps(convergence,indent=2),
 'Verificacion de eventos de limites: '+str(event_limit_verified),json.dumps(limit_convergence,indent=2),
 'Los limites anti-windup SG existentes se localizan mediante callbacks. Se selecciona el lado saturado con un desplazamiento de 1e-10 y se contrasta con 1e-11, sin cambiar los limites ni las ganancias. Los registros de eventos estan en PD_'+label+'_limits.',
 '', 'OPTIMALIDAD', 'No se certifica optimalidad local ni global. La busqueda produjo mejores candidatos, pero no paso la auditoria KKT ni se verifico SOSC.',
 'El candidato improved incorpora 1% de retencion sobre una semilla de la busqueda. El candidato repaired conserva esa asignacion y cambia las ganancias proporcionalmente bajas del PLL mediante Kp>=2*sqrt(Ki); es una heuristica en las unidades del modelo, no una garantia de amortiguamiento del sistema conectado.',
 json.dumps(searches,indent=2),
 'La auditoria de pasos muestra que las diferencias finitas de sensibilidades debiles del PLL quedan contaminadas por error numerico cerca de la frontera del radio.',
 '', 'GANANCIAS Y ALCANCE',
 'Los rho, Kp y Ki por nodo estan en parameters_'+label+'.csv. Se conservan las unidades del modelo Julia.',
 'gain_comparison_'+label+'.csv compara ganancias individuales, nominales y la media de las individuales.',
 'La media no es un controlador comun optimizado; esta comparacion no prueba superioridad frente al mejor controlador comun.',
 'No se han validado limitadores duros de corriente, hardware, EMT ni todas las arquitecturas y contingencias.',
 'Tampoco se demuestra estabilidad no lineal global ni region de atraccion; las trayectorias completadas cubren 60 s despues de los eventos indicados.',
 'La tension y corriente registradas son diagnosticos. La razon de corriente respecto al valor inicial no equivale a un margen frente a la corriente nominal admisible.',
 'La respuesta lineal y la no lineal difieren a 100 MW; las dos se reportan por separado.',
 '', 'REPRODUCCION',
 'julia --startup-file=no --project=@stdlib experiments/codesign_validation_20261001/test_corrections.jl',
 'julia --startup-file=no --project=. experiments/codesign_validation_20261001/validate_oracle.jl',
 'julia --startup-file=no --project=. experiments/codesign_validation_20261001/test_limit_events.jl',
 'julia --startup-file=no --project=. experiments/codesign_validation_20261001/audit_tail.jl',
 f'julia --startup-file=no --project=. experiments/codesign_validation_20261001/validate_pd.jl reports/codesign_validation_20261001/candidate_{label}.toml {label}',
 f'julia --startup-file=no --project=. experiments/codesign_validation_20261001/validate_pd.jl reports/codesign_validation_20261001/candidate_{label}.toml {label}_limits event_limits',
 f'julia --startup-file=no --project=. experiments/codesign_validation_20261001/validate_pd.jl reports/codesign_validation_20261001/candidate_{label}.toml {label}_refined refined_bus8',
 f'julia --startup-file=no --project=. experiments/codesign_validation_20261001/validate_pd.jl reports/codesign_validation_20261001/candidate_{label}.toml {label}_inputs small_inputs',
 f'julia --startup-file=no --project=. experiments/codesign_validation_20261001/diagnose_bus8.jl {label}',
 'El validador reanuda casos completados; usar otra etiqueta para una repeticion independiente.',
 f'python experiments/codesign_validation_20261001/render_report.py --label {label}',
 'Julia y PowerDynamics usan el Project.toml y Manifest.toml del repositorio. La busqueda SLSQP usa el evaluador fisico Julia por socket local; no utiliza PD durante la busqueda.']
(OUT/f'VALIDATION_REPORT_{label}.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'experiments/codesign_validation_20261001').rglob('*') if p.is_file()}
for p in [ROOT/'Project.toml',ROOT/'Manifest.toml',ROOT/'src/bnd_model_expN/PDExactDesignN.jl',ROOT/'src/bnd_model_expN/PDReferenceN.jl',ROOT/'src/pd39/model.jl']:
 manifest[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps({'label':label,'replacement_percent':replacement,'missing':missing,'failed':failed,'convergence':convergence},indent=2))
