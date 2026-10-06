"""Assemble a Spanish, source-linked research record from completed audits."""
import argparse,csv,hashlib,json,platform,tomllib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'reports/nonlinear_codesign_20261001'
p=argparse.ArgumentParser();p.add_argument('candidate_json');p.add_argument('candidate_toml');a=p.parse_args()
jpath=Path(a.candidate_json);tpath=Path(a.candidate_toml)
d=json.loads(jpath.read_text());candidate=tomllib.loads(tpath.read_text());sha=hashlib.sha256(tpath.read_bytes()).hexdigest()
def readcsv(path):return list(csv.DictReader(path.open()))
pd=readcsv(BASE/'PD_physical_final/events.csv');long=readcsv(BASE/'long_horizon_final/events.csv')
eq=readcsv(BASE/'equilibrium_final.csv')
tail=readcsv(BASE/'weak_tail_final/events.csv')
for path in [BASE/'PD_physical_final/provenance.toml',BASE/'long_horizon_final/provenance.toml',BASE/'equilibrium_final_provenance.toml']:
    if tomllib.loads(path.read_text())['candidate_sha256']!=sha:raise RuntimeError(f'Candidate provenance mismatch: {path}')
assert len(pd)==len(long)==len(eq)==6
assert all(r['complete']=='true' for r in pd+long)
assert all(r['F_pass']=='true' and r['R_pass']=='true' and r['V_pass']=='true' for r in pd)
for r in long:
    assert float(r['Fpeak_Hz'])<=.5 and float(r['Rpeak_Hz_s'])<=.5
    assert float(r['Vmin'])>=.9 and float(r['Vmax'])<=1.1 and float(r['limiter_fraction'])>=.002
assert all(float(r['physical_alpha'])<0 and float(r['residual'])<1e-6 for r in eq)
assert len(tail)==1 and tail[0]['complete']=='true' and float(tail[0]['last_time'])>=10000.
assert float(tail[0]['Fpeak_Hz'])<=.5 and float(tail[0]['Rpeak_Hz_s'])<=.5
assert tomllib.loads((BASE/'weak_tail_final/provenance.toml').read_text())['candidate_sha256']==sha
kkt=json.loads((BASE/'stationarity_final/stationarity_report.json').read_text())
uniform=readcsv(BASE/'uniform_baseline/summary.csv')
u=max(float(r['rho']) for r in uniform if r['feasible']=='true')
cross=readcsv(BASE/'exchange_checkpoint/finite_exchange_results.csv')
fixed=next(float(r['actual_signature_change_norm']) for r in cross if r['id']=='31_32_1.0_fixed')
tuned=next(float(r['actual_signature_change_norm']) for r in cross if r['id']=='31_32_1.0_tuned')
fmax=max(float(r['Fpeak_Hz']) for r in pd+long);rmax=max(float(r['Rpeak_Hz_s']) for r in pd+long)
vmin=min(float(r['Vmin']) for r in pd+long);vmax=max(float(r['Vmax']) for r in pd+long)
worst_alpha=max(float(r['physical_alpha']) for r in eq)
summary={**d,'candidate_toml_sha256':sha,'PD_completed_cases':len(pd),'long_horizon_s':600.,
    'weakest_case_supplementary_horizon_s':10000.,
    'observed_frequency_max_Hz':fmax,'observed_RoCoF_max_Hz_s':rmax,'observed_voltage_range_pu':[vmin,vmax],
    'post_disturbance_max_real_pole':worst_alpha,'uniform_nominal_gain_feasible_percent':100*u,
    'additional_replacement_percentage_points':d['replacement_percent']-100*u,
    'exchange_31_32_1MW_signature_improvement_factor':fixed/tuned,
    'global_optimality_certified':False,'continuous_disturbance_robustness_certified':False,
    'KKT_audit':kkt}
(BASE/'validated_result.json').write_text(json.dumps(summary,indent=2))
text=f'''RESULTADO DE LA INVESTIGACION — 1 OCTUBRE 2026

Teoria general y experimento particular
La red, puntos de operacion, incertidumbres, perturbaciones y limites son entradas
de la formulacion. IEEE 39, los seis escalones nominales de carga Z de +/-100 MW,
la ventana de 0.5 s y los limites 0.5 Hz / 0.5 Hz/s son una instancia experimental.
No son constantes de la teoria. Las amplitudes representan cambios de consigna
de una carga Z; su potencia efectiva depende de la tension.

Resultado numerico validado
Sustitucion ponderada por despacho: {d['replacement_percent']:.8f} %.
Generacion sincrona conservada: {d['retained_MW']:.8f} MW.
Los seis escenarios se completaron en PowerDynamics compilado durante 60 s
posteriores al evento y en el modelo exacto reducido durante 600 s; se observaron
los 39 nodos. Maximos observados: {fmax:.9f} Hz y {rmax:.9f} Hz/s.
Rango de tension combinado: [{vmin:.9f}, {vmax:.9f}] pu.
Los seis equilibrios posteriores tienen polos fisicos con parte real negativa;
el menos amortiguado tiene parte real {worst_alpha:.9f} 1/s. Esto no certifica
una region de atraccion ni todas las perturbaciones posibles.
Ese margen posterior es muy pequeno: el punto es una frontera nominal de
investigacion, no una recomendacion de capacidad robusta para operar una red.
El caso mas debil (29, -100 MW) tambien se completo a 10000 s, con salida cada
0.1 s, sin nuevos picos mayores. Es una comprobacion suplementaria; los picos
tempranos se verifican por separado con salida cada 0.005/0.01 s.

La referencia uniforme, con ganancias PLL nominales, alcanzo {100*u:.6f} % en
el barrido/refinamiento local de 60 s: mejora de {d['replacement_percent']-100*u:.6f}
puntos porcentuales al co-disenar fracciones y ganancias. La referencia no es un
maximo escalar global demostrado. No atribuir toda la mejora exclusivamente al
analisis modal: falta la ablacion del optimizador sin esa guia.

Estado de optimalidad
El estado correcto es mejor candidato factible encontrado. No se ha calculado
una cota superior no trivial de sustitucion ni un certificado global. El informe
stationarity_final/stationarity_report.json contiene el residuo KKT, tolerancias
y complementariedad calculados con gradientes nuevos en este candidato.
'''
for r in kkt['KKT_audits']:
    text+=f"Tolerancia activa {r['active_tolerance']}: residuo KKT relativo {r['relative_stationarity_inf']:.8g}; violacion primal {r['primal_violation']:.8g}.\n"
text+=f'''
Metodo matematico
Se minimiza sum_i P0_i(1-rho_i), sujeto a las trayectorias DAE no lineales y a
supremos de restricciones sobre la familia declarada de perturbaciones. Se
calculan sensibilidades mediante eliminacion algebraica exacta y derivadas de
las trayectorias; un QP acotado propone cambios conjuntos en rho, log Kp, log Ki.
Los autovalores y sus sensibilidades guian el paso. La aceptacion exige simulacion
no lineal completa. El caso con pico tardio activo amplio el horizonte del gradiente.
La variante con correccion de factibilidad sigue comprobando cada propuesta en
la simulacion. El tratamiento por varios modos cercanos esta disponible; los
polos muy mal condicionados requieren un analisis de subespacio espectral.

El resultado analitico central es v_rho_i=-g_v^(-1)g_rho_i. Para el modelo afín
en tension, g_rho_i=S_i(i_GFL_i-i_SG_i). El adjunto de red pondera ese desajuste
local segun la perturbacion y la restriccion activa. Esto permite derivar la
compensacion de ganancias y una geometria local de intercambio entre nodos.
En la prueba de 1 MW entre 31/32, el reajuste redujo {fixed/tuned:.4f} veces el
cambio observado en el vector seleccionado de restricciones. La compensacion
acotada probada entre 35/36 no recupero el margen modal declarado; no demuestra
que cualquier otro reajuste no lineal sea imposible. Esas pruebas usan el
candidato intermedio del 90.50 %, identificado por sus parametros propios.

Auditoria fisica y numerica
Se corrigio, en una variante aislada, la convencion del balance de energia DC:
C*vdc*vdcdot=Pdc-Pac, junto con el signo correspondiente del regulador DC.
No se modificaron los archivos compartidos del modelo ni la biblioteca instalada.
La igualdad de polos nominales no detectaba esa diferencia no lineal.
La auditoria de gradientes, balance energetico, cierre de red y marco movil exacto
esta en los archivos TOML del directorio. El contraejemplo previo de 93.2213 %
sigue violando RoCoF en el modelo fisico corregido: aproximadamente 0.89145 Hz/s.

Limites del resultado y validacion que falta
1. Cerrar estacionariedad con tolerancia declarada, multiples inicializaciones y
   un limite superior valido si se quiere afirmar maximo global.
2. Buscar adversarialmente amplitudes, ubicaciones, topologias e incertidumbres;
   seis casos no certifican una familia continua.
3. Incorporar limites reales de corriente, energia DC, fuente primaria y retardos.
   Los cocientes de corriente respecto al valor inicial son diagnosticos.
4. Certificar una region terminal/invariante y cubrir cambios de regimen de
   saturacion. El optimizador actual trabaja con margen de no saturacion SG.
5. Probar ramas con retirada fisica completa de dispositivos. La busqueda actual
   usa fracciones interiores de agregados co-localizados; no simula la maniobra
   fisica de retirar una maquina en operacion.
6. Completar ablaciones y una comparacion entre simuladores con controladores
   equivalentes. ANDES REGCP1/REECA1 no se considera equivalente automaticamente.

Propuesta de poster
Beyond Nodal Damping: How Much Synchronous Generation Can We Safely Replace?
La historia combina un contraejemplo del diseno lineal, una sensibilidad de
sustitucion acoplada por la red, un mapa de SG conservados/ganancias y la envolvente
no lineal validada. El poster existente de hipergrafos es una motivacion previa
con otro alcance; no mezclar sus conteos o bloqueadores con estos resultados.
La novedad a defender es esa contribucion concreta, no el uso de gradientes,
Taylor, Fourier del grafo o SQP por si solos. No se promete un premio.
Las simulaciones diferenciables completas ya son antecedente directo, por
ejemplo Bossart y Hodge (2026), doi:10.1016/j.ijepes.2026.111825. La contribucion
a defender debe centrarse en la sustitucion y su compensabilidad por control.

Archivos principales
METHOD_SPEC.txt: formulacion, derivaciones, condiciones y referencias primarias.
candidate_final_physical.toml: candidato congelado.
physical_candidate_parameters.csv: rho, Kp, Ki y MW por nodo.
nonlinear_codesign_result.pdf: figura de resultados.
PD_physical_final/: validacion compilada independiente y trazas.
long_horizon_final/: seis pruebas de 600 s.
stationarity_final/: gradientes nuevos y diagnostico de optimalidad.
exchange_checkpoint/: pruebas finitas de intercambio y figura.
experiments/nonlinear_codesign_20261001/: codigo Julia y coordinacion Python.
'''
(BASE/'RESULTADOS_ES.txt').write_text(text,encoding='utf-8')
files=list((ROOT/'experiments/nonlinear_codesign_20261001').glob('*.jl'))+list((ROOT/'experiments/nonlinear_codesign_20261001').glob('*.py'))
files+=[ROOT/'Project.toml',ROOT/'Manifest.toml',ROOT/'src/pd39/model.jl',ROOT/'src/bnd_model_expN/PDExactDesignN.jl']
files+=list((ROOT/'src/bnd_design_e').glob('*.jl'))+list((ROOT/'reports/experiment_D/inputs').glob('*.csv'))
files+=[ROOT/'reports/experiment_N/TABLE_N01_original_operating_point.csv',ROOT/'experiments/nonlinear_codesign_20261001/problem.toml']
manifest={'python':platform.python_version(),'candidate_sha256':sha,'files':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.is_file()}}
(BASE/'reproducibility_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({k:summary[k] for k in ['replacement_percent','retained_MW','observed_frequency_max_Hz','observed_RoCoF_max_Hz_s','post_disturbance_max_real_pole','additional_replacement_percentage_points']},indent=2))
