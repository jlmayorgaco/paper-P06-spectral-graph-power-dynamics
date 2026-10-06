from pathlib import Path
import json, hashlib, sys, platform
import pandas as pd
import numpy as np
import scipy

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def texnum(x):
    m,e=f'{x:.3e}'.split('e')
    return m+r'\times10^{'+str(int(e))+'}'
r=json.loads((OUT/'RESULTS.json').read_text());p=r['parity']
ev=pd.read_csv(OUT/'EVENT_RESULTS.csv');ev=ev[ev.model=='physical']
ba=pd.read_csv(OUT/'ENERGY_AND_FREQUENCY_AREA.csv');ba=ba[ba.model=='physical']
g=pd.read_csv(OUT/'GOVERNOR_BRANCHES.csv');active=g[(g.model=='physical')&(g.limited_samples>0)]
pair=pd.read_csv(OUT/'PAIRED_COMPARISON.csv');pair=pair[pair.model=='physical']
eq='sum_g d_g integral(nu_g) - Pm39(0) integral(nu39) = Bgov + Aaw + integral(Delta Pdc) - integral(Delta load_and_losses) - Delta(K+Ef+Edc)'
lines=[]
for x in ev.itertuples():
    verdict='PASS' if x.frequency_pass and x.rocof_pass else ('FAIL RoCoF' if not x.rocof_pass else 'FAIL frecuencia')
    lines.append(f'{x.event:17s}  {x.Fpeak:.9f} Hz  {x.Rpeak:.9f} Hz/s  {verdict}')
report=f'''PUENTE FÍSICO Y LEY INTEGRAL DE FRECUENCIA — IEEE39
Campaña limitada, preregistrada y ejecutada el 4 de octubre de 2026 (Bogotá)

RESULTADO PRINCIPAL
Se cierra el defecto de signo del almacenamiento DC en una variante aislada,
y se demuestra/verifica una identidad de frecuencia–energía para el modelo
IEEE39 con gobernadores, cargas ZIP y pérdidas. La corrección conserva el
espectro y prácticamente las respuestas AC de este candidato. No mejora
los dos eventos antes fallidos y no produce un máximo global de sustitución.

1. Qué se cambió y qué se preservó
PhysicalBridge.jl cambia sólo C vdot=(Pdc-Pac)/v y los dos errores PI a v-V.
No se modificó ninguna fuente ni resultado congelado. Los hashes de todos
los inputs preregistrados siguen iguales: {r['all_frozen_inputs_unchanged']}.
Los parámetros rho, Kp, Ki, ratings, P/Q y demás ganancias permanecieron fijos.
Se usó el candidato Q2B: {r['fixed_GFL_MW']:.9f} MW GFL
({r['fixed_GFL_percent']:.9f}%), {r['retained_SG_MW']:.9f} MW SG.
Ese porcentaje procede del candidato anterior; NO es una optimización nueva
ni prueba de que sea robusto o máximo para esta variante física.

2. Matemática cerrada
EXACT_IDENTITY: E=K+Ef+Edc satisface dot(E)=Pm+Pdc-carga-pérdidas
bajo las ecuaciones declaradas; la versión heredada usaba K+Ef-Edc.
EXACT_IDENTITY: en trim, reflexión de sólo delta(vdc) produce Ap=S Ao S,
Ep=S Eo S. Los puertos AC lineales se conservan con las hipótesis escritas.
La reflexión NO es una equivalencia no lineal general; su defecto comienza
en segundo orden. Estos resultados usan álgebra clásica, no una nueva teoría
universal de similitud ni una demostración de pasividad.
EXACT_IDENTITY: para cualquier tiempo finito desde el trim,
{eq}
Los términos terminales del gobernador y su antiwindup no se omiten.
La salida de TGOV1 es torque por velocidad; SG39 conserva torque fijo.
La prueba completa, las excepciones y una ley infinita CONDICIONAL están en
PHYSICAL_BRIDGE_SECTION.tex y MATH_REVIEW.txt.

3. Paridad y comprobaciones independientes
Residual de equilibrio: {p['physical_residual']:.6g}.
Error relativo de similitud del Jacobiano: {p['jacobian_similarity_relative']:.6g}.
Error relativo de masa descriptor: {p['mass_similarity_relative']:.6g}.
{p['finite_pole_count']} polos finitos (incluido un gauge); matching biyectivo
con error máximo {p['bijective_pole_matching_absolute']:.6g}; alpha física={p['physical_alpha']:.12f} 1/s.
La reproducción NUEVA legacy +1MW coincide punto a punto con su traza archivada.
El error de balance positivo gradiente(E) por RHS versus potencias independientes
es como máximo {r['max_energy_RHS_error_MW']:.6g} MW en los cuatro casos físicos.
Gate numérico de energía de campaña: {r['positive_energy_gate_pass']}.
Esto evalúa una identidad del modelo; no es certificación intervalar ni medida experimental.

4. Cinco integraciones, sin barridos ni retuning
Julia/PowerDynamics: una reproducción legacy y cuatro casos físicos.
Rodas5P, atol=rtol=1e-10, 0–61s, muestras cada 5ms.
El evento a t=1s cambia Pset ZIP, con Qset fijo: el consumo realizado depende
del voltaje y no equivale exactamente al rótulo nominal de MW.
Frecuencia y RoCoF se miden por diferencias causales de fase en buses 30–39
con W=0.5s. NO son máximos de los 39 buses ni derivadas instantáneas.
Límites: 0.5 Hz y 0.5 Hz/s.

{chr(10).join(lines)}

Sólo {r['physical_frequency_and_rocof_pass_count']}/4 cumplen ambos límites,
incluido el ensayo pequeño de 1MW. De los tres casos de 100MW, sólo bus16 pasa.
Las diferencias máximas AC respecto a trazas heredadas son
{pair.max_grid_difference_Hz.max():.6g} Hz y {pair.max_RoCoF_difference_Hz_s.max():.6g} Hz/s.
La corrección física no explica ni elimina los incumplimientos previos.
NEGATIVE_RESULT: este diseño fijo de {r['fixed_GFL_percent']:.6f}% no es seguro para la colección ensayada.

5. Balance finito y limitador observado
La identidad se verifica también integrando las trazas, sin usar un endpoint
supuesto a infinito. A 5ms, error máximo de cierre de energía:
{ba.energy_integral_max_error_MJ.max():.6g} MJ. Error máximo de área:
{ba.frequency_area_max_error_MJ.max():.6g} MJ.
Se conservan los errores a 5/10/20ms en QUADRATURE_SENSITIVITY.csv; la capa
rápida del evento no queda resuelta exactamente por el muestreo de 5ms.
No confundir estos errores de cuadratura con el residual punto a punto de RHS.
En bus29/+100MW se observa directamente el gobernador de SG35 limitado:
{int(active.limited_samples.sum())} muestras, duración estimada {active.limited_sample_duration_s.sum():.3f}s,
con corrección antiwindup {active.antiwindup_integral_MJ.sum():.9f} MJ.
El exceso máximo de válvula es ~2.37e-8 pu por el limitador numérico con desigualdad
estricta; no presentarlo como un margen de actuador positivo ni como gran violación.

6. Qué falta exactamente para la contribución de reemplazo por grafo
El balance físico positivo YA no bloquea esta variante. Todavía bloquean un
piso universal IEEE39 las contribuciones dependientes de la trayectoria:
consumo ZIP, pérdidas, limitadores, almacenamiento terminal y fases finales.
La ley infinita requiere además convergencia/integrabilidad y coherencia de
la salida medida. No se demostraron con un horizonte de 61s.
La teoría lossless de resistencia secante permanece válida en su clase.
No se transfiere simplemente reemplazando su Laplaciana por un Jacobiano AC.

Próximo gate concreto: fijar una familia de diseños y acotar uniformemente
el término residual y los endpoints de la ley de área. Un camino comprobable
es combinar una cota de alcanzabilidad finita con una cota de cola y regiones
de gobernador declaradas; si resulta vacuo, reportarlo. Medir un residual
pequeño después de optimizar no crea una cota predictiva ni global.
Sólo después corresponde co-diseñar rho/Kp/Ki y validar todos los eventos.

7. Alcance y decisión de póster
Este trabajo es una corrección física documentada y un puente verificable,
no evidencia suficiente de novedad competitiva, superioridad de control,
máximo de SG→GFL, estabilidad no lineal global o seguridad con retrasos.
La campaña usa tau_PLL=0; W=0.5s es medida, no retraso PLL.
No hay current limiter ni límite energético de la fuente DC demostrados.
POSTER MAIN CLAIM READY: NO para un máximo IEEE39 demostrado por la ley del grafo.
Gate faltante: cota uniforme de los defectos del modelo completo + diseño que
cumpla todos los eventos del contrato. El documento puede mostrar los resultados
reducidos con su alcance y estos ensayos como auditoría, sin venderlos como ese máximo.
'''
(OUT/'REPORT.txt').write_text(report,encoding='utf-8')
readme='''Leer primero REPORT.txt. Resultados regenerables: RESULTS.json, EVENT_RESULTS.csv,
ENERGY_AND_FREQUENCY_AREA.csv, PAIRED_COMPARISON.csv y GOVERNOR_BRANCHES.csv.

Teoría: PHYSICAL_BRIDGE_SECTION.tex (fragmento integrado en el THEORY.tex existente).
Revisiones independientes: MATH_REVIEW.txt, EVIDENCE_REVIEW.txt.
Código: PhysicalBridge.jl, run_validation.jl, analyze.py, build_report.py.
Freeze: PREREGISTRATION.json + hash; RUN_ENVIRONMENT.toml; MANIFEST.json.
Figuras: FIG_PHYSICAL_BRIDGE y FIG_FREQUENCY_AREA_BRIDGE, PNG y SVG.
Trazas completas: *_SENSORS.csv. Balances derivados: *_BALANCE.csv.

REPRODUCCIÓN SIN SOBREESCRIBIR ESTA CAMPAÑA
Desde la raíz del repo, copiar sólo los siguientes archivos a un nuevo directorio
experiments/physical_frequency_bridge_replay_FECHA_HORA/:
PhysicalBridge.jl, run_validation.jl, analyze.py, PREREGISTRATION.json,
PREREGISTRATION.sha256. Mantener Project/Manifest y las fuentes/input hashes.
Luego ejecutar desde la raíz:
julia --project=. experiments/physical_frequency_bridge_replay_FECHA_HORA/run_validation.jl
python experiments/physical_frequency_bridge_replay_FECHA_HORA/analyze.py
No ejecutar prepare.py otra vez: rechaza una preregistración existente.
El replay usa los mismos cinco casos y no ajusta parámetros por sus resultados.
El análisis requiere Python, numpy, pandas, scipy y matplotlib; versiones en MANIFEST.
No hay commits ni pushes en esta campaña.
'''
(OUT/'README.txt').write_text(readme,encoding='utf-8')
claims=[
 ('C1','EXACT_IDENTITY','Positive DC/filter/rotor storage closes under the new declared model','PHYSICAL_BRIDGE_SECTION.tex','Not global passivity or all physical device energy'),
 ('C2','EXACT_IDENTITY','DC-state reflection preserves tangent descriptor/AC transfer at trim','PHYSICAL_BRIDGE_SECTION.tex','Fixed V; no additional vdc modulation dependence; not nonlinear conjugacy'),
 ('C3','NUMERICALLY_VALIDATED','144 finite poles match bijectively; equilibrium and mass/Jacobian parity pass','PARITY.toml;POLE_MATCHING.csv','One gauge excluded for physical alpha'),
 ('C4','EXACT_IDENTITY','Finite rotor-frequency area balance includes governor endpoints and antiwindup','PHYSICAL_BRIDGE_SECTION.tex','Fixed references, trim, actual load and losses'),
 ('C5','NUMERICALLY_VALIDATED','Positive storage RHS gate passes in all four physical runs','EVENT_RESULTS.csv','Five-ms quadrature residuals retained separately'),
 ('C6','NEGATIVE_RESULT','Fixed 94.215489% candidate fails bus8 RoCoF and bus29 frequency','EVENT_RESULTS.csv','Windowed generator-bus outputs, selected events, not a new optimum'),
 ('C7','NUMERICALLY_VALIDATED','SG35 limiter directly active in bus29 case; antiwindup area is nonzero','GOVERNOR_BRANCHES.csv','Sampling duration; not an actuator guarantee'),
 ('C8','BLOCKED','Predictive graph replacement floor for full IEEE39','MATH_REVIEW.txt','Uniform residual/endpoint/tail bounds not established')]
pd.DataFrame(claims,columns=['claim_id','status','claim','evidence','scope']).to_csv(OUT/'CLAIM_LEDGER.csv',index=False)
rows='\n'.join(f"{x.event.replace('_',r'\_')} & {x.Fpeak:.6f} & {x.Rpeak:.6f} & {'pass' if x.frequency_pass and x.rocof_pass else 'fail'} \\\\" for x in ev.itertuples())
numerical=rf'''
\subsection{{Preregistered zero-delay numerical check}}
Five Julia/PowerDynamics integrations comprise one independent replay of the
legacy 1 MW event and four physical-variant events. The fixed candidate has
{r['fixed_GFL_percent']:.6f}\% GFL and {r['retained_SG_MW']:.6f} MW retained SG;
no design variable is retuned. The trim residual is
${texnum(p['physical_residual'])}$, and all {p['finite_pole_count']} finite poles,
including one gauge, match bijectively with numerical error
${p['bijective_pole_matching_absolute']:.4g}$. The physical rightmost pole is
${p['physical_alpha']:.9f}\,\mathrm{{s}}^{{-1}}$.

\begin{{center}}\begin{{tabular}}{{lrrl}}\toprule
ZIP event & $\max|\Delta f|$ [Hz] & RoCoF [Hz/s] & Both limits\\\midrule
{rows}
\bottomrule\end{{tabular}}\end{{center}}
These are sampled maxima at generator buses 30--39, measured by a causal
0.5 s phase window over 0--61 s, with limits 0.5 Hz and 0.5 Hz/s.
The event is a nominal ZIP-setpoint increment, not constant realized demand.
The two failed 100 MW events remain failed after correcting the energy port.

The positive-energy gradient-times-RHS identity has maximum residual
${texnum(r['max_energy_RHS_error_MW'])}$ MW. Independently integrating stored
traces at 5 ms leaves maximum energy and frequency-area residuals of
{ba.energy_integral_max_error_MJ.max():.6f} MJ and
{ba.frequency_area_max_error_MJ.max():.6f} MJ; these are exposed quadrature
errors, not interval certificates. The SG35 governor is directly observed
on its limiting branch in {int(active.limited_samples.sum())} samples of the
bus29 event, giving approximately {active.limited_sample_duration_s.sum():.3f} s
of sampled limiting and an antiwindup integral of
{active.antiwindup_integral_MJ.sum():.6f} MJ. Ignoring it changes the area balance.
The isolated campaign is reproducible from
\texttt{{experiments/physical\_frequency\_bridge\_20261004/}}.
All frozen source and candidate hashes were checked unchanged.
This closes the physical-port audit for the new variant. It does not close
the full-model replacement-floor or global-optimality gates.
'''
section=(OUT/'PHYSICAL_BRIDGE_SECTION.tex').read_text(encoding='utf-8').replace('% NUMERICAL_BRIDGE_RESULTS',numerical)
(OUT/'NUMERICAL_SECTION.tex').write_text(numerical,encoding='utf-8')
target=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
backup=OUT/'THEORY_before_bridge.tex'
if not backup.exists():backup.write_bytes(target.read_bytes())
text=target.read_text(encoding='utf-8');start='% BEGIN PHYSICAL_FREQUENCY_BRIDGE_20261004';end='% END PHYSICAL_FREQUENCY_BRIDGE_20261004'
block=start+'\n'+section+'\n'+end+'\n\n'
if start in text:
    a=text.index(start);b=text.index(end,a)+len(end);text=text[:a]+block.rstrip()+text[b:]
else:
    anchor=r'\section*{Reproducibility and research integrity}';assert text.count(anchor)==1;text=text.replace(anchor,block+anchor)
eol='\r\n' if b'\r\n' in backup.read_bytes() else '\n'
target.write_bytes(text.replace('\r\n','\n').replace('\n',eol).encode('utf-8'))
(OUT/'MANUSCRIPT_INTEGRATION.json').write_text(json.dumps({'target':str(target),'before_sha256':sha(backup),'after_sha256':sha(target),'compiled':False,'scope':'Append physical-energy derivation and executed numerical bridge; retain prior results'},indent=2),encoding='utf-8')
manifest={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'inputs':json.loads((OUT/'PREREGISTRATION.json').read_text())['inputs'],'outputs':{p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='MANIFEST.json'}}
(OUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Saved report, claim ledger and in-place manuscript extension:',target)
