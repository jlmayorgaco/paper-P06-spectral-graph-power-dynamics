"""Consolidate completed evidence; no incomplete trajectory counts as a pass."""
import csv,json
from pathlib import Path

root=Path(__file__).resolve().parents[2]
out=root/'reports/codesign_validation_20261001'
data={k:json.loads((out/f'validation_summary_{k}.json').read_text(encoding='utf-8')) for k in ('improved','repaired')}
a=data['improved'];b=data['repaired']
rows=[]
for name,d in data.items():
    for r in d['events']:
        rows.append({'candidate':name,'event':r['event'],'horizon_after_event_s':60,
                     'Fpeak_Hz':r['Fpeak_Hz'],'Rpeak_Hz_s':r['Rpeak_Hz_s'],
                     'completed':True,'F_pass':r['F_pass'],'R_pass':r['R_pass']})
    for r in d['bus8_diagnostic_prefix']:
        rows.append({'candidate':name,'event':'bus8_prefix_'+r['label'],'horizon_after_event_s':5,
                     'Fpeak_Hz':r['Fpeak_Hz'],'Rpeak_Hz_s':r['Rpeak_Hz_s'],
                     'completed':True,'F_pass':float(r['Fpeak_Hz'])<=.5,'R_pass':float(r['Rpeak_Hz_s'])<=.5})
with (out/'campaign_event_comparison.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

lines=['RESULTADO DE LA CAMPAÑA SG–GFL — 2026-10-01','',
 'La campaña NO demuestra reemplazo máximo ni un ajuste aprobado para todas las perturbaciones requeridas.',
 'Sí encuentra un candidato con menos SG y confirma su estabilidad lineal local. Una prueba no lineal requerida lo rechaza.',
 '', 'DISEÑO ENSAYADO',
 f'SG retenido: {a["candidate"]["retained_SG_MW"]:.6f} MW; sustitución GFL: {a["replacement_percent"]:.6f}% del despacho base.',
 'El candidato anterior retenía 438.446282 MW. Los ensayados aquí retienen 366.238767 MW y eliminan completamente SG en 31, 32, 34 y 36.',
 'Las dos variantes usan la misma asignación rho; repaired aumenta algunas Kp del PLL. Su nombre es una etiqueta experimental, no una reparación certificada.',
 '', 'EVIDENCIA POR VARIANTE']
for name,d in data.items():
    s=d['pd_spectrum'];p=d['bus8_diagnostic_prefix']
    lines += [f'{name}: alpha_PD={s["pd_alpha"]:.12g} s^-1; {s["physical_poles"]} polos físicos; error máximo de reconstrucción={s["max_matched_pole_error"]:.3g}.',
              f'Suite completa aprobada: {d["declared_suite_pass"]}. Casos con violación demostrada: {d["required_failed"]}.',
              f'Casos requeridos sin trayectoria completa: {d["required_missing"]}.']
    for r in p:
        lines.append(f'  Bus 8 +100 MW, prefijo completo de 5 s, {r["label"]}: F={float(r["Fpeak_Hz"]):.9f} Hz, RoCoF={float(r["Rpeak_Hz_s"]):.9f} Hz/s; tolerancia={r["tol"]}.')
    for r in d['incomplete']:
        lines.append(f'  Integración incompleta {r["event"]}: {r["retcode"]}, t absoluto={r["last_time_s"]} s (evento en t=1).')
    if d.get('bus8_prefix_incomplete'):
        lines.append('  El intento adicional de prefijo de 5 s tampoco completó la integración; no se le atribuyen picos ni convergencia.')
    for r in d['events']:
        if r['event'].startswith('small_'):
            lines.append(f'  Contraste a 1 MW, {r["event"]}: error relativo máximo de forma de onda = {100*float(r["linear_relative_waveform_error"]):.4f}%.')
lines += ['', 'QUÉ ESTÁ COMPROBADO Y QUÉ NO',
 'Los límites declarados son 0.5 Hz y 0.5 Hz/s con ventana causal de fase de 0.5 s. Se observan los diez buses generadores.',
 'La violación de RoCoF en un prefijo completado basta para rechazar ese candidato: no depende del posterior fallo del integrador.',
 'Un fallo MaxIters no demuestra por sí solo inestabilidad física, pérdida de sincronismo ni colapso de tensión.',
 'El escalón de 150 MW en bus 29 es estrés adicional fuera de la especificación de diseño de 100 MW; se informa separado.',
 'El radio de robustez usa perturbaciones complejas no estructuradas de la matriz en coordenadas fijas. Su cota de toda la banda es numérica Float64, no una prueba formal ni una garantía ante eventos grandes.',
 'Las tres búsquedas SLSQP sumaron 14 770 evaluaciones. No convergieron con certificación KKT; el menor residuo de estacionariedad normalizado observado fue aproximadamente 0.04747. No se verificó SOSC ni se obtuvo una cota global del objetivo.',
 'La comprobación del evaluador incluye unidades físicas, derivadas, exponencial matricial independiente, respuesta lineal refinada y cota de su cola. Estas comprobaciones pasan y no convierten la aproximación lineal en una garantía no lineal.',
 'El contraste PowerDynamics–modelo lineal a 1 MW deja diferencias de forma de onda de hasta 2.69%. No se ha completado un barrido de amplitud tendiendo a cero que separe el error de linealización del de medición o implementación de la entrada/salida.',
 '', 'CONSECUENCIA MATEMÁTICA PARA LA SIGUIENTE BÚSQUEDA',
 'Mantener min J=sum_i P0_i*(1-rho_i), con restricciones espectrales y del radio, y añadir seguridad de la DAE no lineal para cada escenario.',
 'Para cada bus de perturbación b y amplitud d: sup_t max_i |f_T_i(t;rho,Kp,Ki,b,d)| <= 0.5 y sup_t max_i |R_T_i(t;rho,Kp,Ki,b,d)| <= 0.5.',
 'f_T=[theta(t)-theta(t-T)]/(2*pi*T); R_T=[theta(t)-2*theta(t-T)+theta(t-2*T)]/(2*pi*T^2).',
 'Añadir límites de tensión, corriente y estados del controlador con valores físicos documentados. La corriente relativa al punto inicial no sustituye una corriente nominal admisible.',
 'Usar los casos fallidos como contraejemplos al seleccionar nuevos diseños; reservar perturbaciones y puntos de operación para validación que no haya intervenido en la selección.',
 'Completar un barrido pequeño de amplitudes y localizar las diferencias de entrada/salida antes de atribuir toda discrepancia a un mecanismo no lineal concreto.',
 'Para afirmar óptimo local: gradientes fiables, conjunto activo completo, KKT y condiciones de segundo orden pertinentes. Para máximo global: una cota válida que cierre la brecha entre el mejor diseño factible y la relajación, cubriendo también las arquitecturas SG presentes/ausentes.',
 '', 'LECTURA CIENTÍFICA',
 'El resultado defendible es que estabilidad modal y robustez lineal no bastan para certificar seguridad ante grandes perturbaciones del co-diseño SG–GFL. Esta campaña contiene un contraejemplo numérico reproducible.',
 'No se ha establecido aquí que esa observación sea nueva en la literatura. El póster no debería afirmar máximo reemplazo seguro ni ganancias óptimas mientras falten estas garantías.',
 '', 'ARCHIVOS',
 'VALIDATION_REPORT_improved.txt y VALIDATION_REPORT_repaired.txt: métodos, métricas, fallos y comandos.',
 'parameters_improved.csv y parameters_repaired.csv: rho, potencia SG, Kp y Ki por nodo; son parámetros ensayados, no recomendaciones aprobadas.',
 'campaign_event_comparison.csv: resultados completos y prefijos, con su horizonte explícito.',
 'validation_improved.png/.pdf y validation_repaired.png/.pdf: figuras científicas en formato vertical.',
 'source_manifest.json y los SHA256 de candidatos: trazabilidad. Los scripts están en experiments/codesign_validation_20261001.']
(out/'CAMPAIGN_FINDINGS.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(out/'CAMPAIGN_FINDINGS.txt')
