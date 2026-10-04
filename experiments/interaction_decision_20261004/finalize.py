"""Evidence audit and human-readable final reports; never edits historical data."""
from model import *
import scipy,matplotlib,zipfile

def main():
    load=lambda n:json.loads((OUT/n).read_text())
    originals=load('INPUT_MANIFEST.json')['inputs']
    for path,digest in originals.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    anchor=load('CERT_anchor.json');single=[load('CERT_single30.json'),load('CERT_single37.json')]
    joint=load('CERT_joint.json');corrected=load('CERT_corrected.json');cp=load('CERT_complex_pair.json')
    assert all(load('CERT_'+name+'.json')['count_one_proved'] for name in ['anchor','single30','single37','joint','corrected','complex_pair','policy_below','policy_above'])
    assert len(load('COMMON_CONTOUR_CERTIFICATES.json'))==4 and all(x['one_root'] for x in load('COMMON_CONTOUR_CERTIFICATES.json'))
    gN=single[0]['center'][0]+single[1]['center'][0]-anchor['center'][0]+.05
    interaction=joint['center'][0]-single[0]['center'][0]-single[1]['center'][0]+anchor['center'][0]
    assert gN+3e-7<0 and joint['center'][0]-1e-7>-.05
    assert corrected['center'][0]+1e-7<-.05
    events=pd.read_csv(OUT/'TABLE_12_NONLINEAR_EVENTS.csv');spectra=pd.read_csv(OUT/'TABLE_11_FULL_SPECTRUM.csv')
    parity=pd.read_csv(OUT/'TABLE_10_JULIA_PARITY.csv');design=load('DESIGN_SUMMARY.json');budget=load('BUDGET_RESULT.json');cross=load('POLICY_CROSSING.json')
    assert len(events)==5 and events['pass'].all() and parity['pass'].all()
    assert int(spectra[spectra.design=='joint']['count'].iloc[0])==2
    assert spectra[spectra.design=='corrected'].status.iloc[0]=='PASS_NUMERICAL'
    # Independent determinant identity at three generic nonroot frequencies.
    a=Model(np.array(anchor['p']));pp=np.array(joint['p']);p1=np.array(single[0]['p']);p2=np.array(single[1]['p'])
    ix=np.r_[GROUPS[0],GROUPS[7]];d=action(pp-a.p)[ix];identities=[]
    for s in [-.2+29j,-.1+30.4j,.1+31.2j]:
        H=d[:,None]*a.T(s)[np.ix_(ix,ix)]
        K12=la.solve(np.eye(4)+H[:4,:4],H[:4,4:]);K21=la.solve(np.eye(4)+H[4:,4:],H[4:,:4])
        rhs=la.det(np.eye(4)-K12@K21)
        vals=[np.linalg.slogdet(m.descriptor(s)) for m in [Model(pp),a,Model(p1),Model(p2)]]
        lhs=vals[0][0]*vals[1][0]/(vals[2][0]*vals[3][0])*np.exp(vals[0][1]+vals[1][1]-vals[2][1]-vals[3][1])
        err=abs(lhs-rhs)/max(abs(lhs),abs(rhs));assert err<1e-8
        identities.append(dict(s_real=s.real,s_imag=s.imag,relative_error=float(err)))
    pd.DataFrame(identities).to_csv(OUT/'TABLE_16_DETERMINANT_IDENTITY.csv',index=False)
    status={'claim_status':'NUMERICALLY_VALIDATED','experiment_status':'FINITE_DECISION_AND_CORRECTED_DESIGN_CLOSED',
        'rigorous_endpoint_root_certificates':8,'rigorous_common_contour_counts':4,
        'certified_joint_real':joint['center'][0],'coordinate_radius':1e-7,
        'interaction_real':interaction,'interaction_radius':4e-7,'corrected_real':corrected['center'][0],
        'all_root_spectrum':'NUMERICAL_CONTOUR_NOT_INTERVAL_CERTIFIED','five_events_passed':int(events['pass'].sum()),
        'F_Hz_max':float(events.F.max()),'RoCoF_Hz_s_max':float(events.R.max()),'SG_slack_min':float(events.slack.min()),
        'GFL_MW':design['total_GFL_MW'],'GFL_percent':design['GFL_percent'],'SG_MW':design['retained_SG_MW'],
        'new_tuning_rule':'equal modal-decay-budget fixed point at fixed rho and Kp',
        'global_optimality':'NOT_ESTABLISHED','uniform_contraction_certificate':'NOT_ESTABLISHED',
        'journal_submission_ready':False,'frozen_inputs_preserved':True}
    save('STATUS.json',status)
    report=f'''# Resultado y reglas de diseño — 4 de octubre de 2026

**Resultado cerrado:** hay un contraejemplo finito certificado en IEEE-39 y una
regla iterativa que corrige su margen sin reducir los MW reemplazados.
La aportación es una decisión concreta; no se ha demostrado máxima penetración.

## Qué se demostró

- Sustitución adicional: {design['added_GFL_MW']:.6f} MW entre buses 30 y 37, con
  rho=.885 en esos dos y .875 en los otros ocho. Retardo fijo de 40 ms en los diez.
- Cada ajuste individual conserva el polo cerca de {anchor['center'][0]:.9f}/s.
- Juntos: {joint['center'][0]:.9f}/s, que viola el margen exigido de -.05/s.
- Interacción real: {interaction:.9f} ±4.01e-7 /s, incluida holgura de redondeo.
- Los cuatro modelos tienen exactamente una raíz en el mismo contorno,
  certificado por intervalos. Cada raíz se encierra en una caja de radio 1e-7.
- Julia reproduce la linealización con error relativo máximo
  {parity.jacobian.max():.3e}; su contorno de espectro completo encuentra 0 raíces
  a la derecha del margen para ambos sitios por separado, y 2 para el conjunto.
  Este conteo global es numérico, no certificado por intervalos.

## Regla derivada y ejecutada

Separar exactamente la respuesta individual y la interacción:

    alpha_joint(t) = alpha_anchor - t + I(t)
    t_next = max(0, I(t) - b),  b = -alpha_anchor - sigma

Cada sitio recibe t/2 de amortiguamiento modal adicional; rho y Kp permanecen
fijos. I(t) se evalúa con el modelo DDE completo. El teorema de contracción
requiere una cota uniforme |I'|<1; aquí solo se comprobó empíricamente su pendiente.
El resultado final, en cambio, tiene su raíz encerrada rigurosamente.

Con sigma=.06/s, {budget['iterations']-1} actualizaciones llevan a
t={budget['t']:.10f}/s y alpha={corrected['center'][0]:.10f}/s.
Ki30={budget['p'][20]:.10f}; Ki37={budget['p'][27]:.10f}.
Kp30=Kp37={budget['p'][10]:.10f}. Todos los demás Ki={budget['p'][21]:.10f}.
El vector completo se encuentra en TABLE_14_DESIGN_VECTORS.csv.

**Consejo de diseño sustentado:** no sumar certificados de polos individuales
para aceptar acciones simultáneas. Calcular la interacción y asignar un margen
adicional, o diseñar el eigenpar colectivo con un patrón compatible. Preservar
un polo complejo local tampoco garantiza composición universal; en este caso
fue una alternativa eficaz, con deriva de solo
{cp['center'][0]-anchor['center'][0]:.3e}/s.

## Frontera de una política, no máximo global

La política registrada de ajuste individual encuentra su primer cruce numérico
entre {cross['added_MW_lower']:.6f} y {cross['added_MW_upper']:.6f} MW adicionales.
Los signos a ambos extremos están certificados. No se probó monotonía uniforme
ni exclusión de todas las demás ganancias. El diseño corregido conserva 7.9 MW.
La estimación cuadrática local predijo h={cross['curvature_prediction_h']:.6f},
frente al cruce cercano a h={(cross['h_lower']+cross['h_upper'])/2:.6f}; no debe
usarse como cota sin controlar su resto.

## Validación no lineal

El diseño corregido pasa 5/5 eventos congelados con método de pasos DDE en Julia:
Fmax={events.F.max():.8f} Hz; RoCoFmax={events.R.max():.8f} Hz/s;
V en [{events.Vmin.min():.8f}, {events.Vmax.max():.8f}];
holgura mínima SG={events.slack.min():.8f}. Frecuencia y RoCoF usan ventana de 0.5 s.
GFL={design['total_GFL_MW']:.6f} MW ({design['GFL_percent']:.8f}%);
SG retenida={design['retained_SG_MW']:.6f} MW.
Es un diseño factible para ese contrato numérico, no un óptimo.

## Qué falló o sigue abierto

- El primer barrido no dio falsas aceptaciones válidas; hubo 5 falsos rechazos.
- 230 combinaciones con conservación del polo complejo no violaron el margen.
- En la tercera campaña hubo 88 casos resueltos y 317 no resueltos; 2 violaciones.
- El criterio de convergencia de la serie por recorridos falla en el contorno
  del caso seleccionado. Se usa la identidad racional exacta y certificados
  directos de las raíces, no una cota de serie inválida.
- El predictor cuadrático ordinario también detecta este caso. La corrección
  de dos pasos mejora la precisión, pero no tiene exclusividad sobre la decisión.
- Sin máximo global, cota superior de reemplazo, garantía sobre todo el rango
  de ganancias, comparación contra optimizadores publicados o incertidumbre
  física robusta. No se modela seguridad de limitador de corriente.

## Juicio editorial

Esto cierra una brecha real respecto de la versión anterior: ahora hay una
decisión que cambia, un certificado y una corrección ejecutada. Es material
para una contribución enfocada, pero no basta para anunciar un journal fuerte
cerrado o un premio. La prioridad bibliográfica no está probada; el siguiente
paso editorial requiere validación independiente más amplia y comparación
directa contra métodos existentes, conservando exactamente este contrato.
'''
    (OUT/'REPORT_ES.md').write_text(report,encoding='utf-8')
    claims=f'''# Claim ledger

| Claim | Status | Evidence and boundary |
|---|---|---|
| Exact pair determinant and rational interaction | EXACT_IDENTITY | Section9; fixed chart, common counts and boundary inverses required; classical tools specialized to physical interventions |
| Scalar interaction-budget convergence | PROVED_REDUCED_MODEL | Conditional uniform derivative and admissible-map hypotheses; these are not all certified for IEEE39 |
| Singleton-additive modal decision reversal | NUMERICALLY_VALIDATED | Rigorous exported-model endpoint/common-contour certificates; I={interaction:.9f} ±4.01e-7/s |
| Corrected enclosed pole | NUMERICALLY_VALIDATED | Rigorous exported-model box centered {corrected['center'][0]:.9f}/s, radius1e-7 |
| Complete spectrum margin and nonlinear events | NUMERICALLY_VALIDATED | Floating-point Julia contour0roots beyondmargin;5/5fixed events; no robust uncertainty claim |
| Iterative controller correction | SUPPORTED_LOCAL | {budget['iterations']-1}updates, fixedrho/Kp, no generic convergence or gain-effort optimum |
| Policy crossing around6.385MW | SUPPORTED_LOCAL | Numerical tracked first crossing with certified endpoint signs; no interval-wide monotonicity proof |
| Complex-pole singleton rule materially fails in tested family | NEGATIVE_RESULT | 0/230violations; maximum drift2.91e-5/s |
| Walk-series certificate applies to selected common contour | NEGATIVE_RESULT | Pair spectralradius>1; numerical agreement ofC2 is not a certificate |
| Global or gain-uniform maximum replacement | BLOCKED | No informative upperbound, no optimalitygap |
| Publication priority, GSP superiority, communication necessity | BLOCKED | Not established by these experiments |
'''
    (OUT/'POSTER_CLAIMS.md').write_text(claims,encoding='utf-8')
    (OUT/'LATEX_COMPILATION.json').write_text(json.dumps({'native':'unavailable: Unable to find standard directories for platform','fallback':'installed MiKTeX pdflatex success twice','pages':38,'overfull_warnings':0,'undefined_references':0,'visual_review':'all changed pages rendered, inspected; no out-of-page text'},indent=2))
    # Post-execution record. Do not imply code was hashed before discovery.
    exclude={'.aux','.log','.out','.pyc','.zip'}
    files=[p for p in OUT.rglob('*') if p.is_file() and p.suffix not in exclude and '__pycache__' not in p.parts and 'qa' not in p.parts and not p.name.endswith('.zip.sha256') and p.name not in ['VALIDATION_MANIFEST.json','BUNDLE_MANIFEST.json']]
    source=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex';files.append(source)
    manifest={'type':'POST_EXECUTION_VALIDATION_MANIFEST','head_before_work':load('INPUT_MANIFEST.json')['head'],
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'python_flint':'0.8.0','arb_bits':160},
        'files':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        'frozen_inputs_preserved':True,'scope':'not retroactive preregistration; hypotheses evolved through timestamped addenda'}
    save('VALIDATION_MANIFEST.json',manifest)
    print(json.dumps(status,indent=2))
if __name__=='__main__':main()
