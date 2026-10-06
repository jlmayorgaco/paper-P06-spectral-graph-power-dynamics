from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports/nonlinear_codesign_20261001/coupled_network_contract"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda n:json.loads((OUT/n).read_text(encoding="utf-8"))
s=read("synthesis_modal.json");v=read("verification.json");a=read("mapping_audit.json")
c=read("refined_center.json");d=read("wrapping_diagnostic.json");g=read("graph_structure_audit.json")
ar=read("arithmetic_audit.json")
assert v["status"]=="NONLINEAR_COUPLED_ENCLOSURE_RECHECK_PASSED_MICROSCOPIC_REGION"
assert all(sha(OUT/n)==h for n,h in v["source_hashes"].items())
assert v["implementation_sha256"]==a["implementation_sha256"]==ar["implementation_sha256"]
npz=np.load(OUT/"certificate_modal.npz")
rows=s["attempts"]
fig,ax=plt.subplots(1,2,figsize=(11.4,4.7))
fig.subplots_adjust(left=.075,right=.98,bottom=.22,top=.79,wspace=.28)
fig.suptitle("Coupled nonlinear certificate: measuring the loss from conservatism",fontsize=14,fontweight="bold",x=.075,ha="left",y=.97)
fig.text(.075,.865,"281 states · all 39 bus instruments · 88 disturbance channels · fixed SG/GFL design",fontsize=10,color="#455764")
r=np.array([q["radius"] for q in rows]);eps=np.array([q["epsilon"] for q in rows]);dec=np.array([q["decay"] for q in rows])
ax[0].plot(r[eps>0],eps[eps>0],"o-",color="#197c80",lw=1.8)
ax[0].plot(v["radius"],v["input_radius"],"o",color="#d75b32",ms=8)
ax[0].set(xscale="log",yscale="log",xlabel="Ellipsoid radius in the modal metric",ylabel="Admissible joint-input radius")
ax[0].set_title("A  Positive, but practically negligible envelope",loc="left",fontsize=10,pad=12)
ax[1].plot(r,dec,"o-",color="#197c80",label="Regional enclosure")
ax[1].axhline(0,color="#555",ls="--",lw=1)
for radius in sorted(set(q["radius"] for q in d["rows"])):
    lo=min(q["point_decay_lower_bound"] for q in d["rows"] if q["radius"]==radius)
    ax[1].plot(radius,lo,"x",color="#d75b32",ms=8,mew=2,label="Sampled points only" if radius==min(q["radius"] for q in d["rows"]) else None)
ax[1].set(xscale="log",xlabel="Ellipsoid radius in the modal metric",ylabel=r"Decay lower bound ($s^{-1}$)")
ax[1].set_title("B  Point checks do not replace a regional bound",loc="left",fontsize=10,pad=12)
ax[1].legend(frameon=False,fontsize=9,loc="lower left")
for aa in ax:
    aa.spines[["right","top"]].set_visible(False);aa.grid(alpha=.15);aa.set_axisbelow(True)
fig.text(.075,.10,"The full nonlinear DAE is enclosed throughout the region. A failed bound does not demonstrate instability.",fontsize=9,color="#455764")
fig.text(.075,.052,"No practical replacement capacity, globally optimal gains, or hardware safety is established by this microscopic certificate.",fontsize=9,color="#455764")
for ext in ("png","pdf","svg"):fig.savefig(OUT/f"coupled_certificate_conservatism.{ext}",dpi=180,facecolor="white")

report=f"""CERTIFICADO NO LINEAL ACOPLADO SG/GFL/RED: PROTOTIPO Y CONSERVADURISMO
2026-10-01. Avance implementado; el objetivo de capacidad útil/óptima sigue abierto.

RESULTADO Y LÍMITE PRINCIPAL
Se implementó y volvió a comprobar un certificado conjunto de invariancia para
281 estados, con todas las SG y GFL presentes, frecuencia común física y medidas
de frecuencia/RoCoF en 39 buses. La clase de entradas admite cualquier señal
medible de 88 canales: dos componentes de corriente por bus y diez entradas de
potencia DC. No se selecciona un evento, nodo perturbado ni frecuencia temporal.

La envolvente demostrada por esta implementación es sólo
  epsilon = {v['input_radius']:.12g},
en las unidades normalizadas de los canales declarados. Es MICROSCÓPICA y no
constituye una capacidad útil frente a perturbaciones de operación. El resultado
comprueba la ruta matemática acoplada y revela su conservadurismo. No se optimizó
rho ni Kp,Ki y no se presenta el porcentaje del candidato como máximo robusto.

1. MODELO CONSERVADO Y CLASE DE ENTRADAS
Partimos de la DAE no lineal física de ReducedDAE.jl: Sauer-Pai con AVR/TGOV,
PLL filtrado, PI de corriente, filtro AC y DC con Edot=Pdc-Pac. En esta instancia,
Rs=0 y Xd''=Xq''=X'' en las diez SG. Esa igualdad se comprueba al exportar; las
fórmulas especializadas siguientes no se atribuyen a máquinas que no la cumplen.

Para una SG con cd,cq lineales en sus flujos internos, velocidad omega y ángulo
delta, el puerto EXACTO es
  i_SG = D(omega) v + c(delta,flujos),
  D(omega) = (Sn/Sbase)/(omega X'') * [[0,-1],[1,0]],
  c = (Sn/Sbase)/X'' * [sin(delta) cd+cos(delta) cq,
                       -cos(delta) cd+sin(delta) cq].
La red mantiene su admitancia AC realificada y cargas Z congeladas:
  G(x,rho) v + h(x,rho) + w_I = 0,
  G = Y + sum_i (1-rho_i) S_i D_i(omega_i) S_i^T,
  v = -G^-1 (h+w_I).
Los GFL aportan sus corrientes dinámicas a h. Las entradas DC entran en
  Cdc vdc vdcdot = Pdc + w_DC - Pac.
El peso rho multiplica el puerto del agregado; no interpola arbitrariamente
ecuaciones internas. Esta entrega mantiene el soporte interior del candidato.

Escribimos w=S_w eta, ||eta(t)||_2<=epsilon casi en todo t. El experimento usa
la normalización identidad sobre los canales en pu indicados; corrientes y
potencias DC no son una única amplitud común en MW. La teoría permite cambiar
S_w/R_w. Esta familia incluye entradas simultáneas y cambiantes sin regularidad
de su derivada, pero no engloba automáticamente fallos/topologías o cambios de
admitancia. Esos requieren sus propios canales o conjuntos de incertidumbre.

Se elimina sólo una posición angular común, usando la última SG como referencia
de coordenadas. Su velocidad permanece, nu=omega_base(omega_refSG-omega_frame).
Se restan los términos exactos de rotación de ángulos y corrientes GFL. El rótulo
de esa referencia es una elección de carta, no un nodo especial de perturbación.
203 estados físicos/de control permanecen después de retirar la simetría angular.

2. FRECUENCIA Y ROCOF EN TODOS LOS BUSES
Para una fase relativa regional phi_b(v) y nu como arriba:
  eta_b_dot = (phi_b-eta_b)/tau_f - nu,
  fhat_b = (phi_b-eta_b)/(2 pi tau_f),
  chi_b_dot = (fhat_b-chi_b)/tau_r,
  rhat_b = (fhat_b-chi_b)/tau_r.
Se conservan esos 78 estados. En esta instancia tau_f=tau_r=0.1 s; son parámetros
del instrumento, no constantes de la teoría. Son medidas filtradas finitas,
diferentes de las ventanas exactas de 0.5 s usadas en los ensayos históricos.
La región verifica una rama de fase con parte real rotada positiva en cada bus.
Así, saltos algebraicos de tensión no obligan a diferenciar w ni producen un
RoCoF ideal impulsivo en el instrumento declarado.

3. CERTIFICADO Y DEMOSTRACIÓN
Sea z el estado aumentado respecto a un centro refinado x_c, H invertible,
y=Hz y Omega_r={{z:||Hz||_2<=r}}. No se almacena energía en variables algebraicas
que pudieran saltar con w. Se obtiene una caja exterior a Omega_r mediante
  |z_j| <= r ||fila_j(H^-1)||_2.
Sobre esa caja y |w_k|<=epsilon se verifican:
  G invertible, dominio de denominadores y régimen interior de limitadores;
  sym(H F_z(z,w) H^-1) <= -a I,          a>0;
  ||H F_w(z,w)||_2 <= b;
  ||H F(0,0)||_2 <= d.
Aquí sym(M)=(M+M^T)/2. F es la ELIMINACIÓN EXACTA de la DAE, incluidos sensores.
El Jacobiano nominal inicializa H; no sustituye F por una planta lineal.
Las derivadas se encierran sobre toda la región mediante diferenciación automática
intervalar de senos, cocientes, saturación AVR activa y solución algebraica.

Para la trayectoria radial (tz,tw), el teorema fundamental del cálculo da
  F(z,w)=F(0,0)+integral_0^1 [F_z(tz,tw) z+F_w(tz,tw) w] dt.
La convexidad de las cajas y la bola de entrada permite aplicar las cotas anteriores.
Con s=||y|| y V=s^2 resulta, para s>0,
  s_dot <= -a s + b epsilon + d,
  Vdot <= -2a V + 2 sqrt(V)(b epsilon+d).
Por tanto,
  b epsilon+d < a r
garantiza orientación estricta hacia el interior en s=r. La rama algebraica única,
el campo localmente Lipschitz y la compacidad permiten continuidad e invariancia
para toda entrada medible declarada. Además,
  s(t) <= exp(-at)s(0)+(b epsilon+d)/a*(1-exp(-at)).
Las cotas de salidas se evalúan sobre la misma región completa. El residuo d no
se descarta: con entrada cero la fórmula da una bola residual alrededor del centro
aproximado, no afirma que ese centro numérico sea exactamente un equilibrio.

4. CÓMO SE CONSTRUYÓ H Y CÓMO PARTICIPA EL GRAFO
Una primera métrica de Lyapunov produjo condicionamiento 1.44e7 y márgenes muy
conservadores. Una métrica de alcanzabilidad resultó peor condicionada y falló.
La métrica seleccionada usa una base real de modos del bloque físico, más una
transformación de Sylvester para la cascada EXACTA de los instrumentos. Se
comprueban después invertibilidad y desigualdades del campo no lineal completo.
No se elimina ningún modo dinámico ni se descarta la frecuencia coherente.

P=H^T H contiene términos cruzados entre buses. En una ablación numérica,
retirar sólo esos bloques de ESTE P cambia la mayor tasa de su forma cuadrática
nominal a +{g['dropped_cross_bus_storage_logarithmic_norm']:.6f} 1/s; pierde la propiedad de Lyapunov.
Esto no demuestra que todos los certificados nodales posibles sean inviables.

Se exportó también un Laplaciano estructural de 39 nodos y {g['graph_edges']} aristas,
con pesos simétricos derivados de las normas de bloques de Y, y su reducción
estructural a diez puertos. Sirve como coordenada de interpretación y futura
restricción de estructura de P. La seguridad calculada usa Y AC completo. No se
identifican sus autovalores espaciales con frecuencias temporales ni se pretende
que Fourier del grafo diagonalice el sistema no lineal heterogéneo.

5. RESULTADO NUMÉRICO VOLVIÓ A COMPROBARSE SIN OPTIMIZAR
  r = {v['radius']:.12g}
  epsilon = {v['input_radius']:.12g}
  a >= {v['certified_decay_lower_bound']:.12g} 1/s
  b <= {v['input_gain_upper_bound']:.12g}
  d <= {v['trim_drift_upper_bound']:.12g}
  a*r-b*epsilon-d >= {v['inward_margin_float']:.12g} > 0
  ||H(x_original-x_c)|| <= {v['initial_original_trim_metric_norm_bound']:.12g} < r
  |fhat_b| <= {v['frequency_bound']:.12g} Hz, todos los buses
  |rhat_b| <= {v['rocof_bound']:.12g} Hz/s, todos los buses
  {v['voltage_min']:.12g} <= |v_b| <= {v['voltage_max']:.12g} pu
La suma exacta de fracciones binarias verifica el signo del margen final.
El centro cambió como máximo {c['maximum_coordinate_offset']:.12g} en una coordenada.
Se representa como valor original más offset separado; redondearlos a una sola
Float64 anularía parte de la mejora. El residuo pasó de 5.11e-9 a 4.18e-25.
Se verificó que el trim original también está dentro de la región calculada.

6. AUDITORÍA Y NIVEL DE GARANTÍA
Treinta puntos contra Julia/las ecuaciones originales dieron errores máximos:
  RHS {a['max_absolute_mapping_errors'][0]:.9g}; tensión {a['max_absolute_mapping_errors'][1]:.9g};
  frecuencia {a['max_absolute_mapping_errors'][2]:.9g}; RoCoF {a['max_absolute_mapping_errors'][3]:.9g}.
El Jacobiano nominal ForwardDiff quedó dentro de la envolvente recomputada.
Se hicieron {ar['checks']} comprobaciones racionales o de 80 dígitos de primitivas de
intervalos. Son tests de implementación, no una prueba independiente universal.

La inversión intervalar usa un inverso aproximado K y una cota E>=|I-KG|,
||E||_infty<1. Una matriz R>=E(|K|+R) encierra |G^-1-K| por la serie de Neumann.
La comprobación espectral usa congruencia por autovectores aproximados y cotas
intervalares de Gershgorin, verificando también su Gram; no supone ortogonalidad
exacta de los autovectores calculados. El verificador recalcula las envolventes,
usa directamente H F_z H^-1 y no reejecuta el optimizador.

Nivel: certificado numérico con envolventes de redondeo bajo el análisis de las
primitivas implementadas, no una prueba mecanizada ni una auditoría independiente
de toda la biblioteca intervalar. El verificador comparte ese evaluador de modelo.
mpmath documenta que su soporte iv es experimental; aquí se usan funciones
elementales y álgebra comprobadas, sin atribuir garantía a iv.findroot.
Documentación primaria: https://mpmath.org/doc/current/contexts.html

7. DÓNDE SE PIERDE UTILIDAD Y QUÉ CAMBIAR
La entrada permitida es del orden de 1e-13 en pu normalizadas: demasiado pequeña
para el problema de sustitución práctica. No se transforma en una cifra de MW,
no se usa para afirmar un máximo robusto del 90.8% y no justifica un póster ganador.
Las cotas de corriente/DC generadas son requisitos del modelo, no ratings de
hardware comprobados. Tampoco hay aquí validación ANDES del nuevo certificado.

La envolvente regional pierde decaimiento cerca de r=4.64e-6. En contraste,
12 puntos de una esfera r=1e-4 conservaron una cota puntual de decaimiento cercana
a 0.05003/s. Este contraste diagnostica conservadurismo, no certifica la esfera
mayor ni excluye estados adversos no muestreados. En particular, una caja de
Jacobianos rompe correlaciones entre entradas de matrices originadas por los
mismos ángulos/corrientes. El condicionamiento de H amplifica esa pérdida.

La siguiente mejora necesaria debe conservar esas dependencias: coordenadas dq
locales exactas y residuos factorizados por dispositivo/puerto; restricciones
sectoriales/multiplicadores de la DAE y estructura de grafo para los términos
cruzados. No basta aumentar muestras ni optimizar rho sobre esta envolvente
microscópica: eso cambiaría el significado práctico del objetivo.

Sólo cuando el certificado cubra una envolvente de operación declarada y útil
se debe ejecutar el co-diseño de rho,Kp,Ki con verificación de cada paso, límites
físicos y comparación contra el ajuste modal. La optimalidad global y la novedad
frente a literatura siguen abiertas. El marco general de disipatividad DAE tiene
antecedentes: https://arxiv.org/abs/2308.08471 ; las herramientas clásicas no se
presentan como teoremas nuevos.

8. REPRODUCCIÓN Y ARCHIVOS
  julia --project=. experiments/nonlinear_codesign_20261001/export_coupled_contract_model.jl
  python experiments/nonlinear_codesign_20261001/audit_coupled_interval_model.py
  python experiments/nonlinear_codesign_20261001/refine_coupled_center.py
  python experiments/nonlinear_codesign_20261001/synthesize_coupled_contract.py --metric modal
  python experiments/nonlinear_codesign_20261001/verify_coupled_contract.py
  python experiments/nonlinear_codesign_20261001/audit_coupled_interval_arithmetic.py
  python experiments/nonlinear_codesign_20261001/diagnose_coupled_wrapping.py
  python experiments/nonlinear_codesign_20261001/audit_coupled_graph_structure.py
  python experiments/nonlinear_codesign_20261001/write_coupled_contract_report.py

Autoritativos: synthesis_modal.json, verification.json, certificate_modal.npz,
model.toml, refined_center.json y sus hashes. synthesis.json y
synthesis_unrecentered_diagnostic.json son etapas anteriores/diagnósticas, no el
resultado final. La tentativa sin recentrar se ejecutó mientras se preparaban
ediciones del código; sus hashes de script no reconstruyen aquella versión y
no se usan como evidencia verificada. Los probes de métricas sólo comparan
opciones descartadas; ningún archivo con status FAIL se promueve a certificado.

Figura: coupled_certificate_conservatism.png/pdf/svg. No se alteró el póster ni
el candidato compartido para presentar este prototipo como capacidad óptima.
"""
(OUT/"CERTIFICADO_ACOPLADO_Y_LIMITES_ES.txt").write_text(report,encoding="utf-8")
manifest={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!="manifest.json"}
(OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
print(json.dumps(dict(report=str(OUT/"CERTIFICADO_ACOPLADO_Y_LIMITES_ES.txt"),status=v["status"],input_radius=v["input_radius"]),indent=2))
