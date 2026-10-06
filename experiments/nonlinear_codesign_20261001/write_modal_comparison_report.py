"""Write the model-agnostic derivation and evidence, without upgrading claims."""
from pathlib import Path
import hashlib
import json
import tomllib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"reports/nonlinear_codesign_20261001"
OUT=BASE/"coupled_dq_contract"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding="utf-8"))
dq=read(OUT/"block_modal_verification.json")
ri=read(BASE/"coupled_network_contract/block_modal_verification.json")
old=read(BASE/"coupled_network_contract/verification.json")
lp=read(OUT/"optimized_comparison_LP_verification.json")
for record,directory in ((dq,OUT),(ri,BASE/"coupled_network_contract"),(lp,OUT)):
    assert all(sha(directory/n)==h for n,h in record["source_hashes"].items())
trajectories=[]
for tag in ("standard","refined"):
    path=OUT/f"nonlinear_validation_julia_{tag}.toml"
    if path.exists():
        data=tomllib.loads(path.read_text(encoding="utf-8"))
        assert data["certificate_sha256"]==sha(OUT/"block_modal_certificate.npz")
        trajectories.append(data)
validation="No hay todavía una campaña temporal completa registrada."
refinement="La auditoría final de sensibilidad a tolerancia está pendiente."
refinement_path=OUT/"trajectory_refinement_audit.json"
if refinement_path.exists():
    audit=read(refinement_path)
    assert all(sha(OUT/name)==value for name,value in audit["source_hashes"].items())
    worst=max(max(row["maximum_scaled_difference"]) for row in audit["cases"])
    refinement=f"Se verificaron seis trayectorias completas; peor diferencia normalizada entre tolerancias = {worst:.10g} (<1e-3). La comprobación se limita a las salidas muestreadas y V; no es una cota rigurosa del error de integración."
if trajectories:
    validation="\n".join(
        f"  {data['tag']}: {data['status']}; tolerancia normalizada {data['normalized_coordinate_abstol_reltol']:.1e}.\n"+
        "\n".join(f"    {r['name']}: completo={r['complete']}, t={r['end_time']:g} s, uso máximo región={r['max_region_usage']:.8g}." for r in data["cases"])
        for data in trajectories)

report=f"""CONTROL NO LINEAL SOBRE GRAFOS: COMPARACIÓN MODAL Y CO-DISEÑO SG -> GFL
2026-10-01. Formulación general y estado verificable de la implementación.

RESULTADO ACTUAL
La garantía se formula para TODAS las señales medibles de una clase acotada,
sin elegir nodos, forma temporal ni frecuencia de perturbación. Se implementó
una región invariante de 184 bloques que conserva los 281 estados de la red.
El radio de entrada pasó de {old['input_radius']:.10g} a {lp['epsilon']:.10g}.
La mejora es de {lp['epsilon']/old['input_radius']:.6g} veces, pero la envolvente
sigue siendo demasiado pequeña para una conclusión de capacidad operativa.
El porcentaje SG/GFL y las ganancias permanecen fijos en este avance.

1. OBJETO GENERAL: UN PROBLEMA DE CONTROL, NO UNA LISTA DE EVENTOS
Sea p=(rho,kappa) el vector de sustitución y control. La DAE index-1 es
  xdot = f(x,v,p,w),       0 = g(x,v,p,w; Y_G).
Su versión matricial es E zdot = F_DAE(z,p,w), E=diag(I,0), z=(x,v).
Se requiere una rama algebraica regular: g_v invertible en el dominio elegido.
La eliminación v=psi(x,p,w) conserva exactamente la DAE en esa rama:
  xdot = F(x,p,w) = f(x,psi(x,p,w),p,w).
No se sustituye F por su linealización. Los instrumentos y controles dinámicos
se añaden a x. Se retira sólo la simetría de ángulo común, no la frecuencia común.

Definimos w=S_w d y W(epsilon)={{d medible: ||d(t)||_2<=epsilon c.t.p.}}.
S_w, epsilon, el dominio inicial y los límites de salida son datos del problema.
La teoría permite cualquier número de buses y canales. No puede garantizarse
un límite físico finito contra entradas arbitrariamente grandes o fuera del modelo.
El conjunto incluye entradas simultáneas, escalones, pulsos y señales cambiantes.
Fallas de líneas, cargas P/Q o incertidumbre paramétrica requieren sus operadores
y conjuntos correspondientes; no se convierten automáticamente en corriente
aditiva exógena. Para incertidumbre theta, las cotas siguientes deben ser uniformes
en theta. Topologías conmutadas requieren regularidad y un certificado común,
incluidas las transiciones, o un análisis separado de tiempos de permanencia.

2. QUÉ APORTAN EL LAPLACIANO Y FOURIER
El operador estructural es L=B diag(a_e) B^T, a_e>0. La red eléctrica usa su
admitancia AC Y_G, con pérdidas, shunts y transformadores cuando correspondan.
L no reemplaza Y_G. L=U Lambda U^T permite coordenadas de Fourier sobre el grafo,
agrupación espacial y parametrizaciones de control. Conservar todos los modos
es un cambio de coordenadas; truncarlos exige acotar el residuo.

En sistemas heterogéneos, U^T no diagonaliza los controles ni las no linealidades.
Por eso usamos y=H(x-x_c), H invertible, y conservamos TODOS los bloques acoplados.
H puede construirse con modos de L, subredes, modos dinámicos o combinaciones.
La implementación actual usa modos reales del Jacobiano dinámico completo y
una transformación exacta de la cascada de instrumentos; no es exclusivamente
una base de Fourier del Laplaciano. La región se verifica después con F no lineal.
El carácter agnóstico se refiere a la formulación; H, Y_G y el diseño dependen
del modelo concreto. Una renumeración de buses debe transformar todos los datos
consistentemente; la selección numérica de bases puede variar en modos repetidos.

3. IDENTIDAD NO LINEAL EXACTA Y TAYLOR CON RESIDUO
Tras centrar las coordenadas escribimos F_c(z,d)=F(x_c+z,p,S_w d).
La identidad integral exacta, en un dominio convexo regular, es
  F_c(z,d)=F_c(0,0)+ integral_0^1 [F_z(tz,td) z + F_d(tz,td) d] dt.
Evaluar/cotar los Jacobianos SOBRE TODO el dominio mantiene la no linealidad.
El Jacobiano nominal sólo selecciona H; no define la garantía regional.
Una expansión de Taylor finita es alternativa si se incluye una cota rigurosa
del residuo de orden superior. Añadir términos sin ese residuo no prueba seguridad.

4. REGIÓN POR BLOQUES Y MATRIZ DE INTERACCIÓN
Particionamos y=Hz en bloques y_k, k=1,...,m. Definimos
  s_k=||y_k||_2,  Omega(r)={{z: s_k<=r_k para todo k}}, r_k>0,
  V(z)=max_k s_k/r_k.
Para encerrar Omega(r), U_jk >= ||fila_j(H^-1[:,bloque k])||_2 da |z|<=U r.
Aquí U es una matriz de cotas, distinta de la base de Fourier citada arriba.

Sobre |z|<=U r y ||d||<=epsilon se verifican uniformemente
  sym[(H F_z H^-1)_kk] <= -a_k I,
  ||(H F_z H^-1)_kj||_2 <= M_kj,  j!=k,  M_kk=0,
  ||(H F_d)_k||_2 <= b_k,   ||(H F_c(0,0))_k||_2 <= delta_k.
Sea D=diag(a_k), M>=0. De la identidad integral y la derivada superior de Dini:
  D^+s <= (M-D)s + b epsilon + delta.                         (A)

CONDICIÓN SUFICIENTE DE INVARIANCIA ROBUSTA:
  (D-M)r > b epsilon + delta,    s(0)<=r,                     (B)
junto con regularidad algebraica y cotas de salidas dentro de sus límites.
En una cara s_k=r_k, (A) apunta estrictamente hacia dentro. Ello prueba que
Omega(r) es positivamente invariante para toda señal de W(epsilon), no sólo
para las simulaciones elegidas. La comparación también proporciona
  s(t) <= exp((M-D)t)s(0) + integral_0^t exp((M-D)(t-tau))(b epsilon+delta)dt.

Si a_k>0, Q=D^-1 M es no negativa y
  varrho(Q) <= max_k (M r)_k/(a_k r_k) < 1.                 (C)
Es una cota de Perron mediante un vector positivo; no depende de aceptar un
autovalor calculado en coma flotante. varrho denota radio espectral, no sustitución.

Además, para dos soluciones con LA MISMA entrada dentro de Omega(r), la misma
envolvente de Jacobianos da contracción incremental en ||H Delta z||_(r,infinito):
  ||Delta z(t)||_H,r <= exp(-gamma t)||Delta z(0)||_H,r,
  gamma = min_k [a_k - sum_j M_kj r_j/r_k] > 0.              (D)
La convexidad de Omega permite aplicar la identidad integral al segmento entre
ambos estados. Para entrada constante, invariancia compacta y contracción dan un
único equilibrio en la región y convergencia exponencial hacia él. El residuo
delta no se elimina artificialmente: x_c puede no ser ese equilibrio exacto.
Para entrada variable se certifica seguimiento incremental acotado, no convergencia
a un equilibrio fijo. La tasa y el dominio son locales/regionales.

Éste es el paso beyond nodal damping: a_k por separado no basta. Importan M_kj,
la orientación de las entradas b_k, los límites de salida y los recursos físicos.
Los bloques actuales son modales, no contratos independientes por nodo.
Comparación, small gain e ISS son herramientas establecidas [1], no novedad por sí solas.

5. ALGORITMO INTERNO: PUNTO FIJO Y UN LP CON COTA DUAL
Un iterador propone r <- c max((D-M)^-1(b epsilon+delta),s(0)), c>1,
recalcula las envolventes en el nuevo dominio y acepta sólo tras verificar (B).
No se presupone que el iterador converja. Su fallo no demuestra inviabilidad física.

Fijados p, H, una caja Xbar y un límite exterior epsilon_bar para calcular a,M,b,
el siguiente problema en (r,epsilon) es LINEAL:
  maximizar epsilon
  sujeto a (D-M)r >= b epsilon+delta+margen,
           U r <= Xbar, r>=s(0), 0<=epsilon<=epsilon_bar.    (E)
Las salidas se acotan sobre Xbar x W(epsilon_bar). Una solución de (E) hereda
esas cotas y la regularidad. Cambiar Xbar, H o p obliga a reconstruir/verificar.
Un dual y sus residuos permiten una cota superior del óptimo de ESTE LP.
No es una cota superior de la capacidad física SG->GFL ni de todos los certificados.

6. CO-DISEÑO EXTERIOR Y CAPACIDAD PARAMETRIZADA
Para una envolvente REQUERIDA epsilon_req, declarada antes de buscar sustitución:
  C_cert(epsilon_req, limites, recursos) =
    sup sum_i P_i^base rho_i / sum_i P_i^base
    sobre rho,kappa,H,r,Xbar sujetos a (B), salidas y restricciones físicas. (F)
El certificado es suficiente: C_cert <= C_fisica para idénticos supuestos,
pero una búsqueda local sólo produce una cota alcanzada <= C_cert.
Optimizar simultáneamente epsilon hasta volverlo microscópico no resuelve (F).
Una frontera C_cert(epsilon_req) muestra explícitamente el compromiso robustez/capacidad.

Para g=0, las sensibilidades exactas de la rama son
  psi_p = -g_v^-1 g_p,   F_p = f_p - f_v g_v^-1 g_p.       (G)
En esta interconexión, a estados fijos,
  g_rho_i=S_i(i_GFL_i-i_SG_i),
  v_rho_i=-G^-1 S_i(i_GFL_i-i_SG_i).
Los términos directos f_p deben mantenerse cuando existan. Para el PLL instalado,
omega_PLL_dot=(xi+Kp e_q-omega_PLL)/tau_PLL, xi_dot=Ki e_q;
sus derivadas directas son e_q/tau_PLL y e_q en sus filas respectivas.
La derivada del centro/equilibrio y de H se incorpora si se actualizan con p.

Una búsqueda exterior razonable usa SQP/región de confianza con H y dominio
congelados durante cada propuesta, seguida de verificación no lineal completa.
Para m_k=a_k r_k-sum_j M_kj r_j-b_k epsilon-delta_k,
  dm_k=r_k da_k+a_k dr_k-sum_j(r_j dM_kj+M_kj dr_j)
       -epsilon db_k-b_k d epsilon-d delta_k.               (H)
Las cotas intervalares/normas/máximos pueden ser no diferenciables; se requieren
derivadas direccionales, subgradientes o un modelo local conservador. No se debe
confundir el gradiente del Jacobiano nominal con el de una cota regional.
Si el Perron dominante es simple, d varrho(Q)=u^T(dQ)v/(u^Tv) ayuda a dirigir
el ajuste; no reemplaza (B). Degeneraciones exigen subespacios/bundles.

Esta búsqueda exterior robusta TODAVÍA NO está implementada/validada en esta etapa.
No se propone EM: no hay aquí variables latentes ni una verosimilitud cuyo EM
garantice mejora. La estructura natural es comparación positiva + LP + SQP.
Los extremos rho=0/1 requieren eliminar estados de dispositivos retirados y
revisar la carta angular. El backend actual usa 0<rho_i<1 en los diez agregados.
No debe extrapolarse hasta una red 100% GFL sin verificar su referencia de tensión.

7. CONTROL LOCAL, VECINOS, FRECUENCIA Y ROCOF
El diseño puede admitir un control dinámico local/distribuido, por ejemplo
  P_i^cmd=P_i^0-k_fi fhat_i-k_ri rhat_i-zeta_i,
  zeta_i_dot=k_ii fhat_i + sum_j a_ij c_ij(fhat_i-fhat_j),
con filtros, retardos, saturación, antiwindup y reservas explícitos cuando existan.
Estos nuevos estados y parámetros se añaden a F y a kappa. Un filtro polinómico
h(L) de grado q usa hasta q saltos de comunicación; su realización y retardos
deben respetar la arquitectura disponible. H se usa para certificar y no obliga
a implementar un controlador que mida todos los modos globalmente.

Los Kp,Ki actuales son los del PLL; NO equivalen a droop, inercia sintética o PI
de soporte de frecuencia. La ley anterior es una extensión propuesta, no probada.
Un GFL necesita potencia/energía DC y margen de corriente para prestar ese soporte.
Más ganancias no crean esos recursos.

Para señales de entrada sólo medibles, una derivada ideal de fase puede ser
ilimitada o no existir. Hay que declarar instrumentos causales de f y RoCoF,
o restringir regularidad/banda de la entrada. Este ensayo usa filtros de 0.1 s,
no una ventana de RoCoF de 0.5 s ni una derivada ideal. Cambiar la instrumentación
cambia el problema y exige otra verificación. Las cotas de 0.5 son ejemplos
del benchmark, no constantes de la teoría.

8. CAMBIO dq EXACTO Y ABLACIÓN
En cada GFL, i_r=cos(theta)i_d-sin(theta)i_q,
                 i_i=sin(theta)i_d+cos(theta)i_q.
Los estados de corriente locales satisfacen exactamente
  iqdot=omega_b/Xf (v_iq-u_q-Rf iq)-(omega_b omega_frame+omega_PLL) id,
  iddot=omega_b/Xf (v_id-u_d-Rf id)+(omega_b omega_frame+omega_PLL) iq.
La velocidad de la carta común se cancela en estas dos ecuaciones; permanece
en theta_dot=omega_PLL-nu y en los instrumentos. No se eliminó la dinámica PLL.
El push-forward de Julia coincide con el modelo RI a 8.74e-12 en los ensayos.

Mismo candidato físico; misma clase y normalización de entradas:
  Una métrica escalar RI: epsilon={old['input_radius']:.12g}.
  Comparación por bloques RI: epsilon={ri['epsilon']:.12g}.
  Comparación por bloques dq: epsilon={dq['epsilon']:.12g}.
  LP con caja dq fija: epsilon={lp['epsilon']:.12g}.
El cambio dq aporta {100*(dq['epsilon']/ri['epsilon']-1):.3f}% respecto del mismo
método por bloques RI. La mayor parte de la mejora procede de separar escalas
modales; no se atribuye al cambio dq ni a una supuesta mejora física de la red.
La comparación escalar combina la tasa más lenta y la máxima sensibilidad de
controles rápidos, causando una pérdida muy grande en este ejemplo.

9. EVIDENCIA VERIFICADA Y LO QUE AÚN FALTA
LP: 184 desigualdades de frontera recomprobadas con aritmética racional sobre
coeficientes encerrados por intervalos, sin violaciones.
  epsilon={lp['epsilon']:.12g}; margen mínimo={lp['min_boundary_margin']:.8g}.
  Cota Perron={lp['comparison_perron_upper']:.12g}<1.
  Tasa regional de contracción >= {lp['contraction_rate_lower']:.12g} /s.
  Brecha relativa del LP normalizado <= {lp['LP_relative_gap_upper']:.12g}.
  |fhat|<={lp['frequency_bound']:.12g} Hz,
  |rhat|<={lp['rocof_bound']:.12g} Hz/s,
  V en [{lp['voltage_min']:.12g},{lp['voltage_max']:.12g}] pu.
La cota dual corresponde a la tabla binaria normalizada guardada. Además se
reverificaron las desigualdades físicas sin normalizar y la contención del dominio.
Las dos verificaciones comparten el backend intervalar; no son una prueba formal
en un asistente de teoremas ni una implementación matemática completamente ajena.

Julia: campo no lineal independiente, Jacobiano por Schur contra ForwardDiff,
error máximo de suma por fila 5.83e-11 en el centro. Las trayectorias usan el
certificado dq ANTERIOR al LP, epsilon={dq['epsilon']:.12g}; no se atribuyen al LP.
{validation}
{refinement}
La campaña inicial a tolerancia 2e-8 agotó 180 s de cálculo por caso alrededor
de t=10--11 s. Sus trazas se conservan como incompletas, no como inestabilidad.
Las simulaciones posteriores no modifican ni sustituyen la prueba de (B).

Pendiente para una afirmación operativa: envolvente de entrada útil y calibrada,
co-diseño exterior rho/ganancias, ratings reales de corriente y tensión DC,
energía/reserva, saturación, incertidumbres y contraste de modelos (incluido ANDES).
Una prueba de óptimo GLOBAL físico requiere una cota superior válida de ese
problema; el dual de (E) sólo acota el subproblema de certificados congelado.
El candidato con 90.785% GFL ya ensayado no es el máximo robusto demostrado.

10. CONTRIBUCIÓN Y PÓSTER
Hipótesis defendible: explicar y utilizar el presupuesto de interacción entre
modos para co-diseñar cuánto SG retirar y cómo ajustar cada GFL, conservando una
región segura no lineal para toda una clase de perturbaciones.
La figura clave sería una frontera capacidad--robustez con ganancias y restricción
limitante identificadas, comparada con damping nominal y contratos independientes.
Todavía no tenemos esa frontera útil: los resultados actuales prueban el mecanismo
matemático y cuantifican conservadurismo. No se promete novedad exhaustiva ni premio.
Título de trabajo: Beyond Nodal Damping: Nonlinear Interaction Budgets for SG-to-GFL Replacement.

ANTECEDENTES PRIMARIOS
[1] Dashkovskiy, Rüffer, Wirth. Small gain theorems for large scale systems and
construction of ISS Lyapunov functions. SIAM JCO 48(6), 4089--4118, 2010.
https://arxiv.org/abs/0901.1842 ; https://doi.org/10.1137/090746483
Respalda composición ISS y small gain; no nuestra capacidad SG/GFL ni prioridad.
[2] Jensen et al. Certifying Stability and Performance of Uncertain Differential-
Algebraic Systems: A Dissipativity Framework. https://arxiv.org/abs/2308.08471
Antecedente de certificados para DAE e incertidumbre; L2 no equivale a límite pico.
[3] Shuman et al. The Emerging Field of Signal Processing on Graphs.
https://arxiv.org/abs/1211.0053
Antecedente Fourier/filtros de grafo; no garantiza desacoplamiento no lineal.

REPRODUCIBILIDAD
Desde la raíz, con Python numpy/scipy/mpmath/threadpoolctl y Julia del proyecto:
  python experiments/nonlinear_codesign_20261001/verify_modal_comparison.py --coordinates dq
  python experiments/nonlinear_codesign_20261001/verify_modal_comparison.py --coordinates ri
  python experiments/nonlinear_codesign_20261001/optimize_comparison_envelope.py
  python experiments/nonlinear_codesign_20261001/verify_comparison_envelope.py
  julia --project=. experiments/nonlinear_codesign_20261001/validate_modal_comparison_julia.jl
  # Segunda ejecución: MODAL_TOL=2e-7, MODAL_TAG=refined (variables de entorno).
  python experiments/nonlinear_codesign_20261001/audit_modal_trajectories.py
  python experiments/nonlinear_codesign_20261001/write_modal_comparison_report.py
El LP conserva el modelo y certificado iniciales. La verificación guarda hashes.
Los scripts de síntesis y exportación completos están en la misma carpeta.
"""
(OUT/"CONTROL_MODAL_COMPARACION_Y_CODESIGN_ES.txt").write_text(report,encoding="utf-8")

values=[old["input_radius"],ri["epsilon"],dq["epsilon"],lp["epsilon"]]
labels=["Single metric\nRI","Modal blocks\nRI","Modal blocks\ndq","Frozen-domain LP\ndq"]
fig,ax=plt.subplots(1,2,figsize=(12,4.8),gridspec_kw={"width_ratios":[1.2,1]})
fig.subplots_adjust(left=.08,right=.98,bottom=.22,top=.78,wspace=.34)
fig.suptitle("Nonlinear interaction bounds: separate the modes, retain the couplings",x=.08,ha="left",y=.965,fontsize=14,fontweight="bold")
fig.text(.08,.865,"Same SG/GFL design · 281 states · 88 input channels · all bounded measurable inputs",fontsize=10,color="#435461")
ax[0].scatter(range(4),values,c=["#768491","#197c80","#197c80","#d75b32"],s=58,zorder=3)
ax[0].plot(range(4),values,color="#aab9bf",lw=1,zorder=1)
ax[0].set(yscale="log",ylabel="Certified normalized input radius",xticks=range(4),xticklabels=labels,ylim=(2e-14,1e-5),xlim=(-.35,3.35))
for i,val in enumerate(values):ax[0].annotate(f"{val:.2e}",(i,val),xytext=(0,10),textcoords="offset points",ha="center",fontsize=9)
ax[0].set_title("A  Most improvement comes from modal comparison",loc="left",fontsize=10,pad=13)
completed=[data for data in trajectories if data["status"]=="THREE_TRAJECTORY_CHECKS_PASSED" and len(data["cases"])==3]
if completed:
    tag=completed[-1]["tag"]
    for name,color in (("origin_waves","#197c80"),("origin_switches","#e2a445"),("boundary_waves","#d75b32")):
        tr=np.loadtxt(OUT/f"trajectory_{tag}_{name}.csv",delimiter=",")
        ax[1].plot(tr[:,0],tr[:,1],label=name.replace("_"," "),color=color,lw=1.5)
    ax[1].axhline(1,color="#555",ls="--",lw=1,label="Region boundary")
    ax[1].set(xlabel="Time (s)",ylabel=r"$V=\max_k\|y_k\|/r_k$",ylim=(0,1.08))
    ax[1].legend(frameon=False,fontsize=8,loc="upper right")
    ax[1].set_title("B  Finite Julia checks of the pre-LP region",loc="left",fontsize=10,pad=13)
else:
    cert=np.load(OUT/"optimized_comparison_LP.npz")
    usage=(cert["M"]@cert["radii"])/(cert["a"]*cert["radii"])
    ax[1].plot(np.sort(usage),color="#197c80",lw=1.8)
    ax[1].axhline(1,color="#555",ls="--",lw=1)
    ax[1].set(xlabel="Modal block (sorted)",ylabel="Internal coupling / decay budget",ylim=(0,1.08))
    ax[1].set_title("B  Couplings remain below each decay budget",loc="left",fontsize=10,pad=13)
for axis in ax:
    axis.spines[["right","top"]].set_visible(False);axis.grid(alpha=.15);axis.set_axisbelow(True)
fig.text(.08,.092,"The envelope remains too small for an operating-capacity claim. The LP gap concerns one fixed certificate problem.",fontsize=9,color="#435461")
fig.text(.08,.043,"Regional nonlinear bounds provide the guarantee; a failed certificate or an incomplete simulation does not prove instability.",fontsize=9,color="#435461")
for ext in ("png","pdf","svg"):fig.savefig(OUT/f"modal_interaction_certificate.{ext}",dpi=180,facecolor="white")
manifest={str(p.relative_to(ROOT)):sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!="manifest.json"}
for name in ("THEORY_GRAPH_ROBUST_CONTROL.txt","METHOD_SPEC.txt","graph_control_identity_audit.json"):
    p=BASE/name;manifest[str(p.relative_to(ROOT))]=sha(p)
scripts=["coupled_dq_model.py","modal_comparison_contract.py","verify_modal_comparison.py","optimize_comparison_envelope.py","verify_comparison_envelope.py","validate_modal_comparison_julia.jl","validate_modal_comparison_julia_timeout.jl","CoupledDQReference.jl","CoupledDQJacobian.jl","write_modal_comparison_report.py","audit_modal_trajectories.py"]
for name in scripts:
    p=ROOT/"experiments/nonlinear_codesign_20261001"/name;manifest[str(p.relative_to(ROOT))]=sha(p)
(OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
print(json.dumps(dict(report=str(OUT/"CONTROL_MODAL_COMPARACION_Y_CODESIGN_ES.txt"),artifacts=len(manifest),trajectory_runs=len(trajectories)),indent=2))
