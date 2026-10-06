from pathlib import Path
import hashlib
import json
import tomllib

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/pll_sector_contract"
synthesis_path, verification_path = OUT/"synthesis.json", OUT/"verification.json"
s = json.loads(synthesis_path.read_text(encoding="utf-8"))
v = json.loads(verification_path.read_text(encoding="utf-8"))
n = json.loads((OUT/"nonlinear_validation.json").read_text(encoding="utf-8"))
j = tomllib.loads((OUT/"julia_mapping_audit.toml").read_text(encoding="utf-8"))
sha = lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
assert v["synthesis_sha256"] == n["synthesis_sha256"] == sha(synthesis_path)
assert n["verification_sha256"] == sha(verification_path)
assert j["ReducedDAE_sha256"] == sha(ROOT/"experiments/nonlinear_codesign_20261001/ReducedDAE.jl")
assert s["specification"]["pll_tau"] == j["pll_tau"]
for filename, recorded in (
    ("pll_sector_contract.py", s["script_sha256"]),
    ("verify_pll_sector_contract.py", v["verification_script_sha256"]),
    ("validate_pll_sector_contract.py", n["validation_script_sha256"]),
    ("audit_pll_contract_mapping.jl", j["script_sha256"]),
):
    assert sha(ROOT/"experiments/nonlinear_codesign_20261001"/filename) == recorded
b, c = s["best"], v["certificates"][-1]
base = v["certificates"][0]
p = s["specification"]
improvement = 100*(c["epsilon_verified"]/base["epsilon_verified"]-1)
report = f"""PRIMER CONTRATO NO LINEAL VERIFICADO: PLL + MEDIDA DE ROCOF
2026-10-01. Resultado parcial del método de co-diseño SG -> GFL.

RESULTADO Y ALCANCE

Se construyó un certificado de invariancia para el PLL de tercer orden realmente
instalado en nuestro GFL, aumentado con un instrumento de RoCoF de primer orden.
La prueba cubre todas las señales medibles de una envolvente de entrada, no una
lista de eventos o ubicaciones. Las desigualdades matriciales se comprobaron con
aritmética racional exacta, reconstruyéndolas desde los parámetros escalares.

El certificado es CONDICIONAL a la envolvente de tensión/frecuencia de referencia.
No demuestra que la red mantenga esas entradas, ni certifica el control de
corriente, DC, energía, SG, RoCoF de tensión de bus o porcentaje de sustitución.
El objetivo completo continúa abierto. Tampoco hay demostración de óptimo global
respecto a Kp,Ki; sí una cota rigurosa del subproblema SDP a ganancias fijas.

MODELO EXACTO EN UNA REFERENCIA GIRATORIA

Sea theta_ref una referencia angular absolutamente continua y nu=theta_ref_dot.
En esa referencia, u_dq=[V0,0]+Delta u, con V0 positivo. Definimos
  delta=theta_PLL-theta_ref,
  z=[delta, omega, xi, chi]^T.
omega es desviación de frecuencia física en rad/s: no se elimina la frecuencia
común al cambiar de referencia. El error de fase del PLL es exactamente
  e=-V0 sin(delta)-sin(delta) Delta u_d+cos(delta) Delta u_q.
Las ecuaciones son
  delta_dot = omega-nu,
  tau omega_dot = xi+Kp e-omega,
  xi_dot = Ki e,
  tau_m chi_dot = omega-chi.
theta_ref es una coordenada de análisis; no introduce una señal remota nueva en
la implementación del PLL ni modifica su ley de control instalada.
Las salidas de este contrato son
  f_PLL=omega/(2 pi), r_PLL=(omega-chi)/(2 pi tau_m).
r_PLL es derivada filtrada de la estimación del PLL; no es RoCoF ideal ni la
ventana móvil usada en los ensayos anteriores de red.

ENTRADAS Y SECTOR NO LINEAL

El conjunto de señales se define por la desigualdad puntual
  ||Delta u||^2/epsilon^2 + nu^2/(kappa epsilon)^2 <= 1,
sin restricción de forma temporal. Delta u puede ser discontinua; nu es medible
y acotada, por lo que theta_ref se obtiene por integración. Este contrato no
incluye impulsos en nu o saltos instantáneos de la referencia escogida.

Sea eta=[(-sin(delta)Delta u_d+cos(delta)Delta u_q)/epsilon,
         nu/(kappa epsilon)]^T. Entonces ||eta||<=1.
Para |delta|<=delta_bar<sqrt(6), la desigualdad rigurosa
  1-delta_bar^2/6 <= sin(delta)/delta <= 1
da a=V0 sin(delta)/delta en [a_min,a_max], con el valor continuo en delta=0,
  a_min=Vmin(1-delta_bar^2/6), a_max=Vmax.
No se reemplazó sin(delta) por delta. Se utiliza una inclusión exacta y una
cota del resto: cos(t)>=1-t^2/2 implica, al integrar, sin(t)>=t-t^3/6.

Las ecuaciones completas del subsistema quedan dentro de
  zdot=A(a)z+epsilon B eta,
  A(a)=[[0,1,0,0],[-Kp a/tau,-1/tau,1/tau,0],
        [-Ki a,0,0,0],[0,1/tau_m,0,-1/tau_m]],
  B=[[0,-kappa],[Kp/tau,0],[Ki,0],[0,0]].
La inclusión permite variación arbitraria de a dentro del intervalo, lo cual
puede ser más conservador que la relación seno real, pero contiene esa relación.

CERTIFICADO Y PROBLEMA CONVEXO INTERIOR

Usamos y=S^(-1)z, Atilde=t0 S^(-1)AS, Btilde=t0 S^(-1)B,
alphatilde=t0 alpha, y V=y^T X^(-1)y, con X positivo definido. Las cotas se
normalizan mediante filas C_j, tales que |C_j y|<=1 significa cumplir el límite.
Para ganancias y alpha fijos, se resuelve el SDP
  maximizar epsilon,
  X>=g I, C_j X C_j^T<=1-g, trace(X)<=T, 0<=epsilon<=E,
  [[Atilde X+X Atilde^T+2 alphatilde X, epsilon Btilde],
   [epsilon Btilde^T, -2 alphatilde I]] <= -g I,
para a=a_min y a=a_max.
El test es afín en a: los extremos prueban todo el intervalo, no son muestreo.
S, g, T y E definen la familia numérica del certificado. Las cotas trace(X) y
epsilon facilitan además una cota dual rigurosa; no se atribuyen a un equipo.

PRUEBA. La congruencia diag(X^(-1),I), evaluada en [y;eta], da en tiempo físico
  Vdot <= -2 alpha V+2 alpha||eta||^2 <= 2 alpha(1-V).
Por ello V<=1 es invariante y, sin entradas, V decae exponencialmente. La identidad
max_(y^T X^-1 y<=1)|C_j y|=sqrt(C_j X C_j^T) garantiza las salidas y mantiene
|delta|<=delta_bar, cerrando la validez del sector. En particular,
  V(t)<=exp(-2 alpha t)V(0)+1-exp(-2 alpha t).
Esta cota y las salidas se cumplen para toda señal de la envolvente y todo
V0 en el intervalo, siempre que el estado inicial esté en la elipsoide.

TUNING MEDIANTE DERIVADA DEL CERTIFICADO

Se diferenciaron las matrices del SDP respecto a Kp,Ki. Si el valor es
diferenciable y se cumplen las condiciones de regularidad de sensibilidad,
  d epsilon_star / d k = -sum_a trace(Z_a dM_a/dk),
con Z_a los multiplicadores duales y X,epsilon evaluados en la solución.
Esta derivada alimenta un descenso sobre -epsilon en log(Kp),log(Ki).
No requiere derivar trayectorias de eventos escogidos. Se mantiene una búsqueda
local no convexa en ganancias, aunque el subproblema de X,epsilon sea convexo.

INSTANCIA NUMÉRICA REPRODUCIBLE

Parámetros, no constantes de la teoría:
  V0 en [{p['voltage_min']}, {p['voltage_max']}] pu;
  delta_bar={p['angle_max']} rad;
  tau_PLL={p['pll_tau']:.17g} s; tau_m={p['measurement_tau']} s;
  alpha={p['decay_rate']} s^-1;
  fmax={p['frequency_max_hz']} Hz; rmax={p['rocof_max_hz_s']} Hz/s;
  kappa={p['frame_rate_per_voltage']:.17g} rad/s por pu;
  t0={p['time_scale']}; g={p['matrix_guard']}; T={p['certificate_trace_bound']}; E={p['input_radius_upper_bound']}.

Ganancias nominales:
  Kp={base['gain_pair'][0]:.12g}, Ki={base['gain_pair'][1]:.12g};
  epsilon certificado={base['epsilon_verified']:.12g} pu.
Mejor par encontrado por la búsqueda local:
  Kp={b['Kp']:.12g}, Ki={b['Ki']:.12g};
  epsilon certificado={c['epsilon_verified']:.12g} pu;
  mejora de radio respecto al nominal={improvement:.6f}%.
El radio de frecuencia de referencia asociado es
  kappa epsilon/(2pi)={c['reference_frequency_envelope_hz']:.12g} Hz.
Los radios no son máximos simultáneos independientes: pertenecen a la elipsoide
conjunta especificada arriba.

TAMAÑO REAL DE LA REGIÓN

Aunque el sector se justificó hasta delta_bar={p['angle_max']} rad, la elipsoide
hallada implica |delta|<={p['angle_max']*b['normalized_output_bounds'][0]:.12g} rad
y |f_PLL|<={p['frequency_max_hz']*b['normalized_output_bounds'][1]:.12g} Hz.
El límite de RoCoF filtrado es el restrictivo. La envolvente de tensión es
pequeña ({100*c['epsilon_verified']:.6f}% de 1 pu), por lo que esto NO constituye
todavía un certificado práctico de reemplazo SG->GFL frente a grandes eventos.
La modestia del dominio debe mostrarse, no ocultarse detrás del porcentaje de mejora.

VERIFICACIÓN Y OPTIMALIDAD: TRES AFIRMACIONES DIFERENTES

1. FACTIBILIDAD DEL CONTRATO. Se reconstruyeron las matrices con fracciones
exactas a partir de parámetros y X almacenados. LDL racional verifica X-gI>0 y
las dos matrices negativas con la guarda. También se verifican exactamente
contención de salidas, trace(X) y cota de epsilon. Una cota racional de pi se
comprueba mediante la identidad de Machin y series alternantes de arctan.
El sector demuestra que el certificado cubre el seno original, no sólo dos modelos.

2. ÓPTIMO DEL SUBPROBLEMA A GANANCIAS FIJAS. Multiplicadores duales reparados
para PSD se verifican por LDL racional. Sus residuos de estacionariedad no se
suponen cero: se acotan usando trace(X)<=T y epsilon<=E. Si R y s son esos
residuos, una cota superior válida es
  U=C+T max_i sum_j |R_ij|+E max(s,0).
C es el término constante de la función dual. Así, para el par elegido,
  {c['epsilon_verified']:.17g} <= epsilon_SDP_star <= {c['fixed_gain_guarded_SDP_upper_bound']:.17g},
  brecha relativa <= {c['fixed_gain_guarded_SDP_relative_gap']:.9g}.
Las fracciones exactas se guardan en verification.json. Es optimalidad del
certificado cuadrático y sus restricciones explícitas, no del alcance físico real.

3. GANANCIAS Y CAPACIDAD DE RED. La búsqueda de Kp,Ki es local. No se ha
certificado su máximo global ni se ha optimizado rho con estos contratos.
Los parámetros nuevos NO se aplicaron a los candidatos de red previos.

CONTRASTES INDEPENDIENTES

Julia comparó {j['cases']} evaluaciones de las ecuaciones de este contrato con el
PLL del GFL completo, variando otros estados del convertidor. Error máximo de RHS:
  {j['max_absolute_rhs_error']:.9g}.
Las tensiones nominales de los diez agregados están dentro del intervalo declarado:
  [{j['nominal_voltage_min']:.12g}, {j['nominal_voltage_max']:.12g}] pu.
Las derivadas duales respecto a ganancias concuerdan con diferencias centrales;
error relativo máximo: {max(a['relative_error'] for a in v['gradient_audit']):.9g}.
Se comprobaron {n['random_points']} puntos de estado/entrada y {len(n['trajectories'])}
trayectorias no lineales, incluyendo una señal que maximiza instantáneamente
Vdot y estados iniciales sobre la frontera de RoCoF. Todas permanecieron dentro
del contrato con las tolerancias declaradas. Estos ensayos auditan la implementación;
la cobertura de todas las señales proviene del certificado, no de esos ensayos.

CÓMO ENCAJA EN BEYOND NODAL DAMPING

Este resultado proporciona una primera pieza local exacta y un ejemplo ejecutable
de tuning guiado por un certificado. Aún hace falta derivar contratos de corriente,
DC y SG, y verificar que la interconexión portuaria suministre entradas compatibles.
Imponer cotas independientes por dispositivo puede ser demasiado conservador;
la siguiente formulación debe conservar correlaciones de red mediante suministros
cruzados, bloques de subred o multiplicadores y los operadores eléctricos completos.
La frecuencia común debe seguir presente al cerrar esos puertos.

No debe multiplicarse este radio por una potencia base para anunciar MW reemplazables.
No hay equivalencia automática entre ruido de tensión, déficit de potencia y retiro
de generación síncrona. Tampoco la prueba local establece que el 100% GFL sea viable.

La estabilidad no lineal del PLL mediante Lyapunov y sectores tiene antecedentes.
Esto es una construcción verificada para nuestra arquitectura y una pieza del
método; no se presenta como nueva teoría general de PLL:
  https://arxiv.org/abs/2105.10957
  https://www.monash.edu/__data/assets/pdf_file/0007/3103891/Milad_Zarif_Mansour_JESTPE_2021_.pdf

REPRODUCCIÓN

Desde la raíz del repositorio:
  python experiments/nonlinear_codesign_20261001/pll_sector_contract.py --optimize
  python experiments/nonlinear_codesign_20261001/verify_pll_sector_contract.py
  julia --project=. experiments/nonlinear_codesign_20261001/audit_pll_contract_mapping.jl
  python experiments/nonlinear_codesign_20261001/validate_pll_sector_contract.py
  python experiments/nonlinear_codesign_20261001/write_pll_contract_report.py

Los hashes se verifican antes de generar este reporte. Las pequeñas variaciones
de un solver no se aceptan automáticamente: deben pasar de nuevo el verificador.
"""
(OUT/"CONTRATO_PLL_ES.txt").write_text(report, encoding="utf-8")
manifest = {path.name:sha(path) for path in OUT.iterdir() if path.is_file() and path.name != "manifest.json"}
(OUT/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
print(json.dumps(dict(report=str(OUT/"CONTRATO_PLL_ES.txt"), improvement_percent=improvement,
                      fixed_gain_SDP_gap=c["fixed_gain_guarded_SDP_relative_gap"],
                      status="PLL_ONLY_CERTIFICATE_WITH_EXPLICIT_SCOPE"), indent=2))
