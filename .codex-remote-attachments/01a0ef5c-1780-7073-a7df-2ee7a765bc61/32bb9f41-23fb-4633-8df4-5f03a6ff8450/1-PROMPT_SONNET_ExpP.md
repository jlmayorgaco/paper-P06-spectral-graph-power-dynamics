# PROMPT PARA SONNET — ExpP incremental
## Co-diseño algebraico SG→GFL sobre el modelo PD-exact de ExpN

### Misión y criterio de trabajo

Trabaja en el repositorio `paper-P06-spectral-graph-power-dynamics` disponible en mi PC. Lee primero este documento, `METODO_COMPLETO_SG_GFL_ExpP.md`, y los artefactos y código de ExpN. Implementa y ejecuta una ampliación incremental de ExpN, NO una reconstrucción desde cero.

Objetivo: calcular y validar candidatos Z=(rho,Kp,Ki) que maximicen SG→GFL para el objetivo declarado. Queremos reemplazar búsquedas gruesas por fronteras algebraicas exactas, corrección continua conjunta, robustez explícita y reducción GSP comprobada. No prometas óptimo global: certifícalo solo si realmente cierras cotas.

No modifiques los artefactos congelados ni candidatos de A–N. No hagas push ni commit sin autorización. Preserva cambios preexistentes. Crea una rama o worktree aislado `research/expP-algebraic-frontiers-robust-codesign` sin perder cambios locales; si ya existe, continúa de forma segura y no sobreescribas archivos congelados. No ejecutes git reset/clean.

Crea, sin duplicar el modelo validado:
- `experiments/bnd_expP/`
- `src/bnd_design_p/` como adaptador y nuevas herramientas, importando funciones de ExpN;
- `test/bnd_expP/`
- `reports/experiment_P/` con subdirectorios por fase;
- un único ejecutor con `--stage P0`, ..., `--stage P6`, `--resume` y `--through P5`.

Ejecuta secuencialmente los módulos obligatorios. Una identidad falsa bloquea las ramas que dependan de ella, NO autoriza cambiar la planta o esconder el fallo. Un módulo opcional GSP/global inconcluso no borra resultados ya válidos. Guarda progreso tras cada fase y respeta presupuestos de tiempo/evaluaciones registrados antes de ejecutarla.

## 1. Punto de partida obligatorio

ExpN reporta:
- Julia 1.11.9, PowerDynamics 5.0.0, NetworkDynamics 1.3.0;
- model SHA `e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a`;
- IEEE39 de la versión instalada: 100 MVA, 60 Hz;
- generación original 5402.761089978847 MW;
- soporte SG {38}; rho38=0.9986453680992127, resto rho=1;
- Kp_i=7.853981633974483, Ki_i=986.9604401089358;
- SG retenida 1.1243444776535036 MW;
- alpha analítica −0.0500000011050941 s^-1;
- alpha PD −0.050000000138500054 s^-1;
- certificado nominal local dentro de soporte {38}, NO global;
- seis pulsos pequeños validados, NO robustez ni capacidad de gran disturbio validada.

Estos números son REGRESIONES, no constantes que debas escribir como resultados calculados. No impongas bus 38 ni ganancias en la esquina para los nuevos problemas. Sí usa ExpN abiertamente como incumbent y semilla: no hay prueba ciega respecto de resultados que ya conocemos.

El informe completo está en `reports/experiment_N/REPORT_EXP_N.md`. Inspecciona APIs reales de `src/bnd_model_expN/`, `src/bnd_opt_expN/` y scripts existentes antes de crear envoltorios. Si las rutas difieren, descubre los archivos y registra las rutas efectivas.

## 2. Contrato físico innegociable

Por generador, no por inyección neta de barra:

    eps_i = 1-rho_i
    S_SG_i = eps_i * S_gen_i0
    S_GFL_i = rho_i * S_gen_i0
    V_i = V_i0 en el nominal

Conservar por separado las cargas inicializadas de 31 y 39. Conservar el trim cerrado de ExpN para `iset_q`, `P_dc`, PLL, filtro, integradores CC/DC. No permitir que initialize_from_pf redistribuya P/Q entre dispositivos.

SG: Sn=eps*Sn0, H en segundos fijo; HSn se escala UNA sola vez.
GFL: agregado paralelo de ExpN, corriente externa=rho*corriente interna equivalente, con estados normalizados invariantes bajo ese contrato. No volver a los modelos mixtos E/H/K anteriores a la corrección P/Q.

No cambies signos dq/DC para que parezcan convencionales: coteja ecuaciones fuente congeladas. En SimpleGFLDC v5 el código alimenta corriente activa en d; no copies un comentario de documentación que diga q si contradice el código.

No corrijas supuestas bases 50/60: ExpM ya descartó esa causa. Verifica bases y Vbase por regresión sin reconstruir la teoría a partir de esa sospecha descartada.

## 3. Objetivos y dominios que debes congelar

Primario nominal, comparable con ExpN:

    minimize J = sum(P_gen_i0 * eps_i)
    all physical finite poles: Re(lambda)<=-0.05 s^-1
    Kp/Kp0, Ki/Ki0 ∈ [0.25,4], independientes.

No llamar a J capacidad instalada. Para un estudio de MW instalados se necesitan ratings activos C_i_MW verificados; crea un objetivo separado solo si existen. No convertir Sn MVA a MW por suposición.

Normaliza variables y restricciones antes de evaluar KKT. Distingue los dos problemas:
- frontera matemática con sigma_req=0.05;
- candidato ejecutable con guard numérico explícito adicional.

Propuesta de guard para NUEVOS candidatos: delta_num=max(1e-6 s^-1,10*error_abscisa_estimado). Congélalo y reporta su efecto en MW. No reinterpretar ese guard como robustez física. ExpN se reproduce con su definición histórica; no falsificar la regresión porque ahora pedimos un guard mayor.

No introduzcas 100 MW, 0.5 Hz o beta como normas universales. Los requisitos heredados se leen de sus configs y se etiquetan como escenarios del proyecto.

## P0 — Regresión breve, no repetir toda la auditoría M

1. Verifica hashes, versión y árbol de código de ExpN; no actualices paquetes. Si falta una biblioteca opcional de certificación usa un entorno separado, nunca mutar el Manifest congelado.
2. Reproduce el candidato ExpN con el modelo analítico. Conserva las tablas de identidad ya existentes y reejecuta las comprobaciones pertinentes.
3. Realiza una construcción PD independiente del MISMO trim para el candidato y uno de los casos mixtos withheld de ExpN, no una malla masiva.
4. Verifica P/Q por componente, signos, pu, mismo espectro físico, misma clasificación de gauge. No borrar polos pequeños por tolerancia.
5. Confirma que el diseñador no evalúa PowerDynamics por candidato: usar contadores/trazas, no solo buscar cadenas `using PowerDynamics`.

Gates orientativos para identidad:
- P/Q <1e-8 pu;
- puerto GFL <1e-8 relativo;
- A reducida <1e-8 relativo;
- matching espectral completo <1e-6 s^-1 con condicionamiento documentado;
- reproducción del modo y del coste de ExpN sin hardcode.

Mantener tolerancias más estrictas de ExpN cuando resulten reproducibles. Si P0 falla: detener todas las etapas de diseño, identificar el motivo y no cambiar ganancias para ocultarlo.

Entregable: `P0_BASELINE.json`, inventario de funciones reutilizadas y tabla corta de regresión.

## P1 — Validar las identidades algebraicas antes de usarlas

Lee las ecuaciones (26)–(33) del documento matemático.

### P1A. Dependencia de un par de ganancias

En una arquitectura fija y con equilibrio ExpN, demuestra desde las ecuaciones:

    P(s,k) = P_ref(s) - (dKp_i*h_p_i+dKi_i*h_I_i)*q_i'

si la estructura realmente lo permite. Incluye los cambios en A y B del puerto y las transformaciones de ecuaciones/estados. E debe ser fijo respecto de ese par o deben incluirse sus términos. No confundir rank-one por dispositivo con rank-one conjunto.

Deriva:

    p = p_ref + dKp_i*p_p_i + dKi_i*p_I_i.

Usa referencia admisible y determinantes escalados coherentemente, LU/bordered solves. No computar adjugadas densas o coeficientes monomiales enormes por defecto.

Pruebas reales:
- una prueba simbólica pequeña, de dimensiones controladas;
- pruebas de identidad en los diez GFL a parámetros no usados para extraer coeficientes;
- dos puntos de frecuencia cerca de modos lentos y otros alejados;
- probar explícitamente que el término cruzado MISMO dispositivo dKp_i*dKi_i es cero bajo la hipótesis;
- mostrar que términos entre dispositivos sí pueden existir;
- medir rank residual y diferencias de pencil, no solo coincidencia de un polo.

Si no es rank-one, conserva un update de rango real demostrado y etiqueta AFFINE_PAIR_UNSUPPORTED. No regreses a un barrido para fingir la fórmula.

### P1B. Retención de un SG

Convención de corriente: T=Ynet+YL−ΣYinj. Define Di=YF_i−YSG_i y

    T=T_F+eps_i*Pi_i*Di*Pi_i'
    Ni=Di*Pi_i'*(T_F\Pi_i)
    det(I2+eps_i*Ni)=1+trace(Ni)*eps_i+det(Ni)*eps_i^2.

Verifica operador/descriptor completo contra identidad y puertos ExpN. Contabiliza polos internos y cancelaciones; `C=I+Delta*G` no muestra los polos del sistema base cuando Delta=0.

Prueba dos puertos con determinante 4x4 y términos cruzados. No sumar dos soluciones de un puerto.

### P1C. Fronteras PI con geometría correcta

Para omega>0 y Q invertible, resuelve el sistema real 2x2 de la ecuación (31).
Para omega=0 usa una ecuación real/recta; NO inviertas Q singular.
Conserva degeneraciones, bordes del gain box y pérdida de regularidad/grado. Una malla de omega construye una aproximación geométrica, no demuestra exhaustividad.

Gates: identidades de operador aproximadamente 1e-10 o mejores donde el condicionamiento lo permite, más error absoluto para cantidades casi nulas. Muestra frecuencia/condición y casos excluidos; no excluir todos los casos difíciles.

Entrega `P1_IDENTITIES.md`, pruebas pequeñas, coeficientes/representaciones evaluables y tabla de errores. Esta fase no declara ningún nuevo óptimo.

## P2 — Recuperar analíticamente la frontera conocida

Con soporte SG {38} y ganancias de ExpN FIJAS, evalúa Ni(s) en s=−0.05. Obtén todas las raíces físicas mediante el cuadrático/eigenvalues de Ni. Reproduce eps38 de ExpN y el coste, dentro de tolerancias condicionadas y del guard histórico.

Verifica cada raíz contra TODOS los polos físicos, no solo el que satisface la ecuación. También calcula la raíz para sigma_req+delta_num y reporta el pequeño coste del guard por separado.

Luego aplica la misma ecuación a los otros nueve soportes de un SG con esas mismas ganancias; son diez evaluaciones algebraicas, no una nueva rejilla de rho. Clasifica raíces inválidas, otros polos limitantes, polos de bloques y falta de raíces.

Prueba las curvas PI condicionadas para un dispositivo con otros parámetros fijos; cruza suavemente la frontera y comprueba el conteo de polos, incluidos cruces reales. No afirmar que los diez pares pueden optimizarse separadamente sin interacción.

Resultado buscado: `EXACT_CONDITIONAL_RETENTION_ROOT`, no `GLOBAL_ZSTAR`.
Si la raíz no reproduce ExpN, diagnostica signo/rango/cancelación antes de seguir. No modificar el punto ExpN para hacerla coincidir.

## P3 — Co-diseño nominal continuo real

Implementa un algoritmo que mueva eps y ganancias simultáneamente. No sustituirlo por 0.02-step grids ni únicamente barridos de gain patterns. Tampoco presentar la proyección de mínima norma de un residual como una minimización de J.

Usa los updates de rango pequeño de P1 para acelerar evaluaciones y P2/curvas PI para semillas o correcciones. La factibilidad la decide el espectro completo.

Dentro de una arquitectura fija:
- derivadas Tz totales, incluidos trim si cambian y Schur/network;
- todas las restricciones espectrales cercanas a activarse;
- sensibilidad de polo simple, no de un eigenvalor ordenado discontinuamente;
- subespacios/ramas para colisiones, y eventos de cambios de soporte;
- normalización de eps, Kp, Ki y multiplicadores;
- paso QP/SQP de active set propio, curvatura positiva regularizada y control de confianza;
- sistema KKT/merit line search implementado y probado;
- KKT final con multiplicadores de bounds, dualidad, complementariedad, LICQ y critical cone/SOSC.

Límite el ámbito y el coste explícitamente: primero soporte {38} con ganancias libres (incluye el incumbent), después vecinos/máscaras capaces de mejorar el coste. Registrar las ramas abiertas. No afirmar que un punto local por cada máscara agota todas las soluciones. No esperar 1–3 iteraciones por decreto.

Usa ExpN como mejor punto nominal conocido. Si termina peor bajo exactamente el mismo dominio, conserva ExpN como incumbent y reporta fallo de mejora; no declarar un nuevo máximo inferior.

### GSP como comparación incremental

Familias:
1. dos ganancias comunes (contiene exactamente el candidato de ExpN);
2. features de un operador;
3. features multigrafo no conmutativas;
4. veinte ganancias independientes como referencia.

Construye operadores sobre las 39 barras físicas y luego restringe filas a generadores si prometes localidad física. El grafo Kron denso no representa comunicación 1-hop.
Usa features exógenas congeladas y orden creciente. Reporta rango/coeficientes: una base de rango 10 no es reducción sobre diez nodos. Presupuesto sugerido 1,2,3,4 features nodales por canal; no aumentar hasta obtener trivialmente 100% de captura sin reportarlo.

No introduzcas señales remotas e_j(t) en los PLL. Son ganancias diseñadas offline. Una arquitectura con comunicación es OTRO experimento.

Comprueba gradiente y direcciones excluidas con el modelo libre. Si hay dirección factible de mejora omitida, amplía la base o declara pérdida de optimalidad estructural. No transferir un certificado local restringido al espacio libre.

Entrega `P3_NOMINAL_CANDIDATES`, KKT ledger, trayectoria J/residuos por iteración, llamadas, tiempo, memoria, reducciones evaluadas y comparaciones equivalentes. No usar simulación PD en el algoritmo.

## P4 — Robustez y seguridad transitoria sobre el mismo modelo

### P4A. Definir la incertidumbre

Leer las configs heredadas y documentar Bw,Cw,Dw, unidades, escalas y tipo de Delta. Reusar el beta normalizado de ExpG solo con la MISMA definición y normalización; no es porcentaje físico ni comparación automática de modelos diferentes.

Para el núcleo usar incertidumbre estática A+BDeltaC o LFT exacta declarada. No certificar incertidumbre de punto operativo congelando arbitrariamente su equilibrio. Variaciones de cargas/voltajes que reequilibran se tratan con ecuaciones consistentes; muestras son escenarios, no certificado continuo.

Si no hay un conjunto de incertidumbre inequívoco, marcar la rama robusta BLOCKED_UNCERTAINTY_DEFINITION y conservar el nominal. No inventar rangos físicos.

### P4B. Certificar H-infinity de forma correcta

Shift A_sigma=A+sigma_req*I en el quotient físico.
No usar el bound modal conservador de ExpG como restricción principal.

Para un Z fijo obtener gamma_L <= ||G_sigma||inf <= gamma_U:
- gamma_L: pico encontrado numéricamente;
- gamma_U: certificado bounded-real/Hamiltoniano con residual verificable o inclusión validada.

Un máximo en una malla es una cota inferior de la norma. Su inverso NO es beta certificado. Reportar beta_cert=1/gamma_U. Manejar feedthrough y well-posedness si Dw no es cero. Para incertidumbre dinámica y decay shift se requieren hipótesis ponderadas adicionales; no transferir automáticamente el resultado estático.

Se permite solver SDP para un certificado de Z fijo; no presentarlo como el algoritmo de diseño ni prohibir herramientas científicas por ser bibliotecas. La optimización primaria sigue usando nuestras ecuaciones/sensibilidades.

### P4C. Co-diseño robusto y transitorio

Agregar la restricción robusta al KKT, con derivadas singulares y múltiples frecuencias activas si es necesario. Leer beta_req existente; presentar varios beta como estudios separados, sin cambiar el requisito para conseguir PASS.

Agregar RoCoF/frecuencia desde el comienzo cuando el escenario de diseño los exige. Reutilizar las definiciones existentes de limites y eventos, etiquetadas como ilustrativas del proyecto. Si se exige 100 MW, una solución que viola ese evento NO es factible aunque sus polos pasen.

Definir entrada en MW y perfil (escalón o pulso con duración), lugar y salida:
- frecuencia COI de SG;
- frecuencia/PLL local cuando esté disponible;
- ventana o filtro RoCoF;
- feedthrough y posibles saltos.

Calcular transitorios lineales por exponencial aumentada o realización equivalente; picos por búsqueda de raíces y endpoints con tail bound. No equiparar un H-infinity L2 con un límite L-infinity de pico de un escalón.

No imponer la fórmula inicial f0 DeltaP/(2 sum eps HSn) sin comprobar que coincide con Cf B en este modelo detallado. No extrapolar capacidad de pulso a escalón sostenido.

Entrega tres resultados claramente separados si son factibles:
- NOMINAL;
- ROBUST_LINEAR;
- ROBUST_WITH_DECLARED_TRANSIENT_CONSTRAINTS.

Conserva Kp/Ki libres como referencia. Compara la familia uniforme y GSP bajo exactamente los mismos requisitos. Un ahorro de parámetros sin pérdida demostrada de MW es un resultado, aunque no mejore el porcentaje nominal.

## P5 — Congelación y validación independiente final

Para cada candidato aceptado:
1. Congelar Z, modelo, config de incertidumbre, perturbaciones, referencias de trim y SHA.
2. Construir PowerDynamics independientemente con el MISMO P/Q por dispositivo. No permitir re-asignación libre ni retuning posterior.
3. Comparar equilibrio, puertos, matrices, TODOS los polos, gauge físico y la rama limitante; llevar errores a presupuestos numéricos explícitos.
4. Repetir los seis pulsos pequeños de ExpN, con sus amplitudes y perfiles exactos, y comparar directamente lineal versus no lineal. El escalado 2:1 solo no basta.
5. Para candidatos robustos/transitorios, verificar el escenario que realmente forma parte del problema, escalando amplitud de forma controlada y guardando fallos. TDS en amplitud finita es evidencia empírica, no una prueba universal de estabilidad no lineal.
6. Usar horizonte suficiente y un criterio de tail/settling. No reclamar máximo para t>=0 a partir de 12 o 60 s sin cota posterior. En ExpN el modo lento no asentó al 2% en 60 s.
7. No afirmar seguridad de límites de corriente no modelados. Exportar corrientes y headroom solo si rating/limite está definido, y señalar ausencia de limitador duro validado.

Si PD falla, candidato rechazado y diagnóstico conservado. Una revisión posterior requiere otro archivo/hash; no sobrescribir el candidato fallido. Un candidato analítico estable a guard numérico menor que el error de validación no obtiene garantía por redondeo.

## P6 — Resultados fuertes adicionales: certificados, no otra búsqueda ciega

Estos módulos pueden ejecutarse tras P3/P4 y no bloquean la validación de un candidato ya válido.

### P6A. GSP/self-energy con resto controlado

Usa ecuaciones (42)–(45). Mantén pivote declarado y repite al menos un pivote alternativo para no presentar la descomposición como un porcentaje físico invariante.

Si ||D^-1 W||<1, evalúa truncación y bound. Si no, retén inversa exacta o agranda bloques/ciclos: no pongas q artificialmente <1. Un bound sobre un contorno debe ser uniforme y controlar analiticidad/polos internos. Si solo se muestrea, status EMPIRICAL, no CERTIFIED.

Comparar coste, error en polos/gradientes, y MW obtenido con autoridad nodal, self-energy exacta y aproximación certificada. No confundir una correlación de graph strength con causalidad ni afirmar que el pivote es donde vive el modo.

### P6B. Cota global dentro del presupuesto del incumbent

Para nominal, cualquier mejora de U debe satisfacer sum w eps<U y eps_i<U/w_i. Usar eso para reducir dominio; no reutilizar U nominal para robustez si es infactible para beta_req.

Si intentas branch-and-bound:
- cajas por arquitectura con límites correctos;
- cotas inferiores del objetivo válidas;
- exclusión de estabilidad por inclusión/conteo uniforme de TODOS los polos;
- tratamiento explícito de singularidades, extremos eps=0 y ramas no suaves;
- outward rounding/aritmética validada para etiqueta FORMAL;
- lista persistente de cajas abiertas y L,U.

En la región nominal J<min w_i todos los GFL están presentes, de modo que hay 2^10 máscaras SG. Fuera de esa región, SG-only/mixed/GFL-only pueden requerir hasta 3^10 arquitecturas; no esconder esa diferencia.

Ni fracaso de Newton, ni ausencia de muestras factibles, ni una LMI suficiente incumplida prueban infactibilidad del problema original. Una condición de estabilidad conservadora interna no da por sí misma un lower bound de nuestra minimización.

Propuesta de tolerancia de globalidad práctica: 0.01 MW, claramente definida como tolerancia del estudio. Solo declarar GLOBAL_EPS_CERTIFIED si U−L<=0.01 con cotas válidas. Si L=0, etiquetar TRIVIAL_PHYSICAL_LOWER_BOUND. No llamar certificado nuevo al simple intervalo [0,U]. Presupuesto por defecto de este módulo: 30 min o 10000 cajas; si se agota, guardar y reportar INCOMPLETE, nunca “todo agotado”.

No es obligatorio demostrar monotonicidad de gains ni dominancia de bus 38. Si se intenta, debe valer en todo el dominio declarado y sobre todos los modos relevantes, respetando cambios de soporte. No promover signos locales a leyes globales.

## 4. Rendimiento, tests y entregables

Los tests deben recalcular propiedades, no solo leer PASS en un JSON. Incluir tests que fallen si:
- se cambia deliberadamente un P/Q compartido;
- se cambia signo del update de puerto;
- se ignora B(K) en el Woodbury;
- se invierte la frontera 2x2 real singular;
- se omite un polo interno;
- se evalúa gradiente congelando un trim que depende de parámetros;
- se declara certificado H-infinity a partir de muestras;
- se cambia un modelo tras congelación.

Por fase guardar:
- código fuente/config SHA y versión;
- residual e hipótesis de cada identidad;
- status EXACT_IDENTITY / CONDITIONAL_FORMULA / LOCAL_NUMERICAL_CERTIFICATE / GLOBAL_BOUND / EMPIRICAL_TDS;
- número de ecuaciones/evaluaciones/eigensolves/LU/llamadas PD;
- tiempo de compilación separado de ejecución y warm runs;
- memoria máxima si medible;
- casos rechazados y motivos.

Figuras científicas necesarias, solo si hay datos:
1. frontera de retención analítica versus polo completo;
2. regiones PI con componente estable etiquetada, incluyendo frontera real;
3. convergencia J, factibilidad y stationarity por iteración;
4. MW versus beta con bandas/cotas claras;
5. uniform/GSP/free: parámetros, MW y coste;
6. espectros analítico/PD y modos activos;
7. frecuencia y RoCoF lineal/no lineal;
8. error de expansión self-energy versus orden, y bound cuando existe.

No rellenar curvas o tablas faltantes con ExpE/G/H. Comparar solo métodos bajo el mismo trim, gain box, objetivo y requisitos.

Escribir:
- `REPORT_EXP_P.md`
- `EQUATIONS_SOURCE_MAP.md` con ecuación, archivo/función/hash y convención;
- `CLAIM_LEDGER_EXP_P.csv`
- `STAGE_STATUS.json`
- candidatos versionados + SHA;
- `README_REPRODUCE.md` con comandos instalados reales;
- `FINAL_SUMMARY_EXP_P.md`.

## 5. Resumen terminal obligatorio

EXP_P_STATUS:
LAST_COMPLETED_STAGE:
EXPN_MODEL_REUSED:
EXPN_REGRESSION:
PQ_CONTRACT_PASSED:
MODEL_SHA:
DESIGN_PD_CALLS:
PLL_SINGLE_DEVICE_RANK:
PLL_PAIR_AFFINE_IDENTITY:
REAL_BOUNDARY_HANDLED_SEPARATELY:
SINGLE_SG_QUADRATIC_IDENTITY:
EXPN_RETENTION_ROOT_REPRODUCED_MW:
CONDITIONAL_ROOT_ALL_POLES_PASS:
NOMINAL_OBJECTIVE_TYPE:
NOMINAL_RHO:
NOMINAL_KP:
NOMINAL_KI:
NOMINAL_RETAINED_SG_MW:
NOMINAL_ALPHA:
NOMINAL_KKT_SCOPE:
NOMINAL_KKT_RESIDUALS:
NOMINAL_NUMERICAL_GUARD:
ROBUST_UNCERTAINTY_DEFINITION:
ROBUST_BETA_REQUIRED:
HINF_GAMMA_LOWER:
HINF_GAMMA_UPPER:
ROBUST_BETA_CERTIFIED:
ROBUST_CERTIFICATE_TYPE:
ROBUST_RHO:
ROBUST_KP:
ROBUST_KI:
ROBUST_RETAINED_SG_MW:
TRANSIENT_PROFILE_AND_LOCATIONS:
TRANSIENT_CONSTRAINTS_ENFORCED:
MAX_ROCOF:
MAX_FREQUENCY_EXCURSION:
TRANSIENT_TAIL_STATUS:
GSP_PARAMETER_COUNT:
GSP_REPLACEMENT_LOSS_MW:
GSP_CERTIFIED_ERROR_OR_EMPIRICAL:
GLOBAL_LOWER_BOUND:
GLOBAL_UPPER_BOUND:
GLOBAL_GAP_MW:
GLOBAL_BOUND_TYPE:
OPEN_ARCHITECTURE_REGIONS:
PD_FULL_SPECTRUM_VALIDATION:
TDS_SMALL_SIGNAL_VALIDATION:
TDS_DESIGN_EVENT_VALIDATION:
NONLINEAR_SECURITY_THEOREM: NOT_CLAIMED salvo demostración independiente
PRIMARY_DESIGN_ITERATIONS:
PRIMARY_DESIGN_TIME:
PRIMARY_DESIGN_EIGENSOLVES:
MAIN_VERIFIED_RESULT:
MAIN_UNRESOLVED_RESULT:
PUSH: NO

## 6. Lo que debe quedar al final

No necesito otra promesa de un experimento siguiente. Necesito ejecución incremental reproducible hasta el último gate que permita la evidencia.

Los resultados más fuertes admisibles son, en orden:
1. Recuperar una frontera de retención conocida mediante el determinante de rango dos, sin rejilla de rho.
2. Derivar y verificar fronteras PI condicionadas usando la estructura de rango uno, sin confundirlas con un óptimo conjunto global.
3. Producir un diseño continuo local validado con menor coste computacional o mejores requisitos de robustez/transitorio.
4. Mostrar cuánta heterogeneidad GSP hace falta realmente, manteniendo PLL locales.
5. Añadir una cota global no trivial o un gap global verificable si el módulo logra demostrarlo.

Si alguna hipótesis falla, escribe exactamente cuál y qué resultados anteriores siguen siendo válidos. No modificar el modelo correcto de ExpN para fabricar un resultado nuevo.
