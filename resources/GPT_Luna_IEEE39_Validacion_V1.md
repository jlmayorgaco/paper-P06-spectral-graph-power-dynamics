# GPT Luna — Ejecutar una validación reproducible de PMUs virtuales en IEEE-39

## 0. Encargo y resultado que debes entregar

Implementa y EJECUTA un experimento completo. No respondas solamente con un plan, pseudocódigo o fórmulas. El objetivo es comprobar si nuestro estimador reconstruye las señales eléctricas de los buses no instrumentados a partir de ocho PMUs durante un transitorio de carga.

Usa Julia + PowerDynamics.jl para generar la planta y sus transitorios. Usa el núcleo Python/NumPy/SciPy adjunto para la inferencia, ampliándolo sólo en los adaptadores que faltan. No lo reescribas en Julia o C++ durante esta validación. La comunicación entre generador, estimador y evaluador debe ser por archivos, principalmente CSV y metadatos JSON. El hecho de usar Julia para la planta no obliga a que el estimador esté en Julia.

El paquete `inputs/pmu_nonlinear_em_reference.zip` es la referencia más reciente. Su README declara correctamente que NO incluye una validación IEEE-39, datos reales, un parser RAW ni un banco completo de eventos. No lo presentes como una solución end-to-end ya probada.

Esta instrucción sustituye los prompts antiguos del experimento, incluidos los que proponían reconstruir las señales simplemente estimando una amplitud escalar y volviendo a simular. Aquí el estimador debe reconstruir trayectorias de estados latentes desde las mediciones.

### Pregunta experimental principal

Con topología, parámetros eléctricos nominales y ubicación admisible de la perturbación conocidos, pero SIN conocer su amplitud, inicio, final, perfil temporal ni estados dinámicos verdaderos, ¿puede el estimador reconstruir V, I de rama, frecuencia y ROCOF en los 31 buses no-PMU?

Esta primera prueba es de **reconstrucción condicionada a un soporte de carga conocido (bus 7)**. NO es una prueba de localización ciega. No vamos a desarrollar todavía el clasificador ni un banco completo H1/H2.

### Experimento fijado

- Red: IEEE New England 39-bus del ejemplo oficial de PowerDynamics.jl.
- Planta: red AC con los generadores dinámicos completos del ejemplo, AVR y gobernadores. No sustituirla por swing de segundo orden.
- PMUs observadas: `[2, 5, 6, 10, 19, 22, 29, 39]`.
- Evento principal: incremento reversible de **10% de la escala de carga de bus 7**, conservando su relación P/Q.
- Intervalo de prueba: `0 <= t <= 30 s`.
- Perturbación: empieza en `5 s`, alcanza 10% en `5.5 s`, permanece hasta `9.5 s` y vuelve EXACTAMENTE a escala nominal en `10 s`.
- Rampas: polinomio quíntico C2, especificado abajo. No escalón discontinuo ni escalera de microcallbacks.
- Exportación PMU: `30 frames/s`, `901` frames, índices `k=0,...,900`.
- Inferencia: batch/full-window, iterativa sobre la trayectoria; puede usar el futuro. No hacer afirmaciones de causalidad ni tiempo real.
- Objetivo primario: reconstrucción eléctrica durante y después del transitorio, no sólo en estado estacionario.

## 1. Seguridad de datos y separación del proyecto

Trabaja en una carpeta de experimento aislada, por ejemplo:

`experiments/IEEE39_VIRTUAL_PMU_VALIDATION_V1/`

Si estás en el repositorio SGSMA, verifica su identidad antes de escribir. La ubicación histórica fue `C:\Users\walla\Documents\Github\SGSMA26-Physics-Informed-AD`; no asumas que ésa es la ruta de esta máquina.

No modificar la planta, artefactos, calibraciones ni estimadores congelados de producción. No hacer push. No inspeccionar, listar, abrir ni usar `RAW0002` ni `Competition_Testing Data Set 2.zip`. No hacer búsquedas recursivas que recorran esos conjuntos.

Los CSV nuevos son **SINTÉTICOS**, aunque tengan el esquema de las PMUs de la competición. No llamarlos RAW0001/RAW0002 ni mezclarlos con los datos reales. Usa identificadores `SIM_PMU39_V1_*`.

Crea `TEST_DATA_GUARD.txt` con las acciones realmente verificadas. No declarar comprobaciones que no ejecutaste.

## 2. Fuente única de verdad de la planta y entorno

Antes de programar, lee la documentación correspondiente a la versión instalada. Fuentes oficiales:

- https://juliaenergy.github.io/PowerDynamics.jl/stable/generated/ieee39_part1/
- https://juliaenergy.github.io/PowerDynamics.jl/stable/generated/ieee39_part2/
- https://juliaenergy.github.io/PowerDynamics.jl/stable/generated/ieee39_part3/
- https://docs.sciml.ai/SciMLBase/stable/interfaces/Callbacks/

Obtén de UNA versión identificada los archivos `bus.csv`, `branch.csv`, `load.csv`, `machine.csv`, `avr.csv` y `gov.csv`. Conserva copias, origen y SHA-256. Usa un entorno Julia local con `Project.toml` y `Manifest.toml`; registra versión de Julia, paquetes y commit/tag cuando aplique. No mezcles parámetros de MATPOWER, ANDES, papers y PowerDynamics para completar la planta.

Comprueba 39 buses, 46 ramas, 10 generadores y los registros de carga efectivos. No fijes por memoria el número de estados diferenciales/algebraicos; obtén ese número del modelo compilado.

### Modelo de cargas de la prueba principal

La validación principal será `MATCHED_Z`: cargas nominales de impedancia constante, para que el modelo algebraico nominal del estimador coincida con la clase de carga de la planta. Conserva los valores nominales P/Q y las bases de los datos públicos.

Inspecciona las fracciones ZIP. Si ya son Z puras, no cambies nada. Si no, configura explícitamente, únicamente para este experimento controlado:

`KpZ=KqZ=1; KpI=KqI=KpC=KqC=0`.

Registra los valores originales y los usados en `load_model_diff.csv`. No presentes esta elección como una propiedad universal del IEEE-39 ni como el modelo verdadero de SGSMA. No cambies generadores, controladores, líneas ni taps para conseguir buen ajuste.

Una prueba secundaria con las fracciones ZIP originales puede medir error de modelo después de cerrar MATCHED_Z; sus resultados deben quedar separados.

### Inicialización y solver

Utiliza la inicialización oficial desde flujo de potencia. En versiones que lo necesiten, aplica las fórmulas de inicialización de `ZIPLoad.Vset` en buses 31 y 39 tal como indica el tutorial. Verifica los nombres reales de símbolos; no pegues APIs inventadas.

Usa un integrador implícito apropiado, por ejemplo `Rodas5P` si es el utilizado por la versión elegida. Punto de partida: `reltol=1e-9`, `abstol=1e-10`, y un `dtmax` que resuelva las rampas. Confirma convergencia endureciendo tolerancias y reduciendo el paso máximo.

Registra residual de KCL, condiciones algebraicas, equilibrio inicial y código de terminación. Conserva una simulación nominal sin evento como control.

Cada simulación debe recibir COPIAS nuevas de estados y parámetros. Los callbacks y los objetos del problema pueden compartir parámetros mutables: compruébalo, no lo supongas.

## 3. Perfil exacto del evento principal

Define:

`S(u) = 6*u^5 - 15*u^4 + 10*u^3`, para `0 <= u <= 1`.

Define el perfil adimensional:

```
w(t) = 0                           t < 5
       S((t-5)/0.5)                5 <= t < 5.5
       1                           5.5 <= t < 9.5
       1 - S((t-9.5)/0.5)          9.5 <= t < 10
       0                           t >= 10
```

Entonces, sólo en bus 7:

`Pset7(t) = Pset7_base * (1 + 0.10*w(t))`

`Qset7(t) = Qset7_base * (1 + 0.10*w(t))`.

Si PowerDynamics usa inyección negativa para consumo, multiplicar ambos setpoints negativos por 1.10 aumenta el consumo. Verifica el signo a partir de las potencias efectivas.

A t>=10 usa los valores originales guardados; no multipliques por 0.90 para revertir.

El 10% es sobre la escala de una carga Z. No garantiza un aumento de 10% en MW instantáneos, porque éstos dependen de |V7|². Exporta tanto setpoints como P/Q realmente consumidos para evaluación privada.

Implementa el perfil como dependencia temporal continua en la carga o en un wrapper mínimo equivalente al componente ZIP nativo. No lo aproximes mediante saltos cada 1/30 s: eso contaminaría la derivada de fase. Si necesitas puntos de parada, utiliza 5, 5.5, 9.5 y 10 s sin reiniciar los estados físicos.

No reinicialices el sistema en un nuevo equilibrio durante el pulso. No fuerces que los ángulos, velocidades ni tensiones vuelvan exactamente al valor inicial en t=30. El parámetro vuelve a normalidad; la respuesta puede conservar un transitorio o desplazamiento angular.

## 4. Corrientes medidas: contrato sintético sin ambigüedades

Para este experimento elegimos explícitamente corrientes que salen del bus PMU hacia la rama indicada:

| Bus PMU | Rama | Terminal de medida | Sentido positivo |
|---|---|---|---|
| 2 | 2–3 | bus 2 | 2 hacia 3 |
| 5 | 5–6 | bus 5 | 5 hacia 6 |
| 6 | 6–7 | bus 6 | 6 hacia 7 |
| 10 | 10–11 | bus 10 | 10 hacia 11 |
| 19 | 19–16 | bus 19 | 19 hacia 16 |
| 22 | 22–21 | bus 22 | 22 hacia 21 |
| 29 | 29–28 | bus 29 | 29 hacia 28 |
| 39 | 39–1 | bus 39 | 39 hacia 1 |

Los endpoints deben existir en la red seleccionada. Determina FROM/TO a partir del catálogo real, no del orden de esta tabla.

Estos son **canales sintéticos elegidos**, no una certificación de los extremos de corriente de los archivos SGSMA reales, especialmente PMUs 2,29,39.

Guarda `pmu_current_map.csv`: pmu_bus, branch_id, from_bus, to_bus, measured_terminal_bus, terminal, sign, voltage_base_kV, current_base_A.

Extrae corrientes nativas de la planta y verifica independientemente que coinciden con los sellos de rama y KCL. El componente nativo podría reportar corriente entrando al bus, mientras aquí usamos corriente saliendo del bus. Corrige ese signo una sola vez y documéntalo.

NO interpretar directamente `r_src` de PowerDynamics como `tap_complex` de MATPOWER sin conversión demostrada. La función `network_matrices` del paquete de referencia espera sellos canónicos ya normalizados, no los CSV nativos de PowerDynamics.

Además exporta todos los terminales de las 46 ramas. Un bus puede tener varias corrientes; no llamar `Ybus*V` a la corriente de una línea.

## 5. Verdad física de V, I, frecuencia y ROCOF

Obtén V positivo-secuencia e I_from/I_to desde la planta dinámica completa. No construyas la verdad con el mismo `Readout` del estimador: necesitamos comparar contra el simulador, no contra una copia de la ecuación inversa.

Genera una trayectoria de evaluación de alta resolución, inicialmente 600 muestras/s. De ella calcula, por bus, la fase temporal continua y:

`f_i(t) = f_frame + Im(Vdot_i/V_i)/(2*pi)`

`rho_i(t) = Im(Vddot_i/V_i - (Vdot_i/V_i)^2)/(2*pi)`.

En la simulación de referencia a 60 Hz, `f_frame=60 Hz`. Si eliges otro marco, registra la transformación y comprueba la identidad.

No uses la velocidad de rotor como frecuencia de bus. No asignes la misma frecuencia a todos los buses. No pongas ROCOF39=0 para imitar un artefacto del RAW real.

Preferencia: derivadas consistentes de la solución/DAE. Si el solver no proporciona derivadas algebraicas fiables, usa diferenciación centrada de orden alto o polinomio local sobre la trayectoria determinista de alta resolución. Especifica coeficientes, ancho de ventana y orden. No diferencies señales con ruido para construir la verdad.

Repite la extracción a 1200 muestras/s o con un paso de diferenciación reducido. Reporta la convergencia de f/rho. Declara una máscara común en los bordes donde la derivada no sea fiable; no ocultes segmentos malos según el error del estimador.

La elección C2 del evento evita saltos ideales de fase. Aun así, verifica regularidad y convergencia en las rampas.

Finalmente muestrea/exporta a `t_k=k/30`, k=0,...,900. Los timestamps decimales del CSV pueden redondearse a 3 decimales para compatibilidad, pero conserva `sample_index` y el tiempo racional exacto en el manifiesto. No infieras tiempos 0.033/0.034 como si el muestreo físico fuera irregular.

## 6. CSV tipo RAW, pero inequívocamente sintéticos

Por bus produce las 17 columnas:

`TIMESTAMP, VA_mag, VA_ang, VB_mag, VB_ang, VC_mag, VC_ang, IA_mag, IA_ang, IB_mag, IB_ang, IC_mag, IC_ang, Frequency, ROCOF, DATA_PRESENT, Event`.

Convenciones:

- Magnitudes V: voltios RMS fase-neutro.
- Magnitudes I: amperios RMS del terminal elegido.
- Ángulos: grados, convención documentada.
- Frecuencia: Hz; ROCOF: Hz/s.
- `DATA_PRESENT=1` en todo este experimento. No estudiar ciber/faltantes.
- `Event` en los archivos que puede leer el estimador: VACÍO, no contiene etiquetas verdaderas.

Para expandir secuencia positiva, con `a=exp(j*2*pi/3)`:

`[VA,VB,VC] = [V1, a^2*V1, a*V1]`.

Lo mismo para cada corriente. Esto representa fasores balanceados, no tres fases independientes ni formas de onda EMT.

Bases:

`Vbase_phase = Vbase_LL/sqrt(3)`

`Ibase = Sbase/(sqrt(3)*Vbase_LL)`.

Utiliza las bases del bus donde se mide cada corriente, particularmente en transformadores.

### Archivos y separación

```
public_model/             # parámetros estáticos permitidos, bases, mapping
calibration_observed/     # sólo ocho PMUs de calibración
observed/                 # sólo ocho CSV de prueba
truth_private/            # 39 buses, todas las ramas, perfil real y etiquetas
reconstruction/           # estimaciones y covarianzas, sin leer truth_private
metrics/                  # evaluador, se ejecuta después
plots/
logs/
```

Para los buses ocultos, un CSV tipo PMU puede elegir una corriente terminal representativa mediante una regla determinista: vecino de menor ID (con desempate por branch_id). Publica `virtual_pmu_current_map.csv`. Para los buses observados, esa regla se sustituye por la tabla de §4.

No basta con ese canal representativo: exporta también todas las corrientes terminales en formato largo y evalúalas por rama/extremo.

Crea nombres como `SIM_PMU39_V1_Bus7_Truth.csv` y `SIM_PMU39_V1_Bus2_Observed.csv`. No reutilices nombres de los RAW reales.

## 7. Cortafuegos entre generador, estimador y evaluador

El generador puede conocer toda la simulación. El estimador no.

### El estimador puede leer

- Y_N / sellos, bases, topología y cargas nominales del experimento.
- Flujo de potencia nominal completo V0. Es metadata de la planta, no verdad transitoria.
- Las ocho PMUs observadas y sus timestamps.
- Mapping de corrientes sintéticas.
- Configuración de calibración que quedó congelada antes del test.
- Que el modelo de reconstrucción admite una perturbación de carga en bus 7.

### El estimador NO puede leer

- Trayectorias verdaderas de los 31 buses ocultos o generadores.
- `alpha_true`, los tiempos 5/10, duración o función w(t) del evento de prueba.
- Potencias/corrientes verdaderas de la carga durante el transitorio.
- Etiquetas Event o el manifiesto privado del evento.
- Ruido de medición realizado, semillas de prueba o valores sin ruido de las PMUs de test.
- Métricas del evaluador para volver a ajustar el modelo.

Si todos estos procesos corren en el mismo equipo, ejecútalos separados, con CLI y lista de rutas de entrada explícita. Registra los archivos realmente abiertos por el estimador. El evaluador se ejecuta sólo después de guardar las predicciones y sus hashes.

No uses `alpha_hat -> nueva simulación completa` como reconstrucción principal. Eso sería otro experimento más débil y no el estimador del paquete.

## 8. La ampliación mínima del estimador que NECESITA este transitorio

El modelo H0 con sólo generadores no contiene un cambio de carga en bus 7. No lo fuerces para que produzca señales ocultas correctas durante la perturbación.

Construye DOS modelos de observación sobre la misma planta y datos:

### M0: control de error de modelo H0

`A0 = Y_N + C_L*diag(yL0)*C_L'`

`A0*M0 = C_G`.

Aquí las cargas nominales están fijas. Es un control; puede fallar durante el evento y eso debe mostrarse.

### M7: modelo principal, soporte de carga conocido

Con corriente de carga residual `d7` definida positiva como INYECCIÓN a la red:

`A0*V = C_G*I_G + e7*d7`.

En desviaciones:

`A0*(V-V0) = C_G*delta_I_G + e7*d7`.

Construye:

`C7 = [C_G, e7]`  (39 x 11 compleja)

`A0*M7 = C7`.

Entonces:

`V = V0 + M7*[delta_I_G; d7]`.

No suministres `d7(t)` al estimador. Es parte de la trayectoria desconocida que debe inferir.

Para MATCHED_Z, la identidad que el EVALUADOR puede comprobar es:

`d7_true(t) = -0.10*w(t)*yL7_base*V7_true(t)`.

Eso justifica que el modelo M7 pueda representar exactamente la intervención elegida en su parte algebraica sin conocer su perfil. Las dinámicas reales de generadores siguen siendo las del simulador, no un random-walk usado para fabricar la verdad.

Esta modificación es sólo añadir una columna de fuente; NO reemplaza el núcleo nonlinear_smoother ni introduce una nueva familia de estimadores.

## 9. Bases observables y nulas

Para cada modelo M0/M7:

`G_complex = [S_P; L_I] * M`

`G_real = realify(G_complex)`.

Ejecuta SVD con varias tolerancias y registra todos los valores singulares, dimensiones y escala utilizada. No copies rango18 o 36/39 de SGSMA a esta red sin medirlos.

Con `W_o` y `N` de la SVD real:

`B = M * (W_o_realpart + j*W_o_imagpart)`

`B_null = M * (N_realpart + j*N_imagpart)`.

Por cada bus/terminal calcula sensibilidad al espacio nulo. Publica los mapas:

`bus_identifiability.csv`

`branch_terminal_identifiability.csv`.

No abortar por rango deficiente. Continuar con coordenadas observables y un prior explícito para lo no observable.

Conserva `B`, su orden y sus signos fijados para todas las ablaciones del mismo modelo. Al retirar mediciones, NO reduzcas silenciosamente el estado ni cambies su prior para favorecer una variante. Registra la información adicional que se pierde.

El prior nulo no se aprende de la misma verosimilitud PMU. Usa media nominal y una covarianza documentada, con sensibilidad de escala `0.1, 1, 10, 100`. Separa métricas de salidas identificables de las asistidas por prior. No fuerces el grupo {20,33,34}: calcula qué salidas están afectadas en esta red y este M7.

## 10. Núcleo matemático a utilizar, sin redefinirlo

Lee `README.md`, `pmu_core.py`, `test_core.py`. Ejecuta las pruebas ANTES de cualquier integración.

Estado:

`s = [u; udot; uddot]`, con cada bloque real de longitud r, y r obtenido de la SVD.

Voltajes y derivadas:

`V = V0 + B*u + B_null*mean_null`

`Vdot = B*udot`

`Vddot = B*uddot`

para un componente nulo constante en la ventana. En las PMUs ese componente debe ser invisible hasta tolerancia.

Corrientes medidas:

`I_P = L_I*V`.

Frecuencia y ROCOF:

`f_i = f_ref + imag(Vdot_i/V_i)/(2*pi)`

`rho_i = imag(Vddot_i/V_i - (Vdot_i/V_i)^2)/(2*pi)`.

Temporal prior:

```
F(h) = [[I, hI, h²I/2],
        [0, I,  hI   ],
        [0, 0,  I    ]]

Q(h) = [[h⁵/20, h⁴/8, h³/6],
        [h⁴/8,  h³/3, h²/2],
        [h³/6,  h²/2, h    ]] kron diag(D).
```

Este prior es jerk blanco integrado, NO swing. La planta generadora sí tiene dinámica física; la inferencia utiliza restricciones algebraicas exactas y un prior temporal estadístico.

Funciones reales disponibles:

- `Readout.evaluate`
- `observable_coordinates`
- `discretize`
- `linear_smoother`
- `nonlinear_smoother`
- `process_moments`
- `update_jerk_density`
- `update_R_local`

Usa la inicialización desde V/I observados y sus covarianzas o desde el nominal con prior propio; nunca desde corrientes/derivadas verdaderas ocultas. La trayectoria inicial puede ser plana, y una inicialización mejor sólo puede usar las mismas ocho PMUs.

El readout está en orden por PMU: `[Vr,Vi,Ir,Ii,f,rho]`. `keep` contiene filas de la medición original. Al aplicar una ablación, `ys`, R, grupos de covarianza y `keep` deben estar en la misma convención; los grupos del M-step se reindexan al vector reducido.

El smoother itera sobre la ventana completa. Cada iteración reinicia el subproblema con el MISMO m0/P0. No reutiliza el posterior como otro prior independiente para los mismos datos.

No etiquetes el resultado como MMSE no lineal exacto. Es suavizado no lineal iterado y EM aproximado con momentos gaussianos locales.

## 11. Ruido y calibración: evitar que EM suprima el transitorio

### N0: prueba numérica sin ruido

La verdad y las ocho PMUs coinciden. El solver del estimador sigue necesitando covarianzas positivas; usa pisos declarados como ponderaciones/tolerancias numéricas, no una R=0 invertida. No interpretes esas covarianzas como incertidumbre de sensor real.

### N1: prueba sintética con ruido conocido

En positivo-secuencia, por PMU y frame, añade ruido gaussiano rectangular independiente con estas desviaciones iniciales de diseño:

- Re(V), Im(V): `0.001 pu` cada una.
- Re(I), Im(I): `0.002 pu` cada una.
- Frecuencia: `0.003 Hz`.
- ROCOF: `0.02 Hz/s`.

Después convierte V/I ruidosos a magnitud/ángulo y ABC para los CSV. Al expandir ABC, las tres fases comparten la perturbación positiva-secuencia; no deben tratarse como tres sensores independientes.

Estas cifras son supuestos sintéticos de benchmark, NO especificaciones de clase M ni ruido verificado de SGSMA. Este es un modelo de PMU fasorial ideal con ruido, no una emulación IEEE C37.118 completa.

Genera todos los canales una vez por semilla. Las ablaciones usan las mismas realizaciones de ruido en los canales comunes.

### Calibración separada

Genera tres conjuntos independientes, que sólo entregan ocho PMUs al entrenamiento:

1. `CAL_NORMAL_TRAIN`: nominal, 60 s, semillas separadas.
2. `CAL_NORMAL_VALID`: nominal, 30 s, sin solapamiento de ruido.
3. `CAL_EXCITED`: 30 s, una perturbación suave de +2% en bus 7 entre 5–10 s y una de -2% entre 18–23 s, con rampas de 0.5 s; no el pulso +10% de test.

Aprende R en el conjunto normal, incluyendo la covarianza posterior del estado en el M-step. Contrasta R aprendido con R inyectado sólo como verificación de calibración.

Para la densidad temporal D, NO pretendas que un equilibrio perfectamente quieto identifica la variabilidad de un transitorio. Conserva un modelo `D_normal` de H0 como control, pero para el reconstructor M7 calibra `D_reconstruction` en CAL_EXCITED usando exclusivamente las ocho PMUs, sin amplitudes, tiempos ni estados verdaderos del generador.

Congela R al estimar D en CAL_EXCITED para que EM no convierta la señal física en ruido. Documenta esta decisión: es calibración de un reconstructor de transitorios, no aprendizaje de un detector H0 puro.

El prior de derivadas m0/P0, escalas y límites de D deben declararse antes del test. No recalibrar D, R, P0, ventanas o suavizado usando errores de los buses ocultos de prueba. Si exploras hiperparámetros, utiliza sólo datos de calibración y criterios predictivos de las PMUs, guarda el espacio explorado y la regla de selección.

Ejecuta primero el pipeline con R inyectado conocido como verificación técnica. Después la variante con R aprendido. No confundas esos resultados.

En el M-step de R, comprueba que los grupos cubren exactamente las filas retenidas y que R sigue definida positiva. `update_R_local` NO implementa ruido temporal AR. Si la generación añade correlación temporal en una prueba secundaria, requiere un modelo correspondiente o debe reportarse como mismatch, no como likelihood exacto.

## 12. Ablaciones que responden la pregunta sobre frecuencia/ROCOF

Sobre el MISMO evento, topología, matriz B, prior nulo y semillas, ejecutar M7 con:

A: sólo V (16 filas reales).

B: V + I (32 filas reales).

C: V + I + f (40 filas reales).

D: V + I + f + ROCOF (48 filas reales).

En la comparación principal conserva la misma calibración/prior temporal predefinido para aislar el efecto de retirar mediciones. Etiqueta esa política claramente. Una variante con recalibración específica por ablación es secundaria y no se mezcla con la anterior.

Ejecuta además:

- `H0_fixed_loads`: mismo núcleo con M0, como control de modelo inadecuado durante la carga extra.
- `Nominal_constant`: voltajes/corrientes del flujo de potencia fijo.
- `Estimator_oracle_sensor_noise`: M7 con R verdadero, sólo como comprobación del ruido.
- `Forward_truth`: la salida del simulador. No es un estimador ni una comparación competitiva.

No presupongas que f/ROCOF van a mejorar todo. Reporta mejoras, empates o empeoramientos.

## 13. Repeticiones y estadística

Primero: todos los controles de modelo e integración en una realización sin ruido.

Después: al menos 20 semillas independientes de ruido de medición para N1 y las cuatro ablaciones. Reutiliza la misma trayectoria física determinista y los mismos parámetros congelados; esto evalúa variabilidad por ruido, NO generalización entre puntos operativos.

Semilla maestra sugerida: `2026092501`. Usa subflujos separados para calibración, test, intervalos de incertidumbre y cualquier aleatoriedad numérica. El estimador no recibe la semilla con la que se generó el ruido del test.

Registra cada semilla, estado del solver, número de iteraciones, costo inicial/final y tiempo. Si alguna realización falla, inclúyela en la tabla y reporta la tasa de fallo; no promediar sólo las exitosas sin decirlo.

Calcula intervalos por semilla/realización, no tratando 901 frames autocorrelacionados como 901 ensayos independientes.

## 14. Métricas obligatorias

Comparar las reconstrucciones contra verdad eléctrica sin ruido, sólo DESPUÉS de finalizar las predicciones.

### Voltajes

Por bus:

`RMSE_complex = sqrt(mean(abs(Vhat - Vtrue)^2))`

`RMSE_mag = sqrt(mean((abs(Vhat)-abs(Vtrue))^2))`

Error angular:

`angle_error = angle(Vhat*conj(Vtrue))`.

Reportar su RMSE en grados y máximo absoluto. Alinear únicamente el marco definido antes de estimar; NO ajustar un offset por bus usando la verdad oculta para reducir el error.

### Corrientes

Las mismas métricas para cada terminal de cada rama. Para ángulos de corrientes casi nulas usa una máscara de magnitud declarada antes de evaluar, y reporta cuántos puntos excluye.

### Frecuencia y ROCOF

RMSE, MAE, p95 de error absoluto y máximo. Utiliza las ecuaciones del readout/estado suavizado para reconstruir, no un postfiltro distinto elegido para favorecer cada curva.

### Error del transitorio, no sólo del nivel de 1 pu

Con una simulación nominal independiente sin evento:

`eta_V = norm(Vhat_hidden - Vtrue_hidden) / norm(Vtrue_hidden - Vnominal_hidden)`

sobre el intervalo 5–15 s, con manejo explícito del denominador casi cero. Calcula también por bus cuando tenga amplitud suficiente.

`improvement = 1 - RMSE_estimator / RMSE_nominal`.

No imponer un porcentaje universal de éxito después de ver resultados. Un valor eta<1 indica mejora frente a ignorar el evento, no precisión industrial certificada.

### Ventanas de evaluación

- preevento: [0,5)
- perturbación: [5,10]
- recuperación temprana: (10,15]
- recuperación tardía: (15,30]
- total: [0,30]

La información sobre estas ventanas es sólo del EVALUADOR.

### Subconjuntos

- 8 buses PMU.
- 31 buses sin PMU.
- Sólo salidas identificables según B_null.
- Salidas asistidas por prior.
- Grupo 20/33/34 por separado, sin asumir a priori su resultado.
- Bus 7 por separado.

No presentar buen error de bus 7 como prueba de reconstrucción global: la PMU 6 y su corriente 6→7 lo hacen especialmente favorable algebraicamente.

## 15. Incertidumbre y espacio nulo

Para V/I, propaga la covarianza observable y suma la contribución del prior nulo. No interpretes la pseudoinversa de una información singular como varianza cero en su núcleo.

Para f/rho, usa el Jacobiano de la salida respecto de u,udot,uddot y también respecto del componente nulo que aparece en el denominador. Alternativamente usa muestras del posterior local y del prior nulo, sin seleccionar muestras con la verdad.

Reporta cobertura puntual de intervalos 50%,90%,95% y anchura media, separadas por variable, bus y fase temporal. Aclara que son intervalos puntuales, no bandas simultáneas con esa misma cobertura.

En comparación contra verdad sin ruido usa covarianza de la variable latente. En leave-one-PMU-out contra medición ruidosa hay que añadir ruido del canal retenido fuera; no mezclar ambas evaluaciones.

El prior nulo puede dar intervalos anchos o falta de cobertura. Ese resultado debe conservarse; no ajustar su media con las trayectorias internas verdaderas.

## 16. Figuras que deben salir

Exporta PNG y SVG, con unidades, leyendas legibles y texto que identifique el escenario:

1. Perfil real de escala de carga y P/Q realmente consumidos en bus 7 (figura del evaluador).
2. Tensiones, frecuencia y ROCOF medidos en las ocho PMUs durante el pulso.
3. Verdad vs reconstrucción vs nominal en buses ocultos 7,18,27,20,33,34; muestra tensión compleja mediante magnitud/ángulo y al menos frecuencia/ROCOF en figuras separadas.
4. Corrientes de ramas 6–7,16–19,19–20,19–33,20–34 cuando existan, para mostrar tanto la región del evento como el pocket.
5. RMSE por cada uno de los 31 buses no-PMU, con identificación de salidas asistidas por prior.
6. Mapa bus-tiempo del error de magnitud y otro de frecuencia/ROCOF.
7. Comparación pareada de las cuatro ablaciones por semilla, no sólo una curva elegida.
8. Costo del smoother por iteración, historia de EM y tasas de fallos.
9. Valores singulares y sensibilidades al espacio nulo.
10. Error normalizado vs incertidumbre/cobertura, con el pocket separado.

No uses una única escala que oculte el transitorio sobre el offset de 1 pu. No recortes curvas problemáticas para que parezcan coincidir.

## 17. Pruebas de integridad numérica

Obligatorias antes de aceptar resultados:

- `test_core.py` original pasa; guardar salida y hash del código base.
- KCL y potencias nominales consistentes.
- Modelo algebraico reproduce la condición nominal de la planta.
- Corrientes de las 8 PMUs coinciden con terminales nativos y sus signos/bases.
- No confundir corriente de línea con corriente neta de bus.
- Ramas y mapping por bus ID, no confundir índices base 0/1 entre Python y Julia.
- Reconstrucción ABC→012 recupera positivo-secuencia y no introduce 3× información.
- La perturbación es realmente +10% de escala y retorna a EXACTAMENTE 1.0.
- Simulación nominal permanece estable.
- Nueva ejecución no hereda parámetros de la anterior.
- Tolerancias y derivadas de alta resolución convergen.
- Jacobianos analíticos vs diferencias finitas también sobre al menos tres estados estimados de esta red, no sólo matrices aleatorias.
- Covarianzas SPD; fallos numéricos se reportan, no se corrigen con jitter ilimitado silencioso.
- Ningún incremento pequeño causado por una búsqueda de paso bloqueada se etiqueta como convergencia.
- Cuando el solver diga LOCAL_STEP_CONVERGED, comprobar además estabilidad de costo, gradiente o residual estacionario en variables escaladas y ausencia de paso truncado patológico.
- Los datos originales se usan una vez estadísticamente en cada solve, aunque haya múltiples iteraciones.
- No usar el futuro oculto: batch permite futuro de las OCHO PMUs, no de los 31 buses ocultos.

## 18. Controles negativos

Ejecuta al menos:

- M0 con cargas fijas sobre el evento: cuantifica el error de no modelar la intervención.
- Estimador sin evento sobre una simulación nominal: no debe crear un pulso espurio.
- Una corrupción diagnóstica deliberada de signo de UNA corriente, fuera del experimento principal: verificar que los controles de física/residual la detectan. No usarla para reescoger el mapping ya fijado.

Estos controles son diagnósticos. No deben alterar las semillas ni parámetros del resultado principal.

## 19. Entregables

Scripts sugeridos:

```
00_preflight.py
01_build_and_simulate.jl
02_export_synthetic_raw.jl
03_build_estimator_maps.py
04_calibrate.py
05_reconstruct.py
06_evaluate.py
07_plot.py
run_all.py
```

Puedes organizar distinto si mejora el código, pero debe existir un comando documentado que ejecute la cadena completa y reanude sin repetir simulaciones ya verificadas.

Entrega:

- código Julia/Python completo;
- Project.toml, Manifest.toml, requirements/lock;
- datos estáticos originales y utilizados con licencia/procedencia;
- config pública del estimador y config privada del evento separadas;
- mapping de canales y bases;
- CSV sintéticos de ocho PMUs, verdad privada y reconstrucciones;
- tablas por bus, rama, ventana, variante y semilla;
- plots;
- logs íntegros de solvers/calibración, tiempos/hardware;
- `REPORT.md`, `LIMITATIONS.md`, `RESULTS_SUMMARY.json`;
- `TEST_DATA_GUARD.txt`, `ESTIMATOR_READS.json`, manifiesto SHA-256;
- ZIP autocontenido de la carpeta del experimento para revisión externa.

Si el volumen es grande, incluye todos los resultados de la realización representativa, todas las métricas de semillas y un manifiesto exacto para regenerar los CSV adicionales. No omitas los casos fallidos.

## 20. Informe y veredicto final

No afirmar de antemano RMSE=0.004 ni rango=18 ni que el método es publicable. Calcula los resultados.

Responde con algo como:

```
PLANT_SOURCE_AND_VERSION:
PLANT_MODEL_AND_LOAD_ASSUMPTION:
SIMULATION_VALID:
NOMINAL_INITIALIZATION_VALID:
EVENT_BUS: 7
EVENT_SCALE_TRUE: 0.10
EVENT_PROFILE_VALID:
PMU_BUSES: [2,5,6,10,19,22,29,39]
CURRENT_MAPPING_VERIFIED:
ESTIMATOR_KNOWS_SUPPORT: YES — LOAD@7 (conditional reconstruction)
ESTIMATOR_KNOWS_AMPLITUDE: NO
ESTIMATOR_KNOWS_ONSET_DURATION_PROFILE: NO
HIDDEN_TRUTH_ACCESSED_BY_ESTIMATOR: NO/FAILED_AUDIT
H0_REAL_RANK:
M7_REAL_RANK:
NON_PMU_OUTPUTS_IDENTIFIABLE_FROM_DATA:
NON_PMU_OUTPUTS_PRIOR_ASSISTED:
VALID_NOISE_REPETITIONS:
FAILED_REPETITIONS:
NON_PMU_VMAG_RMSE:
NON_PMU_COMPLEX_V_RMSE:
NON_PMU_ANGLE_RMSE_DEG:
BRANCH_TERMINAL_I_RMSE:
NON_PMU_FREQUENCY_RMSE_HZ:
NON_PMU_ROCOF_RMSE_HZ_S:
TRANSIENT_ETA_V:
IMPROVEMENT_VI_OVER_V:
IMPROVEMENT_VIF_OVER_VI:
IMPROVEMENT_VIFR_OVER_VIF:
95_PERCENT_COVERAGE_AND_WIDTH:
POCKET_20_33_34_RESULT:
EM_RESULT_AND_CALIBRATION_LIMITS:
NUMERICAL_FAILURES_OR_MODEL_MISMATCH:
WHAT_IS_DEMONSTRATED:
WHAT_IS_NOT_DEMONSTRATED:
FINAL_STATUS: PASS / PARTIAL / FAIL
ZIP_PATH:
```

PASS significa que la cadena está correctamente ejecutada y satisface criterios predeclarados, no observabilidad universal. Puedes tener un experimento ejecutado correctamente y un resultado científico negativo. No conviertas FAIL en PASS retocando la planta o leyendo la verdad.

La conclusión buscada es precisa:

«Bajo este modelo de red/cargas, estas ocho PMUs y una perturbación admisible en bus 7, la reconstrucción de señales ocultas consigue estos errores e incertidumbres; frecuencia y ROCOF aportan —o no— esta mejora medida».

No demuestra todavía detección/localización ciegas, validez de todos los estados de máquina, ni generalización a RAW0002. No añadas esas tareas para retrasar esta validación.

## 21. Instrucción de ejecución

Empieza por inspeccionar el paquete adjunto y el entorno. Ejecuta los tests del núcleo. Implementa la integración completa, corre primero nominal + evento sin ruido, resuelve los fallos verificables y después ejecuta la matriz de ruido/ablaciones predefinida.

No vuelvas a pedir que se elija bus, evento, duración, PMUs, lenguaje o primera pieza: están definidos aquí. No sustituyas el trabajo por una nueva propuesta de arquitectura. Entrega resultados ejecutados o, si existe un bloqueo real de entorno/datos, entrega el código y la evidencia exacta del bloqueo sin inventar resultados.
