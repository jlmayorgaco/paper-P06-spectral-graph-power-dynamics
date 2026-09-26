# Claude Code — cerrar la teoría mediante certificados, búsqueda y diseño

Continúa el repositorio actual, sin reiniciar el programa ni cambiar los resultados congelados. Lee primero `TEORIA_Y_DEMOSTRACIONES.md` y ejecuta `theory_checks.py` en un entorno compatible. Ese script verifica ejemplos matemáticos sintéticos; NO proporciona resultados IEEE-39, Kundur ni IEEE-68 y no debe presentarse como tal.

## 0. Objetivo y límites

El objetivo es cerrar tres contribuciones comprobables:

1. Conteo físico de estabilidad sin esconder modos aperiódicos lentos bajo un cutoff alrededor de cero.
2. Certificados para familias de decisiones binarias y búsqueda exacta de testigos mínimos sin asumir monotonía.
3. Diseño de portfolios seguros ante cualquier implementación parcial, y separación de ese requisito respecto a endpoint estable, secuencia estacionaria segura y transición dinámica segura.

No reclamar novedad de Schur, generalized Nyquist, el punto −1, sum-of-squares, Lyapunov o hypergraphs. No reclamar complejidad polinómica de un problema cuya salida puede ser exponencial. No vender un certificado muestreal como garantía continua.

No reabrir Track B como programa aparte. Kundur sirve aquí para el problema de los cruces aperiódicos del programa ya existente.

## 1. Fuentes y control de cambios

La evidencia existente autorizada es la tabla final TPWRS, G5, G1/G2/G3 y los ledgers actuales. No importar afirmaciones superadas de `METHODS.md` o de los resúmenes antiguos. En particular:

- G1 todavía excluye `|s| <= 1e-3`: debe auditarse, no tratarse como definición física universal.
- ANDES sólo armonizó exactamente la banda electromecánica en el caso de inyecciones estáticas y PSS desactivado. No convertir eso en validación independiente del GFL/PSS completo.
- F11 es afín en parámetros continuos para cada portfolio; no demuestra afinidad en los indicadores binarios de reemplazo.
- La prueba principal de cuatro candidatos IEEE-68 fue negativa. El resultado sobre doce plantas fue una verificación secundaria en tres puntos, no un mapa completo de puertos.
- Un condensador sin amortiguamiento puede introducir su propia inestabilidad. No reutilizar claims anteriores a G1.

Crea una rama nueva desde los últimos commits pertinentes, sin hacer checkout destructivo de cambios ajenos. Snapshot y hash de la evidencia. No sobrescribir F/G ni tocar el póster. Crea una subcarpeta de experimentos `binary_certification_v1` y archivos identificados `BC00`...`BC12`.

Los artefactos grandes deben tener checksum, manifiesto y ruta versionada. Commit de código/configs/docs/tests; datos comprimidos sin pérdida según las reglas del repo. No hacer push ni alterar repositorios remotos.

## 2. Primer bloque obligatorio: BC00 — integridad física y numérica

### BC00-A: MW, MVA y bases

Reconcilia esta alerta antes de producir un nuevo gráfico de beneficio industrial:

`PHASE_E_REPORT.md` llama 4270.7 “MW de PV”, pero la tabla de máquinas original tiene ratings Sn de 1040, 1174.8, 1085.7 y 970.2 MVA, cuya suma es 4270.7.

No supongas que es un bug confirmado: rastrea cada columna desde el dato fuente.

Por máquina, portfolio y reparación guarda:

- S_base_MVA del sistema;
- Sn_machine_MVA;
- inverter_rating_MVA;
- active_nameplate_MW, sólo si existe una fuente válida;
- P_dispatch_pu y P_dispatch_MW;
- Q_dispatch_Mvar;
- available_PV_MW;
- retired_active_MW;
- condenser_MVA;
- P/Q/current capability utilization.

No convertir MVA a MW implícitamente ni suponer factor de potencia uno cuando se usa Q-matched. Si hay error de unidades, genera una tabla corregida nueva y downgrade de los claims afectados; no alteres los espectros anteriores sin causa.

### BC00-B: origen y referencias

Enumera todos los modos próximos a cero, sin borrarlos:

`1e-2, 1e-3, 1e-4, 1e-6, 1e-8` son ventanas de diagnóstico, no reglas de eliminación.

Deriva el generador de la simetría angular a partir de las ecuaciones y las coordenadas del modelo. Identifica cuáles estados/voltajes rotan bajo la misma transformación. Verifica `A R = 0` y, en la DAE, las identidades diferenciales y algebraicas correspondientes.

Si hay dos ceros y sólo una simetría:

- identifica si el segundo es un integrador físico de frecuencia, un estado muerto o una cadena de Jordan;
- NO retires ambos por tamaño;
- el proyecto debe reportar un modo físico marginal aunque nunca sea positivo.

Mantén cuatro estados de clasificación:

`STABLE`, `UNSTABLE`, `BOUNDARY_OR_UNRESOLVED`, `INFEASIBLE`.

Congela criterios usando residuos de autovectores y errores hacia atrás. No uses sólo un umbral de Re(lambda) elegido para no alterar H.

### BC00-C: normalización

En las comparaciones de normas, fija unidades y métricas físicas. Bajo cambios de base transforma también la métrica. No llames al espectro-distancia a −1 un radio de robustez no normal.

**Salida BC00:** tabla de unidades, clasificación de ceros y lista explícita de resultados anteriores potencialmente afectados. Detener los claims afectados, pero conservar los experimentos independientes que sigan siendo válidos.

## 3. BC01 — simetría y los 194 cruces aperiódicos de Kundur

Implementa T1 usando generadores de simetría conocidos, no un eigenvalue cutoff. Empezar en el A_red ya validado; luego generalizar al descriptor si es necesario.

Para cada caso:

1. Construir base ortonormal U de la simetría y Z de su complemento.
2. Construir A_q = Z.T A Z.
3. Verificar la estructura triangular por bloques y la identidad característica.
4. Retener cualquier cero físico adicional en A_q.
5. Implementar opcionalmente la reubicación A_sharp = A - beta U U.T, que desplaza sólo la simetría; comprobar que el espectro físico no depende de beta.
6. Derivar el nuevo puerto DESPUÉS del cociente, no eliminando una fila y columna arbitrarias del puerto viejo.
7. Llevar una contabilidad explícita de modos ocultos y factores de Schur.

Pruebas obligatorias:

- ejemplo Jordan del script;
- un polo físico +1e-7 debe seguir siendo inestable;
- uno exactamente cero físico debe ser BOUNDARY, no STABLE;
- varios valores beta y cambio de referencia angular no cambian espectro físico;
- modo interno no observable al puerto debe seguir apareciendo en el full-model ledger.

Revisar los 194 cruces Kundur ya localizados. Para cada uno, evaluar ambos lados, localizar cero físico y comparar full descriptor, A_q y nueva representación de puerto. Reportar cuáles pueden ahora analizarse en s=0 y cuáles siguen requiriendo el descriptor.

No imponer resultado 194/194. Un polo oculto o una singularidad física de referencia es un límite real.

Regenerar H y kappa sin el disco excluido para puntos auditados de IEEE-39/68. Cambios de etiqueta no son regresiones a ocultar; son correcciones del criterio físico.

## 4. BC02 — puerta de representación de la familia binaria

Separa estrictamente:

A. afinidad en parámetros de control para S fijo;
B. localidad de DeltaY_i en el puerto a equilibrio común;
C. afinidad de una realización dinámica común en los indicadores delta_i.

B puede ser cierta aunque C no lo sea.

Implementa un adapter de investigación con contratos explícitos:

- `equilibrate(portfolio, policy, uncertainty)`;
- `linearize_full(...)` devuelve E,A o bloques fx,fz,gx,gz, nombres y bases;
- `symmetry_generators(...)` derivados físicamente;
- `quotient_model(...)` con espectro preservado;
- `port_model(...)` con factor oculto h_S y contabilidad de polos;
- `affine_binary_model(...)` devuelve una representación comprobada O `NOT_APPLICABLE`;
- `operator_error_bound(cell, ...)` devuelve una cota válida O `UNKNOWN`.

No hacer que el adapter invente una aproximación cuando falta estructura. No igualar dimensiones con estados muertos o integradores. Un relleno estable debe probar equivalencia de vértices y su contribución conocida al determinante.

Para los 16 portfolios IEEE-39 y Kundur, y para los 4096 IEEE-68 ya disponibles, valida qué representación existe. Si sólo puedes obtener una interpolación multilineal a partir de todos los vértices, registra su costo: sirve para BC04 como prueba de principio, no como escalabilidad.

Contabilidad obligatoria:

`det P_S = h_S det T_S`.

En comparaciones de S frente al baseline, el winding de T_S/T_0 debe corregirse por h_S/h_0. No asumir que DeltaY tiene polos estables.

## 5. BC03 — certificados de ganancia y de cardinalidad

Implementa T3 y T3a de la nota.

Primero audita si los bloques normalizados Q_ij son propios y estables, y si la relación Q_S = Q[indices_S,indices_S] es válida. Si falla, no aplicar pequeña ganancia directamente; clasificar `NOT_APPLICABLE_UNSTABLE_PORT` y evaluar una representación estable alternativa cuando sea realizable.

Baselines científicos mínimos:

- pequeña ganancia sin escalado;
- pequeña ganancia escalada/Perron;
- condición mixta ganancia–fase de Huang et al., si sus hipótesis se cumplen;
- condición plug-and-play reciente de Hallinan–Lestas como comparación de alcance, implementando sólo una variante compatible y documentando las restantes;
- Nyquist exhaustivo optimizado del propio repositorio como verdad nominal.

No hace falta reproducir todos los papers completos para obtener el primer resultado, pero no comparar sólo contra SCR.

Por cada caso obtener rbar_ij con:

- cálculo H_infinity/KYP/hamiltoniano adecuado; o
- intervalos de frecuencia y cota Lipschitz + cola a alta frecuencia.

Un raster solo se guarda como `SAMPLED_SCREEN`.

Calcular rho(Rbar), los pesos d y c_r(d) para r=1..m. Reportar la mayor cardinalidad certificada y la cota inferior de kappa. Verificar cada certificado frente al oráculo exhaustivo para m=4,9,12.

Seleccionar antes de mirar resultados puntos base-estables de cada política: uno seguro, uno cerca de frontera y uno con incompatibilidades. Es esperado que algunos no sean certificables.

Métrica primaria: CERO falsos certificados seguros. Métricas secundarias: cobertura, slack, tiempo total y costo de cotas H_inf.

## 6. BC04 — Lyapunov por subcubos y certificados booleanos

### BC04-A: T4

Implementa primero el test barato:

- P de Lyapunov en el ancla estable;
- Xi = parte positiva de Ai.T P + P Ai;
- prueba de la desigualdad final con error residual incluido.

Después la versión SDP con P y Xi como incógnitas. No modificar el venv congelado. Si necesitas CVXPY y un solver SDP, crea un entorno separado y bloqueado; importa sólo matrices/configs. No depender de licencia comercial como único camino.

### BC04-B: T5 / ideal booleano

Reproduce el contraejemplo exacto 2x2 del script: dos vértices estables y promedio inestable. Debe pasar Lyapunov dependiente de delta y fallar cualquier claim de estabilidad de toda la caja continua.

En la familia de cuatro acciones prueba grados d=1,2,3, con un máximo de costo fijado antes. Las igualdades delta_i^2-delta_i=0 son restricciones exactas. Los términos SOS deben estar verificados, no tratados como sugerencias de optimización.

Compara:

1. Lyapunov común;
2. subcubo T4;
3. Lyapunov booleano;
4. caja continua mu/IQC cuando aplicable;
5. oráculo exhaustivo.

Objetivo falsable: el certificado booleano acepta alguna familia segura de decisiones físicas que el certificado de caja no acepta. Si no ocurre, no reclamar menor conservadurismo en redes eléctricas.

Si C no es afín, usar la representación multilineal exacta disponible sólo para prueba de principio, o trabajar en un descriptor local polinómico correctamente derivado. No forzar A0 + sum delta_i Ai.

Salida: matrices P/Gram, residuos, mínimos de desigualdades en todos los vértices pequeños y estado `VERIFIED_NUMERICAL_CERTIFICATE` o `INCONCLUSIVE`.

## 7. BC05 — búsqueda exacta sin monotonía

Implementa un árbol de decisiones sobre (F1,F0,J), con caché por hash del modelo y policy.

Reglas:

- una rama se poda como SAFE sólo con certificado válido T3/T4/T5/contorno;
- un test fallido no implica UNSAFE;
- si un conjunto conocido inestable está incluido en todos los miembros de una rama, esa rama puede excluirse de la búsqueda de testigos MÍNIMOS, sin etiquetar sus vértices como inestables;
- no declarar minimalidad por greedy de una eliminación; usar comprobación exhaustiva de inferiores, certificados o búsqueda por cardinalidad;
- no podar superconjuntos porque un subconjunto sea estable.

Pruebas negativas obligatorias:

- unsafe = { {1}, {1,2,3} }; el triple NO es mínimo aunque todos sus pares sean estables;
- estabilidad recuperada al añadir una acción;
- endpoint estable sin camino de pasos estables;
- hipergrafo exponencial de todos los conjuntos de tamaño r.

Comparar en cuatro, nueve y doce candidatos contra enumeración optimizada M2 y puerto/Nyquist M3. Contar builds, eigensolves, factorizaciones, SDPs, Lyapunov, tiempo CPU, wall-clock, memoria, setup y amortización. Los labels del oráculo no pueden guiar la búsqueda; sólo validar el resultado después.

Salida exacta sólo cuando cada rama restante esté certificada/excluida/evaluada. Si agota presupuesto, devolver [kappa_lower,kappa_upper], testigos encontrados y volumen desconocido, sin proclamar H completo.

No prometer subexponencialidad. Es una comparación de esfuerzo con garantía, no una conjetura de complejidad universal.

## 8. BC06 — robustez con cuantificadores correctos

Usa un conjunto físico de incertidumbre fijado: carga, disponibilidad PV, K/T de excitación, parámetros SG. Mantén corredores de despacho factibles; toda inviabilidad queda separada.

Distingue:

- kappa_nominal;
- mínimo observado en escenarios;
- cota inferior certificada de kappa_worst_case;
- testigo de kappa_worst_case obtenido con una realización concreta.

Para una garantía continua usar intervalos/cotas uniformes/SOS, no “pasó 2000 Monte Carlo”. Empieza en celdas pequeñas alrededor de puntos congelados y refina sólo donde el bound no cierre.

Para un conjunto FINITO de escenarios, construir H_exists = mínimos de la unión de H por escenario. Eso certifica sólo esos escenarios. H_forall NO es simplemente la intersección de hyperedges mínimas: evaluar la cuantificación sobre la inseguridad antes de minimizar.

Estadística: si se reporta Clopper–Pearson sobre tasa de éxito, utilizar muestreo iid o justificar el diseño; un único Latin Hypercube no son Bernoulli independientes. Los puntos de una malla tampoco lo son.

Preguntas principales:

- ¿kappa_nominal=4 pero kappa_wc menor bajo incertidumbre certificada?
- ¿qué testigos persisten y cuáles aparecen sólo en algunas celdas?
- ¿el diseño seguro para cualquier escenario conserva más/menos PV que el nominal y por cuánto?

No recalibrar la mitigación usando el holdout y llamarla una reparación congelada.

## 9. BC07 — diseño endpoint / orden específico / cualquier orden

Implementa T6.

A. Optimización de endpoint estable, mediante enumeración exacta para referencia pequeña.

B. Búsqueda de alguna secuencia estacionaria segura en el lattice.

C. Optimización independiente del orden con restricciones:

`sum(x_i for i in H) <= |H|-1` para cada hyperedge del H COMPLETO.

Si hay hyperedges desconocidas, usar un oráculo de separación y repetir hasta certificar; no afirmar seguridad sólo porque no contiene las ya encontradas.

Mantener restricciones de capacidad, Q/current, despacho y costo explícitas. Usar MW activos rastreados, no Sn en MVA.

Comprobar las implicaciones A/B/C en los resultados. No afirmar que un hitting set de nodos automáticamente dimensiona un condensador: ese soporte cambia las dinámicas y requiere un H nuevo o una certificación.

La salida útil es una decisión: qué reemplazos elegir, cuáles pueden retrasarse sin riesgo, cuándo se necesita orden específico y cuánto soporte físico debe conservarse. Validar los planes obtenidos contra todos sus prefijos/subconjuntos correspondientes.

## 10. BC08 — calibración local de closure

Implementar T7, con derivadas del Q reducido y autovectores izquierdos/derechos.

No buscar la frecuencia del puerto usando el eigenvalue full que después usarás como verdad de validación. La búsqueda de m_port debe ser independiente; el espectro full sólo se consulta al comparar.

En cruces simples y transversales medir:

- m_local;
- |f_s|;
- m_local/|f_s| contra |alpha|;
- error de frecuencia;
- diferencia entre mínimo local y global;
- lado estable/inestable obtenido mediante conteo, no por el signo de m.

Probar ambos lados de la frontera. Seleccionar antes casos IEEE-39, Kundur y, si hay fronteras localizadas, IEEE-68. Incluir un fold donde la hipótesis falla y no aplicar la ley transversal allí.

No llamar al margen eigenvalue-distance un radio de robustez. Validar el contraejemplo no normal del script.

## 11. BC09 — certificación de reutilización reequilibrada

F11 mostró que reutilizar ingenuamente el operador bajo PF unitario puede equivocarse. Mantener ese resultado.

Para cada celda pequeña:

- construir un ancla exacta;
- estimar y, cuando posible, certificar la cota de error de la familia A o del pencil;
- aplicar (25) o (29);
- permitir reutilización sólo si cierra la desigualdad con margen y errores numéricos incluidos;
- rechazar y resolver directamente cuando no cierra.

Distinguir bound analítico/intervalar de máximo medido en muestras. El segundo no es certificado.

Métrica primaria: cero errores en H entre celdas autorizadas. Métrica secundaria: tiempo neto contando cálculo de cotas. Si las cotas cuestan más que relinearizar, no reivindicar velocidad.

## 12. BC10–BC11 — puente no lineal y transición real, sólo tras cerrar lo anterior

Aplicar T8 en pocos casos claros: baseline, un portfolio certificado, uno reparado y un punto cerca de limitador. Derivar cotas de la no linealidad sobre la variedad algebraica regular. Incluir distancia a límites de corriente, capacidad Q, cambios de active set y singularidad algebraica.

Una prueba aleatoria de Hessianos es sólo un screen. La etiqueta de certificado requiere cota justificada.

Evaluar TDS desde condiciones dentro de una elipsoide certificada. Una simulación fuera que se mantenga estable no invalida una cota suficiente conservadora.

Para secuencias de commissioning:

- definir rampas y resets físicos;
- resolver consistencia algebraica tras cada cambio;
- comprobar inclusión de elipsoides (42) o un bound con dwell time y desplazamiento de equilibrio;
- simular al menos una ruta que pasa todos los endpoints pero viola la transición si se ejecuta demasiado rápido.

No llamar a un régimen slow commissioning “arbitrary switching”. No atribuir potencia activa sostenida a PV sin reserva/energía disponible. Respetar P_available, current limits y prioridades P/Q.

## 13. BC12 — auditoría final de novedad y artículo

La literatura mínima a revisar está en [L1]–[L10] de la nota; especialmente Hallinan–Lestas Parts I/II (arXiv:2609.07315/2609.07331), Huang et al. (2309.08037), el radio de estabilidad disperso (1810.10578), Lyapunov/SDP y certificados semialgebraicos de Parrilo.

Comparar lo que realmente se ejecutó, no afirmar que Nyquist no puede dar H: al evaluarlo por subset sí puede y ya está demostrado en F10.

Posibles contribuciones, sólo si los gates pasan:

1. conteo completo respetando la simetría y cruces por cero;
2. certificados binarios/cardinalidad y poda válida para una familia no monótona;
3. diseño seguro ante cualquier orden parcial, comparado con endpoint y secuencia específica;
4. extensión local no lineal como resultado secundario.

El teorema de constancia por regiones, el small-gain clásico y la existencia de Lyapunov polinómico en un conjunto finito son fundamentos, no claims de prioridad.

## 14. Ejecución por fases y presupuesto

No lanzar todos los bloques costosos simultáneamente. Usar recursos detectados y threadpool BLAS limitado. Nunca detener el job ajeno que ya exista en la máquina.

**Fase I:** BC00, BC01, BC02 y todos los tests sintéticos. Reportar discrepancias de datos/modos antes de continuar con claims afectados.

**Fase II:** BC03/BC04 sobre cuatro acciones y puntos preregistrados. Si no hay representación válida, no forzar la matemática; usar el camino alternativo explícito. Congelar resultados aunque los certificados no acepten nada.

**Fase III:** BC05 sobre 9 y 12 acciones, BC07 pequeño y BC08. Tiempo máximo por solver/caso, número de reinicios y tolerancias congelados antes.

**Fase IV:** BC06/BC09 y sólo después BC10–BC11. No condicionar el muestreo a conseguir un claim favorable.

Primera respuesta solicitada: tabla de disponibilidad de datos/APIs, status de BC00, ecuación exacta de la simetría y representación binaria aplicable. Después empezar Fase I. No pedir al usuario que elija buses nuevos para fabricar un fenómeno.

## 15. Entregables

- `BC00_DATA_AND_ORIGIN_AUDIT.md`;
- pruebas y código de quotient/relocation, certificates, search y planning;
- configs preregistradas, raw tables y manifests;
- matrices/cotas/certificados reutilizables;
- `BC_THEOREM_APPLICABILITY.md` con una fila por hipótesis y benchmark;
- `BC_FAILED_AND_UNKNOWN.md`;
- `BC_FINAL_EVIDENCE_TABLE.md`;
- `BC_WHAT_IS_NEW_AND_WHAT_IS_KNOWN.md`;
- lista de claims del artículo que siguen válidos, los corregidos y los aún no probados.

En cada resultado usar sólo uno de:

`PROVED_UNDER_ASSUMPTIONS`, `VERIFIED_ON_MODEL`, `CERTIFIED_FAMILY`,
`EMPIRICAL_ONLY`, `UNKNOWN`, `NOT_APPLICABLE`, `FALSIFIED`.

No editar el póster ni el manuscrito automáticamente. Detenerse al entregar la auditoría.
