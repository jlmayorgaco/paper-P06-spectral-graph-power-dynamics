# Óptimo local numéricamente certificado — cierre de Q2B

## Resultado y alcance

**PASS_LOCAL_SECURE — certificación numérica.** Se ejecutó co-diseño conjunto de las cinco retenciones presentes y las 20 ganancias PLL, con SQP explícito, múltiples picos activos y corrección de segundo orden. El modelo y el evento permanecen congelados. El óptimo es local para el **problema analítico de seguridad especificado**; la factibilidad no lineal se valida independientemente. Esto no demuestra optimalidad del problema con restricciones no lineales de trayectoria ni globalidad.

- Retención: **312.523288435 MW**; GFL: **5090.237801544 MW**, **94.215489%** del despacho original.
- Soporte: **30, 33, 35, 37, 39**. Mejora frente al candidato congelado Q2B anterior: **125.922994 MW** menos de SG.
- El punto KKT antes del resguardo de redondeo cuesta 312.523286071441 MW. El representante factible añade solamente 2.36369192e-06 MW.
- No se modificó el candidato anterior ni el modelo ExpN; 17/17 hashes y ambos archivos de dependencias coinciden. Sin commit, push o cambio de rama.

## Certificado local

| Prueba | Punto KKT | Representante congelado |
|---|---:|---:|
| Residuo primal | 0 | 0 |
| Estacionariedad | 1.40622e-07 | 8.27784e-07 |
| Complementariedad | 1.45101e-11 | 1.30214e-09 |
| Rango activo / restricciones | 23 / 23 | 23 / 23 |
| Condición del Jacobiano normalizado | 4.89299 | 4.89298 |

Todos los multiplicadores activos son positivos. Hay dos direcciones tangentes libres. Los autovalores de la Hessiana reducida son 9.052908 y 374.393104; la variación entre tres pasos de diferenciación es 0.009032. La menor curvatura positiva supera ampliamente esa variación. Objetivo dividido por 1000; epsilon y ganancias normalizadas a sus intervalos. Estos residuos no equivalen a una distancia en MW al óptimo global.

Es un **certificado numérico a tolerancias declaradas**, no una prueba con intervalos del cero exacto de las ecuaciones KKT. El pequeño resguardo factible deja residuos no nulos documentados. Se verifican LICQ, signos, complementariedad y SOSC; KKT por sí solo no habría bastado.

Las restricciones activas son beta, dos picos de frecuencia del bus 34 (aproximadamente 7.29445 y 14.62293 s) y las 20 cotas de ganancias. La restricción modal está holgada; su multiplicador es cero. No se presenta un modo modal activo ficticio ni una explicación BND de una restricción que no limita este óptimo.

## Factibilidad analítica completa

| Magnitud | Resultado |
|---|---:|
| Alpha, 143 polos físicos | -0.050057324502 1/s |
| Beta observada | 1.69912086867e-06 |
| Beta mínima cubierta por certificado | 1.69912069992e-06 |
| Pico de frecuencia | 0.499999994973 Hz |
| Frecuencia estacionaria | 0.481297135318 Hz |
| Pico de RoCoF | 0.167030273074 Hz/s |

La robustez se cubre con una desigualdad estricta bounded-real/Riccati a 256 bits: X positivo y residual negativo definido. Residuo de Riccati: 3.9982e-24; margen mínimo de la desigualdad: 3.3982e-14. Cubre toda frecuencia y la cola infinita. No depende de que una malla haya encontrado el máximo. Redondeo al más cercano: **NUMERICALLY_CERTIFIED**, no FORMALLY_CERTIFIED. El intento previo de cubrir la banda por subdivisiones se interrumpió por costo; no se usa como evidencia de certificación.

Para tiempo se combinaron raíces de las derivadas, cotas de Taylor sobre 1844 intervalos y una cota modal de la cola desde 60 s. La comprobación es numérica; el margen de redondeo se declara en el código. Las mallas densas son contraste y visualización.

## Validación independiente PowerDynamics

- Alpha PD: **-0.050057324250 1/s**; error con diseño: 2.527e-10 1/s.
- Emparejamiento biyectivo de todos los polos: error máximo 2.337e-09 1/s; error relativo de A en las mismas coordenadas físicas: 1.085e-15.
- Trim: residuo 6.574e-12; error P 3.979e-15 pu; error Q 2.723e-15 pu. Cargas 31/39 separadas e invariantes.
- Se obtuvo además el certificado Riccati directamente sobre la matriz reconstruida por PD.
- Evento no lineal bus16/100 MW: **frecuencia 0.46875546 Hz; RoCoF 0.17177144 Hz/s**. Ambos límites de 0.5 pasan. Frecuencia media final: 0.42021662 Hz.
- Error máximo de balance: 2.145e-10 MW; cambio Pdc: 0. Valores y trazas por evento en `PD/EVENTS.csv` y `PD/*SENSORS.csv`.
- No se declara protección por límite de corriente GFL ni seguridad para perturbaciones arbitrarias. El evento cambia Pset del ZIP exactamente como en el proyecto; el consumo efectivo y pérdidas se registran durante el evento. Ventana de medida T=0.5 s, sin escogerla después de ver el resultado.

El contraste local de 1 MW tiene aproximadamente 0.2705% de discrepancia máxima relativa de frecuencia. A 100 MW la discrepancia alcanza 13.7542%: la respuesta no lineal pasa los límites, pero no se afirma identidad lineal para una perturbación grande. Véase `PD/LINEAR_NONLINEAR_SCALING.csv`.

**Fuera del caso de diseño:** bus8_100MW: F=0.410549 Hz, R=0.512750 Hz/s; bus29_100MW: F=0.588571 Hz, R=0.281889 Hz/s. Los eventos de buses 8/29 son validación externa; sus fallos se conservan y limitan el alcance. No se incorporaron como restricciones ocultas ni se retocó el candidato. En el instante exacto del escalón se reconstruye el límite algebraico derecho manteniendo todos los estados diferenciales fijos: la interpolación de PD había combinado inicialmente parámetros posteriores con voltajes anteriores. La corrección afecta sólo la exportación de esa muestra, no la trayectoria integrada.

## Parámetros congelados

| Bus | epsilon | rho | SG MW | Kp | Ki |
|---|---:|---:|---:|---:|---:|
| 30 | 0.604447568 | 0.395552432 | 151.111892 | 7.853981634 | 986.960440109 |
| 31 | 0.000000000 | 1.000000000 | 0.000000 | 125.663706144 | 61.685027507 |
| 32 | 0.000000000 | 1.000000000 | 0.000000 | 125.663706144 | 986.960440109 |
| 33 | 0.078433012 | 0.921566988 | 49.569663 | 7.853981634 | 986.960440109 |
| 34 | 0.000000000 | 1.000000000 | 0.000000 | 125.663706144 | 986.960440109 |
| 35 | 0.077635538 | 0.922364462 | 50.463100 | 7.853981634 | 986.960440109 |
| 36 | 0.000000000 | 1.000000000 | 0.000000 | 7.853981634 | 986.960440109 |
| 37 | 0.059888862 | 0.940111138 | 32.339986 | 7.853981634 | 986.960440109 |
| 38 | 0.000000000 | 1.000000000 | 0.000000 | 7.853981634 | 986.960440109 |
| 39 | 0.107136992 | 0.892863008 | 29.038648 | 125.663706144 | 61.685027507 |

SHA-256 del candidato: `66513a8d3a1b4cc37c9d6f2b1a5824371fed2d0c5621d2ccdb7540490fe59124`.

## Cambios de arquitectura y globalidad

Se construyeron límites exactos unilaterales para insertar SG en 31, 32, 34, 36 y 38. Los cinco violan estrictamente beta en omega=0; la evaluación del inverso fue a 256 bits. Las matrices con estados adicionales se proyectan mediante una coisometría sobre el sistema anterior, de modo que su norma resolvente no puede ser menor. Esto excluye también inserciones simultáneas infinitesimales por continuidad, bajo el modelo de incertidumbre congelado. Los residuos se conservan en `ADJACENT_EXACT_LIMITS.csv`. La afirmación es local y numérica; no excluye inserciones finitas ni ramas desconectadas.

Para el problema global, la única cota inferior admitida aquí es **L=0 MW**. La cota superior factible es **U=312.523288435 MW**. Por tanto:

    0 <= J_global* <= 312.523288435 MW
    0 <= J_candidato - J_global* <= 312.523288435 MW

La mejora global todavía posible está acotada por **5.784511 puntos porcentuales del despacho original**. Esto es una cota máxima, no una estimación de cuánto falta. No se puede calcular un error porcentual respecto de J_global* usando una cota inferior cero.

El intento de cota mediante eliminación DC exacta y polinomios de Bernstein no dio una cota útil y carece de envolventes con redondeo exterior. No se usa el heurístico afín de 57.2 MW. La identidad estática se comprobó en 100 casos (error máximo 1.64e-11 Hz), pero eso no convierte sus coeficientes Float64 en una prueba global. Los 1024 soportes siguen abiertos en cuanto a completitud global; uno contiene este punto local certificado. Véase `GLOBAL_SUPPORT_LEDGER.csv`.

## Relevancia para el paper

La contribución que estos datos sostienen es un co-diseño físico de retención SG y PLL locales, con restricciones espectrales, robustas y temporales, evidencia de optimalidad local y validación independiente. SQP, la corrección de segundo orden y el lema bounded-real son herramientas conocidas; no se reclaman como teoremas nuevos. El resultado nuevo del caso es la asignación de seguridad y su comprobación reproducible. La novedad frente a toda la literatura requiere una comparación bibliográfica específica.

Los valores marginales ponderados son aproximadamente uno en los cinco buses interiores. Buses 30 y 39: domina la frecuencia; buses 33,35,37: domina beta. La explicación está en `MARGINAL_SECURITY_VALUES.csv`, con signo coherente con minimizar SG retenido.

## Tiempo, obstáculos y reproducción

La ejecución SQP registrada tardó 574.291 s sin compilación, en 68 iteraciones externas. La certificación robusta tardó 10.533 s y la temporal 1.128 s. Los contadores completos de evaluaciones del SQP no se guardaron en esa ejecución; el código ahora los registra al reproducirla. No se inventa ese dato. Las primeras pruebas fallidas de gradientes, corrección Newton y cobertura por subdivisión se conservan para trazabilidad.

Se conserva la secuencia previa fallida en `LOCAL_CERTIFICATE`; este directorio la supera con evidencia nueva, sin reescribirla. Reproducción: véase `REPRODUCE.md`. Falta una cota global fuerte, completar ramas y una prueba formal con intervalos. Esas limitaciones no anulan el resultado local.

Referencias metodológicas: [Boyd et al., Linear Matrix Inequalities in System and Control Theory](https://web.stanford.edu/~boyd/lmibook/), [Boyd, Balakrishnan y Kabamba, norma H-infinity](https://stanford.edu/~boyd/papers/bisection_hinfty.html), [Gould y Robinson, SQP y corrección de segundo orden](https://optimization-online.org/2008/12/2192/).
