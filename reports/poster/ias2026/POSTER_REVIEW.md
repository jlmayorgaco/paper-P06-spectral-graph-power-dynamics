# Revisión de jurado · Póster IAS 2026

## Alcance

Revisión del PDF de 36 × 48 pulgadas, de las quince vistas modulares, del código LaTeX, de las cinco referencias del pie y de los archivos fuente congelados ahora incluidos en la rama. Se evaluó una lectura a distancia (título y resultados), una lectura de un minuto (problema, caso y consecuencia) y una lectura técnica de cerca (ecuaciones, denominadores y límites). `generate_poster_assets.py` regenera los valores y figuras; `audit_poster_evidence.py` comprueba hashes y recuentos. La revisión reproduce las exportaciones del póster, no las simulaciones originales.

## Contrato científico que comunica el póster

- **Problema:** reemplazos SG→IBR que parecen aceptables por separado pueden formar una coalición inestable.
- **Resultado central:** en el caso IEEE-39 mostrado, los 15 subconjuntos propios de los cuatro candidatos son estables y la coalición completa no lo es.
- **Mecanismo expuesto:** la factorización separa factores locales y retorno colectivo; el barrido mostrado aproxima la singularidad del término colectivo mientras los factores locales permanecen regulares.
- **Consecuencia:** evaluar solo el destino final no garantiza una secuencia de implementación segura; debe comprobarse la estabilidad de cada prefijo. La condición de no contener ninguna hiperarista caracteriza la garantía más fuerte de que **todos** los órdenes son seguros.
- **Límite:** “seguro” significa estable tras reequilibrar el modelo fasorial DAE. El póster no modela maniobras físicas, EMT ni límites de corriente o enlace DC.

## Lectura como visitante y como jurado

| Tiempo de lectura | Lo que llega con claridad | Riesgo restante |
| --- | --- | --- |
| 10 segundos | Título, cuatro buses rojos, 15/15 frente a la coalición roja y conclusión del pie. | La novedad exacta de “hipergrafo mínimo” exige acercarse. |
| 1 minuto | Problema → caso IEEE-39 → diferencia entre destino y ruta → orden de reemplazo. | Los paneles de mecanismo y políticas requieren conocimientos previos de dinámica de potencia. |
| 5 minutos | Factorización, barrido de frontera, sensibilidad a política, límites del conjunto sintético y comprobaciones del modelo. | La validación independiente de otro modelo dinámico sigue fuera del alcance de este póster. |

## Evaluación por pieza

| Pieza | Diagnóstico y decisión editorial |
| --- | --- |
| Header blanco y verde | Buena identidad y un título dominante. El subtítulo define el campo; las siglas SG/IBR se explican visualmente en el módulo 1. Se conservó. |
| 01 Motivación | Pregunta concreta y diagrama inmediato. Se conservó porque funciona como entrada para un asistente no especialista. |
| 02 Resultado principal | Los cuatro buses destacados y el 15/15 forman la evidencia principal. Se eliminaron el código interno de los ajustes de control y la frecuencia de 0,6223 Hz, cuya identidad modal entre análisis no está resuelta. Los 39 rótulos grises son información para lectura cercana. |
| 03 Destino y ruta | La secuencia A/B/C y la implicación están bien jerarquizadas. Se sustituyó “F10” por un rótulo que describe qué miden las barras y mantiene el denominador 52. “Las recíprocas pueden fallar” expresa el resultado general sin atribuir ambos contraejemplos al caso IEEE-39. |
| Franja de resultados | Resume tres rasgos del caso nominal y, en la cuarta celda, los 3/327 destinos dependientes del orden del censo IEEE-39 con política matched; el alcance de ambos resultados queda distinguido. Se retiró 0,6223 Hz. |
| 04A Factorización | Las ecuaciones tienen espacio y una secuencia visible. La nota delimita el resultado a la factorización de orden completo; no se muestra una equivalencia reducida que aún no supera el criterio estricto. |
| 04B Frontera | El gráfico contrasta retorno local y colectivo con colores y formas. La escala logarítmica y el valor de frontera sustentan la lectura; se conservó. |
| 05 Política | El atlas A muestra cambio estructural en el plano $g\times k$ con $t=1.5,h=1$. Los recuentos A/B/C de hipergrafos de espectro completo indican sus planos y parámetros fijos: B varía $g\times t$ con $k=h=1$ y C varía $g\times h$ con $k=t=1$. La secuencia de tamaños corresponde a una línea de B con $t\approx0.852$. |
| 06 Conjunto sintético | El panel ahora muestra la partición exacta de los 967 casos válidos: 412 sin bloqueador, 190 con $H_4$ mínimo y 365 con triples mínimos. Los últimos dos suman 555 con bloqueador colectivo. Así cada porcentaje tiene una interpretación exclusiva y el subtotal queda explícito. |
| 07 Reajuste | Muestra beneficio y deterioro, una fortaleza de honestidad científica. Los 594 casos que mejoran y 373 que empeoran particionan los 967 válidos; los 186 rescates son un subconjunto de las mejoras. La nota distingue consistencia con archivos congelados de validación independiente. |
| 08 Rutas | La formulación exige añadir una unidad por paso y comprobar $\alpha_\perp(S_t,\theta)<0$ en cada prefijo. Evitar todas las hiperaristas en cada $S_t$ sería la garantía C y no distinguiría órdenes. El panel identifica el censo IEEE-39 de política matched, distinto del caso nominal P4: 3 de 327 destinos estables dependen del orden. Las rutas roja y verde se describen también con palabras. |
| 09A Comprobaciones | Los dos denominadores son visibles. Se sustituyó TDS por “time-domain” y se explicitan los 30 fallos entre 90 ejecuciones. |
| 09B Alcance | Es la defensa necesaria ante preguntas sobre maniobras físicas y EMT. Se conservó sin prometer validación que el póster no demuestra. |
| Footer | Conclusión recordable, contacto y referencias identificables. El margen blanco inferior se conserva. |

## Comprobación de cifras y procedencia

| Afirmación | Comprobación disponible aquí | Estado para un jurado |
| --- | --- | --- |
| 15/15 subconjuntos estables; $H_4=\{30,33,35,37\}$ inestable; $+0.127\,\mathrm{s}^{-1}$ | `generated/results.tex`, `POSTER_NUMBERS_FINAL.csv`, `20260911_CLAIM_MATRIX.csv` y `evidence/claim_audit.md`. El generador y las entradas están incluidos en esta rama. | Caso nominal P4, auditado a nivel de exportaciones. La frecuencia de 0,6223 Hz se retiró del póster por identidad modal no resuelta. |
| Frontera en $g\approx0.207681$ | `generated/tables/f2_local_collective_closure.csv` contiene el punto con $\sigma_{\min}$ local 0.489235 y colectivo $4.39367\times10^{-8}$; la gráfica coincide. | Trazable al CSV incluido. |
| 22/52, 31/52 y 38/52 coincidencias | Se recalculan de las filas B3/B4/B5 de `research/results/20260911_BASELINE_COMPARISON.csv` incluido. | Tabla y generador versionados. |
| 30/36/16 hipergrafos distintos y mapa de políticas | La rama incluye el CSV de 120 390 puntos del mapa A y la matriz auditada que registra A/B/C. La definición de los tres planos se muestra en el panel 5 y en `EVIDENCE.md`. | El mapa A se regenera; los recuentos B/C se cruzan con la matriz congelada, sin recalcular esas dos mallas completas. |
| 967/1000 válidos; partición 412/190/365; subtotal 555/967 | `SCENARIO_METRICS.csv` y `MC_EVENT_RATES.csv` incluidos se auditan fila por fila. 412+190+365=967 y 190+365=555; los porcentajes de la partición suman 100 % salvo redondeo. | Recuento interno reproducible de un conjunto sintético congelado; interpretación científica pendiente de aprobación final. |
| 594 mejoran, 373 empeoran, 186 rescates | `SCENARIO_METRICS.csv` y `MC_EVENT_RATES.csv` incluidos confirman 594+373=967 y que los 186 rescates están entre las mejoras. Los CSV editables de las gráficas conservan 594 y 373 filas. | Recuento interno reproducible; no es una validación de otro modelo. |
| 100/100 entre códigos; 56/60 temporal; 30/90 fallos | Las macros coinciden con `IAS26-FINAL_EVIDENCE_STRIP.csv` y `IAS2026_POSTER_READINESS.json` incluidos. Los cuatro desacuerdos a frecuencia cero y los 30 fallos de solver constan en la auditoría versionada. | Alcance y limitaciones explícitos; paridad entre códigos del mismo modelo. |
| 3/327 destinos estables dependen del orden | `FC10_census_lattice.csv` y `FC10_summary.json` incluidos documentan 512 carteras bajo la política reactiva matched. Para $T=\{30,32,33,34,35\}$, convertir 34 al final deja el prefijo inestable $\{30,32,33,35\}$; otro orden llega al mismo destino estable. | Recuento reproducible; ningún destino estable del censo carece de ruta segura. El contraejemplo A sin B pertenece a la teoría general, no a este censo. |

La rama incluye ahora el generador, 21 fuentes congeladas con sus SHA-256 y una auditoría ejecutable de las cifras visibles. Esto cierra la procedencia **del artefacto póster**. Sigue abierto un trabajo científico diferente: reproducir los cálculos dinámicos originales y completar una validación independiente del modelo. Los resultados del conjunto sintético no se presentan como probabilidades ni como evidencia de campo, y la identidad modal de 0,6223 Hz no se afirma.

## Referencias

Las referencias del pie son identificables y sus datos bibliográficos principales concuerdan con [el registro IEEE de Hatziargyriou et al.](https://ieeexplore.ieee.org/abstract/document/9286772), [el registro institucional del artículo de Kundur et al.](https://orbi.uliege.be/handle/2268/3467), [el registro bibliográfico del artículo de Dörfler et al.](https://pubmed.ncbi.nlm.nih.gov/23319658/) y [la ficha editorial del libro de Milano](https://link.springer.com/book/10.1007/978-3-642-13669-6). Estas obras aportan contexto; no son la procedencia de los resultados numéricos originales del póster.

## Verificación de cierre

Las quince piezas se compilaron con `-StrictLayout`; el PDF final es de una página de 36 × 48 pulgadas. Se inspeccionaron el PNG completo y las vistas de las piezas modificadas, sin texto cortado ni elementos superpuestos. El generador y `audit_poster_evidence.py` se ejecutaron desde la rama. La aprobación visual y la auditoría de las exportaciones no sustituyen la reproducción de las simulaciones originales.
