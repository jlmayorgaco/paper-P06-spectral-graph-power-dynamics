# Revisión de jurado · Póster IAS 2026

## Alcance

Revisión del PDF de 36 × 48 pulgadas, de las quince vistas modulares, del código LaTeX, de `generated/results.tex`, de los seis CSV incluidos y de las cinco referencias del pie. También se consultaron **en modo lectura** el generador y la auditoría de afirmaciones conservados en el checkout principal (`reports/poster/ias2026/generate_poster_assets.py` y `build/claim_audit.md`), junto con las tablas originales de F10, políticas y campaña. Se evaluó una lectura a distancia (título y resultados), una lectura de un minuto (problema, caso y consecuencia) y una lectura técnica de cerca (ecuaciones, denominadores y límites). Esta revisión comprueba consistencia editorial y visual; no reproduce las simulaciones originales.

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
| 5 minutos | Factorización, barrido de frontera, sensibilidad a política, límites del conjunto sintético y comprobaciones del modelo. | Faltan archivos de procedencia para auditar independientemente varios recuentos. |

## Evaluación por pieza

| Pieza | Diagnóstico y decisión editorial |
| --- | --- |
| Header blanco y verde | Buena identidad y un título dominante. El subtítulo define el campo; las siglas SG/IBR se explican visualmente en el módulo 1. Se conservó. |
| 01 Motivación | Pregunta concreta y diagrama inmediato. Se conservó porque funciona como entrada para un asistente no especialista. |
| 02 Resultado principal | Los cuatro buses destacados y el 15/15 forman la evidencia principal. Se eliminó la repetición de números en el pie de la red y el código interno de los ajustes de control. Los 39 rótulos grises son información para lectura cercana, no para la lectura a distancia. |
| 03 Destino y ruta | La secuencia A/B/C y la implicación están bien jerarquizadas. Se sustituyó “F10” por un rótulo que describe qué miden las barras y mantiene el denominador 52. “Las recíprocas pueden fallar” expresa el resultado general sin atribuir ambos contraejemplos al caso IEEE-39. |
| Franja de resultados | Útil como resumen a distancia, aunque repite valores del módulo 2. Se cambió “critical frequency” por “target-band frequency”: la auditoría local indica que la identidad modal sigue sin resolver. |
| 04A Factorización | Las ecuaciones tienen espacio y una secuencia visible. Se cambió la nota “M1” por un límite comprensible: la reducción de un modo aún no supera el criterio estricto. |
| 04B Frontera | El gráfico contrasta retorno local y colectivo con colores y formas. La escala logarítmica y el valor de frontera sustentan la lectura; se conservó. |
| 05 Política | El atlas muestra cambio estructural, pero las etiquetas F7A/F7B/F7C eran códigos internos. Ahora los tres recuentos se alinean con etiquetas de barrido A/B/C. Falta documentar qué distingue esos barridos. |
| 06 Conjunto sintético | Los tres porcentajes se solapan y suman más de 100 %. Se indica expresamente que las categorías se superponen y que los totales son provisionales; el denominador 967 permanece visible. |
| 07 Reajuste | Muestra beneficio y deterioro, una fortaleza de honestidad científica. Se identifica el carácter sintético de los recuentos, se amplían las siglas de la comprobación temporal y se muestra 30/90 para los fallos. La cifra de rescates mide otro resultado y no es una tercera parte de la partición 594/373. |
| 08 Rutas | La formulación ahora exige añadir una unidad por paso y comprobar $\alpha_\perp(S_t,\theta)<0$ en cada prefijo. La condición previa de evitar todas las hiperaristas en cada $S_t$ era demasiado fuerte: con una política fija y conjuntos anidados equivalía a la garantía C sobre el destino, por lo que no distinguía órdenes. El panel incluye el resultado observado: 3 de 327 destinos estables dependen del orden. Las rutas roja y verde se describen también con palabras. |
| 09A Comprobaciones | Los dos denominadores son visibles. Se sustituyó TDS por “time-domain” y se explicitan los 30 fallos entre 90 ejecuciones. |
| 09B Alcance | Es la defensa necesaria ante preguntas sobre maniobras físicas y EMT. Se conservó sin prometer validación que el póster no demuestra. |
| Footer | Conclusión recordable, contacto y referencias identificables. El margen blanco inferior se conserva. |

## Comprobación de cifras y procedencia

| Afirmación | Comprobación disponible aquí | Estado para un jurado |
| --- | --- | --- |
| 15/15 subconjuntos estables; $H_4=\{30,33,35,37\}$ inestable; $+0.127\,\mathrm{s}^{-1}$ y 0.6223 Hz | Los valores están en `generated/results.tex`. La auditoría del checkout principal los atribuye a entradas congeladas y aclara que 0.6223 Hz es una frecuencia de banda objetivo, sin identidad modal resuelta. | Caso nominal auditado; el paquete de reproducción no está incluido en esta rama. |
| Frontera en $g\approx0.207681$ | `generated/tables/f2_local_collective_closure.csv` contiene el punto con $\sigma_{\min}$ local 0.489235 y colectivo $4.39367\times10^{-8}$; la gráfica coincide. | Trazable al CSV incluido. |
| 22/52, 31/52 y 38/52 coincidencias | Coinciden con las filas B3/B4/B5 de `research/results/20260911_BASELINE_COMPARISON.csv` en el checkout principal. | Recuento cruzado con la tabla fuente local; esa tabla no está en esta rama. |
| 30/36/16 hipergrafos distintos y mapa de políticas | La auditoría del checkout principal los vincula a F7A/F7B/F7C y a `generated/tables/policy_atlas_source.csv`; la rama solo incluye el PNG y las macros. | Procedencia localizada; falta empaquetar la definición de los tres barridos con el póster. |
| 967/1000 válidos; 190/967, 555/967 y 365/967 | Se comprobaron los porcentajes a una décima: 96.7 %, 19.6 %, 57.4 % y 37.7 %. `MC_EVENT_RATES.csv` del checkout principal confirma 190 y 555; las tres categorías se superponen. | Aritmética consistente; totales expresamente provisionales. |
| 594 mejoran, 373 empeoran, 186 rescates | Los CSV de mejora y deterioro de la rama contienen 594 y 373 filas; suman los 967 casos válidos. `MC_EVENT_RATES.csv` del checkout principal confirma el recuento de 186. | Recuentos cruzados con fuentes locales; revisión de campaña aún pendiente. |
| 100/100 entre códigos; 56/60 temporal; 30/90 fallos | Las macros coinciden con `build/claim_audit.md` del checkout principal: acuerdo entre códigos sobre el mismo modelo, cuatro desacuerdos a frecuencia cero y 30 fallos de solver. Las trazas incluidas tienen 1500 filas cada una. | Alcance y limitaciones explícitos; registro completo no empaquetado en esta rama. |
| 3/327 destinos estables dependen del orden | `research/docs/FINAL_TRANSACTION_THEORY_AND_EVIDENCE.md` y `research/experiments/final_closure/FC10_monotone_and_paths.py` del checkout principal documentan el censo de 512 carteras y el caso $T=\{30,32,33,34,35\}$: al convertir 34 al final aparece el prefijo inestable $\{30,32,33,35\}$; otro orden llega al mismo destino estable. Ningún destino estable del censo carece de ruta segura. | El panel usa el recuento observado; el contraejemplo A sin B pertenece a la teoría general, no a este censo. |

La rama del póster contiene exportaciones suficientes para compilar, pero no el generador, la auditoría original ni la mayoría de las tablas que lo alimentaron. Esos archivos sí existen como trabajo local no integrado en el checkout principal, incluidos `build/claim_audit.md`, `generate_poster_assets.py` y `research/configs/ias2026/ieee39_network.json`. **Prioridad antes de entregar un paquete reproducible:** seleccionar y versionar el manifiesto, las tablas por caso y el generador con rutas reproducibles. No conviene convertir “provisional” en “validado” antes de cerrar la revisión de campaña ni afirmar identidad modal para 0.6223 Hz.

## Referencias

Las referencias del pie son identificables y sus datos bibliográficos principales concuerdan con [el registro IEEE de Hatziargyriou et al.](https://ieeexplore.ieee.org/abstract/document/9286772), [el registro institucional del artículo de Kundur et al.](https://orbi.uliege.be/handle/2268/3467), [el registro bibliográfico del artículo de Dörfler et al.](https://pubmed.ncbi.nlm.nih.gov/23319658/) y [la ficha editorial del libro de Milano](https://link.springer.com/book/10.1007/978-3-642-13669-6). Estas obras aportan contexto; no son la procedencia de los resultados numéricos originales del póster.

## Verificación de cierre

Las quince piezas se compilaron con `-StrictLayout`; el PDF final es de una página de 36 × 48 pulgadas. Se inspeccionaron el PNG completo y las vistas de las piezas modificadas, sin texto cortado ni elementos superpuestos. La aprobación visual no sustituye la reproducción de los resultados originales señalados arriba.
