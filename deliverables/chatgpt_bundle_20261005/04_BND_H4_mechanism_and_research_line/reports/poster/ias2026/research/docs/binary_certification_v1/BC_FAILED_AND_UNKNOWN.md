# BC — Lo que falló y lo que sigue desconocido (Fase I)

Los fallos se reportan con la misma prioridad que los éxitos. Nada de esta
lista se reintentó con otro criterio para cambiar el resultado.

## Falló (`FALSIFIED` o negativo)

| # | qué | evidencia | consecuencia |
|---|---|---|---|
| F1 | Puerto **original** evaluado en `s = 0` con el criterio preregistrado | 24/194 evaluables. El winding coincide en 134; el signo cambia en 77; el signo medio dice «estable» en 138 aunque el lado alto es inestable en los 194 (`BC01_kundur_crossings.csv`) | El cociente `det T_S / det T_0` no sirve en el origen. El cero doble estructural (rotación y modo neutro) se parte con el ruido de diferencias finitas |
| F2 | Conteo crudo del pencil del descriptor en el lado estable de los cruces | 166/194 dan 1 o 2 autovalores «positivos» donde el conteo físico es 0 | Sin cociente, un conteo sin disco es falso; con disco, oculta el cruce |
| F3 | Estabilidad asintótica estricta de cualquier portfolio sin gobernador | `A_q` tiene un cero exacto (modo neutro) en el 100 % de los casos con `D = 0` | «Estable» significa estable módulo el modo neutro. Ver BC_THEOREM_APPLICABILITY |
| F4 | Ranker «MW reemplazados» (N16/E39) | AUC de 0.77–0.86 a 0.57–0.73 con MW activos | El claim se detiene |
| F5 | Afinidad de la `A` reducida en `delta` | no afinidad relativa de 0.36 en IEEE-39 y de 4.1 en Kundur (BC02); solo el jacobiano del descriptor es afín | T4 no es aplicable tal como está enunciado |
| F6 | Padding C1 (fantasmas) como realización estable en IEEE-68 | `Re` máxima de los fantasmas entre +3e-15 y +7e-14: cero marginal | en IEEE-68 solo se admite el relleno C2 |
| F7 | Tamaño de soporte del incremento de puerto en `BC02_summary.json` (umbral 1e-8) | informa entre 4 y 54 canales en IEEE-39 y 12 en Kundur; la comprobación directa muestra que el exceso es la tolerancia del flujo de carga en la slack | la métrica congelada no mide el soporte físico; se reporta junto con la comprobación directa |

## Desconocido (`UNKNOWN`)

| # | qué | por qué | qué lo resolvería |
|---|---|---|---|
| U1 | 212 etiquetas G1 de IEEE-39 (camino rápido, `g <= 0.005`) | un par oscilatorio (2.4–4.4 rad/s) con `d_axis` menor que `10 eps` del ensamblaje afín | evaluar esos puntos por el camino directo (jacobiano exacto), BC05 |
| U2 | 428 puntos IEEE-39 con algún subconjunto no resuelto (424 con `g <= 0.005` y 4 aleatorios) | igual que U1 | igual que U1 |
| U3 | 1 subconjunto no resuelto en Kundur | `d_axis` por debajo del umbral | camino directo con un paso menor o aritmética extendida |
| U4 | `d_axis` como certificado continuo | es malla más refinamiento: criterio numérico, no intervalos | cota de Lipschitz de `sigma_min` entre muestras, o aritmética de intervalos |
| U5 | Cota del error del operador (`operator_error_bound`) | no hay cota analítica ni por intervalos de `A(h) - A_exacto` | derivadas analíticas o automáticas (NL) o intervalos |
| U6 | Puerto reubicado en IEEE-39 e IEEE-68 | solo se ejecutó en los cruces de Kundur | BC03 |
| U7 | Familia de doce plantas del IEEE-68 (4096 vértices × 3) | aplazada a BC05 (preregistrado) | BC05 |
| U8 | Carácter de novedad del puerto reubicado | la reubicación de rango uno (Brauer) es clásica; lo aplicado es su uso dentro del bloque de dispositivos para ver cruces aperiódicos por el puerto | revisión bibliográfica dirigida; hasta entonces no se reclama |
| U9 | T3, T3a, T4 (versión descriptor), T5, (29), búsqueda sin monotonía, diseño seguro (T6 y (32)), T7, T8, transiciones | Fases II–IV no ejecutadas | BC03–BC12 |
| U11 | Ceros marginales en los bloques de dispositivo de IEEE-68 (`Re` máx. 2e-13) | no se identificó qué estado los produce | identificarlo antes del contorno de (16) en BC03 |
| U12 | La no afinidad de 5e-7 en IEEE-68 | norma relativa dominada por las entradas rígidas; no prueba afinidad | medirla en coordenadas escaladas |
| U10 | Que el relleno C2 preserve la simetría | la rotación no debe actuar sobre los estados de relleno | comprobar `A R = 0` con `R` nulo en el relleno (BC04) |

## Desviaciones del prerregistro (declaradas)

1. **Puerto reubicado.** Se añadió después de congelar `BC00_config.yaml`,
   cuando el puerto original falló (F1). El criterio original se reporta tal
   como estaba congelado; el reubicado es un análisis adicional, no un
   sustituto.
2. **Clasificador.** El diseño se fijó en la prueba de humo de 8 casos, antes
   de recalcular cualquier `H`: balanceo, retirada exacta de estados muertos,
   `known_zero`. No se cambió después de ver resultados de `H`.
