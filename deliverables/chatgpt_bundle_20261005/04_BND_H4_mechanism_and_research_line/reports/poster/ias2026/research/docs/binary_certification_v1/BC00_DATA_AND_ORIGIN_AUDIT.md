# BC00–BC02 — Auditoría de datos, origen y representación (Fase I)

Rama `ias2026/binary-certification-v1`, creada desde `9a216145` (snapshot de
G1–G6). Los criterios se congelaron en `7a77da30`
(`configs/binary_certification_v1/BC00_config.yaml`) antes de recalcular
cualquier hipergrafo. No se sobrescribió ningún resultado de F1–F12 ni de G1–G6
y no se tocó el póster. El entorno congelado (`.venv/tx3-analysis`) no se
modificó. Las comprobaciones simbólicas se ejecutaron en un entorno aparte,
`.venv/bc-cert` (`configs/binary_certification_v1/ENV_bc-cert.txt`).

Etiquetas: `PROVED_UNDER_ASSUMPTIONS`, `VERIFIED_ON_MODEL`,
`CERTIFIED_FAMILY`, `EMPIRICAL_ONLY`, `UNKNOWN`, `NOT_APPLICABLE`,
`FALSIFIED`.

## 0. Nota teórica y script sintético

`spectral_portfolio_theory_v1/theory_checks.py` se volvió a ejecutar sin
cambios (`results/BC/BC00/theory_checks/theory_check_results_rerun.json`). Sus
404 comprobaciones son **sintéticas**: matrices pequeñas construidas a mano. No
son ensayos IEEE ni validación estadística. Reproducen T1/V1, la identidad
`det P = h det T`, V2 (`kappa >= 5`), V3, T5, el no-go T2, el ejemplo diagonal
de (33) y la conmutación inestable de §15.1.

## 1. Datos y API disponibles

| elemento | disponible | dónde | límite declarado |
|---|---|---|---|
| DAE IEEE-39 (L0): `solve_case`, `Ieee39Dae.f/g`, jacobianos por diferencias centrales | sí | `src/ibr_cycles/models/ieee39_*.py` | AVR IEEEX1 fusionado (KA y TE); sin gobernadores ni límites |
| Mapas congelados F7A/F7B/F7C y puntos F8 | sí | `results/F7`, `results/F8/F8_points.json` | camino rápido afín en `(g, k)` fuera de los nodos exactos |
| H, kappa y fronteras de G1 (con disco `abs(s) <= 1e-3`) | sí | `results/G1` | es el criterio auditado aquí |
| Kundur (F12): K12A, K12B | sí | `experiments/ias2026/F12_kundur.py`, `configs/kundur` | SEXS sin límites; `D = 0` |
| IEEE-68 (G3) | sí | `experiments/ias2026/G3_ieee68.py`, `configs/ieee68` | 4 candidatos: resultado negativo (`H` vacío). Doce plantas: solo 3 puntos |
| ANDES | parcial | reconciliación F1 | solo estática armonizada con el PSS desactivado. No valida GFL ni PSS dinámicos |
| Afinidad en `theta` (F11) | sí | F7/F11 | en los parámetros de control con `S` fijo, **no** en `delta` |
| EMT, límites, anti-windup, enlace DC, gobernadores | no | — | fuera del modelo congelado |
| Cotas rigurosas del operador (intervalos, `H_inf` verificada) | no | — | `operator_error_bound` devuelve `UNKNOWN` |
| sympy | solo en `.venv/bc-cert` | — | el venv congelado no se tocó |

Adaptadores: `src/ibr_cycles/certification/adapter.py` (`Ieee39Adapter`,
`KundurAdapter`, `Ieee68Adapter`). Contrato: `equilibrate`, `linearize_full`,
`symmetry_generators`, `quotient_model`, `port_model`, `affine_binary_model`,
`operator_error_bound`. Pruebas: `tests/test_bc_adapter.py`.

## 2. BC00-A — Unidades

`replaced_mw` suma la potencia nominal de los convertidores (`Sn`, en MVA). No
es potencia activa. Detalle en `results/BC/BC00/BC00_units_*.csv` y
`BC00_units_summary.json`.

| magnitud (flagship 30+33+35+37) | valor |
|---|---|
| columna `replaced_mw` (etiquetada MW) | 4270.7 (son **MVA**) |
| potencia activa desplazada | 2096.6 MW |
| potencia reactiva desplazada | 420.6 Mvar |
| nominal activa de las máquinas (`pmax`) | 2983 MW |
| razón reportado / activo | 2.04 |

Efecto en el ranker «MW reemplazados» de E39 (AUC; más alto significa más
inestable):

| tamaño | n | inestables | columna antigua (MVA) | MW activos reales |
|---|---|---|---|---|
| 4 | 126 | 10 | 0.813 | 0.572 |
| 5 | 126 | 51 | 0.773 | 0.636 |
| 6 | 84 | 78 | 0.861 | 0.726 |

Los espectros **no** cambian: la etiqueta se asigna después del cálculo
dinámico.

## 3. BC00-B — Origen, ceros y simetría

### 3.1 Ecuación exacta de simetría

Las ecuaciones son invariantes a una rotación común de todos los ángulos. El
generador se deriva de las ecuaciones, no del espectro:

- `R_x = 1` en cada ángulo de rotor (`delta`) y en cada ángulo de PLL
  (`theta_pll`), y 0 en el resto;
- `R_z = j V` en cada tensión de barra, es decir `(vx, vy) -> (-vy, vx)`.

Se cumple

    f_x R_x + f_z R_z = 0,    g_x R_x + g_z R_z = 0    =>    A R_x = 0.

Si **todas** las máquinas tienen `D = 0` y no hay gobernadores, existe además un
compañero de Jordan `w`, el desplazamiento uniforme de frecuencia: velocidad
o deslizamiento `+1`, integrador del PLL `+omega_B`, filtro de inercia
sintética `+1` y, en el 68 barras, el lavado del PSS `+K`. Cumple

    A w = omega_B R_x.

Con `U` y `Z` bases ortonormales de `im R_x` y de su complemento, y
`A_q = Z^T A Z`:

    det(sI - A) = s det(sI - A_q),
    A_q (Z^T w) = 0      (exacto: A Z Z^T w = A w = omega_B R_x, y Z^T R_x = 0),
    det(sI - A_q) = s det(sI - A_qq).

`A_qq` es `A_q` deflactada por el autovector exacto `Z^T w`
(`certification/symmetry.py: deflate_eigenvector`). El origen tiene un
bloque de Jordan de tamaño 2. T1 retira solo la rotación. El segundo cero es
**físico**: la frecuencia deriva sin fuerza restauradora porque no hay
gobernador ni amortiguamiento. No es un error numérico ni una gauge, y se
conserva en el ledger.

Verificación en los tres benchmarks (`results/BC/BC01/*_summary.json`):

| benchmark | casos | `max abs(A R_x)` | `max abs(A w - omega_B R_x)` | compañero verificado |
|---|---|---|---|---|
| Kundur | 47 336 | 2.4e-10 | 3.9e-12 | 100 % |
| IEEE-39 | 220 288 | 2.5e-10 | 3.9e-12 | 100 % de los casos con `D = 0` |
| IEEE-68 | 1 936 | 5.6e-16 | 3.7e-15 | 100 % |

En el DAE (antes de eliminar `z`) la identidad se cumple a ~1e-12
(`tests/test_bc_adapter.py`). La rotación finita no lineal (NL00) conserva los
residuos a ~3e-13.

**Consecuencia.** Ningún portfolio sin gobernador de estos modelos es
asintóticamente estable en sentido estricto. Su estado físico es
`BOUNDARY_OR_UNRESOLVED` por el modo neutro de frecuencia. `STABLE` y
`UNSTABLE` se reportan para `A_qq`, es decir **módulo el modo neutro
declarado**. Así debe leerse toda «estabilidad» de F, G y BC.

### 3.2 Estados muertos

Con `blend = 0` hay estados congelados (filas nulas de `A`), por ejemplo en los
condensadores de `D = 2`. Producen ceros exactos. Se retiran exactamente
(estructura triangular por bloques) y se anotan en el ledger. No se eliminan
por magnitud.

### 3.3 Clasificador de cuatro estados

(`certification/classify.py`, congelado en `BC00_config.yaml`.)

1. Balanceo diagonal (`matrix_balance`), aplicando la misma similitud a la
   matriz de error.
2. `eps = abs(A(h) - A(2h)) + error hacia atrás`.
3. `d_axis = min_w sigma_min(A - i w I)`. Si `d_axis > 10 eps`, el conteo del
   semiplano derecho es exacto: ningún autovalor puede cruzar el eje con una
   perturbación de norma `eps`.
4. Si no, un autovalor con `Re(lambda) > 10 kappa(lambda) eps` cuenta como
   positivo resuelto (`UNSTABLE`; el conteo es una cota inferior). En otro
   caso el estado es `BOUNDARY_OR_UNRESOLVED`.
5. Ningún autovalor se elimina por magnitud. Las ventanas 1e-2 … 1e-8 son
   solo diagnósticas.

`d_axis` es un criterio **numérico** (malla más refinamiento acotado), no un
certificado por intervalos.

Pruebas obligatorias (`tests/test_bc_quotient.py`, 10/10):

- el ejemplo de Jordan V1;
- un polo físico `+1e-7` sigue `UNSTABLE`;
- un cero físico exacto es `BOUNDARY`;
- la reubicación con varios `beta` y el cambio de referencia angular no
  cambian el espectro físico;
- un modo interno no observable en el puerto sigue en el ledger completo.

## 4. BC00-C — Normalización

- Las comparaciones de normas se hacen en coordenadas balanceadas, con la
  matriz de error transformada por la misma similitud. Una norma euclídea en
  coordenadas crudas no es invariante. En el 68 barras daba `d_axis`
  artificialmente pequeño antes del balanceo.
- Con un cambio de base de potencia, la métrica se transforma por congruencia.
  NL00 verifica que el cambio `Sbase` de 100 a 1000 MVA deja los residuos
  invariantes a ~1e-15.
- La distancia espectral a −1 **no** es un radio de robustez no normal.
  `theory_checks` lo reproduce: `Q_L` tiene distancia 1 para todo `L`,
  mientras que `sigma_min(I + Q_L)` es aproximadamente `1/L`. El «margen de
  cierre» de F10 se lee como descriptor local, no como radio.

## 5. BC01 — H físico sin el disco y los 194 cruces de Kundur

### 5.1 H y kappa regenerados

(`experiments/binary_certification_v1/BC01_physical_H.py`)

| benchmark | puntos | casos de subconjunto | puntos exactos | puntos con subconjunto no resuelto | H igual a G1 |
|---|---|---|---|---|---|
| Kundur (K12A + K12B, todos los nodos) | 5 917 | 47 336 | 5 916 | 1 | 5 917 / 5 917 |
| IEEE-68 (122 nodos auditados) | 121 | 1 936 | 121 | 0 | 121 / 121 |
| IEEE-39 (F8, lengua G2, F7 `g <= 0.005` y 600 aleatorios) | 13 735 | 220 288 | 13 307 | 428 | 13 513 / 13 725 |

- En IEEE-39 los 10 puntos F8 y de la lengua no tienen etiqueta G1
  comparable. Sus `H` coinciden con F8 y G2: P1 `kappa = 1`, P2 2, P3 3,
  P4 4, `P_fold` 4, `P_inf` vacío. En la línea `k = 1.30`: 4 a `g = 0.020`,
  vacío a 0.042, 4 a 0.102 y vacío a 0.173.
- Los **212** cambios de etiqueta respecto de G1 están todos en puntos no
  exactos del camino rápido con `g <= 0.005` (etiqueta `G_LE_0.005`). En
  los 212 la diferencia es exactamente un subconjunto
  `BOUNDARY_OR_UNRESOLVED`.
- **El disco no interviene en IEEE-39.** En los 433 subconjuntos no resueltos
  ningún autovalor de `A_qq` tiene `abs(lambda) < 1e-2`. El modo cercano al
  eje es un **par oscilatorio** de 2.4 a 4.4 rad/s, con `abs(Re)` de 5e-6 a
  2e-2. El error del ensamblaje afín del camino rápido (umbral
  `10 eps`, mediana 1.8e-4) supera `d_axis`: G1 decidía el signo de `Re`
  sin criterio de resolución. Son puntos cerca de fronteras oscilatorias,
  no cruces aperiódicos.
- **Ningún cambio** ocurre en puntos exactos. Esas 212 etiquetas de G1 pasan a
  `UNKNOWN`. No son regresiones: son la corrección del criterio. Resolverlas
  exige evaluar esos puntos por el camino directo (no ejecutado en la Fase I).
- Kundur tiene 3 nodos donde un subconjunto inestable tiene abscisa menor que
  1e-3. G1 ya lo contaba, porque el subconjunto no era mínimo, así que `H` no
  cambia.

### 5.2 Los 194 cruces aperiódicos de Kundur

(`BC01_kundur_crossings.py`; `results/BC/BC01/BC01_kundur_crossings.csv`)

| resultado | cuenta |
|---|---|
| cruces (K12A 122, K12B 72), todos `APERIODIC_THROUGH_ORIGIN` | 194 |
| lado bajo: `A_qq` `STABLE`; lado alto: `UNSTABLE` (1 autovalor positivo) | 194 / 194 |
| compañero neutro verificado en ambos lados | 194 / 194 |
| `g0_physical - g*` (G1) | de −2.0e-5 a −3.6e-6 |
| descriptor completo frente a `A` (distancia máxima de espectros finitos) | 2.4e-7 |
| puerto **original** evaluable en `s = 0` (criterio preregistrado) | 24 / 194 |
| puerto **reubicado** evaluable en `s = 0`, con cambio de signo a través del cruce | **194 / 194** |
| resultado independiente de `beta` | 194 / 194 |
| el bloque de dispositivos cambia de signo | 0 / 194 |

**Mecanismo.** El modo que cruza es el integrador con fuga del regulador Q/V.
A `g = 0` está en `−0.05` (la fuga). Cruza el origen en `g` del orden de
3e-4, y en el lado alto su valor real está entre 0.0115 y 0.0134. El disco de
G1 desplazaba la frontera entre 4e-6 y 2e-5 en `g`, siempre hacia arriba.

**Puerto original.** `det T(s)` tiene un cero doble estructural en `s = 0`
(la rotación y el modo neutro). El ruido de las diferencias finitas lo parte.
El criterio preregistrado (igual winding de `det T_S` y `det T_0` en
`abs(s) = 1e-5`, estabilidad de la media y cambio de signo) se cumple en
24/194. Además:

- el winding coincide solo en 134;
- el signo cambia solo en 77;
- el signo medio indica «estable» en 138, aunque el lado alto es inestable en
  los 194.

El cociente de determinantes original **no sirve** en `s = 0`.

**Puerto reubicado** (`certification/port_origin.py`). Dentro del bloque
`f_x` de dispositivos se aplican dos reubicaciones de rango uno, una sobre la
rotación `R_x` y otra sobre el compañero `w`, con `beta` y `beta2`. Luego se
toma `T##(0) = g_z - g_x f_x##^-1 f_z`. Resultados:

- la identidad de determinantes se cumple con residuo de hasta 8.7e-13;
- el signo de `det T##(0)` cambia a través del cruce en 194/194;
- el signo de `det A` (reubicada) también cambia en 194/194, y el del bloque
  de dispositivos en ninguno.

Es decir, el cruce lo ve la red a través del puerto. Esta representación se
añadió **después** del prerregistro, al fallar el puerto original. Se declara
como desviación: el criterio original se reporta tal como estaba congelado.

**Conteo crudo del descriptor.** En el lado bajo el conteo crudo del pencil da
1 o 2 autovalores «positivos» en 166/194 casos, cuando el conteo físico es 0:
son los ceros estructurales partidos. Sin el cociente, un conteo sin disco
sería **falso**; con el disco, ocultaría el cruce. El cociente resuelve ambos
problemas.

## 6. BC02 — Representación de la familia binaria

(`BC02_binary_representation.py`; `results/BC/BC02/BC02_vertices.csv` y
`BC02_summary.json`; 13 familias, 245 vértices. Cada celda de la tabla da el máximo por familia; los
rangos van de la familia más baja a la más alta.) Las tres afirmaciones se
separan así:

- **A. Afinidad en los parámetros de control con `S` fijo** (F11/G3). Es de
  `theta = (g, k, …)`, no de `delta`. No se reabre aquí: se verificó en F7 y
  G3.
- **B. Localidad de puertos con equilibrio común.**
  `T_S = T_0 + sum_i (T_i - T_0)`, con incrementos soportados en las barras
  del candidato.
- **C. Realización común afín en `delta`.** Cada candidato lleva ambos
  dispositivos con pesos `2(1 - delta_i)` y `2 delta_i` sobre el caso
  `rho = 0.5`. El ausente queda en su equilibrio intensivo, que es el mismo en
  todos los vértices con despacho emparejado. Hay dos variantes:
  - C1 mantiene la dinámica del ausente (fantasma), alimentada por la tensión
    de barra y sin inyectar;
  - C2 la reemplaza por el relleno `-(x - x0)`.

| comprobación | IEEE-39 (6 puntos F8 × 16 vértices) | Kundur (3 nodos × 8) | IEEE-68, 4 candidatos (3 × 16) | IEEE-68, 12 plantas (77 vértices) |
|---|---|---|---|---|
| C: jacobiano del descriptor afín en `delta` (error relativo) | C1 0; C2 1.9e-19 | 1.2e-19 y 2.2e-19 | 0 y 1.7e-23 | 0 y 6.7e-23 |
| C: equilibrio común, `f = g = 0` en todos los vértices | 1.5e-12 | 9.9e-13 | 8.9e-13 | 1.0e-12 |
| C: equivalencia de vértices (espectro = portfolio real + fantasma o relleno) | 4.4e-5 a 8.0e-5 | 5.2e-5 a 6.7e-5 | 1.4e-5 a 2.1e-5 | 1.9e-5 (25 vértices) |
| C1: `Re` máxima de los fantasmas | −0.05 (Hurwitz) | −0.05 (Hurwitz) | +3e-15 a +2.5e-14: **cero marginal** | +7e-14: **cero marginal** |
| `A` reducida: no afinidad (segunda diferencia relativa) | 0.36–0.37 | 4.1 | 5e-7 | 2.2e-6 |
| B: `T_S = T_0 + sum (T_i - T_0)` en tres valores de `s` | 2.6e-8 a 8.0e-8 | 2.8e-8 a 5.8e-8 | 2.4e-11 a 4.1e-11 | 1.3e-10 |
| (15): `det P = h det T` (error de log) | 5.7e-13 | 8.5e-14 | 1.9e-12 | 2.0e-12 |
| ledger: bloques de dispositivo con `Re > 1e-9` | 0 | 0 | 0 (hay ceros marginales, `Re` máx. 2e-13) | 0 |

Lectura:

1. **C existe y es exacta en los tres benchmarks.** El jacobiano del
   descriptor es afín en `delta` a redondeo, con un equilibrio común.
   - La equivalencia de vértices se cumple a 1e-5–8e-5. Es del orden de la
     raíz del error del jacobiano (1e-9), lo esperable si el máximo lo aporta
     el cero doble defectivo; no se aisló qué autovalor lo da.
   - En IEEE-68, el padding C1 **no** es estable: algún bloque fantasma tiene
     un cero marginal. Allí solo C2 cumple la regla «padding estable, nunca
     estados muertos marginales». En IEEE-39 y Kundur ambas variantes son
     Hurwitz.
2. **La `A` reducida no es afín en `delta`.** La no afinidad relativa es de
   0.36 en IEEE-39 y de 4.1 en Kundur. En IEEE-68 da 5e-7–2e-6, pero es una
   norma de Frobenius relativa dominada por las entradas rígidas del modelo
   (`abs(f)` hasta 2e5). **No** es evidencia de afinidad. T4 no se aplica a
   `A(delta)` tal como está enunciado.
3. **B se cumple** a la precisión del equilibrio.
   - Una comprobación directa en Kundur (barras 2 y 3) y en IEEE-39 P4
     (barras 30 y 33) da estos órdenes para el incremento `T_i - T_0` fuera de
     los dos canales del candidato, relativo al incremento dentro:

     | benchmark | fuera / dentro |
     |---|---|
     | Kundur | 8e-8 |
     | IEEE-39 P4 | 8e-7 a 1.5e-6 |

   - El resto está en la barra de la máquina slack y coincide con la
     tolerancia del flujo de carga: `abs(z_S - z_0)` entre 2e-8 y 1.5e-7.
   - El tamaño de soporte de `BC02_summary.json` (umbral 1e-8 del máximo) está
     dominado por esa tolerancia. El soporte físico es de 2 canales.
4. **Ledger (prerrequisito de BC03).** Ningún bloque de dispositivo de ningún
   portfolio real tiene un autovalor en el semiplano derecho: `N_h = 0`. En
   IEEE-68 hay bloques de dispositivo con un cero marginal, así que el
   contorno de (16) debe sortear `s = 0` también por `h_S`.

## 7. Resultados anteriores potencialmente afectados

Se detienen los claims afectados. Los experimentos independientes que siguen
siendo válidos se conservan.

| resultado anterior | efecto | acción |
|---|---|---|
| E21/E34: «PV MW» y 4270.7 MW del flagship | son MVA nominales; la potencia activa desplazada es 2096.6 MW | corregir la etiqueta; el espectro no cambia |
| E14/E18: comparaciones «MW-matched» | el emparejamiento se hizo en MVA | reetiquetar como «emparejado en MVA nominal» o rehacerlo en MW |
| N16/E39: ranker «MW reemplazados» | AUC de 0.77–0.86 a 0.57–0.73 con MW activos | detener el claim de ranker por MW |
| toda «estabilidad» F/G/BC | estable módulo el modo neutro; estrictamente `BOUNDARY` | añadir la calificación |
| etiquetas G1 en IEEE-39 con `g <= 0.005` en el camino rápido | 212 no certificables: par oscilatorio con `d_axis` menor que el error del ensamblaje (no es el disco) | `UNKNOWN`; los puntos exactos no cambian |
| fronteras aperiódicas de Kundur (G1) | desplazadas entre −2e-5 y −4e-6 en `g` | sin efecto a la precisión citada (1e-4) |
| margen de cierre y distancia a −1 (F10) | no es un radio de robustez | leer como descriptor local |
| claims pre-G1 con condensador sin amortiguar | el condensador sin amortiguar es inestable | no reutilizar |
| ANDES | solo estática con el PSS desactivado | no citar como validación dinámica |
| IEEE-68 | 4 candidatos: negativo; doce plantas: 3 puntos | sin generalización |
