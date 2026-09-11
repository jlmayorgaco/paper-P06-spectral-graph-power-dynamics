# BC — Tabla de evidencia (estado al cierre de la Fase I)

Solo se ejecutaron BC00–BC02. Las filas de BC03–BC12 figuran como `UNKNOWN`.
Los claims en formato CSV están en `BC_CLAIMS_PHASE_I.csv`.

| id | claim | estado | método | evidencia | límite |
|---|---|---|---|---|---|
| C01 | `replaced_mw` son MVA nominales; el flagship desplaza 2096.6 MW activos | `VERIFIED_ON_MODEL` | traza de `Sn` y del despacho por máquina | `results/BC/BC00/BC00_units_*` | la etiqueta no altera el espectro |
| C02 | El ranker de E39 por MW activos tiene AUC de 0.57–0.73, no de 0.77–0.86 | `VERIFIED_ON_MODEL` | recálculo con MW activos | `BC00_E39_ranker_corrected.csv` | mismas etiquetas de estabilidad |
| C03 | `A R = 0` con `R` derivado de las ecuaciones, en los tres benchmarks | `VERIFIED_ON_MODEL` | residuo en 269 560 casos | `results/BC/BC01/*_summary.json` | diferencias finitas (2.5e-10) |
| C04 | Con `D = 0`, `A w = omega_B R`: el origen tiene un bloque de Jordan y un cero físico marginal | `PROVED_UNDER_ASSUMPTIONS` (identidad) + `VERIFIED_ON_MODEL` | derivación y residuo (3.9e-12) | BC00 §3.1 | desaparece si `D != 0` o hay gobernador |
| C05 | Ningún portfolio sin gobernador es asintóticamente estable en sentido estricto | `VERIFIED_ON_MODEL` | C04 | BC00 §3.1 | limitación del modelo, no de la red |
| C06 | El H físico (sin disco) es igual al de G1 en Kundur (5917/5917) y en IEEE-68 (121/121) | `VERIFIED_ON_MODEL` | clasificador de cuatro estados sobre `A_qq` | `BC01_*_points.csv` | 1 punto de Kundur con un subconjunto no resuelto |
| C07 | En IEEE-39 hay 0 cambios en puntos exactos; 212 etiquetas del camino rápido con `g <= 0.005` son no certificables | `VERIFIED_ON_MODEL` / `UNKNOWN` | `d_axis` frente a `10 eps` | `BC01_ieee39_points.csv` | par oscilatorio cerca del eje; no es el disco |
| C08 | Los 194 cruces aperiódicos de Kundur son el integrador con fuga Q/V, que cruza el origen | `VERIFIED_ON_MODEL` | bisección sobre `A_qq`, ambos lados clasificados | `BC01_kundur_crossings.csv` | frontera física entre 4e-6 y 2e-5 por debajo de G1 |
| C09 | El puerto original no es evaluable en `s = 0` (24/194) | `FALSIFIED` (para el puerto original) | criterio preregistrado | ídem | — |
| C10 | El puerto reubicado `T##(0)` es regular y cambia de signo en 194/194 cruces; el bloque de dispositivos no cambia | `VERIFIED_ON_MODEL` | dos reubicaciones exactas de rango uno | ídem; `certification/port_origin.py` | añadido tras el prerregistro; novedad `UNKNOWN` |
| C11 | Existe una realización común del descriptor afín en `delta`, con equilibrio común, en los tres benchmarks. La `A` reducida no es afín | `VERIFIED_ON_MODEL` | 245 vértices: afinidad a 1e-19, equilibrio a 1.5e-12, equivalencia de vértices a 8e-5 | `results/BC/BC02`; BC00 §6 | en IEEE-68 solo C2 (relleno) es un padding estable |
| C11b | Localidad de puertos `T_S = T_0 + sum (T_i - T_0)` y `det P = h det T` | `VERIFIED_ON_MODEL` | 3 valores de `s`; comprobación directa del soporte | ídem | localidad a la tolerancia del flujo de carga (1.5e-7) |
| C11c | Ningún bloque de dispositivo tiene polos en el semiplano derecho (`N_h = 0`) | `VERIFIED_ON_MODEL` | ledger de los 245 vértices | ídem | IEEE-68: ceros marginales en `h_S` |
| C12 | Certificados binarios y de cardinalidad (T3, T4, T5, (29)) | `UNKNOWN` | — | — | Fases II–III |
| C13 | Búsqueda exacta de testigos mínimos sin monotonía | `UNKNOWN` | — | — | BC05 |
| C14 | Diseño seguro ante cualquier orden frente a endpoint, secuencia estacionaria y transición | `UNKNOWN` en modelos; T6 `PROVED_UNDER_ASSUMPTIONS` | — | — | BC07 y BC11 |
