# BC — qué es conocido y qué podría ser nuevo (versión de la Fase I)

Estado al cierre de la Fase I (BC00–BC02 y pruebas sintéticas). Las Fases
II–IV (BC03–BC12) **no se han ejecutado**. Toda afirmación de novedad que
dependa de ellas figura como `UNKNOWN`.

## Conocido. No se reclama prioridad

| herramienta u objeto | por qué es conocido | uso en este programa |
|---|---|---|
| Complemento de Schur, `det P = h det T` | álgebra lineal clásica | contabilidad de polos del puerto (BC02: identidad verificada) |
| Nyquist generalizado, punto −1, diferencia de retorno | clásico (N23, F10) | localizador de fronteras oscilatorias |
| Cociente por una simetría (T1) y reubicación (T1a) | reducción de simetría estándar; teorema de Brauer para la reubicación de rango uno | conteo sin disco alrededor de cero (BC01) |
| Radio de estabilidad complejo, `min_w sigma_min(A - iwI)` | Hinrichsen–Pritchard, clásico | criterio numérico de resolución del signo (BC00-B) |
| Hipergrafos, clutters, elementos minimales | combinatoria estándar | objeto `H` (F2C) |
| Pequeña ganancia, Perron–Frobenius, Lyapunov, LMI, SOS | [L2, L3, L7, L8] | certificados de BC03/BC04 (no ejecutados) |
| Amortiguamiento no monótono con soporte de tensión | literatura de interacción de convertidores | observado (N22), no reclamado |
| Radio de estabilidad real disperso | [L1] Katewa–Pasqualetti | distinto de activaciones binarias; no se reclama |
| Certificados plug-and-play, subsistemas inestables | [L3, L4] Hallinan–Lestas | antecedente directo de BC03; no se reclama |

## Hechos nuevos del modelo encontrados en la Fase I

Son hechos verificados sobre estos modelos (`VERIFIED_ON_MODEL`), no
contribuciones metodológicas.

1. **Estructura exacta del origen en los tres benchmarks.**
   - La simetría rotacional `R` y su compañero de Jordan `w` (modo neutro de
     frecuencia) se derivan de las ecuaciones.
   - Las identidades se cumplen a precisión de diferencias finitas.
   - Todo portfolio con `D = 0` y sin gobernador tiene un modo físico
     marginal. En sentido estricto, **ningún portfolio de estos modelos es
     asintóticamente estable**. `STABLE` significa estable salvo ese modo
     declarado.
2. **Unidades.** `replaced_mw` es la potencia nominal del convertidor en
   MVA.
   - El flagship desplaza 2096.6 MW activos, no 4270.7 MW.
   - El ranker «MW reemplazados» de E39, recalculado con MW reales, tiene un
     AUC de 0.57–0.73, no de 0.77–0.86.
3. **Los cruces aperiódicos de Kundur son el integrador del regulador de
   tensión.** Nacen del polo de fuga −0.05 a `g = 0` y cruzan cero a
   `g ≈ 3e-4`.
   - El disco de G1 desplazó la frontera hacia arriba, entre 4e-6 y 2e-5 en
     `g`.
   - El puerto original no puede evaluarlos en `s = 0` (24/194). El puerto
     reubicado sí puede: su signo cambia a través del cruce en 194/194 (ver
     BC01).
   - Los 212 cambios de etiqueta de IEEE-39 **no** se deben al disco. Son
     pares oscilatorios no resueltos a la precisión del camino rápido.
4. **Existe una realización común del descriptor, afín en los
   indicadores.**
   - Se verificó en los tres benchmarks: 245 vértices, afinidad a 1e-19 y
     equilibrio común a 1.5e-12.
   - La matriz reducida `A(delta)` **no** es afín: no afinidad de 0.36 en
     IEEE-39 y de 4.1 en Kundur.
   - En IEEE-68 solo el relleno C2 es un padding estable.

## Contribuciones candidatas: estado

| pieza | qué haría falta | estado |
|---|---|---|
| A. Conteo físico completo con la simetría y los cruces por cero | T1 más el ledger de ceros (hecho); puerto regular en `s = 0` tras el cociente (hecho en Kundur) | `VERIFIED_ON_MODEL` para el conteo. El carácter de novedad del puerto reubicado es `UNKNOWN`: la reubicación de Brauer es clásica, y lo específico es aplicarla dentro del bloque de dispositivos para certificar cruces aperiódicos por puertos. |
| B. Certificados binarios y de cardinalidad, con poda válida | BC03–BC05 | `UNKNOWN` (no ejecutado) |
| C. Diseño seguro ante cualquier orden parcial | BC07 | `UNKNOWN`. T6 está `PROVED_UNDER_ASSUMPTIONS` (sintético). |
| D. Puente no lineal local | BC10–BC11 | `UNKNOWN` |

## Qué no se afirmará

- Complejidad polinómica de la búsqueda de `H`: la salida puede ser
  exponencial (§11.4 de la nota).
- Novedad de Schur, Nyquist, el punto −1, SOS, Lyapunov o hipergrafos.
- Un certificado muestreal como garantía continua: el criterio `d_axis` es
  numérico, no un certificado por intervalos.
- Validación ANDES del GFL/PSS completo. La reconciliación ANDES solo cubre
  inyecciones estáticas con el PSS desactivado.
