# BC — Aplicabilidad de los teoremas de la nota v1 a los modelos (Fase I)

Nota: `spectral_portfolio_theory_v1/TEORIA_Y_DEMOSTRACIONES.md`. Cada fila dice
qué hipótesis se verificaron en los modelos después de BC00–BC02 y cuáles
faltan. Una demostración sintética no convierte un teorema en resultado sobre
IEEE; el estado sobre el modelo se da por separado.

| resultado | estado matemático | hipótesis en los modelos | estado en los modelos |
|---|---|---|---|
| T1, cociente de una simetría declarada | `PROVED_UNDER_ASSUMPTIONS` (álgebra lineal clásica) | `A R = 0` con `R` derivado de las ecuaciones: verificado a 2.5e-10 en los tres benchmarks | `VERIFIED_ON_MODEL` |
| T1, «lo que no dice»: cero físico adicional | — | cadena de Jordan `A w = omega_B R_x` en todo caso con `D = 0`: el cociente deja un cero **físico** (modo neutro de frecuencia) | `VERIFIED_ON_MODEL`; la estabilidad estricta es `BOUNDARY` en todos los portfolios sin gobernador |
| Deflación del modo neutro (`A_qq`) | exacta: `A_q Z^T w = 0` | requiere `w` verificado; no existe si alguna máquina tiene `D != 0` (condensador de `D = 2`) | se usa solo con `w` verificado. Sin `w` se clasifica `A_q` directamente |
| T1a, reubicación `A - beta U U^T` | `PROVED_UNDER_ASSUMPTIONS` | independencia de `beta`: pruebas sintéticas y 194/194 cruces | `VERIFIED_ON_MODEL` |
| (15), `det P_S = h_S det T_S` | clásica (Schur) | sintético 1.8e-15; en los modelos, BC02 hasta 2.0e-12 (245 vértices) | `VERIFIED_ON_MODEL` |
| (16)–(17), contabilidad de polos del puerto | clásica | el término `N_h` es obligatorio. En los cruces de Kundur el puerto **original** en `s = 0` falla (24/194) por el cero doble estructural | `FALSIFIED` para el puerto original en `s = 0` (no para la identidad) |
| §7, puerto después del cociente | procedimiento | hecho con dos reubicaciones exactas dentro del bloque de dispositivos. `det T##(0)` es regular y cambia de signo en 194/194 | `VERIFIED_ON_MODEL` (Kundur). IEEE-39 e IEEE-68: `UNKNOWN` (no ejecutado) |
| (18)–(20), acciones localizadas `T_S = T_0 + sum E_i dY_i E_i^T` | álgebra exacta **si** el equilibrio es común | con despacho emparejado, el equilibrio de red coincide en todos los vértices a la tolerancia del flujo de carga (`abs(z_S - z_0) <= 1.5e-7`). La suma se cumple a 8e-8 o mejor; el soporte físico es de 2 canales por candidato | `VERIFIED_ON_MODEL` (en los 13 conjuntos de BC02) |
| T2, no-go de orden fijo | `PROVED_UNDER_ASSUMPTIONS` (contraejemplo elemental) | no requiere el modelo | se usa solo para impedir leer `kappa` como grado de interacción |
| T3, pequeña ganancia sobre todas las subfamilias | `PROVED_UNDER_ASSUMPTIONS` | (1) baseline y acciones individuales internamente estables tras el cociente: el baseline es estable módulo el modo neutro. Ledger de BC02: ningún bloque de dispositivo tiene `Re > 1e-9`, pero en IEEE-68 hay ceros marginales en `h_S`. (2) `Q` estable: no verificado. (3) cotas rigurosas `H_inf`: no disponibles | `UNKNOWN` (BC03) |
| T3a, cota inferior de `kappa` | `PROVED_UNDER_ASSUMPTIONS` | igual que T3 | `UNKNOWN` |
| T4, mayorante LMI de un subcubo | `PROVED_UNDER_ASSUMPTIONS` | requiere `A(delta) = A_F + sum delta_i A_i`. La `A` reducida **no** es afín en `delta`: su no afinidad relativa es 0.36 en IEEE-39 y 4.1 en Kundur (BC02). Solo el jacobiano del descriptor es afín, a 1e-19 | `NOT_APPLICABLE` tal como está enunciado. La versión para el descriptor (Lyapunov generalizado) es `UNKNOWN` |
| T4a, transferencia al modelo reequilibrado | `PROVED_UNDER_ASSUMPTIONS` | requiere una cota válida de `abs(E)`. `operator_error_bound` devuelve `UNKNOWN` | `NOT_APPLICABLE` hasta tener la cota |
| T5, completitud de Lyapunov booleana | `PROVED_UNDER_ASSUMPTIONS` | dimensión común: sí, en C1/C2. En IEEE-68 C1 tiene fantasmas marginales, así que solo sirve C2. Vértices Hurwitz: **no** con el modo neutro. Debe aplicarse a `A_qq` en la realización común; hay que comprobar que el relleno de C2 preserva la simetría (`R` nulo en los estados de relleno) | `UNKNOWN` (BC04) |
| (29), homotopía de celda | `PROVED_UNDER_ASSUMPTIONS` | pencil no singular en el contorno: el pencil original es singular en `s = 0` por la simetría. Se requiere el cociente | `UNKNOWN` (BC05) |
| §11.2, búsqueda exacta sin monotonía | algoritmo | la poda de superconjuntos solo sirve para minimalidad; no etiqueta estabilidad | `UNKNOWN` (BC05) |
| T6, hipergrafo y seguridad ante cualquier orden | `PROVED_UNDER_ASSUMPTIONS` | (a) todos los portfolios factibles: sí, 0 `INFEASIBLE` en los conjuntos auditados; (b) baseline estable: módulo el modo neutro; (c) `H` completo: enumeración exhaustiva de 16 u 8 subconjuntos en puntos exactos | `VERIFIED_ON_MODEL` en los puntos exactos, con «estable» leído como «estable módulo el modo neutro». En los 428 puntos IEEE-39 con subconjuntos no resueltos, `H` es `UNKNOWN` |
| (32), plan máximo seguro | exacto si `H` es completo y `p_i` está en MW | `p_i` debe ser la potencia activa (BC00-A), no `replaced_mw` | `UNKNOWN` (BC07) |
| (33), tres nociones de planificación | `PROVED_UNDER_ASSUMPTIONS` (ejemplo ejecutado) | — | BC07 |
| T7, ley local del retorno | `PROVED_UNDER_ASSUMPTIONS` | cruce simple y transversal; `f_s != 0`. En los cruces aperiódicos `omega = 0` y el retorno original es singular | `UNKNOWN` (BC08) |
| T8, entorno no lineal local | `PROVED_UNDER_ASSUMPTIONS` | `L` debe venir de cotas válidas del Hessiano. NL02 solo da curvatura muestreada. `A_q` no es Hurwitz (modo neutro): hay que usar coordenadas relativas que también retiren la frecuencia media | `UNKNOWN` (NL/BC10) |
| §15.1, transición híbrida | condición suficiente (42) | requiere el mapa de reset: `reset_map` es la identidad en el active set fijo (NL01); hace falta el cambio de equilibrio | `UNKNOWN` (BC11) |

## Modificación declarada de la noción de seguridad

La definición (8) de la nota exige `alpha_q < 0`. Con `D = 0` y sin gobernador,
`alpha_q = 0` **exactamente** en todos los portfolios: `F_safe` sería vacío y T6
trivial. Por eso en F, G y BC «estable» significa

    alpha(A_qq) < 0,   con el cero de A_q de autovector Z^T w verificado y declarado,

es decir, estabilidad **relativa**: los ángulos y la frecuencia se miden
relativos al centro de inercia. No es estabilidad asintótica del sistema
completo: la frecuencia media es marginal porque el modelo no tiene
gobernadores. Añadir un gobernador sería una nueva versión del modelo (L2), no
una corrección del congelado.
