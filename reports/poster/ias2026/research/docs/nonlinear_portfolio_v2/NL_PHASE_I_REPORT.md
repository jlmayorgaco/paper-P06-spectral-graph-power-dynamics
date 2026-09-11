# NL — Informe de la Fase I (NL00–NL02): primer stop

Rama `ias2026/nonlinear-portfolio-v2`, creada desde la cabeza de la Fase I de
BC. Versión del modelo: `L1-v2.0`. Etiquetas de claims: **D** (definición),
**T** (teorema o identidad demostrada), **V** (verificado numéricamente en el
modelo), **E** (empírico).

Por instrucción, el programa se detiene aquí. No se ejecutó nada de
NL03–NL12.

## 1. Modelo encontrado y hashes

- **L0** es el código que corre hoy:
  - `ieee39_devices.py`: máquina de dos ejes, AVR/PSS y GFL de 10 estados más
    `x_v` (y `y_vi` opcional);
  - `ieee68_devices.py`: máquina subtransitoria, DC4B/ST1A;
  - `ieee39_case.py`: DAE, KCL y portfolios;
  - `ieee39_network.py`: Ybus y flujo de carga;
  - las redes congeladas de IEEE-39, Kundur e IEEE-68.
- Los hashes están en `results/NL/NL00/NL00_hashes.json`, en
  `NL00_MODEL_PROVENANCE.md` §1 y en el manifiesto
  `configs/nonlinear_portfolio_v2/NL_SNAPSHOT_MANIFEST.json`. Ese manifiesto
  contiene los blobs git de `src/` y `configs/` en el commit base, antes de
  los archivos NL.
- **L1** es la reescritura matricial:
  - `nonlinear/network.py`: `Cf`, `Ct`, `Yff`/`Yft`/`Ytf`/`Ytt` con taps
    complejos, flujos y pérdidas por rama, realificación apilada o
    intercalada con permutación explícita;
  - `nonlinear/model.py`: `PhasorModel` con el contrato del v2.

  Se probó igual a L0 (`tests/test_nl01_model.py`, 12/12):

  | comparación | error |
  |---|---|
  | Ybus de L1 frente a L0 | `< 1e-12` relativo |
  | KCL fuera del equilibrio | `< 1e-12` |
  | pérdidas por rama | `< 1e-10` |
  | espectro físico tras reordenar las barras | `< 1e-6` |

- Manifiesto de la corrida (código, salidas en sha256, entorno, semillas):
  `results/NL/NL_PHASE_I_MANIFEST.json`. Reejecutar NL00 y NL02 tras
  formatear el código reprodujo todas las salidas byte a byte, salvo el campo
  `elapsed_s`.
- Ecuaciones completas: `NONLINEAR_POWER_SYSTEM_MODEL.md`. Trazabilidad
  ecuación → código → prueba: `EQUATION_CODE_TRACEABILITY.csv` (33 filas).

## 2. Discrepancias con el v2

Detalle en `NONLINEAR_POWER_SYSTEM_MODEL.md` §9. L0 se mantiene sin cambios.

| punto | L0 | estado |
|---|---|---|
| Excitación de IEEE-39 | primer orden con KA y TE fusionados (armonizado) | `MODEL_CHANGED`, declarado |
| Mecánica | forma de potencia (dos ejes); forma de par (68 barras) | ambas válidas a velocidad nominal; no intercambiables fuera de ella |
| Gobernador (N.27) | ausente en los tres benchmarks | L2 |
| Límites y anti-windup del GFL (N.36–37), retardo de modulación | ausentes | L2 |
| Enlace DC (N.38–40) | ausente (`P_ref` = despacho, DC ilimitado) | L2 |
| Cargas ZIP (N.15) | PQ en IEEE-39 y Kundur; Z en IEEE-68 | sin cambio |
| Masa dependiente del estado | `M = I` | declarado |
| Paquete de código v2 (625 comprobaciones) | **no entregado** con el PDF | no reproducible; solo se reprodujeron las 404 de v1 |

No se ajustó ningún parámetro de IEEE-39 para reproducir autovalores.

## 3. Unidades y simetría

- **Unidades (V).**
  - Cada dispositivo está en su propia base, con `gamma = Sn/100` solo en la
    inyección.
  - El cambio a base del sistema es una identidad de comportamiento: corriente
    y derivadas coinciden a 3e-15 (NL00-A).
  - `replaced_mw` = MVA nominales. El flagship desplaza 2096.6 MW activos, no
    4270.7 (BC00-A).
- **Simetría (T + V).** La rotación exacta es `R_x = 1` en `delta` y
  `theta_pll`, con `R_z = jV`:
  - se cumple en la tangente (`A R = 0`, 2.5e-10);
  - se cumple en el modelo no lineal con rotación finita (3e-13).
- **Modo neutro (V).** Con `D = 0` y sin gobernador, el compañero de Jordan
  `w` cumple `A w = omega_B R_x`. El cociente deja un cero **físico**, el modo
  neutro de frecuencia. En sentido estricto ningún portfolio es
  asintóticamente estable; «estable» significa estable módulo ese modo (BC00
  §3.1).
- **Slack (V).** No hay fuentes fijas. La slack del flujo de carga es una
  máquina dinámica: IEEE-39 barra 39, Kundur máquina 1, IEEE-68 máquina 16.

## 4. Derivadas disponibles y método (NL02)

`nonlinear/manifold.py`. `G_z` se factoriza una vez con LU y nunca se invierte.

| objeto | fórmula | método | verificación en 11 puntos (errores relativos máximos) |
|---|---|---|---|
| `J_psi = -G_z^-1 G_x` | carta algebraica | LU sobre jacobianos por diferencias centrales | 2.8e-8 frente a `psi` resuelto por Newton |
| `A = f_x + f_z J_psi` | campo reducido | igual | 1.4e-8 |
| `H_psi[a,b]`, `H_r[a,b]` | N.47–N.50 | producto Hessiano-vector: diferencia central de 4 puntos de los residuos **reales** (sin paso complejo: el modelo tiene conjugados y módulos) | ver abajo |
| `B`, `C`, `D` | `f_u - f_z G_z^-1 G_u`, etc. | diferencias centrales en `u` | 1.0e-7, 2.8e-8, 8.9e-8 |

El Hessiano reducido se comparó con la segunda diferencia directa de
`r(x) = f(x, psi(x))`, con `psi` resuelto en cada evaluación, para 9 pasos
de 1e-1 a 1e-5 (`NL02_convergence.csv`). Mejor acuerdo independiente
(excluido el paso de la fórmula):

| benchmark | acuerdo |
|---|---|
| IEEE-39 (6 puntos: base, flagship, fuera del equilibrio, lengua, condensador `D = 2`) | 1.3e-6 a 5.7e-6 |
| Kundur (3 puntos) | 1.9e-7 a 5.0e-7 |
| IEEE-68 (2 puntos) | 1.1e-3 y 1.7e-4. **No certificado**: `abs(f)` llega a 2.3e5 (escalado rígido) y la sensibilidad de la fórmula al paso alcanza 0.10 con `h = 1e-5` |

Las 10 fallas de carta ocurren solo con pasos de 0.1 y 0.03: el punto sale de
la carta. Se reportan como límite, no se ocultan.

Estado de la curvatura: **V (SAMPLED_CURVATURE)**. Son derivadas en puntos,
no cotas sobre un dominio. No hay diferenciación automática disponible (ni jax
ni casadi); las cotas de NL03/NL04 requieren derivación analítica o
aritmética de intervalos.

## 5. Qué NL puede correr sin cambiar la física

| bloque | L0 suficiente | condición o reserva |
|---|---|---|
| NL03, existencia algebraica (N4, Krawczyk) | sí | intervalos: `mpmath.iv` en `.venv/bc-cert` o redondeo dirigido propio. La cota de `I_PQ` (N1) aplica a IEEE-39 y Kundur; IEEE-68 tiene carga Z (lineal) |
| NL04, radio no lineal | sí, con reserva | `A_q` **no** es Hurwitz: el modo neutro es una familia de equilibrios relativos. El Lyapunov debe hacerse en coordenadas que retiren también la frecuencia media (`A_qq`); la conclusión es estabilidad del conjunto, no del punto. El flagship inestable es `NOT_APPLICABLE` |
| NL05, admitancia de segundo orden | sí | usa la maquinaria de NL02 |
| NL06, rango de validez con TDS | sí | una perturbación sostenida de generación necesita gobernador (ausente): usar pulsos |
| NL07, cancelación N.34 | la comprobación simbólica y numérica en L0, sí | la parte con filtro, retardo, error de desacoplamiento y límites es **L2** |
| NL08, PV sin batería | no | **L2**: arquitectura DC completa |
| NL09, límites e híbrido | no | **L2** |
| NL10, grafo y energía | N8 sobre una reducción sin pérdidas y de tensión fija (modelo derivado distinto); N9 sobre L0 | N9 puede terminar `UNKNOWN` |
| NL11-A, retuning continuo | sí | — |
| NL11-B, cambio discreto | parcial | el mapa de reset al añadir o retirar un dispositivo es semántica nueva: la persistencia y la inicialización deben definirse y versionarse |
| NL12 + BC03–BC12 | sí | realización común C1/C2 (BC02); T4 no aplica a la `A` reducida |
| Tercer orden y Hopf | opcional | no antes de NL00–NL06 |

## 6. Estimación de costo

Referencias medidas en esta máquina, con 12 workers:

| corrida | tiempo |
|---|---|
| BC01 IEEE-39, 220 288 clasificaciones | 2.0 h |
| BC01 IEEE-68, 1 936 | 44 min |
| BC01 Kundur, 47 336 | 6.7 min |
| NL02, 11 puntos | 6–28 s |

Estimación para la Fase II (NL03–NL06 en 4–6 casos de IEEE-39 más los casos
aperiódicos de Kundur):

| bloque | desarrollo | cómputo |
|---|---|---|
| NL03 | 1–2 días (intervalos y cotas del estator y del GFL) | minutos por caso |
| NL04 | 1 día; depende de NL03 | Lyapunov trivial (n de unos 150). Las TDS de falsación, cientos de corridas, ~1 h |
| NL05 | 1–2 días | unas 200 TDS, ~1 h |
| NL06 | 0.5 día | unas 200 TDS, ~1 h |

La Fase III (L2: límites, DC y transición) tiene un costo dominado por el
desarrollo de modelos nuevos (días) y reejecuta solo los casos congelados. La
Fase IV (IEEE-68, familia de 12 plantas) cuesta del orden de 4096 × 3
clasificaciones directas, estimadas en 1–2 h, más SDP (cvxpy con
Clarabel/SCS, presentes en el venv congelado).

## 7. Primeras pruebas (hechas)

- `tests/test_nl01_model.py` (12/12):
  - Ybus L1 = L0 en las tres redes;
  - pérdidas por rama con taps;
  - ramas sintéticas con desfasadores (40 casos, `Yft != Ytf`);
  - una sola permutación entre las dos realificaciones;
  - equilibrio y potencias de los dispositivos;
  - KCL de L1 = residuo de L0 fuera del equilibrio;
  - invariancia del espectro físico al reordenar las barras;
  - las entradas solo alteran las filas declaradas.
- `tests/test_bc_quotient.py` (10/10) y `tests/test_bc_adapter.py` (5/5):
  cociente, modo neutro y adaptadores.
- NL00 y NL02 como experimentos reproducibles, con semillas y pasos
  declarados.

## 8. Qué no se afirma

- Ninguna cota certificada: NL02 es curvatura muestreada.
- Ninguna validación EMT: L3 no existe, y un modelo RLC no se llamará EMT.
- Ninguna fidelidad v2 añadida al congelado.
