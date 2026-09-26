# NL00 — Procedencia del modelo, unidades y origen

Script: `experiments/nonlinear_portfolio_v2/NL00_audit.py`. Resultados en
`results/NL/NL00/`. No se añadió física: se audita el modelo L0 tal como corre.

## 1. Modelo encontrado y hashes

Base: el snapshot `9a216145` (G1–G6), más la Fase I de BC en
`ias2026/binary-certification-v1`. El programa NL trabaja en
`ias2026/nonlinear-portfolio-v2`, que nace de esa rama. Los hashes SHA-256 de
las fuentes están en `results/NL/NL00/NL00_hashes.json`; el manifiesto completo
de `src/` y `configs/`, en `configs/nonlinear_portfolio_v2/NL_SNAPSHOT_MANIFEST.json`.

| fuente | sha256 (prefijo) | papel |
|---|---|---|
| `configs/ias2026/ieee39_network.json` | `705f5329…` | red, máquinas y controles IEEE-39 congelados |
| `configs/kundur/kundur_network.json` | `d63986cc…` | Kundur (F12) |
| `configs/ieee68/ieee68_network.json` | `0cee6c94…` | IEEE-68 (G3) |
| `data/raw/ieee39_full.xlsx` | `9c2048dc…` | fuente ANDES del IEEE-39 |
| `src/ibr_cycles/models/ieee39_devices.py` | `bcf0c76b…` | máquina de dos ejes, AVR/PSS, GFL |
| `src/ibr_cycles/models/ieee68_devices.py` | `3cb4259c…` | máquina subtransitoria del 68 barras, DC4B/ST1A |
| `src/ibr_cycles/models/ieee39_case.py` | `1f38a8e5…` | DAE, KCL y ensamblaje de portfolios |
| `src/ibr_cycles/models/ieee39_network.py` | `f42144cd…` | Ybus y flujo de carga |

La nota v2 existe solo en PDF (`TEORIA_NO_LINEAL_MATRICIAL_Y_CERTIFICACION.pdf`,
`713fe27b…`). Su paquete de código (`nonlinear_models.py`,
`nonlinear_checks.py`, 625 comprobaciones) **no se entregó**, así que las
comprobaciones v2 no pudieron reproducirse. Las diferencias ecuación por
ecuación están en `NONLINEAR_POWER_SYSTEM_MODEL.md` §9 y en
`EQUATION_CODE_TRACEABILITY.csv`.

## 2. NL00-A — Unidades y `gamma = Sn / Sbase`

- `NL00_units.csv` traza, por máquina:
  - `Sn`, el despacho P/Q, `pmax`, la utilización del convertidor y la
    tensión;
  - la potencia activa realmente retirada al reemplazarla.

  La tabla viene de BC00-A y se verificó de nuevo aquí.
- **Prueba de comportamiento del cambio de base** (`NL00_base_conversion.csv`).
  - Cada máquina se reexpresa en la base del sistema: reactancias y `ra`
    divididas por `gamma`; `M`, `D` y `Pm` multiplicadas por `gamma`; en el
    PSS de entrada de potencia, la ganancia dividida por `gamma`; peso 1.
  - En un punto fuera del equilibrio, la corriente de red y las derivadas
    (con los estados de potencia escalados) coinciden:

| caso | desajuste máximo de corriente | desajuste máximo de derivadas |
|---|---|---|
| IEEE-39 base | 1.8e-15 | 5.9e-16 |
| IEEE-39 30@0.5 | 2.7e-15 | 7.6e-16 |
| IEEE-39 flagship | 9.2e-16 | 5.8e-16 |
| Kundur 2 | 9.2e-16 | 1.3e-15 |
| IEEE-68 3+4+6+9 | 0 | 0 |

`P_sys = gamma P_local`, `I_sys = gamma I_local` y `X_sys = X_local / gamma`
se cumplen como identidades de comportamiento, no solo de etiqueta.

## 3. NL00-B — La cifra 4270.7

Es la suma de `Sn` (MVA) de las cuatro máquinas del flagship
(barras 30, 33, 35, 37: 1040 + 1174.8 + 1085.7 + 970.2). La potencia activa desplazada es 2096.6 MW,
la reactiva 420.6 Mvar y la nominal activa 2983 MW. Ver BC00-A.

## 4. NL00-C — Polos cercanos a cero y T1

`NL00_symmetry_audit.csv`, reutilizando las filas de BC01 (el mismo código):

| benchmark | casos | `max abs(A R)` | `max abs(A w - omega_B R)` | modo neutro verificado | autovalores de `A` en `abs(lambda) < 1e-3` (máx.) | autovalores de `A_qq` en `abs(lambda) < 1e-3` (máx.) |
|---|---|---|---|---|---|---|
| IEEE-39 | 220 288 | 2.5e-10 | 3.9e-12 | 100 % | 2 | 0 |
| Kundur | 47 336 | 2.4e-10 | 3.9e-12 | 100 % | 3 | 1 |
| IEEE-68 | 1 936 | 5.6e-16 | 3.7e-15 | 100 % | 2 | 0 |

- Los dos ceros de `A` son la rotación y el modo neutro de frecuencia
  (bloque de Jordan).
- El tercero en Kundur es el integrador del regulador Q/V cerca de su cruce
  aperiódico (BC01 §5.2). Es físico y se conserva.

**Rotación finita no lineal** (`NL00_rotation_nonlinear.csv`, `phi = 0.37` rad,
punto fuera del equilibrio). Todos los ángulos se desplazan `+phi` y todas las
tensiones se multiplican por `e^{j phi}`:

| caso | `max abs(Δf)` | error de covarianza de KCL | `max abs(f)` |
|---|---|---|---|
| IEEE-39 base | 5.4e-15 | 3.0e-13 | 0.49 |
| IEEE-39 30@0.5 | 7.8e-14 | 1.3e-13 | 2.6 |
| IEEE-39 flagship | 1.6e-13 | 1.4e-13 | 6.1 |
| Kundur 2 | 2.8e-13 | 3.1e-14 | 0.89 |
| IEEE-68 3+4+6+9 | 2.8e-13 | 2.4e-13 | 2.3e5 |

La simetría es exacta en el modelo no lineal, no solo en su tangente.

## 5. NL00-D — Slack frente a fuente infinita

`NL00_slack_audit.csv`:

- En el modelo dinámico, ninguna barra tiene tensión ni ángulo fijos. Todas
  las tensiones son incógnitas algebraicas (78, 20 y 136).
- La barra slack del flujo de carga lleva una máquina dinámica ordinaria:
  IEEE-39 barra 39, Kundur máquina 1, IEEE-68 máquina 16.
- Por eso existe la simetría rotacional, y por eso el modo neutro de
  frecuencia es marginal.
- Declarado: la barra 39 es una inercia grande pero finita, no una barra
  infinita.

## 6. NL00-E — Alcance de la evidencia ANDES

La reconciliación ANDES (F1) cubre solo inyecciones estáticas armonizadas,
con el PSS desactivado. **No valida** la dinámica del GFL, del PSS ni del AVR
fusionado. En IEEE-39 el AVR combina KA y TE en un primer orden
(`MODEL_CHANGED`, declarado). ANDES no puede citarse como validación
dinámica de L0.

## 7. Resumen de estado

| ítem | estado |
|---|---|
| unidades y `gamma` | `V` (identidad de comportamiento a ~1e-15) |
| 4270.7 | `V`: son MVA |
| simetría exacta (tangente y finita) | `V` |
| modo neutro marginal | `V`: físico, no numérico |
| slack | `V`: sin fuentes fijas |
| ANDES | alcance limitado: solo estática |
| fidelidad v2 (gobernador, límites, DC, ZIP) | `ABSENT`; es L2 y no se añade al congelado |
