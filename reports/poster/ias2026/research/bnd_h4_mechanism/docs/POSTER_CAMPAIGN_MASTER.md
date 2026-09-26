# IAS 2026 — plan maestro de evidencia para Beyond Nodal Damping

**Estado:** plan de ejecución, versión 1.0, 2026-09-26. Este archivo organiza trabajo futuro; no registra un experimento nuevo ni altera los gates numéricos existentes.

**Contratos vigentes:** [MODEL_CONTRACT.md](MODEL_CONTRACT.md), [FIGURE_DESIGN_CONTRACT.md](FIGURE_DESIGN_CONTRACT.md), [M1A_RECONCILIATION_CONTRACT.md](M1A_RECONCILIATION_CONTRACT.md) y [H4_MODAL_SCHUR_CONTRACT.md](H4_MODAL_SCHUR_CONTRACT.md).

**Control operativo:** `../tickets/STATUS.json` y un solo ticket activo.

## 1. Objetivo y afirmación central

En el IEEE-39/TX4 registrado, demostrar mediante evidencia trazable que una combinación de reemplazos SG→GFL puede producir una inestabilidad oscilatoria aunque todos sus subconjuntos propios sean estables; que los factores físicos locales permanecen regulares en la frontera mientras el cierre colectivo se aproxima a singularidad; y que una acción de control registrada mejora el margen, con alcance medido en escenarios operativos y en dinámica fasorial no lineal.

La frase del póster debe poder decirse sin símbolos: **cuatro reemplazos seguros en cada combinación menor fallan juntos por una realimentación colectiva de la red; un ajuste de control puede corregir el caso estudiado**. La última cláusula queda condicionada al contraste fuera del caso canónico.

Cinco preguntas gobiernan el programa:

1. ¿Existe el bloqueador H4 y cómo se sitúa en la familia V9?
2. ¿El cierre matricial local/colectivo representa la frontera del DAE completo?
3. ¿Una intervención física elegida antes del holdout mejora el espectro completo?
4. ¿En qué proporción persisten el fenómeno y la mejora bajo envolventes declaradas?
5. ¿La respuesta fasorial no lineal muestra las tasas y frecuencias predichas en pequeña señal?

## 2. Punto de partida verificado

| Evidencia | Estado al crear este plan | Límite práctico |
|---|---|---|
| Baseline H4 Python G0–G11 | PASS, `results/20260925T204017_549c3c07_h4_crossmode_v3` | Modelo TX4/P4 registrado. |
| Paridad H4 Python–Julia GFL11 | PASS en la reconciliación del baseline | Reproducción de un mismo modelo; no validación física independiente de todos los portfolios. |
| F1 V4 | Asset `READY`, 16 portfolios y 1.216 eigenvalues transversales con eigenvectores en `results/20260926T083621_h4_f1_f2_assets_v1` | Auditar `alpha_perp`, `alpha_Omega`, identidad modal y procedencia de cada caso antes de cerrar IAS26-020. |
| F2 | Asset `READY`, 72 factores físicos locales y 17 puntos colectivos en el mismo run | El gráfico actual no prueba todavía la continuación de una raíz reducida general en `s`. |
| M1A | `M1A_COMPLETE`, `BLOCKED_M1_STRICT`; primer fallo `base+deltaY`, residual relativo `2.7385662723e-8`, `safe_to_run_m2=false` en `results/20260926T083325_h4_m1a_reconciliation_v3` | Reparar bloques de carga base/flagship antes de F3/M2. |
| TDS triple/H4/H4 retuneado | Evidencia fasorial histórica en `experiments/tx4_tds_final.py` y sus resultados | Reutilizar sólo tras auditar modelo, perturbación, trayectorias crudas y ajuste temporal. |
| V9 y robustez CDW | Campañas históricas existentes | Reutilización condicionada a igualdad de modelo, política, equilibrio, parámetros, banda y hashes. |

En el NPZ de F1, H4 tiene 86 estados DAE y 84 coordenadas transversales; estas dimensiones no deben confundirse. Una auditoría de lectura del NPZ encontró que `alpha_perp` y el máximo en la banda 0.3–1.5 Hz coinciden en los 16 portfolios V4 nominales. Esto no autoriza extender esa coincidencia a V9, incertidumbre o GFM.

## 3. Alcance y decisiones congeladas

- Modelo central: IEEE-39, H4 = `{30,33,35,37}`, P4 = `(g,k,t,h)=(0.03625,1.425,1.5,1)`, red y GFL fijados por `MODEL_CONTRACT.md`.
- Acción primaria registrada: retuning Q/V `g: 0.03625 → 0.25`. Es un candidato, no una reparación universal ni un óptimo.
- Métrica primaria de estabilidad: abscisa del **espectro transversal completo** `alpha_perp`. Guardar por separado `alpha_Omega` para la familia 0.3–1.5 Hz y `MAC_to_H4`, `eigenvalue_gap`, `mode_family_id` cuando haya continuación o sustitución tecnológica.
- Umbral de decisión `tau_dec`, envolventes, límites físicos, semillas, tolerancias de integración, disturbios y reglas de ajuste se fijan en la configuración versionada de cada ticket **antes** del run correspondiente. Este plan no inventa valores faltantes.
- Cada corrida científica usa un `RUN_ID` nuevo e inmutable. Los outputs históricos se preservan como legado; no se sobrescriben ni se renombran para fingir que pertenecen a una campaña nueva.
- No mezclar resultados de TX4/GFL11 con BND39, GFL10, una configuración de carga distinta, ni ParaEMT en un mismo claim.
- No afirmar «damping local positivo», interacción física irreducible de orden cuatro, robustez universal, causalidad física de una ablación `eta`, TDS como EMT, ni optimalidad de la reparación sin pruebas específicas.

## 4. Secuencia de tickets y dependencias

| Ticket | Pregunta / entregable | Dependencia | Estado inicial |
|---|---|---|---|
| IAS26-010 | Reparar discrepancia `base+deltaY`; decidir M1/F3 | contratos M1A/M1 | **ACTIVO** |
| IAS26-020 | Auditar V4/F1 exacto y separar espectro completo/banda | plan | evidencia existente, pendiente de cierre |
| IAS26-030 | Auditar y completar frontera `g`/F2 con identidad de modo | 020; M1 sólo para raíz reducida general | evidencia existente, pendiente de cierre |
| IAS26-050 | Congelar intervención `g=0.25` y comparador, si lo hay | 030 | pendiente |
| IAS26-060 | MC de punto operativo y contraste emparejado | 050 + config previa completa | pendiente |
| IAS26-080 | TDS fasorial no lineal canónica y holdout estratificado | 050 + 060 | pendiente |
| IAS26-040 | Ablación matemática `eta` / F3 | **M1 PASS** + 030 | **bloqueado** |
| IAS26-070 | MC de parámetros, separado de MC operativo | 060 + envolvente legítima registrada | ampliación |
| IAS26-090 | Lattice GFL/GFM / F4 | núcleo cerrado + contrato GFM | ampliación |
| IAS26-100 | Censo N−1 factible | núcleo cerrado + contrato de contingencia | ampliación |
| IAS26-110 | Integrar cuatro figuras y franja de evidencia | 020, 030, 060, 080; 040/090 sólo si pasan | pendiente |
| IAS26-120 | Auditoría adversarial de claims y freeze del póster | 110 | pendiente |

**Camino crítico:** `010 → 020 → 030 → 050 → 060 → 080 → 110 → 120`. Si M1 no pasa, IAS26-040 se cierra como `BLOCKED`/`REFUTED` según la causa; el camino crítico sigue. Al terminar IAS26-080 se revisa si 070, 090 o 100 aportan evidencia necesaria antes del freeze. Sólo un ticket se ejecuta a la vez.

**V9:** el censo de 512 portfolios y múltiples políticas es un suplemento contextual. IAS26-020 audita su compatibilidad con P4/TX4 antes de emplearlo. Un censo exhaustivo de una familia finita no lleva intervalo binomial.

## 5. Diseño de evidencia por bloque

### Fenómeno — IAS26-020

V4 contiene exactamente 16 portfolios. Conservar espectros y eigenvectores, `alpha_perp`, `alpha_Omega`, frecuencia, dimensiones DAE/transversal, residual del equilibrio y árbol de inclusión. Gate: 15 subconjuntos propios por debajo de `-tau_dec` y H4 por encima de `+tau_dec` para `alpha_perp`; cualquier punto dentro de la banda de decisión es indeterminado. Verificar el NPZ existente antes de resolver otra vez.

### Mecanismo y frontera — IAS26-010/030/040

Primero, el escalón `A_red → T_raw → T_port → T_base+deltaY → T_action → h_k` debe cerrar con los umbrales ya registrados. M1 exige `|h_k(lambda_c)| < 1e-8`, raíz independiente próxima al polo completo y complemento finito. El diagnóstico M1A fija `base+deltaY` y `action-space` a residual relativo `≤1e-12`; registrar también error absoluto y condicionamiento. Una sola corrección focalizada y un nuevo run de verificación; si no cierra, no se ejecuta `eta`.

La frontera física varía la ganancia Q/V `g`, no una ganancia PLL. IAS26-030 alinea `alpha_perp(g)`, la raíz de la familia registrada, `sigma_min(I+M_ii)` físico y `sigma_min(I+Q_H)` en el mismo equilibrio, frecuencia y orden de puertos. La derivada por puertos se compara con una continuación independiente del DAE sólo si la representación ya pasó su gate. F3, de ejecutarse, es una **ablación de operador** con base fijada antes de observar el mapa; su interpretación no es una maniobra de control físico.

### Intervención — IAS26-050/060

Comparación emparejada entre H4 original y `g=0.25` con la misma red, cargas, despacho, ubicaciones, ratings y reglas de factibilidad. La predicción de `Delta alpha` se congela antes de evaluar el DAE de holdout. Guardar acierto de signo, error de magnitud, `alpha_perp` final, frecuencia, identidad modal y casos de empeoramiento. Un comparador adicional requiere selección en desarrollo y holdout independiente; si no queda preespecificado, el contraste primario es solamente `g=0.25`.

### Robustez — IAS26-060/070

MC-OP: **1.000 IDs de escenario pre-generados**, semilla y regla de perturbación operativa pre-registradas; 16 portfolios V4 por ID y un H4 retuneado por ID. Son aproximadamente 16.000 evaluaciones del lattice + 1.000 de retuning. MC-PAR: otro ensamble independiente de 1.000 IDs, parámetros solamente y mismo esquema, después del núcleo. No cambiar tamaño por un intervalo que parezca estrecho.

Por escenario `s`: `x_s=max_{R⊊H4} alpha_perp,s(R)`, `y_s=alpha_perp,s(H4)`, `m_s=min(y_s,-x_s)`, `B_s=[x_s<-tau_dec]∧[y_s>tau_dec]`, `Delta alpha_s=alpha_repaired,s-alpha_original,s`. Distinguir persistencia del **mismo H4** de existencia de **otro bloqueador mínimo**; esto último sólo se afirma dentro de la familia efectivamente censada. Declarar inválidos, límites violados e indeterminados con sus IDs; no reemplazar draws fallidos después de verlos. Reportar cotas conservadoras de la proporción entre todos los IDs, y Wilson 95 % condicional entre los escenarios válidos; la unidad estadística es el escenario, no cada portfolio. Para márgenes, cuantiles y bootstrap por escenario completo. No presentar estas proporciones como frecuencia natural de eventos en redes reales.

### Dinámica temporal — IAS26-080

Propuesta de volumen: 3 configuraciones canónicas (triple estable, H4, H4 retuneado) × 3 amplitudes pequeñas × 2 ubicaciones de pulso = **18 trayectorias**. Selección automática pre-registrada de 18 escenarios MC × 3 configuraciones × 2 amplitudes = **108 trayectorias**. Total propuesto: **126**, más repeticiones de tolerancia. Si una clase de estrato no existe, registrar que está vacía y aplicar una regla de sustitución fijada antes de leer los resultados; no seleccionar trazas a mano.

Perturbación temporal con inicialización algebraica consistente y observables relativos del DAE fasorial. Estimar `alpha_TDS` y `f_TDS` de las trayectorias integradas mediante ventana y criterio de ajuste predefinidos; contrastar con `alpha_eig` y `f_eig`. Registrar el rango de amplitud de pequeña señal y discrepancias; no inferir recuperación de frecuencia o nadir con un modelo sin regulación primaria. Las TDS históricas pueden reducir trabajo sólo si contienen señales crudas y el mismo contrato.

## 6. Figuras y decisión de arquitectura

`FIGURE_DESIGN_CONTRACT.md` permanece congelado: F1 y F2 ya existen; F3 depende de M1, F4 depende de un modelo GFM defendible. Este plan **no renombra ni promociona** nuevas figuras principales. IAS26-110 decidirá, con evidencia completa, si el póster usa F3/F4 o si solicita una **enmienda versionada** para sustituir espacios por robustez y reparación/TDS. La enmienda debe escribirse antes de exportar la versión final; ningún resultado fallido ocupa un espacio mediante una figura provisional.

La historia científica que debe sobrevivir a cualquiera de esas decisiones es: **combinación problemática → cierre colectivo → corrección física → alcance estadístico y temporal**. La paridad Python–Julia es una franja de evidencia; sus números no se interpretan como validación de campo.

## 7. Reglas de ejecución y cierre

1. Leer `tickets/STATUS.json` y ejecutar **sólo** `active_ticket`. Al terminar, registrar gate, `RUN_ID`, evidencia, resultado negativo y próximo ticket autorizado; no iniciar automáticamente el siguiente.
2. Antes de resolver un caso, auditar si ya existe un artefacto con iguales hashes de modelo, equilibrio, parámetros, banda, solver y alcance. La coincidencia parcial no autoriza reutilización silenciosa.
3. Configuración, semillas, tolerancias, selection rules, denominadores y tratamiento de fallos se congelan antes de cada campaña; cambios posteriores requieren enmienda y nuevo `RUN_ID`.
4. Cada ejecución científica escribe `results/<RUN_ID>/` con `config/`, `inputs/`, `raw/`, `derived/`, `figures/`, `claims/`, `report/`, `environment/` y `logs/`. Conservar espectros, eigenvectores necesarios, tracking, estados de fallo y trayectorias. No dejar outputs ni zip en la raíz.
5. El código reutilizable y los entry points futuros deben separar lógica de modelo, ejecución y renderizado. Los scripts heredados permanecen auditables; cualquier refactor respeta los hashes y gates del contrato del modelo, o requiere enmienda explícita.
6. Clasificar conclusiones como `PROVED`, `NUMERICALLY_VALIDATED`, `EMPIRICAL_HOLDOUT`, `BENCHMARK_SPECIFIC`, `CONDITIONAL`, `REFUTED`, `OPEN` o `BLOCKED`. Un gate fallido detiene sólo las dependencias reales y conserva el negativo.
7. Auditoría final: minimalidad, espectro completo versus banda, regularidad local versus estabilidad local, misma física en comparación, selección de reparación previa al holdout, TDS medida de señales, paridad same-model, TX4 versus BND39 y ausencia de claim EMT.
8. Sólo `main` local; no push ni ramas remotas. El plan puede versionarse en un commit local. Los resultados numéricos quedan en sus carpetas de corrida.

**Condición de cierre del programa:** exactamente cuatro figuras principales y una franja, cada claim trazable a una tabla/manifiesto, todos los negativos visibles en el ledger, y un argumento que un experto puede resumir en 15 segundos sin ver el suplemento.
