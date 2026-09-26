# IAS26-030 — frontera y F2

**Pregunta:** ¿coincide el cruce de `alpha_perp(g)` con la pérdida de regularidad colectiva, manteniéndose regulares los factores físicos locales?

**Entrada a reutilizar:** baseline `results/20260925T204017_549c3c07_h4_crossmode_v3`, `TX4_PHYSICAL_LOCAL_FACTORS.csv`, `TX4_CONTEXTUAL_RETURN_SWEEP.csv`, `TX4_CONTEXTUAL_RETURN_AT_BOUNDARY.csv`, y F2 en `results/20260926T083621_h4_f1_f2_assets_v1`. Ya hay 17 puntos de `g` y `g*=0.20768140519037842`.

**Trabajo:** alinear `alpha_perp`, frecuencia, `sigma_min(I+M_ii)` físico y `sigma_min(I+Q_H)` por mismo equilibrio, frecuencia y orden de puertos. Auditar continuidad/identidad del modo, añadir el panel de `alpha(g)` si los artefactos existentes bastan. Si falta precisión de frontera, congelar una malla suplementaria y root solve antes de ejecutarlos. Contrastar una raíz reducida general en `s` sólo tras M1 PASS.

**Gate:** boundary completo y colectivo coinciden dentro del umbral vigente del baseline; factores locales físicos permanecen por encima del mínimo pre-registrado en P4 y frontera; identidades Schur/puerto pasan; toda comparación de polos reducidos queda condicionada a M1.

**Salida:** F2 y tabla de procedencia. Claim permitido: «los factores locales permanecen regulares mientras el cierre colectivo se aproxima a singularidad». No afirmar damping local positivo.
