# IAS26-090 — intervención tecnológica GFL/GFM

**Prioridad:** ampliación. **Dependencia:** núcleo evaluado y contrato GFM aprobado/versionado. **Pregunta:** ¿qué sustituciones de tecnología mueven H4 a estabilidad y conservan la identidad del modo?

**Antes del run:** congelar ecuaciones GFM, ratings, `P/Q`, droop, inercia virtual si existe, voltage loop, límites, despacho y equilibrio; validar mismo network/load y factibilidad. Sin contrato defendible, estado `BLOCKED_MODEL` y ninguna F4.

**Diseño:** las 16 asignaciones GFL/GFM en buses 30/33/35/37, matriz 4×4 en el orden del contrato visual. Guardar espectro completo, `alpha_perp`, raíz/frecuencia, `MAC_to_H4`, gap, familia, singularidades comparables y estados inválidos. Si se selecciona una configuración ganadora, exigir holdout nuevo; propuesta 500 escenarios × 16 asignaciones con comparaciones emparejadas.

**Gate:** 16/16 celdas resueltas o explícitamente inválidas bajo mismo contrato; ningún valor inferido para celda faltante. `r_GFM` sólo se reporta entre asignaciones factibles y evaluadas; no decir que GFM es universalmente mejor.
