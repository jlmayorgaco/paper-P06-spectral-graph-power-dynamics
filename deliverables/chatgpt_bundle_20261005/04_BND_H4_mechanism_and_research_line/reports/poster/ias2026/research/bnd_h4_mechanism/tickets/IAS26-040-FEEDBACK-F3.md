# IAS26-040 — ablación matemática de feedback / F3

**Estado inicial:** BLOQUEADO por M1. **Dependencias:** IAS26-010 con M1 PASS e IAS26-030 cerrado.

**Pregunta:** ¿la variación de los términos cruzados en un operador reducido fijado antes del barrido mueve la misma raíz hacia la región inestable?

**Diseño:** `h_k(s;eta)=h_local(s)+eta h_cross(s)` con `eta∈[0,1]`. Congelar base/partición, malla `g×eta`, solver, tolerancias y regla de tracking antes de ejecutar; propuesta inicial 41×41 sujeta a coste y resolución. Guardar raíz reducida, raíz DAE cuando la comparación tiene el mismo modelo físico, residual, condición del complemento, `MAC_to_H4`, gap, familia, derivada analítica y diferencia finita. Vigilar polos del subsistema abierto al aproximar `eta=0`.

**Gate:** M1 pasa en el punto base y en los puntos requeridos; raíces independientes y derivadas verifican las tolerancias pre-registradas; no hay cambio de familia sin marca explícita. Si el gate falla, conservar diagnóstico y eliminar F3 del póster mediante decisión documentada en IAS26-110.

**Interpretación:** ablación del operador, no implementación física. El retuning `g=0.25` y sus pruebas avanzan por IAS26-050/060 aunque esta ablación fracase.
