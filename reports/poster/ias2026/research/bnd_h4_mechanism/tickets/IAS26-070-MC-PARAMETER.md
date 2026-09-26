# IAS26-070 — Monte Carlo paramétrico

**Prioridad:** ampliación después de IAS26-060. **Pregunta:** ¿cambia el bloqueador o la eficacia del tratamiento cuando varían parámetros registrados, manteniendo fijo el punto operativo?

**Pre-run:** auditar envolventes CDW para igualdad de modelo y significado. Si no son compatibles, declarar una nueva envolvente con límites físicamente motivados antes de producir datos. Separar incertidumbre de controlador de incertidumbre de planta/red cuando el presupuesto lo permita; no mezclarla con MC-OP. Propuesta: 1.000 IDs por ensamble efectivamente ejecutado, semilla nueva.

**Medidas:** mismos 16 portfolios, H4 retuneado, `alpha_perp`, `alpha_Omega`, identidad modal, persistencia de H4 y de cualquier bloqueador mínimo dentro de V4, resultado emparejado y fallos. Reportar cada familia de incertidumbre y denominador por separado; no sumar ensambles para inflar N.

**Gate:** contrato de parámetro y conjunto de IDs completos; estadística reconstruible y negativos visibles. Puede permanecer `OPEN` sin bloquear IAS26-080 ni el núcleo del póster.
