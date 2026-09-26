# IAS26-060 — Monte Carlo operativo

**Pregunta:** ¿con qué frecuencia, bajo una distribución de estrés operativa declarada, persisten H4 y la mejora del retuning?

**Pre-run obligatorio:** config con distribución de `P_L,Q_L,P_G`, correlaciones, límites, redispatch, 1.000 IDs de escenario y semilla, `tau_dec`, fallo/timeout, resolución modal, regla de no reemplazo y criterio de validez. No iniciar mientras un campo esté abierto.

**Por ID:** resolver los 16 portfolios V4 y H4 retuneado en el mismo escenario. Registrar espectros/identidad de modo o evidencias suficientes para `alpha_perp`, `alpha_Omega`, `x=max alpha_perp(proper)`, `y=alpha_perp(H4)`, `m=min(y,-x)`, evento H4, evento «otro bloqueador mínimo» dentro de V4, `Delta alpha`, estado numérico y factibilidad. Aproximadamente 17.000 evaluaciones; 1.000 unidades independientes.

**Análisis:** proporciones con Wilson 95 % entre válidos y cotas conservadoras sobre los 1.000 IDs; proporción inválida explícita. Contraste emparejado de reparación; mediana y q05/q25/q75/q95; bootstrap por escenario. Registrar deterioros, cambios de familia e indeterminados. No transformar una proporción del ensamble sintético en probabilidad natural de falla.

**Gate:** los 1.000 IDs tienen resultado o código de fallo auditable; ninguna exclusión post-hoc; denominadores y CIs reproducibles desde tablas crudas.
