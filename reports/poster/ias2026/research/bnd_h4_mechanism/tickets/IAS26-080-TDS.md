# IAS26-080 — dinámica fasorial no lineal

**Dependencias:** IAS26-050/060. **Pregunta:** ¿las tasas y frecuencias estimadas de señales integradas coinciden localmente con la predicción espectral?

**Reutilización:** auditar `experiments/tx4_tds_final.py`, sus tres trazas y el protocolo G2 antes de repetir. El nuevo contrato debe fijar pulso temporal, dos ubicaciones, tres amplitudes nominales, observable relativo, inicialización consistente, ventana/filtro/ajuste y tolerancias de integración. No introducir una perturbación permanente para afirmar restauración de frecuencia.

**Diseño propuesto:** 18 trayectorias nominales (3 casos × 3 amplitudes × 2 ubicaciones) y 108 trayectorias holdout (18 IDs MC seleccionados por estratos pre-registrados × 3 casos × 2 amplitudes), más repeticiones numéricas. En escenarios con cambio de bloqueador, la tercera configuración será el subconjunto relevante elegido por una regla determinista fijada antes del run. La selección de estratos debe definir qué pasa si una clase está vacía antes de examinar los datos. Conservar toda trayectoria y metadatos.

**Gate:** `alpha_TDS` y `f_TDS` salen de la señal integrada con criterio previo; comparar con `alpha_eig` y `f_eig`, reportando diferencias, fallos de ajuste y rango de amplitudes. Un ajuste fallido no se sustituye por la envolvente teórica. Reportar explícitamente «DAE fasorial no lineal»; no EMT.
