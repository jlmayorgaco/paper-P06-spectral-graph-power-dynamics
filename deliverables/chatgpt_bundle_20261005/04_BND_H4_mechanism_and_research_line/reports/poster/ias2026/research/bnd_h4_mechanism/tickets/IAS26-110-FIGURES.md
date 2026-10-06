# IAS26-110 — integración visual final

**Dependencias:** IAS26-020, 030, 060 y 080 cerrados; decisiones de 040/090 registradas. **Pregunta:** ¿puede el lector entender fenómeno, mecanismo, alcance y acción en 15 segundos sin claims no respaldados?

**Regla de arquitectura:** partir de `docs/FIGURE_DESIGN_CONTRACT.md` y de F1/F2 ya producidas. F3 sólo si M1 y el mapa pasan; F4 sólo si el contrato GFM y sus 16 celdas pasan. Si la evidencia central de robustez/reparación debe ocupar un panel principal, redactar una **enmienda versionada** del contrato visual antes del render final. No cambiar silenciosamente F3/F4 por otros paneles.

**Entrega:** exactamente cuatro figuras principales y franja pequeña de paridad/TDS, más suplemento regenerable. Scripts de figuras consumen datos congelados; no ejecutan ciencia. Cada caption declara modelo, operador, denominador, incertidumbre y limitación relevante. Revisar tamaño de lectura, color accesible y SVG/PDF/PNG en su `RUN_ID`.

**Gate:** cada marca visual tiene fila fuente y estado de claim; ninguna curva, celda o tasa procede de un run fallido o de una fórmula presentada como medición.
