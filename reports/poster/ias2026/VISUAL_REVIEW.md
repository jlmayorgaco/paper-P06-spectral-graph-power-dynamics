# Revisión visual por módulo · IAS 2026

Se revisaron los PNG individuales de las trece piezas a 100 dpi y el póster completo a 75 dpi. La evaluación cubrió alineación con la retícula, jerarquía, contraste, espacio interior, legibilidad y posibles cortes o superposiciones. Los valores científicos y las medidas exteriores no se modificaron.

| Pieza | Evaluación y decisión |
| --- | --- |
| Header blanco | Logos equilibrados en los extremos y aire central limpio. Se copiaron ambos a `assets/` para que el proyecto compile desde GitHub. |
| Header verde | Título dominante, subtítulo y autor con separación y contraste adecuados. Se conserva. |
| 01 Motivación y pregunta | Secuencia clara: contexto, transición SG→IBR, pregunta y definición. La caja amarilla y la ecuación conservan aire suficiente. |
| 02 Resultado IEEE‑39 | Se reforzó el contraste de líneas y buses secundarios sin competir con los cuatro candidatos rojos. Las tres columnas mantienen títulos y pies legibles. |
| 03 Destino y ruta | Las garantías A/B/C, la implicación y las barras de evidencia mantienen un recorrido visual ordenado. Se conserva. |
| Franja de resultados | Cuatro cifras con anchos y separadores consistentes; legibles sobre verde oscuro. Se conserva. |
| 04 Cierre de red | Fórmulas, mecanismo, gráfica y lectura contextual distinguen evidencia de interpretación. Los avisos quedan dentro de la caja. Se conserva. |
| 05 Atlas de políticas | Se acortó la nota de estado y se marcó «CONDITIONAL» con el color semántico dorado. El mapa, cifras y leyenda tienen jerarquía clara. |
| 06 Envolvente sintética | Se aclaró que las fracciones usan 967 IDs válidos y se reunió el aviso de revisión junto al cierre del módulo. Las tres tarjetas quedan equilibradas. |
| 07 Reajuste | Se trasladó la leyenda de la serie temporal a la zona superior sin datos para liberar el centro de la curva. Cifras, dos gráficas y estados quedan separados. |
| 08 Rutas seguras | Formulación, comparación de rutas y validaciones se leen en ese orden; rojo y verde tienen también textos y formas distintivas. Se conserva. |
| 09 Evidencia y alcance | Las seis métricas se reparten de forma regular; la frase de alcance cabe completa y cierra el argumento. Se conserva. |
| Footer | Conclusión, referencias, contacto y QR tienen zonas definidas. El QR está ahora dentro de `assets/` para una compilación portátil. |

Las piezas se recompilaron con `./compile_section.ps1 -Section all -StrictLayout` sin avisos de composición. Tras cada ajuste visual se revisó la pieza afectada y la página completa. Para futuros cambios, sigue el contrato de `SECTIONS.md` y los roles de `DESIGN_SYSTEM.md`.
