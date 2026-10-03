# Revisión visual · Póster IAS 2026

Se inspeccionaron las quince vistas individuales a 100 dpi y el póster completo a 75 dpi. La revisión cubrió jerarquía, alineación, contraste, espacio interior, legibilidad y recortes. El póster conserva una página de 36 × 48 pulgadas.

| Pieza | Resultado de la revisión |
| --- | --- |
| Header blanco | Los dos logos mantienen sus márgenes y el centro queda libre. |
| Header verde | Título, subtítulo y autor conservan la jerarquía. |
| 01 Motivación | Contexto, reemplazo SG→IBR y pregunta se leen en ese orden. |
| 02 IEEE-39 | Los candidatos rojos destacan frente a la red secundaria; el pie deja de repetir los números ya visibles y los ajustes de control se nombran sin código interno. Se retiró la frecuencia de identidad modal no resuelta. |
| 03 Destino y ruta | Las garantías A/B/C, la fórmula y las barras de evidencia siguen una lectura clara. El rótulo de las barras explica la medida sin el código F10. |
| Franja de resultados | Cuatro cifras y sus etiquetas permanecen alineadas; la cuarta resume la dependencia del orden en el censo de política matched. |
| 04A Factorización | Las cuatro ecuaciones ocupan filas amplias; la nota inferior delimita el resultado a la factorización de orden completo. |
| 04B Frontera | La gráfica y los valores local/colectivo comparten una escala visual clara; la leyenda queda dentro del panel. |
| 05 Políticas | El atlas y sus cifras conservan jerarquía. Los tres recuentos están alineados y cada barrido indica parámetros variables y fijos; el rótulo distingue el mapa A de los recuentos de espectro completo. |
| 06 Envolvente | Los 967 casos válidos aparecen como una partición 412/190/365; el subtotal de 555 con bloqueador colectivo se explica sin sumar porcentajes superpuestos. |
| 07 Reajuste | Las gráficas permanecen visibles; los recuentos se identifican como sintéticos, los fallos temporales muestran denominador y se delimita la validación independiente. |
| 08 Rutas | Las rutas insegura y segura se distinguen por color, forma y texto; la ecuación comprueba la estabilidad de cada prefijo y el panel indica la frecuencia observada de órdenes dependientes. |
| 09A Comprobaciones | Dos cifras con sus denominadores y la nota de fallos del solver sustituyen el inventario repetido de seis métricas; “time-domain” reemplaza la sigla TDS. |
| 09B Alcance | Tres frases explican qué modela la ruta segura y qué fenómenos no se simularon. |
| Footer | Conclusión, contacto y referencias ocupan tres bloques. Los cinco trabajos incluyen datos bibliográficos suficientes para identificarlos sin QR. El pie conserva aire interior y una franja blanca más amplia bajo el póster. |

Los paneles 04A/04B y 09A/09B son archivos independientes con caja exterior fija. La nueva composición mantiene márgenes de 23 mm, separación horizontal de 8 mm y separación vertical de 4 mm. Los encabezados de todos los módulos se revisaron tras añadir tracking y corregir la alineación óptica con los círculos; el título más largo de la primera fila permanece en una línea. La compilación final con `-Section all -StrictLayout` no emitió avisos de composición; se inspeccionó el PNG de cada panel y el de la página completa después de corregir un recorte inicial en 04B.

## Auditoría tipográfica

Las dos variantes de TeX Gyre Heros cubren texto, títulos y cifras; Latin Modern Math cubre las fórmulas. Los tamaños de texto de los quince módulos y de las figuras TikZ generadas se definen ahora en `poster_design_system.tex`; no quedan declaraciones `\fontsize` en esos archivos. Se igualaron las leyendas y notas de métricas a 24/27 pt, los nodos SG/IBR a 30/32 pt y las cifras locales a roles consistentes. El texto ordinario más pequeño son las referencias de 20 pt y el mayor es el título de 90 pt; los glifos inferiores a 20 pt en el PDF corresponden a subíndices, exponentes y marcas matemáticas. El módulo 9A se revisó a 100 dpi tras dar aire a la nota inferior sin solaparla con los denominadores.
