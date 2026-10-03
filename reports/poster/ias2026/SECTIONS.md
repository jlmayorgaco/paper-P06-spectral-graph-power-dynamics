# Póster IAS 2026: edición modular

Cada archivo de `sections/` contiene una pieza editable y puede compilarse por separado. `main.tex` solo arma la página. `poster_layout.tex` fija las medidas exteriores y las vistas individuales. `poster_design_system.tex` define colores, fuentes, tipos de texto y componentes; [DESIGN_SYSTEM.md](DESIGN_SYSTEM.md) explica cuándo usar cada uno. `poster_common.tex` carga los paquetes y esos archivos compartidos. Los datos y afirmaciones numéricas viven en `generated/` y se regeneran desde las fuentes congeladas descritas en [EVIDENCE.md](EVIDENCE.md). Los logos usados por el póster están en `assets/`; la evaluación visual de cada pieza está en [VISUAL_REVIEW.md](VISUAL_REVIEW.md) y la revisión de contenido en [POSTER_REVIEW.md](POSTER_REVIEW.md).

## Contrato de medidas

Las medidas siguientes están definidas en `poster_layout.tex` y se usan tanto en la página completa como en los PDF individuales. **Para editar un módulo, conserva su ancho y alto exteriores.** Puedes reorganizar texto, imágenes, diagramas, tipografía y espacios *dentro* de esa caja. Si necesitas otra medida exterior, hay que rediseñar la retícula completa.

| Pieza | Ancho × alto exterior |
| --- | ---: |
| Header blanco | 898,4 × 54 mm |
| Header verde | 868,4 × 86 mm |
| 01 Motivación y pregunta | 224 × 252 mm |
| 02 Resultado IEEE-39 | 412,4 × 252 mm |
| 03 Destino y ruta | 216 × 252 mm |
| Franja de resultados | 868,4 × 44 mm |
| 04A Factorización | 430,2 × 181 mm |
| 04B Frontera colectiva | 430,2 × 181 mm |
| 05 Atlas de políticas | 430,2 × 178 mm |
| 06 Envolvente sintética | 430,2 × 178 mm |
| 07 Reajuste | 430,2 × 210 mm |
| 08 Rutas seguras | 430,2 × 210 mm |
| 09A Comprobaciones | 430,2 × 84 mm |
| 09B Alcance del modelo | 430,2 × 84 mm |
| Footer | 868,4 × 80 mm |

La retícula tiene 868,4 mm de ancho y márgenes laterales de 23 mm, alineados con los bordes exteriores de los logos. La primera fila suma 224 + 8 + 412,4 + 8 + 216 = 868,4 mm. Las filas de dos columnas suman 430,2 + 8 + 430,2 = 868,4 mm. El espacio vertical entre filas es 4 mm.

## Archivos

| Archivo | Contenido |
| --- | --- |
| `sections/header_white.tex` | Franja blanca y logos |
| `sections/header_green.tex` | Franja verde, título y autor |
| `sections/01_motivation_question.tex` | Motivación y pregunta |
| `sections/02_flagship_ieee39.tex` | Resultado IEEE-39 |
| `sections/03_target_vs_path.tex` | Destino estable y ruta segura |
| `sections/result_strip.tex` | Franja de cuatro resultados |
| `sections/04_factorization.tex` | Derivación en cuatro pasos |
| `sections/04_boundary_scan.tex` | Gráfica y valores en la frontera |
| `sections/05_policy_atlas.tex` | Atlas de políticas |
| `sections/06_synthetic_envelope.tex` | Envolvente sintética |
| `sections/07_retuning.tex` | Reajuste de control |
| `sections/08_safe_paths.tex` | Rutas seguras |
| `sections/09_model_checks.tex` | Reproducción cruzada y TDS |
| `sections/09_model_scope.tex` | Alcance de la afirmación de ruta segura |
| `sections/footer.tex` | Conclusión, referencias y contacto |

## Cómo editar un solo módulo

1. Abre únicamente el archivo de `sections/` correspondiente. Edita el contenido dentro de `\PosterSectionContent`; conserva el `minipage` exterior, la altura de `TicketTwoPanel` o `PosterPanel`, y el bloque de compilación individual.
2. Desde esta carpeta, ejecuta `./compile_section.ps1 -Section 05_policy_atlas -StrictLayout` (sustituye el nombre). El PDF y PNG de esa sección aparecerán en `build/sections/`.
3. Revisa visualmente ese PNG: legibilidad, cortes, superposiciones y aire interior. Corrige dentro del módulo.
4. Ejecuta `./compile_section.ps1 -Section full -StrictLayout` y revisa `../../../output/pdf/IAS2026_Poster_Pulido.png` para verificar la integración. `-Section all -StrictLayout` recompila las quince piezas y el póster.

Un aviso de composición detiene la compilación con `-StrictLayout`. La compilación limpia no reemplaza la revisión visual: TikZ puede dibujar fuera de una caja sin emitir un aviso.

## Encargo listo para otro chat

> Edita solo `reports/poster/ias2026/sections/05_policy_atlas.tex`. Lee `SECTIONS.md` y `DESIGN_SYSTEM.md`. Conserva la caja exterior de 430,2 × 178 mm, los márgenes, el sistema de colores y tipografía, y la posición de los otros módulos. Mejora el contenido y el espacio interior sin alterar `main.tex` ni la retícula. Usa los roles de `poster_design_system.tex`. Compila esa sección con `-StrictLayout`, revisa su PNG, compila el póster completo y revisa la integración visual. Informa qué cambió y dónde quedaron los PDF/PNG.

Sustituye el nombre y las medidas por la fila correspondiente de la tabla. Si la solicitud afecta varios módulos o toda la composición, edita la retícula explícitamente y revisa todos los módulos afectados.
