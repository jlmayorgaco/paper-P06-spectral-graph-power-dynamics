# Sistema de diseño · Póster IAS 2026

Este documento describe el estilo **ya aplicado** al póster. Las definiciones ejecutables están en `poster_design_system.tex`; las medidas exteriores de la página y los módulos están en `poster_layout.tex`. `poster_common.tex` carga paquetes, el sistema visual y los datos compartidos. Un cambio a un token visual debe revisarse en el póster completo porque puede afectar a varias secciones.

## Principio visual

Póster científico de 36 × 48 pulgadas, legible por capas: título dominante, bandas de sección verdes, resultados destacados y texto técnico sobre blanco. El verde organiza; el dorado marca navegación y cifras clave; el azul identifica evidencia o estados estables; el rojo identifica inestabilidad o fallo. Los diagramas repiten esas distinciones con texto y formas, además del color.

## Colores

| Token TeX | Hex | Uso |
| --- | --- | --- |
| `IASWhite` | `#FFFFFF` | Fondo de página y de módulos |
| `IASDarkText` | `#17372D` | Párrafos, ecuaciones y etiquetas |
| `IASHeaderGreen` | `#004A3A` | Banda principal del título |
| `IASMidGreen` | `#00664B` | Encabezados de módulos |
| `IASDarkGreen` | `#003B2D` | Títulos internos y números oscuros |
| `IASDeepGreen` | `#00291F` | Franja de resultados y pie |
| `IASGold` | `#C99A20` | Índices, énfasis y cifras clave |
| `IASBlue` | `#176494` | Evidencia y estados estables |
| `IASVermillion` | `#BC3030` | Inestabilidad, bloqueo o deterioro |
| `IASGray` | `#5F6B66` | Estado pendiente y texto secundario |
| `IASPaleGreen` | `#EAF2ED` | Fondos de apoyo y texto claro sobre verde |
| `IASPaleGold` | `#FFF4D6` | Recuadros de pregunta o advertencia |
| `IASRule` | `#CFDDD5` | Bordes y separadores suaves |

`IASIvory` es un alias heredado de `IASWhite`; no introduce otro color. Usa los tokens por su función. En gráficos, acompaña azul/rojo/dorado con círculos, cruces, rombos o etiquetas; no dependas solo del color.

## Fuentes y jerarquía

Fuente de texto: **TeX Gyre Heros**. Títulos y cifras destacados: **TeX Gyre Heros Cn** (`\PosterCondensed`). Matemáticas: **Latin Modern Math**. Los tamaños son del PDF final; la segunda cifra es el interlineado.

| Rol | Comando | Tamaño / interlineado | Uso |
| --- | --- | ---: | --- |
| Título principal | `\PosterHeroTitle` | 90 / 94 pt | Una línea dominante en el header verde |
| Sobrelinea del título | `\PosterHeaderEyebrow` | 25 / 28 pt | Contexto breve en dorado |
| Subtítulo del header | `\PosterHeaderSubtitle` | 36 / 40 pt | Explicación científica breve |
| Autor y afiliación | `\PosterHeaderByline` | 30 / 34 pt | Línea final del header |
| Título de módulo | `\PosterPanelHeadingType` | 43 / 46 pt | Bandas verdes generales |
| Título de primera fila | `\PosterTicketHeadingType` | 44 / 47 pt | Tres módulos superiores |
| Número de módulo | `\PosterPanelBadgeType` | 32 / 32 pt | Círculo dorado |
| Párrafo principal | `\PosterParagraph` | 34,2 / 39 pt | Texto normal, alineado a la izquierda y con final irregular |
| Subtítulo interior | `\Subhead` | 38 / 41 pt | Entrada a figura o bloque técnico |
| Ecuación | `\EquationSize` | 36 / 41 pt | Fórmulas corrientes |
| Ecuación protagonista | `\HeroEquation` | 44 / 48 pt | Fórmula que debe dominar un bloque |
| Cifra interior | `\MetricSize` | 38 / 42 pt | Resultados dentro de módulos |
| Tabla o nota compacta | `\TableSize` | 26 / 29 pt | Datos densos con lectura cercana |
| Pie de figura | `\FigureSize` | 24,7 / 27,2 pt | Leyendas breves y etiquetas |
| Franja: cifra | `\PosterMetricValueType` | 40 / 42 pt | Valor dorado |
| Franja: etiqueta | `\PosterMetricLabelType` | 27 / 30 pt | Descripción blanca |
| Franja: detalle | `\PosterMetricDetailType` | 24,2 / 26,2 pt | Nota verde pálido |
| Pie: rótulo inicial | `\PosterFooterKickerType` | 24 / 27 pt | «TAKEAWAY» en dorado |
| Pie: titular | `\PosterFooterTitleType` | 56 / 59 pt | Conclusión principal |
| Pie: explicación | `\PosterFooterSummaryType` | 29 / 32 pt | Instrucción de lectura |
| Pie: referencias | `\PosterFooterReferenceType` | 19 / 22 pt | Citas densas |
| Pie: rótulo de referencias | `\PosterFooterReferenceHeadingType` | 23 / 26 pt | Etiqueta dorada |
| Pie: rótulo de contacto | `\PosterFooterContactHeadingType` | 25 / 28 pt | Etiqueta dorada |
| Pie: nombre de contacto | `\PosterFooterContactNameType` | 23 / 26 pt | Nombre en blanco |
| Pie: datos de contacto | `\PosterFooterContactLineType` | 21 / 24 pt | Correo, afiliación y ciudad |
| Pie: etiqueta del QR | `\PosterFooterQrLabelType` | 19 / 21 pt | Texto bajo el código |

Hay tamaños locales en diagramas y tablas que resuelven casos concretos. Para texto nuevo, empieza por un rol de esta tabla y ajusta localmente solo si la figura lo exige. No comprimas ni escales fuentes para forzar que quepan.

## Retícula y aire

`poster_layout.tex` fija página de 914,4 × 1219,2 mm, margen útil de 23 mm, ancho útil de 868,4 mm, separación de 8 mm entre columnas y 4 mm entre filas. El header blanco ocupa 898,4 × 54 mm; la banda verde ocupa 868,4 × 86 mm. Las cajas exteriores de cada módulo están enumeradas en `SECTIONS.md`.

La banda blanca contiene solo los logos: IAS a la izquierda, Universidad de los Andes a la derecha. Sus anchos son **137 mm** y **120 mm**, respectivamente (`\PosterIASLogoWidth`, `\PosterUniandesLogoWidth`). Conserva la proporción original de cada archivo y el aire blanco entre ambos; el título comienza en la banda verde inferior.

Dentro de los módulos, `PosterPanel` usa una banda de título de **32 mm** y `TicketTwoPanel` una de **36 mm**. Ambos usan borde de **1,1 pt**, esquinas de **3 mm**, padding horizontal de **3 mm**, superior de **4 mm** e inferior de **2 mm**. El número va en un círculo de **14 mm**. La franja de métricas usa celdas de **40 mm** de alto. Estas medidas están definidas como tokens en `poster_design_system.tex`.

Los espacios locales de 1–6 mm son ajustes ópticos de figuras y fórmulas; mantenlos dentro de cada módulo. Si el contenido no cabe, reordena o edita el interior antes de tocar la caja exterior. Conserva suficiente blanco entre bloques para distinguir pregunta, método, evidencia y conclusión.

## Componentes y estados

| Componente | Función |
| --- | --- |
| `PosterPanel` | Módulo general con banda verde y número dorado |
| `TicketTwoPanel` | Módulo de primera fila con encabezado más alto |
| `MetricCell` | Cifra, etiqueta y detalle de la franja de resultados |
| `\StatusPass`, `\StatusConditional`, `\StatusPending` | Estados de evidencia azul, dorado y gris |
| `\StableMark`, `\UnstableMark`, `\BoundaryMark`, `\UnresolvedMark` | Marcas con forma y color para diagramas |

Uso recomendado en una sección:

```tex
{\PosterParagraph Texto principal del módulo.\par}
{\Subhead\color{IASDarkGreen}Título interior\par}
{\FigureSize\color{IASDarkText}Leyenda breve de la figura.\par}
```

Para una edición independiente: lee `SECTIONS.md`, modifica un solo archivo de `sections/`, compila esa sección con `-StrictLayout`, inspecciona su PNG y después revisa la composición completa. Cambia `poster_design_system.tex` cuando quieras cambiar una regla visual compartida; cambia `poster_layout.tex` cuando quieras mover o redimensionar cajas exteriores.
