# Evidencia congelada del póster IAS 2026

Esta rama incluye los archivos fuente necesarios para **regenerar los valores y las figuras numéricas del póster** a partir de salidas ya calculadas. `EVIDENCE_MANIFEST.json` registra tamaño y SHA-256 de cada fuente; el hash normaliza CRLF a LF para funcionar en distintos sistemas. `audit_poster_evidence.py` comprueba esos hashes y vuelve a contar las afirmaciones visibles. Este paquete no ejecuta de nuevo los flujos de potencia, las DAE ni los problemas de autovalores.

## Reproducción

Desde `reports/poster/ias2026`:

```powershell
python generate_poster_assets.py
python audit_poster_evidence.py
.\compile_section.ps1 -Section all -StrictLayout
```

El generador usa Python 3, NetworkX y Pillow; `poster_requirements.txt` fija las versiones con las que se verificó esta entrega. Los logos aprobados están versionados en `assets/` y no son transformados por el generador. Los resultados editables quedan en `generated/`; la página final queda en `output/pdf/IAS2026_Poster_Pulido.pdf` desde la raíz del repositorio.

## Relación entre afirmaciones y fuentes

| Resultado del póster | Fuente congelada de esta rama | Alcance |
| --- | --- | --- |
| Caso nominal IEEE-39, 15/15, $H_4$, $+0.127\,\mathrm{s}^{-1}$, MW/MVA | `research/bnd_h4_mechanism/results/.../tables/POSTER_NUMBERS_FINAL.csv`, `research/results/20260911_CLAIM_MATRIX.csv`, `research/bnd_h4_mechanism/configs/IAS26-050_INTERVENTION_V1.json` | Caso P4 nominal, no regla universal. |
| Frontera y factores local/colectivo | `research/bnd_h4_mechanism/results/20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1/tables/` | Barrido auditado; equivalencia de modo reducido no demostrada. |
| Atlas de políticas A y recuentos A/B/C | `research/results/FINAL_CLOSURE/figures/FIG1_policy_map_source.csv`, `research/results/20260911_CLAIM_MATRIX.csv` | A: $g\times k$, $t=1.5$, $h=1$; B: $g\times t$, $k=h=1$; C: $g\times h$, $k=t=1$. Los recuentos 30/36/16 corresponden a hipergrafos de espectro completo; la imagen representa A. |
| 52 puntos de comparación | `research/results/20260911_BASELINE_COMPARISON.csv` | B3/B4/B5: 22/52, 31/52 y 38/52. |
| Conjunto sintético de 1000 ID y retuning | `research/bnd_h4_mechanism/results/20260927T144629Z_d0fecb32_ias26_060_operating_v1/derived/SCENARIO_METRICS.csv`, tablas `MC_EVENT_RATES.csv`, `IAS26-FINAL_P4_MC_DATA.csv` y trazas TDS del mismo run | 967 casos físicamente válidos; frecuencias descriptivas del conjunto congelado. |
| Dependencia del orden | `research/results/FINAL_CLOSURE/FC10_census_lattice.csv` y `FC10_summary.json` | Censo IEEE-39 de nueve candidatos con política reactiva *matched*, distinto del caso nominal P4. Hay 327 destinos estables, 3 dependientes del orden y ninguno sin ruta estable en ese censo. |
| 100/100 y 56/60; 30/90 fallos | `IAS26-FINAL_EVIDENCE_STRIP.csv`, `IAS2026_POSTER_READINESS.json` y `evidence/claim_audit.md` | Acuerdo entre códigos sobre el mismo modelo; la comprobación temporal tiene cuatro desacuerdos de frecuencia cero y 30 fallos de solver. |

En el conjunto sintético, los 967 casos válidos se dividen en **412 sin bloqueador en $V_4$, 190 con $H_4$ mínimo y 365 con triples mínimos**. Los dos últimos grupos suman los 555 casos con bloqueador colectivo. Bajo la intervención, 594 mejoran y 373 empeoran; los 186 rescates de $H_4$ son un subconjunto de las mejoras. La auditoría comprueba esas relaciones fila por fila.

## Límites de la evidencia

La auditoría de esta rama verifica procedencia, hashes, aritmética, denominadores y coherencia de las exportaciones; no constituye validación independiente del modelo dinámico. Las cifras de la campaña sintética permanecen sujetas a revisión científica final y no son probabilidades de una red real. El valor de 0,6223 Hz no aparece ya como resultado destacado porque la identidad modal entre análisis sigue sin resolverse. El análisis de rutas se refiere a estados fasoriales DAE reequilibrados, no a maniobras físicas ni a EMT.

Las decisiones de liberación originales están conservadas en `evidence/claim_audit.md`, `evidence/scientific_claim_audit.md` y `research/bnd_h4_mechanism/tickets/STATUS.json`.
