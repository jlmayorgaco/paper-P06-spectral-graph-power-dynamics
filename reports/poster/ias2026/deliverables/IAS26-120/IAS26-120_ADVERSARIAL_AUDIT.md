# IAS26-120 — auditoría adversarial y freeze

**Revisión científica y de claims:** `PASS_WITH_EXPLICIT_LIMITATIONS`  
**Estado del gate IAS26-120:** `INCOMPLETE_HOUSEKEEPING`  
**Alcance:** auditoría de la entrega IAS26-110 contra evidencia congelada. No se ejecutaron nuevos solves, barridos, campañas ni cálculos científicos.  
**Run de ensemble:** `20260927T144629Z_d0fecb32_ias26_060_operating_v1`  
**PDF auditado:** `../IAS26-110/IAS2026_IAS26-110_FINAL.pdf`  
**Decisión:** se congela el paquete de póster como resultado del programa actual, no como validación universal ni como cierre de los gates bloqueados.

## Resultado ejecutivo

El paquete contiene exactamente cuatro paneles principales (P1–P4) y una franja de reproducibilidad. La sustitución de los paneles planeados para F3/F4 está autorizada por el contrato visual versionado `IAS2026_POSTER_ARCHITECTURE_FINAL_V1`: P3 presenta el alcance sintético de operación y P4 la intervención fija más evidencia TDS ya existente. Esto no cambia los protocolos, ecuaciones, parámetros ni umbrales científicos.

Se confirmó que la fuente final `main.tex` compila en dos pasadas desde su carpeta de entrega y reproduce byte por byte el PDF final actualizado. El PDF es de una página, 2592 × 3456 pt (36 × 48 in); no hay advertencias TeX de cajas `Overfull`/`Underfull`. La página reconstruida es pixel-identical a la exportación revisada. Las figuras vectoriales fuente P1–P4 y la franja permanecen en PDF/SVG/PNG bajo el run inmutable IAS26-060. No se generaron ZIP ni resultados científicos nuevos.

La evidencia respalda el testigo nominal H4, la proximidad a singularidad del cierre colectivo en su frontera nominal, frecuencias descriptivas del ensemble sintético V4, un efecto mixto de `g=0.25`, reproducción Python–Julia del mismo modelo y una validación TDS condicionada. No respalda universalidad sobre otros puntos de operación, familias completas de dispositivos, redes, validación EMT/campo, ni el gate M1 estricto.

## Ledger `claim → run → archivo/fila → gate → figura/caption`

| Claim | Run y fuente primaria | Fila/clave y gate | Ubicación y wording permitido |
|---|---|---|---|
| C1 — H4 nominal es un blocker transversal de inclusión mínima | IAS26-020 `20260926T161116Z_d0fecb32_ias26_020_f1_audit_v1`; `tables/F1_FINAL_AUDIT.csv`; `tables/F1_FULL_SPECTRA_AUDITED.csv` | 16 portfolios; `portfolio_id=30+33+35+37` y sus 15 proper subsets; `tau_dec=1e-8 s^-1`; 15 `STABLE`, 1 `UNSTABLE`; PASS nominal | P1: “en el caso nominal auditado”; H4 `alpha_perp=+0.1270064678 s^-1`, `f=0.6222796695 Hz`. No es una identidad universal. |
| C2 — factores locales físicos regulares y cierre colectivo casi singular en la frontera nominal | IAS26-030 `20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1`; `tables/F2_BOUNDARY_PORT_AUDIT.csv`; `report/IAS26-030_SUMMARY.md` | Filas `device_bus=30,33,35,37`; `g*=0.20768140519037842`; orden `30+33+35+37`; mínimo local `sigma_min(I+M_ii)=0.4892354858`; `sigma_min(I+Q_H)=4.3936703687e-8`; PASS para la afirmación estrecha de F2 | P2: factores físicos no singulares frente a cierre colectivo próximo a singularidad. El propio archivo marca `m1_reduced_comparison=BLOCKED_NOT_PERFORMED`; no se afirma amortiguamiento local positivo ni validación de polos reducidos. |
| C3 — persistencia y cambio de testigo en el V4 sintético congelado | IAS26-060 `20260927T144629Z_d0fecb32_ias26_060_operating_v1`; `tables/MC_EVENT_RATES.csv`; `derived/SCENARIO_METRICS.csv`; `derived/MC_CASES.csv`; memo `IAS26-060_V4_BLOCKER_ATLAS_DESCRIPTIVE.md` | E1 `190/967`, Wilson 95% `[17.266%,22.271%]`; E2 `555/967`, `[54.254%,60.475%]`; 33 IDs físicamente inviables retenidos; 0 reemplazos; gate descriptivo condicionado | P3: solo la familia de 16 portfolios V4 y la envolvente sintética declarada. El atlas post hoc muestra 365 escenarios con testigo alternativo; no demuestra invariancia del mecanismo `I+Q_B` para dichos testigos. |
| C4 — el retuning fijo tiene efecto mixto | IAS26-060 y config congelada IAS26-050 | E3 mejora `594/967`, Wilson `[58.320%,64.444%]`; E4 rescata H4 `186/967`, `[16.874%,21.839%]`; E5 deteriora `373/967`, `[35.556%,41.680%]`; gates por definiciones preregistradas | P4: `g=0.03625 → 0.25`; candidato fijado antes del ensemble a partir del barrido canónico preexistente; no es óptimo ni cura universal. |
| C5 — paridad de veredictos entre implementaciones | IAS26-060 `derived/CROSSCODE_SUMMARY.json`; `tables/IAS26-FINAL_EVIDENCE_STRIP.csv` | 50 escenarios, 100 comparaciones, `100/100` veredictos; máximo `|Δalpha|=3.74954e-9 s^-1`, máximo `|Δf|=1.53494e-10 Hz`; 0/100 identidades de modo resueltas | Franja: reproducción Python–Julia del mismo modelo GFL11, no validación física independiente. |
| C6 — concordancia TDS condicional | IAS26-060/080 `tables/IAS26-080_TDS_RESULTS.csv`, `tables/IAS26-FINAL_P4_TDS_TRACES.csv`, trazas en `raw/tds/traces/` | 90 IDs de trayectoria retenidos: 60 `COMPLETED`, 30 `SOLVER_FAILED`; signos concuerdan 56/60; los 4 desacuerdos son de frecuencia TDS cero y R² ausente | P4: solo tres trazas canónicas D2/A1 como ejemplos; los fallos y desacuerdos se declaran. Evidencia de DAE fasorial, no EMT/campo. |
| C7 — exclusiones/gates bloqueados | `tickets/STATUS.json`; `config/IAS2026_POSTER_ARCHITECTURE_FINAL_V1.json`; readiness report del run | IAS26-010 `BLOCKED_M1_STRICT`; IAS26-040 bloqueado; campaña GFL/GFM no ejecutada | F3 y F4 no aparecen como resultados. No claims de eta/M2/M3, GFM/F4 ni pase de M1. |

La configuración IAS26-050 registra `frozen_at_utc=2026-09-26T17:39:53Z` y selecciona `g=0.25` desde un barrido canónico preexistente, “not optimized on development or holdout scenarios”. El pre-run IAS26-060 fija `tau_dec=1e-8`, `alpha_perp` de espectro transversal completo como métrica primaria y `alpha_Omega` de 0.3–1.5 Hz solo como métrica secundaria separada. En F1 se verificó que ambas columnas coinciden en las 16 filas nominales; esto no autoriza sustituir una por otra en general.

## Diez ataques obligatorios

| # | Ataque | Resultado de revisión | Disposición / límite explícito |
|---:|---|---|---|
| 1 | ¿H4 es minimal usando `alpha_perp` y el umbral congelado? | **PASS nominal.** Las 16 filas de F1 dan 15 subsets propios estables y H4 inestable con `tau_dec=1e-8`; H4 tiene dimensión transversal 84 y 84 valores espectrales auditados. | El claim es del H4 nominal IEEE-39/TX4 congelado; no de todos los escenarios. |
| 2 | ¿Se confundió regularidad local con estabilidad local? | **PASS de wording.** La cantidad es el menor valor singular de los factores físicos `I+M_ii`; todos son no singulares en la frontera auditada. | No inferir autovalores locales, amortiguamiento positivo, ni estabilidad autónoma de cada dispositivo. |
| 3 | ¿Coexisten cierre casi singular y crossing en misma frecuencia/orden? | **PASS para la frontera auditada.** F2 usa el mismo run/operating point, H4 y orden de puertos; todas las filas dan `f=0.706424788006417 Hz`, `g*=0.20768140519037842` y el mismo `sigma_min(I+Q_H)=4.39367e-8`. La diferencia de frecuencia local-colectiva máxima es `1.63328e-11 Hz`. | El cruce del retorno y la identidad de Schur de F2 son la evidencia declarada. La comparación estricta de reducción M1 sigue `BLOCKED_NOT_PERFORMED`; no se promociona como validación reducida–DAE. |
| 4 | ¿Los porcentajes MC ocultan denominadores o exclusiones? | **PASS.** 1,000 IDs precongelados; 967 factibles y 33 `GENERATOR_Q_LIMIT`; 17,000 filas de caso; cero reemplazos y `outcome_dependent_filtering=false`. Wilson 95% y denominador 967 aparecen donde aplican. | Frecuencias condicionales de una envolvente sintética, no probabilidades del IEEE-39 real ni de la población de redes. |
| 5 | ¿Se eligió `g=0.25` después de mirar el holdout? | **PASS.** IAS26-050 está congelado antes de IAS26-060 y enlaza la elección a un barrido canónico ya existente. | Solo se prueba este candidato fijo. No se afirma optimalidad ni ajuste adaptativo por escenario. |
| 6 | ¿Los números TDS provienen de trazas medidas/ajustadas y se muestran los fallos? | **PASS con limitación visible.** Se conservaron 90 trayectorias, 60 ajustes finitos completados y 30 fallos de solver; 56/60 signos concuerdan. Los 4 desacuerdos tienen `frequency_tds_hz=0` y `fit_r2_diagnostic` ausente. P4 usa tres trazas canónicas D2/A1, no una selección estadística de las 60. | No borrar fallos, no inferir identidad de modo y no llamar estos ensayos EMT/field. |
| 7 | ¿Python–Julia son validaciones físicas independientes? | **PASS de alcance.** Hay 100/100 acuerdos para el mismo modelo GFL11 y discrepancias numéricas pequeñas declaradas. | Todas las 100 identidades de modo están sin resolver; ambas implementaciones no son modelos físicos independientes. |
| 8 | ¿Se cambió el espectro completo por una banda conveniente? | **PASS.** F1 conserva `F1_FULL_SPECTRA_AUDITED.csv`; para H4 se auditan 84 valores transversales. `alpha_perp` es el criterio primario; `alpha_Omega` se reporta separadamente como secundario. | La igualdad observada alpha-perp/alpha-Omega en las 16 filas nominales no se extrapola a la campaña MC. No se hace claim de tracking modal. |
| 9 | ¿Se mezclaron BND39, TX4, otras tecnologías o modelos? | **PASS.** Los paneles auditados proceden del benchmark IEEE-39 con la configuración TX4 congelada y H4/puertos `30,33,35,37`; el ensemble usa su contrato/configuración IEEE-39 TX4. | No importar claims del censo BND39/V9, GFL10, red/carga distinta o ParaEMT. F4/GFM no se corrió. |
| 10 | ¿Hay un claim EMT/campo? | **PASS.** Caption y alcance de P4 identifican simulación DAE fasorial. | Sin evidencia EMT ni validación de campo; ambas quedan explícitamente fuera. |

## Argumento oral auditado

**15 segundos.** En el IEEE-39/TX4 nominal, H4 es inestable aunque sus 15 subconjuntos propios son estables. En su frontera nominal, los factores locales físicos siguen no singulares mientras el cierre colectivo se acerca a singularidad. Bajo 967 puntos sintéticos factibles, H4 no es universal y `g=0.25` a veces empeora; el póster muestra ambas limitaciones.

**60 segundos.** El primer panel establece el fenómeno en el caso nominal: el conjunto de reemplazos `{30,33,35,37}` cruza la estabilidad transversal, mientras los 15 subsets propios permanecen del lado estable con `tau_dec=1e-8`. El segundo panel ubica la frontera nominal: a `g*=0.2076814`, los cuatro factores físicos `I+M_ii` son no singulares, pero `I+Q_H` tiene valor singular mínimo `4.39×10^-8`; esto no prueba estabilidad local ni supera el gate M1 estricto. El ensemble congelado es sintético: 967 de 1,000 IDs fueron factibles. H4 mínimo aparece en 190; algún blocker mínimo dentro del V4 censado aparece en 555. El retuning fijo, preseleccionado en `g=0.25`, mejora 594 casos pero deteriora 373 y rescata H4 en 186. TDS y Python–Julia aportan chequeos condicionados de DAE fasorial y reproducción same-model, no EMT ni validación independiente.

**3 minutos.** Este trabajo estudia composabilidad de reemplazos en una red IEEE-39 con la configuración TX4 congelada, no una población de redes ni el benchmark BND39. Primero, el censo nominal exacto de 16 portfolios usa la abscisa espectral transversal completa y un umbral de `1e-8 s^-1`: los 15 subsets propios de H4 son estables y H4 es inestable con `alpha_perp=+0.1270 s^-1` cerca de `0.6223 Hz`. La auditoría conserva los espectros; la banda 0.3–1.5 Hz se mantiene como métrica secundaria y no sustituye la decisión primaria. Segundo, en la frontera nominal auditada en el mismo orden de puertos, los factores físicos locales tienen mínimos singulares de aproximadamente 0.489 a 0.727, mientras `I+Q_H` baja a `4.39×10^-8` alrededor de `0.7064 Hz`. Esto respalda una afirmación de regularidad local frente a cierre colectivo casi singular, pero no demuestra amortiguamiento local positivo ni completa la reconciliación de reducción M1, que continúa bloqueada. Tercero, el pre-run congeló 1,000 IDs sintéticos; 33 fueron físicamente inviables y se retuvieron, sin reemplazos. Entre los 967 válidos, H4 sigue siendo un blocker mínimo en 190; alguna incompatibilidad mínima no singleton en el sublattice V4 aparece en 555. El atlas descriptivo encuentra 365 testigos alternativos, pero no recalculó el mecanismo colectivo para cada nuevo testigo, por lo que no afirmamos invariancia del mecanismo. Cuarto, `g=0.25` fue fijado antes del ensemble. Mejora el margen en 594, rescata H4 en 186 y lo deteriora en 373; no es una cura universal ni un óptimo. Finalmente, los signos TDS concuerdan en 56 de 60 fits finitos, con 30 fallos de solver retenidos y cuatro desacuerdos de frecuencia cero/R² ausente; Python y Julia concuerdan en 100 de 100 veredictos del mismo modelo, pero no resuelven identidad modal ni constituyen validación física independiente. El alcance final es por tanto un resultado nominal más evidencia sintética condicionada, con limitaciones y gates bloqueados visibles.

## Claims seguros y no seguros

**Seguros:** (a) H4 es un blocker mínimo transversal en el caso nominal congelado; (b) en la frontera nominal auditada los factores locales físicos son no singulares y el cierre colectivo está próximo a singularidad; (c) las frecuencias de persistencia y mejora son condicionales a los 967 escenarios sintéticos factibles y al V4 censado; (d) el testigo cambia en parte de esos escenarios; (e) el candidato fijo `g=0.25` tiene resultado mixto; (f) la paridad Python–Julia y los chequeos TDS tienen exactamente el alcance declarado.

**No seguros/prohibidos:** H4 universal o global; invariancia de `I+Q_B` para los testigos alternativos; `g=0.25` como óptimo o cura; frecuencia de MC como probabilidad real; estabilidad por mera regularidad local; validación independiente o identidad modal resuelta por Python–Julia; cruce de banda como sustituto general de `alpha_perp`; claims de BND39/V9, GFM/F4, eta/M2/M3, M1 estricto aprobado, EMT o campo.

## Regeneración, hashes y organización

Construcción verificada desde `reports/poster/ias2026/deliverables/IAS26-110/`, en dos pasadas de XeLaTeX. El resultado construido coincide byte a byte con el PDF de entrega:

| Archivo | SHA-256 |
|---|---|
| `../IAS26-110/IAS2026_IAS26-110_FINAL.pdf` | `9F0A4957162BA96F627A19ADEE714347006A8C68BE82853D8B1A39B6EC03930C` |
| `../IAS26-110/main.tex` | `7EABD3E71E2F9861C192508134A633982881C9B68E08CA26BAAF6AABAD613F91` |
| Config preregistrada IAS26-050 | `54F5B18FB6752EFE0BB17C1474ABF6C89C5AE9D0F7064BE9D67487E5695691E1` |
| Config congelada IAS26-060 | `75DB11E64299757BBE170F5ECADCB8C1A3F2EA969583E3C95BCDEBC008BFB4B4` |
| Manifest de 1,000 IDs | `45C425FA4D146AA3A0A3BF26CFBDDA1C99590DA67EA6C4318CF788C1006B6EB3` |
| F1 audit CSV | `4EA8042A80BB1E68BEE1F272B0DAF7BAA4F23A9611FA6B5975CD3239AC11F718` |
| F1 espectros auditados | `9DA3ED87AFB445CF155C815FEF50FE1A2AF0D28D4A753C6B968ADA698E92B919` |
| F2 boundary/port audit | `B45E64A564CECA1D75DD5D548F1D39A6765376930F9D5D8D1D7BD93F62509B9B` |
| MC event rates | `A0AC215148F5AA30D76C6F2ED20A68C9252D5691A4B62965506383F43230FCD2` |
| TDS results | `F374AF001F59450FC8067E2503CB1DDCF67CECF9ABEE9946F9C490A40BBADC8E` |
| Cross-code summary | `51EFC35929E3F9E9B7D495F13F001923598F035481863A29D11766D3A49EF3A9` |
| Versioned poster mapping | `D255237D0E6D239A140217AF4F902484DD2372AB19D7C1A9985CB06B532E49F6` |

La fuente TeX, el PDF, suplemento e informe están bajo `reports/poster/ias2026/deliverables/`. **No pude completar la limpieza de scratch:** `exec_command` rechazó `Remove-Item` por política, incluso después de comprobar los paths absolutos dentro del repo. Por tanto permanecen auxiliares producidos durante esta QA en `reports/poster/ias2026/deliverables/IAS26-110/$qa/` (incluye un PDF recompilado y auxiliares TeX), `reports/poster/ias2026/deliverables/IAS26-110/texput.log`, `reports/poster/ias2026/build/IAS26-110/tmp/pdfs/` (build/preview de QA), `reports/poster/ias2026/.tmp-ias26-110-qa/` (vacío) y `<repo>/$qaDir/` (dos logs TeX). No apunté a runs congelados ni al borrador histórico. El borrador preexistente `reports/poster/ias2026/main.tex` y su antiguo `build/` permanecen intactos. No se creó ningún ZIP. Esta limitación impide cerrar el requisito de organización de IAS26-120 hasta retirar esos scratch.

## Cierre del ticket

**Gate IAS26-120: INCOMPLETE_HOUSEKEEPING.** La auditoría de los cuatro paneles y la franja pasa con limitaciones explícitas; el PDF regenera desde la fuente; los resultados negativos y bloqueados no se ocultan. El ticket 110 queda completo y congelado. IAS26-120 no se marca cerrado hasta retirar los scratch listados arriba. IAS26-010 sigue `BLOCKED_M1_STRICT`; F3/M1, eta/M2/M3 y F4/GFM no están establecidos. No se inició ninguna campaña siguiente.
