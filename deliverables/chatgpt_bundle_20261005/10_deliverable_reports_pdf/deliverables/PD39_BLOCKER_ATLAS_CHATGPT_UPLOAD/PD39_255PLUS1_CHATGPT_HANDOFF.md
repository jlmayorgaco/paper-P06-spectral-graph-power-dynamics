# PD39 255+1 ChatGPT handoff

La campaña nueva se ejecutó en la rama `research/pd39-255plus1-mechanism-validation`, creada desde `902403cf`. No se usó ni se modificó la rama confirmatoria antigua.

## Resultado ejecutivo

`FINAL CASE: C`

La lectura correcta es parcial y negativa en puntos esenciales: Gate B pasa, pero Gate A falla su criterio independiente de alpha; Gate C no generaliza el patrón exacto; Gate D es parcial; Gate E está bloqueado por semántica/API ausente. Por eso la recomendación para el poster IAS es **NO** rediseñarlo alrededor de una afirmación universal 255+1.

## Gates

- Numerical audit: FAIL. 243/243 filas de tolerancia y 81/81 auditorías independientes se ejecutaron. La comparación 1e-10 vs 1e-12 fue estable, y el descriptor generalizado concordó; pero 72/81 diferencias de alpha por Jacobiano FD excedieron 1e-4 s^-1 (máximo 0.138325400017). No hubo cambios de signo, pero sí 67 cambios de clase 0.05.
- High-PLL TDS: PASS. 12/12 casos, bus 39, pulso +1% P/Q de 1.0 a 1.1 s, integración a 20 s, 12/12 con signo consistente. Medias de tasa: V8 0.1119606096, best 7/8 -0.2850800408, worst 7/8 -0.1140069162 s^-1.
- H0 over 24 holdouts: total 42; vector [4, 2, 2, 0, 3, 0, 2, 2, 0, 0, 0, 19, 0, 0, 0, 1, 1, 2, 0, 3, 0, 0, 0, 1].
- H0.05 over 24 holdouts: total 47; vector [4, 2, 2, 0, 5, 0, 2, 3, 0, 0, 0, 19, 0, 0, 0, 2, 1, 2, 0, 4, 0, 0, 0, 1].
- Penetration/composition/mixed: 90 comparaciones aggregate-matched con diferencias MW/MVA agregadas nulas y efectos alpha no nulos; no se ejecutó intervención mixta y no se reclama una ventaja mixta.
- Mode mechanism: critical modes classified electromechanical/control. V8 high-PLL is oscillatory at about 0.32-0.37 Hz with positive alpha; 7/8 modes are portfolio-dependent and often distinct in common bus-voltage MAC. Homotopy continua SG->GFL: BLOCKED.
- Closure: BLOCKED. No existen objetos K,D,Q exactos documentados ni identidad determinant/Q->-1 utilizable; no se fabricó proxy.
- IAS poster redesign: NO.

## Los ocho 7/8 margins de discovery

- missing 34: portfolio 30;32;33;35;36;37;38; converted MW 4112.0; m9 0.0858347015467; nominal alpha -0.0967153322641
- missing 30: portfolio 32;33;34;35;36;37;38; converted MW 4370.0; m9 0.126228983466; nominal alpha -0.130545971751
- missing 35: portfolio 30;32;33;34;36;37;38; converted MW 3970.0; m9 0.136356176224; nominal alpha -0.138251947158
- missing 38: portfolio 30;32;33;34;35;36;37; converted MW 3790.0; m9 0.13801971639; nominal alpha -0.138019716734
- missing 32: portfolio 30;33;34;35;36;37;38; converted MW 3970.0; m9 0.138022315453; nominal alpha -0.138022318748
- missing 36: portfolio 30;32;33;34;35;37;38; converted MW 4060.0; m9 0.138057659073; nominal alpha -0.13805765951
- missing 33: portfolio 30;32;34;35;36;37;38; converted MW 3988.0; m9 0.138295566944; nominal alpha -0.138295569249
- missing 37: portfolio 30;32;33;34;35;36;38; converted MW 4080.0; m9 0.138325419498; nominal alpha -0.138325422033

## Artefactos

- `docs/PD39_255PLUS1_FINAL_REPORT.pdf`
- `docs/PD39_255PLUS1_FINAL_REPORT.md`
- `results/PD39_255PLUS1_HEADLINE.json`
- `results/PD39_255PLUS1_MASTER_CLAIMS.csv`
- `results/PD39_255PLUS1_NUMERICAL_TRUTH_AUDIT.csv`
- `results/PD39_255PLUS1_HIGH_PLL_TDS.csv`
- `results/PD39_255PLUS1_HOLDOUT_CENSUS.csv`
- `results/PD39_255PLUS1_HOLDOUT_CONDITION_SUMMARY.csv`
- `results/PD39_255PLUS1_PENETRATION_COMPOSITION.csv`
- `results/PD39_255PLUS1_MODE_TRACKING.csv`
- `results/PD39_255PLUS1_HOMOTOPY_AUDIT.csv`
- `results/PD39_255PLUS1_CLOSURE_AUDIT.csv`
- `deliverables/PD39_255PLUS1_CHATGPT_UPLOAD.zip`

No push was performed. The excluded campaigns (weak nodes/links, structured radius, repairs, planners, co-design, IEEE-68, EMT, new controller search) were not run.

FINAL CASE: C
