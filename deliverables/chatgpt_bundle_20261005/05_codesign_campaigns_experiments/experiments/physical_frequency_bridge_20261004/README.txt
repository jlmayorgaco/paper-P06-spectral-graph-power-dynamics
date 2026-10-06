Leer primero REPORT.txt. Resultados regenerables: RESULTS.json, EVENT_RESULTS.csv,
ENERGY_AND_FREQUENCY_AREA.csv, PAIRED_COMPARISON.csv y GOVERNOR_BRANCHES.csv.

Teoría: PHYSICAL_BRIDGE_SECTION.tex (fragmento integrado en el THEORY.tex existente).
Revisiones independientes: MATH_REVIEW.txt, EVIDENCE_REVIEW.txt.
Código: PhysicalBridge.jl, run_validation.jl, analyze.py, build_report.py.
Freeze: PREREGISTRATION.json + hash; RUN_ENVIRONMENT.toml; MANIFEST.json.
Figuras: FIG_PHYSICAL_BRIDGE y FIG_FREQUENCY_AREA_BRIDGE, PNG y SVG.
Trazas completas: *_SENSORS.csv. Balances derivados: *_BALANCE.csv.

REPRODUCCIÓN SIN SOBREESCRIBIR ESTA CAMPAÑA
Desde la raíz del repo, copiar sólo los siguientes archivos a un nuevo directorio
experiments/physical_frequency_bridge_replay_FECHA_HORA/:
PhysicalBridge.jl, run_validation.jl, analyze.py, PREREGISTRATION.json,
PREREGISTRATION.sha256. Mantener Project/Manifest y las fuentes/input hashes.
Luego ejecutar desde la raíz:
julia --project=. experiments/physical_frequency_bridge_replay_FECHA_HORA/run_validation.jl
python experiments/physical_frequency_bridge_replay_FECHA_HORA/analyze.py
No ejecutar prepare.py otra vez: rechaza una preregistración existente.
El replay usa los mismos cinco casos y no ajusta parámetros por sus resultados.
El análisis requiere Python, numpy, pandas, scipy y matplotlib; versiones en MANIFEST.
No hay commits ni pushes en esta campaña.
