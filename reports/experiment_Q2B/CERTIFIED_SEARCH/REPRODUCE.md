# Reproducción del cierre local

Ejecutar desde la raíz del repositorio, Julia 1.11.9, Project/Manifest existentes. No actualizar paquetes. El modelo congelado de ExpN se comprueba por hashes. PowerDynamics solamente se carga en el último comando.

Para conservar los resultados publicados, usar un directorio diferente:

```powershell
$env:EXPQ2B_CERT_OUT = Join-Path (Get-Location) 'reports/experiment_Q2B/REPRO_LOCAL'
julia --startup-file=no --project=. experiments/bnd_expQ2B/certified_search/run_search.jl reports/experiment_Q2B/CERTIFIED_SEARCH/MULTIPEAK_CURRENT.toml 120
julia --startup-file=no --project=. experiments/bnd_expQ2B/certified_search/run_certify.jl
julia --startup-file=no --project=. experiments/bnd_expQ2B/certified_search/run_full_audit.jl
julia --startup-file=no --project=. experiments/bnd_expQ2B/certified_search/freeze.jl
julia --startup-file=no --project=. experiments/bnd_expQ2B/certified_search/validate_independent_pd.jl
```

Los scripts que ofrecen consola al terminar aceptan `EXIT`. En ejecución por lotes terminan al cerrar stdin. `freeze.jl` rechaza sobrescribir un candidato ya congelado. Una tolerancia incumplida bloquea el freeze; una discrepancia de validación debe reportarse, nunca retocarse sobre el mismo archivo.

El seed es el punto de la búsqueda anterior donde dos extremos temporales competían. No se impone el resultado final como seed. `SQP.jl` reutiliza únicamente el álgebra explícita QP/NNLS de `local_certificate/solve.jl`; usa retenciones y ganancias conjuntamente. `MultiPeak.jl` separa picos temporales y `retract` corrige simultáneamente las restricciones activas. El punto congelado se evalúa con 143 polos físicos, prueba Riccati de banda completa y cotas temporales.

Los comandos siguientes ensamblan/verifican **el resultado publicado en CERTIFIED_SEARCH** después de que su validación PD haya terminado:

```powershell
python experiments/bnd_expQ2B/certified_search/finalize_results.py
python test/bnd_expQ2B/test_certified_artifacts.py
python experiments/bnd_expQ2B/certified_search/package_certified_results.py
```

`StaticBound.jl` y `run_static_bound.jl` conservan el intento no concluyente de cota global. Sus resultados Float64 no forman parte del certificado global. `test_riccati.jl`, los historiales sin SOC y los diagnósticos del punto estancado son trazabilidad de desarrollo, no pasos necesarios de la cadena aceptada.

El criterio publicado es local **numérico**, no un certificado formal con aritmética de intervalos. La comprobación de Hessiana usa tres tamaños de paso. La robustez usa aritmética de 256 bits al más cercano y una desigualdad estricta. El resultado no certifica todas las arquitecturas ni sustituye la cota global faltante.
