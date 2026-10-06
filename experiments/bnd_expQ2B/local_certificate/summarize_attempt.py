"""Summarize measured local-certificate evidence without altering old artifacts."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports/experiment_Q2B/LOCAL_CERTIFICATE"
FINAL = OUT / "FINITE_EXTREMA"
status = tomllib.loads((FINAL / "LOCAL_ATTEMPT_STATUS.toml").read_text())
audit = tomllib.loads((FINAL / "SECOND_ORDER_AUDIT.toml").read_text())
gate = list(csv.DictReader((FINAL / "TABLE_gradient_gate.csv").open(newline="")))
freeze = json.loads((ROOT / "reports/experiment_N/MODEL_FREEZE.json").read_text())
integrity = []
for name, expected in freeze["inputs"].items():
    actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    integrity.append({"file": name, "expected": expected, "actual": actual, "pass": actual == expected})
assert all(row["pass"] for row in integrity), "Frozen ExpN file changed"
original = ROOT / "reports/experiment_Q2B/CORE/Z_Q2B_SECURE_T05_FINAL.toml"
assert hashlib.sha256(original.read_bytes()).hexdigest() == "f53250383cdb906a2d0cb3130ca9c1e72c63b472717c5c82aea1d5838ac785ce"
provenance = {
    "parent_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
    "model_sha": freeze["MODEL_SHA"], "frozen_model_inputs": integrity,
    "original_Q2B_candidate_unchanged": True,
    "current_sources": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in Path(__file__).parent.glob("*.jl")},
    "dependencies": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                     for name in ("Project.toml", "Manifest.toml")},
}
(OUT / "AUDIT_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
rows = []
for parent in (OUT, OUT / "WITHOUT_38", OUT / "WITHOUT_38_REFINED", FINAL):
    p = parent / "LOCAL_ATTEMPT_STATUS.toml"
    if not p.exists():
        continue
    d = tomllib.loads(p.read_text())
    rows.append({"attempt": parent.name, **{key: d.get(key, "NOT_RECORDED") for key in
        ("status", "retained_MW", "stationarity", "primal", "complementarity", "elapsed_s", "model_evaluations",
         "robust_frequency_factorizations", "local_certified", "global_certified")}})
with (OUT / "TABLE_LOCAL_ATTEMPTS.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
passed = sum(r["pass"].lower() == "true" for r in gate)
text = f"""# Resultado de la búsqueda de certificado local

## Conclusión

**NO SE HA CERTIFICADO UN NUEVO ÓPTIMO LOCAL SEGURO.** Se ejecutó una corrección
conjunta de retenciones y 20 ganancias, con derivadas analíticas y SQP de conjunto
activo. El resultado final es `{status['status']}`. El punto obtenido es
exploratorio, no un candidato congelado, una cota superior segura ni un resultado
validado en PowerDynamics.

Se conserva la rama `{provenance['branch']}`, el modelo ExpN y el candidato Q2B
congelado. Se verificaron {len(integrity)}/{len(integrity)} hashes del modelo. No
se ejecutó PowerDynamics dentro del optimizador. No se hizo commit ni push.

## Problema realmente evaluado

- Retenciones SG y las 20 ganancias PLL libres dentro de los límites congelados.
- Modelo físico y reparto P/Q de ExpN; dimensión reconstruida al retirar SG38.
- Evento bus16, +100 MW sostenidos, ventana de fase causal T=0.5 s; límites
  0.5 Hz y 0.5 Hz/s. Buses8/29 no son restricciones ocultas.
- Espectro físico completo, radio robusto observado, frecuencia estacionaria,
  extremos de frecuencia y RoCoF. La cobertura global de frecuencia/tiempo sigue
  pendiente: esta distinción impide aceptar factibilidad completa.

## Resultado medido del último intento

| Comprobación | Valor |
|---|---:|
| Retención exploratoria | {audit['retained_MW']:.9f} MW |
| Residuo primal normalizado | {audit['primal']:.9g} |
| Estacionariedad normalizada | {audit['stationarity']:.9g} |
| Complementariedad | {audit['complementarity']:.9g} |
| Rango / filas activas | {audit['LICQ_rank']} / {audit['active_count']} |
| Alpha | {audit['alpha']:.12g} s^-1 |
| Beta observado (no cota inferior) | {audit['beta_observed']:.12g} |
| Pico lineal de frecuencia | {audit['Fpeak']:.12g} Hz |
| Pico lineal de RoCoF | {audit['Rpeak']:.12g} Hz/s |
| Condición de autovectores | {audit['eigenvector_condition']:.9g} |
| SOSC | {audit['SOSC']} |
| Tiempo SQP, sin compilación | {status['elapsed_s']:.3f} s |
| Evaluaciones del modelo en esta ejecución | {status['model_evaluations']} |
| Factorizaciones para robustez | {status['robust_frequency_factorizations']} |

El objetivo se dividió por 1000; epsilon y ganancias se normalizaron por sus
intervalos. Los residuos se refieren a esas unidades, no son distancias en MW
al óptimo. Los números anteriores se reevaluaron con una malla temporal más fina
para el contraste; las raíces temporales se refinan con bisección.

## Ecuaciones y pruebas

La auditoría final contrastó {passed}/{len(gate)} combinaciones familia/coordenada
con diferencias centrales en un punto interior cercano. Incluye el espectro,
radio robusto, ganancia DC, extremos de frecuencia y RoCoF. La aceptación usa
error relativo <1e-5 o absoluto físico <1e-7 para derivadas pequeñas, seleccionando
el paso compatible con el condicionamiento. Se conservan todos los errores,
incluidos pasos que fallan, en los CSV. Esto valida los puntos ensayados y no
certifica derivadas uniformemente en todo el dominio.

La discrepancia del pico en la coordenada17 se estudió con cinco pasos y malla
temporal más fina. En h=2e-4, el máximo error físico fue aproximadamente
1.474e-8, con condición de autovectores aproximadamente1.046e5. No se cambió el
umbral de aceptación.

Se corrigieron en este evaluador nuevo el salto de fase en radianes, la derivada
del límite DC y la localización estacionaria de la frecuencia robusta. Se
separaron los extremos transitorios del límite estacionario para evitar contar
una muestra final arbitraria como restricción activa independiente. Los archivos
históricos no se reescribieron. Nueve comprobaciones sintéticas verifican el
refinamiento de la frecuencia robusta y su derivada de envolvente.

## Alcance y bloqueo

No se ha completado la cadena factibilidad -> KKT -> cualificación -> Hessiano
positivo en el cono crítico -> cobertura robusta/temporal -> validación PD.
El intento no autoriza un nuevo mínimo local, ni global, ni una distancia al
óptimo seguro. Un pico robusto puntual sólo da una **cota superior de beta_star**;
no prueba beta_star >= beta_req. Se mantienen abiertas las arquitecturas
competidoras y las ramas no exploradas.

La dificultad actual pertenece a la convergencia y certificación numérica del
problema continuo, no demuestra inexistencia de un diseño seguro. No se sustituyó
esa dificultad por un barrido ni por una afirmación de optimalidad.

## Qué sí sabemos sobre el problema nominal anterior

ExpP/P3 registra un certificado local numérico del problema **sólo modal** en
soporte{{38}}, sujeto a sus salvedades de guarda y precisión. No se extiende al
problema robusto con evento100 MW. Bajo la factibilidad nominal aceptada en
ExpN/P, J>=0 y el candidato de1.124344478 MW dan:

    0 <= J*_modal <= 1.124344478 MW.

Por tanto, aquel candidato está a como máximo1.124344478 MW del óptimo modal
global (0.0208106% del despacho original). Es una cota trivial útil; no es un
certificado de igualdad ni satisface la brecha objetivo0.01 MW. No proporciona
una cota inferior ni superior del problema seguro.

## Reproducción

Desde la raíz del repositorio, en PowerShell:

```powershell
$env:EXPQ2B_LOCAL_LABEL='FINITE_EXTREMA'
julia --startup-file=no --project=. experiments/bnd_expQ2B/local_certificate/audit.jl reports/experiment_Q2B/LOCAL_CERTIFICATE/WITHOUT_38_REFINED/CURRENT_LOCAL_CANDIDATE.toml
julia --startup-file=no --project=. experiments/bnd_expQ2B/local_certificate/run_and_certify.jl reports/experiment_Q2B/LOCAL_CERTIFICATE/WITHOUT_38_REFINED/CURRENT_LOCAL_CANDIDATE.toml 120
python experiments/bnd_expQ2B/local_certificate/summarize_attempt.py
```

Los primeros intentos de TABLE_LOCAL_ATTEMPTS.csv documentan iteraciones de
desarrollo del evaluador. La receta anterior corresponde a su versión final,
cuyos hashes figuran en AUDIT_PROVENANCE.json. Consultar
CERTIFICATE_SCOPE_AND_NOVELTY.md para las condiciones matemáticas y el contraste
bibliográfico preliminar.
"""
(OUT / "RESULTADO_CERTIFICADO_LOCAL.md").write_text(text, encoding="utf-8")
print(text.split("## Problema realmente evaluado")[0])
