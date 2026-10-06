# Resultado de la búsqueda de certificado local

## Conclusión

**NO SE HA CERTIFICADO UN NUEVO ÓPTIMO LOCAL SEGURO.** Se ejecutó una corrección
conjunta de retenciones y 20 ganancias, con derivadas analíticas y SQP de conjunto
activo. El resultado final es `MERIT_STAGNATION`. El punto obtenido es
exploratorio, no un candidato congelado, una cota superior segura ni un resultado
validado en PowerDynamics.

Se conserva la rama `research/expQ2B-secure-optimum`, el modelo ExpN y el candidato Q2B
congelado. Se verificaron 17/17 hashes del modelo. No
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
| Retención exploratoria | 314.545169229 MW |
| Residuo primal normalizado | 5.90233721e-05 |
| Estacionariedad normalizada | 0.238002837 |
| Complementariedad | 1.94705813e-08 |
| Rango / filas activas | 11 / 11 |
| Alpha | -0.0500479390425 s^-1 |
| Beta observado (no cota inferior) | 1.69902041208e-06 |
| Pico lineal de frecuencia | 0.500000083586 Hz |
| Pico lineal de RoCoF | 0.160256397207 Hz/s |
| Condición de autovectores | 102896.471 |
| SOSC | NOT_TESTED_FIRST_ORDER_GATE_FAILED |
| Tiempo SQP, sin compilación | 65.287 s |
| Evaluaciones del modelo en esta ejecución | 377 |
| Factorizaciones para robustez | 754 |

El objetivo se dividió por 1000; epsilon y ganancias se normalizaron por sus
intervalos. Los residuos se refieren a esas unidades, no son distancias en MW
al óptimo. Los números anteriores se reevaluaron con una malla temporal más fina
para el contraste; las raíces temporales se refinan con bisección.

## Ecuaciones y pruebas

La auditoría final contrastó 125/125 combinaciones familia/coordenada
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
soporte{38}, sujeto a sus salvedades de guarda y precisión. No se extiende al
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
