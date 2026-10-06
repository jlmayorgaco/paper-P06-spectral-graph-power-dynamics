# Experimento ejecutado: amortiguamiento colectivo con puertos físicos

**Resultado: sí existe una contribución colectiva física y material en este
modelo; su signo aislado no proporciona una ley suficiente de estabilidad PLL
ni una nueva solución de reemplazo óptimo.**

Se mantuvo el reemplazo en **88.4551408%**, 4779.019928 MW GFL y
623.741162 MW SG. No se optimizó rho, Kp o Ki en esta campaña.
N = ganancias nominales; Z = ajuste anterior sin retraso; T = ajuste anterior
de margen de latencia. Sus vectores completos están en models/*/design.toml.

## Resultado físico — EXACT_IDENTITY + NUMERICALLY_VALIDATED

Se derivó la impedancia de par/velocidad Z(s), exacta para la linealización de
orden completo con retrasos puros, al eliminar estados ocultos a frecuencia
finita. No es una transferencia exacta del sistema no lineal. Su parte hermítica D_H describe trabajo
mecánico incremental para un movimiento armónico prescrito. Se separó la diagonal
nodal del término entre nodos sin llamarlos laplacianos ni matrices constantes.
Las hipótesis y unidades están en THEORY.md.

En T, 45 ms, sobre el movimiento del modo PLL seguido:

- contribución nodal: **38.658333**;
- contribución colectiva: **-737.485824**;
- suma: **-698.827491**.

Son coeficientes de amortiguamiento en la normalización de velocidad fijada
(norma del fasor de velocidad = 1), no MW de generación reemplazable. Aquí
omitir el término colectivo cambia el signo. Su fracción absoluta es
95.019%. No todos los diseños/frecuencias muestran
el mismo efecto; se incluyen todos los casos del protocolo.

## Polos del modelo con exponenciales exactas — SUPPORTED_LOCAL

Máximo real entre las ramas PLL rastreadas (no todo el espectro infinito):

| design | tau_ms | real | frequency_hz |
| --- | --- | --- | --- |
| N | 40 | 0.25266355 | 5.2104607 |
| N | 45 | 2.2815843 | 4.8355693 |
| T | 40 | -1.8101291 | 4.5844044 |
| T | 45 | 0.43900926 | 4.6192868 |
| Z | 40 | 1.2705099 | 5.4022206 |
| Z | 45 | 3.1993408 | 4.9931953 |


Primer cruce Re(s)=0 entre esas ramas:

| design | tau_ms | frequency_hz | mode_damping_nodal | mode_damping_cross |
| --- | --- | --- | --- | --- |
| N | 39.48362 | 5.2518676 | 1652.8743 | -1652.8742 |
| T | 43.916999 | 4.6890663 | 4151.0249 | -4151.0249 |
| Z | 37.477197 | 5.6316245 | -1824.0752 | 1824.0752 |


Estos cruces de cero son distintos de los límites anteriores de seguridad
Re(s)=-0.05. El barrido no certifica ausencia de todas las demás raíces DDE.
En particular, los modos lentos pueden ser más cercanos al eje que la familia
PLL cuando esa familia está bien amortiguada.

## Validación no lineal con retraso puro — SUPPORTED_LOCAL

Se ejecutó el método de pasos con Rodas5P, usando historial de un modo exacto,
dos amplitudes (1e-5 y 1e-4 rad PLL), seis casos N/T a 0, 40 y 45 ms, y 3 s.
La referencia lineal es la solución analítica Re(a v exp(lambda t)) del DDE.

| case | predicted_alpha | fitted_alpha | relative_trajectory_error |
| --- | --- | --- | --- |
| N_0ms | -0.080148382 | N/A | 1.9939193e-06 |
| N_40ms | 0.25266355 | 0.25266353 | 3.4052177e-05 |
| N_45ms | 2.2815843 | 2.2817186 | 0.00061385762 |
| T_0ms | -0.080148351 | N/A | 1.9535041e-06 |
| T_40ms | -1.8101291 | -1.8101288 | 1.278242e-05 |
| T_45ms | 0.43900926 | 0.4390093 | 0.00027159157 |


El ajuste de crecimiento se omite a 0 ms: el modo de ~0.015 Hz no completa un
ciclo en 3 s; la trayectoria sí se compara. La tabla completa mantiene ambas
amplitudes y cualquier refinamiento requerido. Máximo error relativo de
trayectoria entre ejecuciones estándar: 0.00730766.
Si el error crece con amplitud y no con refinamiento, es efecto no lineal.
**No son los cinco eventos de +/-100 MW, ni una demostración de seguridad robusta.**

## Acción colectiva sobre el polo — seguimiento exploratorio

Se derivó la sensibilidad del polo con Z_p y Z_s, conservando las partes reactiva
y disipativa y los modos izquierdo/derecho. Se separó la acción diagonal y
colectiva, validada contra la característica completa y diferencias centradas.

| design | slope_per_ms | collective_fraction |
| --- | --- | --- |
| N | 0.48031959 | 0.088525677 |
| T | 0.52395287 | 0.26953399 |


Para 42 derivadas (retraso uniforme y 20 ganancias por diseño), máximo error
relativo frente a diferencias finitas: 1.78236e-07; máximo
desacuerdo entre formulación completa y de puertos: 4.81583e-08.
La fracción colectiva anterior es la fracción absoluta de las dos contribuciones
a la sensibilidad, no MW/ms de reemplazo. Todos los componentes están en TABLE_08.
Esta parte fue añadida tras el hallazgo del límite del criterio de signo y está
marcada exploratoria. La identidad de sensibilidad es conocida; no se vende
como un nuevo teorema de optimización.

## Resultado negativo principal — NEGATIVE_RESULT

El signo del amortiguamiento armónico proyectado **no clasifica universalmente
la estabilidad**: N a 40 ms tiene un polo seguido con parte real positiva
y amortiguamiento proyectado positivo. La dinámica reactiva, los estados ocultos
y la geometría de los modos importan. Por tanto queda rechazada la sustitución
del problema PLL completo por una condición escalar de signo de D_H.

La reducción anterior en s=0 sigue bloqueada para su partición original; este
experimento no la repara retroactivamente. La presente reducción es dependiente
de frecuencia y solo válida donde el bloque eliminado es invertible. El máximo
condicionamiento observado en los bloques ocultos de las ramas es
2.686e+11: se reporta y no se presenta como certificación
en aritmética validada. Máximo residuo de raíz completo 3.06e-15;
máximo residuo de Schur 8.5e-09. Máximo error de
identidad transferencia/impedancia 1.2e-10.

## Qué falta para la contribución del póster

1. Un diseño predictivo que use la interacción colectiva y supere un comparador
   con las mismas variables, restricciones, tolerancias y presupuesto de cálculo.
2. Volver a optimizar rho junto con Kp/Ki y ejecutar los cinco eventos congelados
   con retraso no lineal; el mecanismo por sí mismo no levanta el límite de reserva.
3. Verificar heterogeneidad espacial y cualquier compresión GSP en datos reservados.
4. Una cota superior válida antes de decir máximo global o casi óptimo.

**Decisión: resultado de mecanismo útil; aún no justifica un póster que anuncie
máximo reemplazo óptimo mediante Beyond Nodal Damping.** El parentesco con
literatura de impedancias/par complejo y control basado en grafos se documenta
en LITERATURE_POSITION.md; no se afirma novedad por renombrar esas herramientas.
