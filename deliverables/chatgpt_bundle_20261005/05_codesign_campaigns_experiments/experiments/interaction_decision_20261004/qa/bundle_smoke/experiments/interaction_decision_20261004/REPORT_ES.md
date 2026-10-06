# Resultado y reglas de diseño — 4 octubre2026

**Resultado cerrado:** hay un contraejemplo finito certificado en IEEE39 y una
regla iterativa que corrige su margen sin reducir los MW reemplazados.
La aportación es una decisión concreta; no se ha demostrado máxima penetración.

## Qué se demostró

- Sustitución adicional: 7.900000MW entre buses30 y37, con
  rho=.885 en esos dos y .875 en los otros ocho. Delay fijo40ms en los diez.
- Cada ajuste individual conserva el polo cerca de -0.065708927/s.
- Juntos: -0.041399606/s, que viola el margen exigido-.05/s.
- Interacción real: 0.024309321 ±4.01e-7 /s, incluida holgura de redondeo.
- Los cuatro modelos tienen exactamente una raíz en el mismo contorno,
  certificado por intervalos. Cada raíz se encierra en una caja de radio1e-7.
- Julia reproduce la linealización con error relativo máximo
  4.752e-18; su contorno de espectro completo encuentra0raíces
  a la derecha del margen para ambos sitios por separado, y2 para el conjunto.
  Este conteo global es numérico, no certificado por intervalos.

## Regla derivada y ejecutada

Separar exactamente la respuesta individual y la interacción:

    alpha_joint(t) = alpha_anchor - t + I(t)
    t_next = max(0, I(t) - b),  b = -alpha_anchor - sigma

Cada sitio recibe t/2 de amortiguamiento modal adicional; rho y Kp permanecen
fijos. I(t) se evalúa con el modelo DDE completo. El teorema de contracción
requiere una cota uniforme |I'|<1; aquí solo se comprobó empíricamente su pendiente.
El resultado final, en cambio, tiene su raíz encerrada rigurosamente.

Con sigma=.06/s, 11 actualizaciones llevan a
t=0.0170049249/s y alpha=-0.0600000000/s.
Ki30=272.4169693558; Ki37=288.7033226487.
Kp30=Kp37=28.8398205600. Todos los demás Ki=281.2837254310.
El vector completo se encuentra en TABLE_14_DESIGN_VECTORS.csv.

**Consejo de diseño sustentado:** no sumar certificados de polos individuales
para aceptar acciones simultáneas. Calcular la interacción y asignar un margen
adicional, o diseñar el eigenpar colectivo con un patrón compatible. Preservar
un polo complejo local tampoco garantiza composición universal; en este caso
fue una alternativa eficaz, con deriva de solo
8.957e-07/s.

## Frontera de una política, no máximo global

La política registrada de ajuste individual encuentra su primer cruce numérico
entre 6.384455 y 6.386041MW adicionales.
Los signos a ambos extremos están certificados. No se probó monotonía uniforme
ni exclusión de todas las demás ganancias. El diseño corregido conserva7.9MW.
La estimación cuadrática local predijo h=0.009412,
frente al cruce cercano a h=0.008083; no debe
usarse como cota sin controlar su resto.

## Validación no lineal

El diseño corregido pasa5/5eventos congelados con método de pasos DDE en Julia:
Fmax=0.46587397Hz; RoCoFmax=0.21906007Hz/s;
V en[0.97204369,1.07877556];
holgura mínimaSG=0.00891844. Frecuencia y RoCoF usan ventana.5s.
GFL=4735.315954MW (87.64622153%);
SG retenida=667.445136MW.
Es un diseño factible para ese contrato numérico, no un óptimo.

## Qué falló o sigue abierto

- El primer barrido no dio falsas aceptaciones válidas; hubo5falsos rechazos.
-230combinaciones con conservación del polo complejo no violaron el margen.
- En la tercera campaña hubo88casos resueltos y317no resueltos;2violaciones.
- El criterio de convergencia de la serie por recorridos falla en el contorno
  del caso seleccionado. Se usa la identidad racional exacta y certificados
  directos de las raíces, no una cota de serie inválida.
- El predictor cuadrático ordinario también detecta este caso. La corrección
  de dos pasos mejora la precisión, pero no tiene exclusividad sobre la decisión.
- Sin máximo global, cota superior de reemplazo, garantía sobre todo el rango
  de ganancias, comparación contra optimizadores publicados o incertidumbre
  física robusta. No se modela seguridad de limitador de corriente.

## Juicio editorial

Esto cierra una brecha real respecto de la versión anterior: ahora hay una
decisión que cambia, un certificado y una corrección ejecutada. Es material
para una contribución enfocada, pero no basta para anunciar un journal fuerte
cerrado o un premio. La prioridad bibliográfica no está probada; el siguiente
paso editorial requiere validación independiente más amplia y comparación
directa contra métodos existentes, conservando exactamente este contrato.
