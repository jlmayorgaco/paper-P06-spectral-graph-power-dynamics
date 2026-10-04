# Leyes de red para el ajuste de PLL — resultado del experimento

## Resultado central

Para este modelo SG/GFL, manteniendo planta, rho y Ki, los cambios finitos de
Kp y retardo cumplen exactamente una ley entre coeficientes:

    H_b(s)-H_a(s) = -s² kappa H0 + O(s³)
    kappa = sum_i w_i (delta_tau_i/Ki_i - delta_Kp_i/Ki_i²)

Los pesos w salen de la singularidad rotacional del retorno de red completo;
no se eligieron como centralidades. Tienen signos distintos. Con estabilidad,
el área firmada de la diferencia de frecuencias es cero y su primer momento
temporal es kappa*f_infinity. No es una cota de nadir, RoCoF ni seguridad.
Cambiar Ki mueve todas las áreas normalizadas en una misma dirección:
sus diferencias entre buses son invariantes dentro de esta arquitectura.

La derivación para control matricial muestra además que
delta_KP=KI*S*KI, S*1=0, preserva H0,H1,H2. Un polinomio de Laplaciano sin
término constante cumple esa condición. Es un resultado teórico condicionado;
ese controlador con comunicación no fue implementado ni comparado.

## Qué se ejecutó en IEEE-39

204 estados, diez PLL filtrados, retardos exponenciales puros en el detector,
red AC con pérdidas, 39 salidas, tres entradas linealizadas. Las ecuaciones
se comprobaron con recurrencia de Laurent, transferencia completa204D y
derivadas por diferencias finitas. Los signos de los pesos y las hipótesis de
regularidad se encerraron con aritmética de bolas de 192 bits para la exportación
binaria con la restauración rotacional explícita. No cubre incertidumbre física.
Los 14 diseños colectivos redondeados satisfacen |kappa|<1e-12 mediante encierro
intervalar; no se atribuye un cero literalmente exacto a los gains numéricos.

Conservamos **4735.315954 MW GFL (87.646222%)** y
**667.445136 MW SG**. No se optimizó ni aumentó ese porcentaje.
Solo cambió Kp; Ki y rho están congelados en los valores previos.

Al pasar de 40 a 41 ms:

| Método | Parte real crítica del catálogo (1/s) | Norma de cambio relativo Kp |
|---|---:|---:|
| Sin ajuste | 0.413951581 | 0.0000000 |
| Cancelación PLL por PLL | 0.490669709 | 0.0313241 |
| Solo momento, mínimo cuadrático | 0.405520336 | 0.0110501 |
| Solo margen modal | -0.060001359 | 0.0910625 |
| Momento colectivo + margen modal | -0.058572898 | 0.0858146 |

El último estado del algoritmo es **MAX_ITERATIONS**: no equivale a optimalidad
global ni, si agotó iteraciones, a convergencia local. El coste cuadrático frente
al otro candidato local modal cambia -11.194%; es una comparación de dos
resultados del algoritmo, no el precio óptimo probado de la restricción.
Además, el candidato colectivo no alcanza el objetivo auxiliar de -0.0600001/s,
mientras que el modal directo sí. Ambos se evalúan con el límite original -0.05/s;
sus costes no permiten inferir superioridad a igual margen modal alcanzado.

Espectro completo: **PASS_NUMERICAL**, conteo numérico por
contorno con cota de región y exponenciales exactas; no certificado intervalar.
Cinco eventos no lineales de 60 s: **5/5 pasan**.
Máximo |delta f|=0.465928588Hz; RoCoF=0.219224385Hz/s;
tensión mínima=0.972043689, máxima=1.078775559;
holgura SG mínima=0.008899887. Ventana heredada de 0.5 s.
Los eventos son cambios de carga Z parametrizados en MW, no potencia constante
exacta durante toda la trayectoria. No hay evidencia de seguridad de limitador
de corriente de convertidor: no está modelado.

## Resultados negativos que se conservan

- Igualar momentos de baja frecuencia no conserva el amortiguamiento modal:
  la regla por PLL hace inestable el caso41ms aunque kappa sea cero.
- Un solo modo protegido falló; el catálogo completo es necesario aquí.
- El SQP sin globalización falló a 41 ms; sus resultados permanecen sin cambios.
- La ley es asintótica. A s=.0025/s, el error uniforme del primer término sigue
  siendo 6.62%.
- Las pantallas SDP y la cota causal del experimento hermano no proporcionaron
  una barrera continua de todos los gains ni un piso SG útil.

## Novedad y publicación

La identidad específica y su restricción de autoridad son candidatas a una
aportación útil. Residuos, momentos, cancelación de consenso y QP son clásicos.
No se encontró una coincidencia exacta en la búsqueda acotada, pero eso no
prueba prioridad. Las revisiones independientes están en REVIEW_*.md.

**No se declara listo para journal ni ganador de póster.** Falta demostrar un
beneficio operativo de preservar este momento frente al ajuste modal directo,
transferir el resultado a otro punto/red, y desarrollar/validar el controlador
de grafo si se quiere vender esa arquitectura. La máxima sustitución SG→GFL,
su cota superior y la ventaja global frente a otros métodos siguen abiertas.

## Reproducción y entregables

THEORY_MOMENT_SECTION.tex contiene pruebas; el documento THEORY.tex abierto
integra la sección y los datos generados. TABLE_09 y globalized/TABLE_14 guardan
rho,Kp,Ki,tau completos. TABLE_11 contiene el espectro independiente y TABLE_12
los cinco eventos. Los protocolos, fallos, fuentes, versiones, hashes y revisión
se preservan junto al código. REPRODUCE.md da los comandos.
