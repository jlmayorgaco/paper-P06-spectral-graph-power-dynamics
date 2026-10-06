# Extensión propuesta: compensación de retardo con dos Kp y todas las Ki fijas

## Estado
Derivación nueva de esta revisión y prueba numérica exploratoria sobre el banco SG-GFL de nueve buses previamente entregado. No es una validación IEEE-39, no hay trayectorias no lineales nuevas y no se ha establecido prioridad bibliográfica.

## Contrato
Planta, equilibrio, reparto SG/GFL y arquitectura fijados. El detector de cada PLL es compartido por sus inyecciones proporcional e integral; las ganancias no entran en otras ecuaciones. Los cambios de delay son escenarios prescritos, no un controlador que altera físicamente el retardo. Se conserva un polo complejo simple; el patrón modal puede cambiar. El control durante la operación permanece local y no se añade comunicación entre convertidores.

## 1. Límite de la compensación en un solo PLL
Sea lambda=alpha+j*omega, omega>0, y a_i=d lambda/d kp_i no nulo. La identidad del PLL da
  d lambda/d tau_i=-(lambda*kp_i+ki_i)*a_i.
Para compensar tau_i -> tau_i+h cambiando solo kp_i y manteniendo ki_i fija sería necesario
  d kp_i/dh=lambda*kp_i+ki_i.
Para kp_i>0 esto no es real. Por tanto no existe una rama diferenciable local de ese ajuste que conserve el polo completo. No excluye mantener solo su parte real ni modificar otros controles.

## 2. Cooperación entre sitios
Si solo se modifican kp_i,kp_j, la tangente que conserva el polo resuelve
  J_ij [kp_i',kp_j']^T = -[Re(lambda_tau*r),Im(lambda_tau*r)]^T,
  J_ij=[[Re a_i,Re a_j],[Im a_i,Im a_j]].
Si Im(conj(a_i)*a_j)!=0, la función implícita proporciona una rama local exacta de ganancias reales para pequeños cambios de retardo. También se requieren raíz simple, modelo suave y límites de ganancia interiores. La condición no garantiza un intervalo finito grande ni el resto del espectro. En caso singular puede existir solución para una dirección particular; la singularidad no prueba imposibilidad global.
Las a_i contienen los residuos de toda la red. Para medir condicionamiento/coste deben normalizarse los cambios de ganancia.

## 3. Solución finita por una cuadrática real
En una frecuencia compleja lambda fija y delays nuevos prescritos, escribir x=Delta kp_i, y=Delta kp_j. Por la estructura de rango uno de cada actualización, el determinante completo tiene exactamente la forma
  f(x,y)=c0+c1*x+c2*y+c12*x*y.
Si c2+c12*x!=0:
  y=-(c0+c1*x)/(c2+c12*x).
La condición y real es una cuadrática real:
  Im(c1*conj(c12))*x^2
  + Im(c0*conj(c12)+c1*conj(c2))*x
  + Im(c0*conj(c2))=0.
Deben conservarse todas las raíces reales, tratar las ramas degeneradas y filtrar por límites físicos. Cada candidata requiere verificación de su determinante completo, la raíz buscada y el espectro. No basta la cuadrática.

Si el operador de referencia a lambda es invertible puede utilizarse R=C_S Delta_ref(lambda)^(-1)P_S, con el signo incluido en P. Entonces c0=1,c1=R11,c2=R22,c12=det R. Si la referencia es singular, este cociente no puede utilizarse: escoger otra referencia regular o formar los coeficientes del determinante completo sin dividir por cero. No usar el cociente en el polo de referencia a h=0.

## 4. Relación con los momentos
En las hipótesis del apéndice J de Network Moment Laws (bloque oculto regular cerca de cero, gauge simple, d0 no nulo, planta fija), mantener todas las Ki fijas conserva H0 y H1 de la transferencia física de frecuencia, aunque cambien Kp y delay. Así H_new-H_old=O(s^2). Cuando ambas respuestas son estables y los momentos convergen, el área firmada de la diferencia de respuestas a escalón es cero. No es igualdad de nadir/RoCoF, ni preservación de toda la respuesta de baja frecuencia, ni prueba de seguridad no lineal.
Se comprobaron en punto flotante las hipótesis algebraicas del germen del banco: ver moment_hypothesis_checks.json. No es un certificado intervalar ni una integral temporal calculada.

## 5. Ejemplo ejecutado
Referencia: rho=0.75, kp_i=44.42840330706686, ki_i=986.9604401089358, tau_i=0.020 s.
Se protege lambda=-0.52786372+j49.58027295 aproximadamente; el cálculo usa la raíz refinada, no estos decimales redondeados.
El delay del detector del PLL 1 sube a 0.021 s; los otros quedan en 0.020 s.

Sin ajuste: alpha=+0.2753720376, dos raíces inestables según conteos numéricos.
Reparación por dos Kp: kp=(44.4284033071,49.3363371092,34.0084500074); todas las Ki siguen iguales a 986.9604401089.
Alpha de catálogo refinado=-0.1282547935. Error en la raíz protegida=8.26e-15. Conteos numéricos refinados: cero raíces a la derecha de -0.1/s y cero raíces inestables.
El PLL cuyo delay aumentó permanece sin retuning: actúan los Kp de los sitios 2 y 3, mediante acoplamiento eléctrico, sin comunicación nueva.
Comparador de transporte local: kp_1=45.3129956606, ki_1=876.6238310208; alpha=-0.1273958870, también pasa los conteos. La nueva reparación NO demuestra menor esfuerzo ni mejor desempeño temporal que este comparador.

## 6. Qué se ejecutó
- probe.py: varias raíces de la referencia, todos los pares, h=0.01/0.1/1 ms, escenarios comunes y locales. Se conservan las filas fallidas y soluciones inadmisibles.
- delay_scan.py: 3 sitios de delay x 10 incrementos x 3 pares; 90 problemas bilineales, hasta dos candidatos regulares por problema. Devolvió 74 registros con ganancias dentro del intervalo [0.25,2] del valor original. El ejemplo reportado fue seleccionado exploratoriamente, no es un test ciego.
- check_remote_repair.py: comparación del ejemplo exitoso, órdenes Padé 6/8/10 solo para semillas, refinamiento con exponentiales exactas, contornos de densidad 1/2 con cota de cola del código previo. No intervalos; no ejecución temporal.
- check_selected.py: conserva otro ejemplo donde la perturbación no requería retuning para ser estable, para no ocultarlo.

El catálogo por sí solo no prueba espectro completo. Los conteos del ejemplo principal usan el procedimiento del paquete previo; no se presenta como implementación independiente del contador.
La exploración incluye candidatos inestables y escenarios sin candidato admisible; el método no es universal.

## 7. Cierre necesario antes de claim principal
Repetir en IEEE-39 con el mismo lambda, reparto, límite de ganancias y delays que el comparador. Evaluar todos los pares: no escoger después solo el ganador. Pasar N_0 y N_sigma, todas las raíces y eventos originales, refinamiento temporal y muestras fuera del punto base. Comparar contra transporte individual, su Taylor correcta y ajuste de Kp con optimizador completo. Reportar también coste normalizado. No convertir el criterio de rango ni el área firmada en garantía de nadir, robustez o máximo delay.

## Fuentes de la idea y frontera de novedad
El manuscrito de proyecto Network Moment Laws ya contiene transporte, derivadas PLL y conservación de momentos a Ki fija, además de un ajuste colectivo linealizado. Esta extensión propone su unión con una selección de dos ganancias proporcionales por diversidad de fase y una solución finita bilineal exacta. La álgebra de la cuadrática, la función implícita y la asignación de control son herramientas conocidas.
Antecedentes cercanos: Wilson/Nelson/Farhang-Boroujeny (2009), Parameter Derivation of Type-2 Discrete-Time PLLs Containing Feedback Delays; Huang et al., arXiv:1903.05489; Ionescu/Iftime/Stefan (ECC 2023), A moment matching-based loop shaping design with closed-loop pole placement; literatura de control allocation.
No se afirma que nadie haya publicado una formulación equivalente. Hace falta auditoría por arquitectura, objetivo y condiciones, no solo búsqueda de palabras.

## Reproducción
Python, NumPy, SciPy y pandas. Ejecutar desde cualquier directorio:
  OPENBLAS_NUM_THREADS=1 python probe.py
  OPENBLAS_NUM_THREADS=1 python delay_scan.py
  OPENBLAS_NUM_THREADS=1 python check_remote_repair.py
Los scripts resuelven nuevamente el caso y guardan nuevos resultados en este directorio. Conservar copias antes de reejecutar si se necesita comparar las tablas entregadas.

Versiones de esta ejecución: {'numpy': '2.3.5', 'scipy': '1.17.0', 'pandas': '2.2.3'}
