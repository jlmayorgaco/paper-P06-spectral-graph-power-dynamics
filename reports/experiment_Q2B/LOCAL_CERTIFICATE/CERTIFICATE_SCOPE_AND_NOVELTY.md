# Alcance de optimalidad y novedad — auditoría del 1 de octubre de 2026

## Problema que se intenta certificar

Se conserva el modelo físico de ExpN, su reparto P/Q, sus límites de PLL y el
evento de diseño de +100 MW en el bus 16. La frecuencia primaria es la diferencia
causal de fase con T=0.5 s en los diez buses generadores. Los eventos de buses 8 y
29 son validación fuera de muestra; no se añaden como restricciones ocultas.

Los archivos congelados y los resultados anteriores permanecen como evidencia.
Este directorio contiene una auditoría y un intento de corrección posterior.
La campaña `codesign_validation_20261001` se lee como antecedente: su candidato
de 366.238767 MW no está certificado y sus eventos adicionales no completaron
una validación no lineal satisfactoria. Además, aquella campaña optimizó tres
buses de perturbación; sus residuos KKT no corresponden al problema de bus 16
solamente que se estudia aquí.

## Garantía local

Para un soporte fijo, con funciones C2 y todas las desigualdades físicas
representadas como g(z)<=0, sea L(z,mu)=J(z)+mu' g(z). Una condición suficiente
de mínimo local estricto es:

1. Factibilidad completa: trim, límites de variables, todos los polos físicos,
   radio robusto en toda la banda, y extremos temporales de las salidas declaradas.
2. KKT: grad L=0, mu>=0 y mu_i*g_i=0.
3. Cualificación de restricciones; se audita LICQ y se identifican redundancias.
4. Positividad de d' Hess(L) d para todo d no nulo del cono crítico.

Con complementariedad estricta y LICQ, se puede usar una base Z del espacio
creado por las igualdades activas y exigir Z' Hess(L) Z positiva definida. Si hay
restricciones activas con multiplicador cero, comprobar únicamente ese espacio
nulo puede omitir direcciones: se necesita el cono crítico. Si hay polos,
frecuencias robustas o picos temporales empatados, se incluyen las ramas
correspondientes. En colisiones defectivas no se aplica el certificado suave.

Residuos pequeños y autovalores positivos en Float64 producen evidencia de
certificación numérica, no una prueba formal. Una garantía matemática validada
requiere acotar errores y encerrar la solución KKT y la positividad mediante
aritmética validada o argumentos analíticos equivalentes.

Un certificado en un soporte no demuestra optimalidad frente a inserciones de
SG ausentes. Para extender su alcance deben analizarse explícitamente los
soportes adyacentes y sus límites unilaterales, sin diferenciar a través de la
desaparición de estados. Tampoco un óptimo del problema linealizado es un
óptimo de un problema distinto con restricciones de trayectorias no lineales.
La validación PowerDynamics posterior verifica el candidato ensayado.

## Garantía global

Se necesitan L<=J_global<=U. U debe proceder de un candidato que satisfaga todas
las restricciones del problema, incluida robustez. L debe ser una cota válida
de una relajación resuelta globalmente o una cota que cubra todos los soportes y
subdominios continuos capaces de mejorar U. La enumeración de soportes seguida
de optimización local no basta. Una solución local de CORE tampoco es una cota
inferior global del problema con RoCoF. El umbral de brecha de proyecto es
U-L<=0.01 MW; no se ha demostrado aquí.

## Cambio numérico verificable

El evaluador nuevo reutiliza las matrices físicas de ExpN y deriva exactamente
la clausura de Schur dentro del soporte. Con Y=Gy^-1 C y L=Gy^-1 b_load:

    dY = Gy^-1 (dC-dGy Y)
    dL = -Gy^-1 dGy L
    dAred = dA-dB Y-B dY
    dBload = -dB L-B dL

Incluye las derivadas de entrada y salida de frecuencia y del salto algebraico
de fase en radianes. Para los picos usa la derivada de Fréchet de phi2(A,t),
contratada con los residuos, y el teorema de la envolvente en extremos simples.
Se contrastan con diferencias centrales y se registran los condicionamientos.

El radio puntual se evalúa mediante 1/||M^-1||_2 y sus vectores singulares,
para mejorar la precisión relativa cuando sigma_min(M) es diminuto frente a
||M||. Se conserva exactamente la métrica original. La búsqueda de unas pocas
frecuencias sigue siendo una estimación observada y requiere la auditoría de
toda la banda antes de declarar factibilidad robusta.

## Novedad: evaluación preliminar, no revisión exhaustiva

La utilización de SQP/KKT, sensibilidades de autovalores o certificados locales
no es novedosa por sí misma. Ya existen métodos SL/QP para optimización de la
abscisa espectral: Kungurtsev, Michiels y Diehl, *An Inequality Constrained
SL/QP Method for Minimizing the Spectral Abscissa* (2014),
https://arxiv.org/abs/1411.2362.

También existen planificación conjunta de recursos para estabilidad con
synchronous condensers y GFM: Chu y Teng, *Coordinated Planning for Stability
Enhancement in High IBR-Penetrated Systems* (2024),
https://arxiv.org/abs/2404.14012; y colocación de GFM usando propiedades de una
Laplaciana reducida: Yang et al., *Placing Grid-Forming Converters to Enhance
Small Signal Stability of PLL-Integrated Power Systems* (2020),
https://arxiv.org/abs/2007.03997.

Los precios sombra de estabilidad nodal también aparecen en trabajos recientes,
por ejemplo Wang y Geng, *Decentralized Stability-Constrained Optimal Power Flow
for Inverter-Based Power Systems* (2026), https://arxiv.org/abs/2604.17603.

La hipótesis de aporte del proyecto debe ser más específica: una clausura
paramétrica física para sustitución fraccionaria SG/GFL convencional, tratamiento
explícito de arquitecturas, identidades analíticas de baja dimensión y una
explicación de sensibilidad local/colectiva contrastada contra el sistema
completo. El co-diseño con seguridad y certificado local sería un resultado
que fortalece ese aporte. El código reproducible y la coincidencia con PD son
evidencia de validez; no demuestran por sí solos novedad en la literatura.

No se sostiene prioridad bibliográfica ni optimalidad con esta revisión
preliminar. La observación general de que estabilidad de pequeña señal no
garantiza seguridad de grandes perturbaciones tampoco debe presentarse como
descubrimiento nuevo.
