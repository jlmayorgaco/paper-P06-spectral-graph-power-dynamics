# Teoría candidata: frontera de compensabilidad finita con estructura de red

**Estado:** formulación de investigación; no es todavía un resultado de novedad ni una garantía validada.

## 1. Problema de diseño

Considérese una DAE de red

\[
\dot x=f(x,z,\rho,\kappa,w),\qquad 0=g(x,z,\rho,\kappa,w),
\]

donde \(x\) agrupa estados dinámicos, \(z\) variables algebraicas, \(\rho_i\in[0,1]\) la fracción de reemplazo síncrono-a-GFL en el sitio \(i\), \(\kappa_i=(\log K_{p,i},\log K_{i,i})\), y \(w\) una perturbación admisible. El modelo puede incluir ramas, cargas, controles y modos de limitador, siempre que cada modo se represente con sus propias ecuaciones y región de validez.

Sea \(h_j(\rho,\kappa,w)\le 0\) una restricción de desempeño evaluada en equilibrio o trayectoria: polos, frecuencia, RoCoF, tensión u otras cantidades. El diseño robusto ideal sería

\[
\max_{\rho,\kappa}\; c^\top\rho
\quad\text{sujeto a}\quad
h_j(\rho,\kappa,w)\le0\quad\forall j,\;\forall w\in\mathcal W.
\]

En la prueba IEEE-39 de este informe se fijó \(\rho\) y se cambió \(\kappa\), para aislar el efecto causal del retuning en seis eventos. No se resolvió este problema robusto ni se optimizó globalmente \(\rho\).

## 2. Reducción DAE y sensibilidades

En una rama regular donde \(g_z\) es no singular, el teorema de la función implícita define localmente \(z=\phi(x,\eta)\), con \(\eta=(\rho,\kappa)\). La dinámica reducida es

\[
\dot x=F(x,\eta,w):=f(x,\phi(x,\eta,w),\eta,w).
\]

Para \(a\) una coordenada de diseño,

\[
g_z\phi_a+g_a=0,\qquad
F_a=f_a+f_z\phi_a.
\]

Las sensibilidades de trayectoria \(s_a=\partial x/\partial\eta_a\) obedecen

\[
\dot s_a=F_xs_a+F_a.
\]

Para un par \(a,b\), la sensibilidad de segundo orden \(S_{ab}=\partial^2x/(\partial\eta_a\partial\eta_b)\) satisface

\[
\dot S_{ab}=F_xS_{ab}
+F_{xx}[s_a,s_b]+F_{xa}s_b+F_{xb}s_a+F_{ab}.
\]

Las derivadas \(F_{ab}\) incluyen la curvatura de la variedad algebraica. En particular,

\[
g_z\phi_{ab}
=-D^2g\big[(e_a,\phi_a),(e_b,\phi_b)\big],
\]

con los términos de estado incluidos en \(e_a\) cuando \(a\) también representa una perturbación de \(x\). Así, un cálculo de curvatura que diferencie solo \(f\) y mantenga fijo \(z\) no corresponde a la DAE acoplada.

## 3. Candidato a certificado de exclusión finita

Apílense las restricciones relevantes en \(h(\eta)\le0\). Para cualquier vector dual \(\lambda\ge0\), defínase \(q(\eta)=\lambda^\top h(\eta)\). Todo punto factible debe cumplir \(q(\eta)\le0\).

En un recuadro de diseño \(\mathcal B=\{\eta_0+\delta:\delta\in\mathcal D\}\), supóngase una cota uniforme verificada

\[
\|\nabla^2q(\eta)\|_2\le L_\lambda
\quad\forall\eta\in\mathcal B.
\]

Entonces el resto de Taylor da, para todo \(\delta\in\mathcal D\),

\[
q(\eta_0+\delta)\ge
q(\eta_0)+\nabla q(\eta_0)^\top\delta
-\frac{L_\lambda}{2}\|\delta\|_2^2.
\]

Si el mínimo global de la derecha sobre \(\mathcal D\) es positivo, ningún punto de ese recuadro satisface todas las restricciones. Esto sería una exclusión finita de todos los ajustes PLL del recuadro para el nivel de reemplazo considerado. Para una caja, la parte lineal se minimiza mediante su función soporte; el término cuadrático conservador puede acotarse por el radio de la caja. Una caja logarítmica hace multiplicativamente simétricos los límites de \(K_p,K_i\).

La desigualdad de Taylor y la separación dual son herramientas conocidas. La contribución potencial no sería renombrarlas, sino construir una cota \(L_\lambda\) útil, rigurosa y escalable que respete la estructura de la red, la DAE, las métricas de trayectoria y los modos de limitador; luego demostrar que el certificado anticipa la frontera observada por integración no lineal.

## 4. Dónde entra “beyond nodal damping”

Escríbase el Jacobiano algebraico disperso como \(g_z=D-C\), donde \(D\) contiene bloques nodales y \(C\) acoplamientos entre nodos. Si se verifica \(\|D^{-1}C\|<1\) en una norma declarada, entonces

\[
g_z^{-1}=(I-D^{-1}C)^{-1}D^{-1}
=\sum_{k=0}^{\infty}(D^{-1}C)^kD^{-1}.
\]

Cada potencia suma caminos de longitud \(k\) en el grafo algebraico; los retornos contienen ciclos. La cola truncada tras \(K\) pasos tiene la cota

\[
\left\|\sum_{k=K+1}^{\infty}(D^{-1}C)^kD^{-1}\right\|
\le
\frac{\|D^{-1}C\|^{K+1}}{1-\|D^{-1}C\|}\|D^{-1}\|.
\]

Esto sugiere un mecanismo concreto: propagar sensibilidades y curvaturas por caminos locales, agrupar retornos por ciclos o bandas espectrales, y acotar la cola para producir la cota Hessiana. La construcción falla si la condición de contracción no se verifica; en redes fuertes o Jacobianos mal condicionados habría que usar una factorización/resolvente con otra cota certificada. La expansión no debe presentarse como identidad útil sin comprobar su convergencia numérica en el punto operativo y en toda la caja.

## 5. Condiciones matemáticas que aún faltan

1. **Caja definida y parametrización:** acotar simultáneamente \(\rho\), \(K_p,K_i\), estados iniciales y la región algebraica.
2. **Regularidad:** demostrar que \(g_z\) sigue no singular y que no se cruza una bifurcación dentro de cada región.
3. **Métricas no suaves:** frecuencia máxima sobre tiempo y nodos es un máximo de funciones y puede cambiar de rama activa. Hay que certificar ramas activas, usar una cota superior diferenciable o tratar el máximo mediante epígrafos; no aplicar Hessianas de una rama como si fueran globales.
4. **Híbridos y saturación:** dividir la caja por modo de limitador y validar transiciones, o incorporar explícitamente las superficies de cambio.
5. **Resto verificado:** obtener \(L_\lambda\) mediante aritmética de intervalos, cotas matriciales rigurosas u otro método con error controlado; diferencias finitas solas no certifican una cota.
6. **Falsación empírica:** comparar el certificado con una continuación de la frontera, búsquedas de factibilidad multistart y simulaciones PowerDynamics fuera de muestra. Un certificado conservador puede ser válido pero inútil si excluye demasiado poco.
7. **Alcance de perturbación:** reemplazar los seis escalones por una clase \(\mathcal W\) explícita antes de llamar robusto al resultado.

## 6. Qué mostró la primera prueba

En el IEEE-39 actual, para el mismo \(\rho\) (90.7854% GFL), las ganancias co-diseñadas redujeron los seis picos de frecuencia; la peor frecuencia bajó de 0.516675 a 0.499003 Hz y las violaciones del umbral exploratorio de 0.5 Hz pasaron de 2/6 a 0/6. El peor RoCoF aumentó de 0.242987 a 0.334563 Hz/s. Eso valida una compensación de desempeño finita en esta instancia, con tradeoff medible.

La envolvente lineal permite solo \(1.6964\times10^{-5}\) de aumento uniforme en \(\rho_i\) (0.001696 puntos porcentuales) dentro de la caja de ganancias indicada. La sensibilidad modal analítica coincide con la diferencia central fresca con error relativo \(2.6\times10^{-6}\); esto valida esa derivada local, no las derivadas de segundo orden. El valor de envolvente es tangente, no una cota física. No se probó que el retuning aumente el máximo factible de reemplazo, que la solución sea global, ni que pase límites de corriente o energía DC.

## 7. Hipótesis falsable para el proyecto

**Hipótesis:** una cota de curvatura estructurada por el grafo para restricciones del DAE puede producir certificados de no factibilidad finitos más estrechos que una cota nodal o una caja intervalar indiferenciada, y predecir el desplazamiento de la frontera de reemplazo causado por el retuning PLL.

**Criterio de éxito:** en un modelo reducido, certificar el error Hessiano y excluir un subconjunto no vacío de diseños PLL para niveles de \(\rho\) donde una búsqueda densa e independiente no encuentra factibilidad; después, en IEEE-39, comparar la frontera predicha y la medida con PowerDynamics, reportar brecha de capacidad, costo computacional y todos los casos adversos. Si la estructura de ciclos/resolvente no mejora el certificado sobre los controles nodales o densos, esa parte de la hipótesis queda refutada.

El nombre “beyond nodal damping” describe la motivación, no la novedad demostrada. La hipótesis solo se sostiene como contribución si la cota de curvatura realmente usa y mejora con información de caminos/ciclos del grafo.
