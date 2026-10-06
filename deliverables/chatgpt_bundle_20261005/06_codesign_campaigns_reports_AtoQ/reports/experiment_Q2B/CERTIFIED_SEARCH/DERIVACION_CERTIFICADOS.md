# Método ejecutado y alcance de los certificados

## 1. Problema finito y restricciones temporales

Dentro del soporte S={30,33,35,37,39}, se minimiza J/1000 con 25 variables normalizadas: cinco epsilon y veinte ganancias. Se mantienen todos los polos físicos y las funciones temporales exactas del modelo ExpN. No se optimiza con PowerDynamics.

Cada máximo temporal regular satisface f_t(z,t)=0. Si la segunda derivada temporal es no nula, el teorema de la función implícita define t(z) localmente. La derivada del valor extremo es sign(f) f_z: el término f_t t_z desaparece. Dos máximos con igual altura son **dos desigualdades**; no se diferencia su máximo como si fuera una función suave única.

La ejecución precedente se estancaba al intercambiar el pico ganador. `MultiPeak.jl` conserva ambos picos y añade otros detectados. SQP resuelve conjuntamente las direcciones epsilon/Kp/Ki. La corrección de segundo orden impone simultáneamente las superficies activas tras el predictor. La evaluación exacta decide aceptar cada paso.

## 2. Optimalidad local

Con g<=0, L=J/1000+mu'g. Se verifican factibilidad, mu>=0, mu'g=0 y grad L=0. Las cotas de ganancias forman parte de g. Se comprueba rango completo del Jacobiano activo (LICQ).

Todos los multiplicadores activos son estrictamente positivos. En este caso el cono crítico se reduce al espacio tangente Z=ker(J_activo). La condición suficiente de segundo orden es Z' Hess(L) Z > 0. La Hessiana se evalúa diferenciando gradientes analíticos con tres pasos y comparando su variación con la menor curvatura. Se incluyen las posiciones variables de los extremos al continuar sus ramas.

En aritmética exacta, estas condiciones dan un mínimo local estricto del problema suave de ese soporte. Aquí se suministra una verificación **numérica**: residuos, rango, condición y estabilidad frente al tamaño de paso. No se demuestra mediante intervalos la existencia de un cero exacto de KKT en una caja. El representante congelado tiene un pequeño resguardo factible; sus residuos propios se reportan aparte del punto KKT.

## 3. Robustez en toda frecuencia

Sea As=A+0.05I en las coordenadas físicas del cociente ortogonal congelado. La incertidumbre completa heredada usa M(s)=(sI-As)^(-1). Para b=beta_req, es suficiente encontrar X=X'>0 tal que

    As' X + X As + b (I + X^2) < 0.

Tomando P=X/b, esto equivale a

    As' P + P As + I + b^2 P^2 < 0.

El complemento de Schur da la desigualdad bounded-real para ||M||_infinity < 1/b. Cubre toda la recta de frecuencias, no sólo los puntos muestreados. Es una herramienta clásica, no una nueva teoría: [Boyd et al., sección 2.7.3](https://web.stanford.edu/~boyd/lmibook/lmibook.pdf).

Se inicia X desde el subespacio estable del Hamiltoniano y se refina mediante Newton/Lyapunov. La ecuación se resuelve con b_solve algo mayor que b; después se verifica la desigualdad estricta a b. La aritmética de los residuos y las pruebas de Cholesky es BigFloat de 256 bits. No se cambian la incertidumbre ni beta requerida. Se verifica nuevamente con A obtenida independientemente de PD.

El redondeo es al más cercano, no exterior. Por tanto la etiqueta es NUMERICALLY_CERTIFIED. El residual se compara con el margen estricto de la desigualdad; una mera casi-igualdad Float64 no sería aceptada.

## 4. Cobertura temporal

La medida congelada es una diferencia causal de fases con T=0.5 s; RoCoF usa la segunda diferencia. Incluye el salto algebraico inicial y los puntos donde cambia la ventana. Se mantiene idéntica en todas las arquitecturas.

Las trayectorias se evalúan con residuos modales. Para tiempo tardío se separa el valor estacionario de los exponentes, evitando la cancelación entre términos que crecen linealmente. En cada intervalo se usa

    |f(t)| <= |f(c)| + h |f'(c)| + h^2 sup |f''| / 2.

Los módulos de los residuos y las partes reales de los polos acotan f''. Se subdivide hasta cubrir el límite. La cola infinita se acota por la envolvente exponencial absoluta o por el signo del modo lento dominante. La cobertura final está documentada para las diez salidas F/R, con margen aritmético numérico declarado. No es una prueba formal intervalar.

## 5. Inserciones SG unilaterales

Para una inserción con epsilon->0+, se mantienen las ecuaciones internas nuevas y se toma exactamente cero su contribución de corriente en el **límite matemático**. No se acepta ese dispositivo como arquitectura física en epsilon=0.

El estado anterior es una proyección del ampliado. Tras los cocientes ortogonales, existe P con PP'=I y P A_nuevo=A_anterior P. Entonces

    R_anterior(s) = P R_nuevo(s) P'
    ||R_anterior(s)||_2 <= ||R_nuevo(s)||_2.

Si una inserción ya viola beta en omega=0, añadir más estados sin corriente en ese límite no restaura el requisito. Los cinco límites individuales violan estrictamente beta; por continuidad cada arquitectura adyacente tiene un entorno inviable alrededor de ese límite. Los residuos de coisometría, entrelazamiento y del inverso se entregan. Esto respalda la exclusión **local numérica**, no las inserciones finitas o las ramas distantes.

## 6. Valor marginal

Para cada epsilon interior, la estacionariedad física implica

    1 = sum_a (1000 mu_a) [- (1/P_i) partial g_a / partial epsilon_i].

El factor 1000 compensa la escala del objetivo. Se verificó esta suma en los cinco buses. Las contribuciones modales y de RoCoF son cero porque esas restricciones están holgadas; no se les atribuye una participación artificial. Una contribución robusta negativa es posible y no se trunca.

## 7. Brecha global

Todo despacho retenido es no negativo: L=0 es una cota inferior válida. El candidato factible da U. Ni el valor de los multiplicadores locales ni la curvatura local producen una cota global para este problema no convexo y de dimensión variable. El intervalo [0,U] sigue abierto.

La eliminación estática exacta se ensayó para producir una cota más fuerte mediante Bernstein. No produjo evidencia admisible y no se utiliza. Tampoco se usa una ley afín de frecuencia que el proyecto había rechazado. La comparación con el candidato anterior mide mejora experimental, no distancia al óptimo global.
