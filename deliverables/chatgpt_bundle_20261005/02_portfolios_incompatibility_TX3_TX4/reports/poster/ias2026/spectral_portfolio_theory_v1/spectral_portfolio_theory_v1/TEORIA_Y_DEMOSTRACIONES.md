---
title: "Certificación binaria de portfolios SG–IBR"
subtitle: "Desarrollo matemático, límites y protocolo para validación independiente"
author: "Nota de investigación para el proyecto Vancouver IAS Annual Meeting Poster"
date: "10 de septiembre de 2026"
lang: es
---

# 1. Alcance y estatus de esta nota

Esta nota desarrolla una extensión del programa existente: pasar de **describir cuáles combinaciones fallan** a **certificar cuáles familias no pueden fallar, encontrar las combinaciones mínimas peligrosas con poda válida y diseñar decisiones robustas**. Se separan cuatro tipos de material:

- **D — Documentado:** resultados que constan en los informes suministrados. No equivalen a una nueva ejecución independiente en esta sesión.
- **T — Derivado:** proposiciones con demostración incluida. Su corrección no implica que la proposición sea nueva en la literatura.
- **V — Verificado aquí:** ejemplos sintéticos ejecutados en `theory_checks.py`; no son experimentos IEEE-39 ni ANDES.
- **E — Por validar:** utilidad, conservadurismo, escalabilidad y alcance físico que debe medir Claude Code sobre los modelos congelados.

La documentación disponible incluye `TPWRS_FINAL_EVIDENCE_TABLE.md`, `G5_NOVELTY_STATEMENT.md`, G1/G2/G3 y los informes F1/F7/F8/F10/F11/F12. Se leyeron directamente del archivo `docs.zip`; no se dispone aquí del repositorio completo con las implementaciones ANDES y los datos de todos los experimentos. Ningún resultado nuevo de esta nota debe incorporarse al ledger como validación de una red eléctrica hasta ejecutar su experimento correspondiente.

## 1.1 Lo que se conserva

Según la tabla final, el proyecto ya tiene mapas de incompatibilidad dependientes de política en IEEE-39 y Kundur; un resultado negativo para cuatro candidatos de IEEE-68; y un resultado secundario de alta penetración sobre doce plantas con órdenes mínimos 6, 9 y 11 en tres políticas. La reducción de puertos coincide con el modelo interno; la reconciliación ANDES se limita a la configuración de inyecciones estáticas, PSS desactivado y equivalencia de banda electromecánica. No valida por sí sola el GFL completo ni todos los modos rápidos.

La novedad declarada en G5 es la **estructura mínima de incompatibilidad dependiente de política**, no Schur, Nyquist, el punto −1 ni los hipergrafos. Esta disciplina se mantiene.

## 1.2 Cuatro correcciones prioritarias

**Origen.** La tabla final llama RHP a

$$\{s:\Re s>0,\ |s|>10^{-3}\}.$$

Eso es un semiplano con un disco excluido, no todo el semiplano derecho. Un polo físico $+10^{-5}$ puede quedar oculto. Eliminar una simetría angular conocida es legítimo; eliminar todo lo pequeño no lo es.

**Métrica.** La distancia $\min_{\mu\in\sigma(Q)}|1+\mu|$ es invariante bajo similitud, pero no es un radio general de robustez ante perturbaciones. Tampoco determina por sí sola de qué lado de una frontera se encuentra el sistema.

**Unidades.** `PHASE_E_REPORT.md` utiliza 4270.7 como MW de PV. En la tabla de generación comunicada durante el proyecto, 1040, 1174.8, 1085.7 y 970.2 son ratings $S_n$ en MVA y suman precisamente 4270.7. Esto es una alerta de trazabilidad, no una corrección numérica definitiva sin los archivos de origen. El protocolo exige reconstruir por separado MW inyectados, MW de capacidad activa, MVA de convertidor y MVA de condensador.

**No monotonicidad.** Encontrar un subconjunto inestable no permite declarar inestable a todo superconjunto. Sí permite excluir sus superconjuntos de la búsqueda de *testigos mínimos*, que es otra tarea.

# 2. Literatura: dónde termina lo conocido

La revisión se centra en fuentes primarias. Las referencias completas están al final.

| Referencia | Aporte existente relevante | Consecuencia para este proyecto |
|---|---|---|
| Katewa y Pasqualetti [L1] | Radio de estabilidad real con patrón de dispersión prescrito | No presentar la perturbación estructurada mínima como concepto nuevo; distinguir activaciones físicas binarias y cardinalidad de amplitud continua. |
| Huang et al. [L2] | Certificados descentralizados de pequeña ganancia y fase para redes multiconvertidor | Un certificado de ganancia por sí solo no basta como novedad. Debe servir para certificar subfamilias de decisiones y medir su utilidad. |
| Hallinan y Lestas, partes I y II [L3–L4] | Certificados plug-and-play, representaciones alternativas y tratamiento de subsistemas inestables | Es imprescindible no suponer que toda admitancia aislada es estable; comparar cualquier certificado nuevo con esta línea. |
| Zhu et al. [L5] | Participación por impedancias, sensibilidad y trazado de causas con modelos de caja gris | La localización de dispositivo/parámetro y el retuning ya tienen antecedentes fuertes. |
| Rodríguez-Ortega et al. [L6] | Equivalencia entre espacio de estados y Nyquist incorporando dinámica de frecuencia | Reproducir autovalores con un retorno de puertos es verificación metodológica, no por sí sola novedad. |
| Boyd et al. [L7] y Parrilo [L8] | LMI, Lyapunov y certificados polinómicos/semialgebraicos | Sustentan las dos familias de certificados propuestas; no se reclama prioridad de las herramientas. |
| Yazdani y Lotfifard [L9] | Discusión reciente de limitaciones de usar Nyquist nominal como análisis robusto MIMO | Obliga a separar identidad nominal, incertidumbre física y certificación robusta. Es un preprint, no una autoridad que invalide Nyquist. |

Los preprints [L3–L4] fueron depositados el 7 de septiembre de 2026. Hacen especialmente débil un claim genérico de «nuevo certificado descentralizado para IBR». La oportunidad más específica es **certificar conjuntos de opciones discretas físicamente definidas sin exigir propiedades de configuraciones fraccionarias que no se van a implementar**.

No se encontró en esta revisión acotada una coincidencia exacta con toda la combinación propuesta. Eso no demuestra prioridad. Antes del envío deben contrastarse también las búsquedas “finite-mode stability”, “Boolean uncertainty”, “minimal destabilizing sets”, “noncoherent fault trees”, “parameter-dependent Lyapunov” y “cardinality-constrained robust control”.

# 3. Modelo físico y tres niveles de decisión

Sea $\mathcal V=\{1,\ldots,m\}$ el conjunto **finito y congelado** de proyectos candidatos. La decisión

$$\delta\in\{0,1\}^{m},\qquad S=\{i:\delta_i=1\}$$

elige qué generadores se reemplazan. $\theta$ recoge políticas y parámetros controlables; $\xi$ recoge incertidumbre. En un modo de control/limitador fijo, la red se describe por

$$\dot x=f(x,z;\delta,\theta,\xi),\qquad 0=g(x,z;\delta,\theta,\xi). \tag{1}$$

Un modelo descriptor $E\dot x=F$ también es válido; las expresiones siguientes se aplican después de una reducción índice uno correctamente justificada. Se exige equilibrio, regularidad, restricciones de capacidad y un régimen de limitadores explícito.

## 3.1 El grafo físico permanece explícito

Con una incidencia de red $B$, su ampliación bifásica es $\mathcal B=B\otimes I_2$. En el caso sin taps, una representación es

$$Y_{\rm net}(s)=\mathcal B\,\operatorname{blkdiag}(Y_e(s))\,\mathcal B^T+Y_{\rm sh}(s). \tag{2}$$

Los transformadores requieren operadores de conexión con tap y rotación, no reemplazarlos a ciegas por una incidencia ±1. Las admitancias incrementales de carga PQ se incorporan en la linealización. Se define $Y_i$ con corriente *absorbida* positiva para fijar el signo de KCL.

El hecho de que la red sea un grafo no hace que su inversa sea dispersa. Los acoplamientos de puertos pueden ser densos. Umbralizar entradas pequeñas es una aproximación y necesita una cota del error descartado.

## 3.2 Tres objetos diferentes

1. **Portfolio estático:** elegir $S$ y comprobar el equilibrio final.
2. **Familia de implementaciones:** garantizar seguridad para múltiples subconjuntos, por ejemplo todos los que puedan quedar operativos durante una transición.
3. **Trayectoria híbrida real:** ejecutar cambios con estados iniciales, rampas, conmutaciones, limitadores y protecciones.

Un certificado del segundo objeto no certifica automáticamente el tercero. La Sección 15 introduce una extensión local que respeta esa diferencia.

## 3.3 No confundir dos afinidades

F11 verificó ensamblaje afín en ciertos **parámetros de control para cada $S$ fijo**. Eso no prueba

$$A(\delta)=A_0+\sum_i\delta_i A_i. \tag{3}$$

La retirada SG y la incorporación GFL pueden cambiar la dimensión y el equilibrio. Todo certificado que use (3) tiene una puerta de aplicabilidad: construir una realización común exacta o trabajar sobre un modelo de puertos/pencil que conserve los términos restantes. Añadir rellenos estables puede igualar dimensiones en los vértices, pero no demuestra afinidad ni localidad por sí solo.

# 4. Estabilidad física y eliminación correcta de la referencia

Al resolver las restricciones algebraicas se obtiene

$$A=f_x-f_zg_z^{-1}g_x. \tag{4}$$

La inversión se implementa como solución lineal. La condición de $g_z$ depende de unidades; se reportan también escalado, valor singular mínimo y error hacia atrás, no sólo un umbral de condición.

## Teorema T1 — Cociente de una simetría declarada

Sea $R\in\mathbb R^{n\times r_g}$ de rango completo un conjunto de generadores de simetría **identificado a partir de las ecuaciones físicas**, con

$$AR=0.$$

Sea $U$ una base ortonormal de $\operatorname{im}R$ y $Z$ una base ortonormal de su complemento. Entonces

$$[U\ Z]^TA[U\ Z]=\begin{bmatrix}0&B_g\\0&A_q\end{bmatrix},\qquad A_q=Z^TAZ, \tag{5}$$

por lo que

$$\det(sI-A)=s^{r_g}\det(sI-A_q). \tag{6}$$

**Demostración.** $AU=0$ anula la primera columna de bloques. El determinante de una matriz triangular de bloques es el producto de sus determinantes diagonales. La variable cociente $y=Z^Tx$ satisface $\dot y=A_qy$ independientemente del representante angular. $\square$

**Lo que NO dice.** No se pueden eliminar otros autovalores en cero. Si el modelo tiene una cadena de Jordan de longitud dos y sólo una simetría, el cociente deja un cero físico. La ausencia de gobernador/amortiguamiento puede producir esa situación. Tampoco se pueden eliminar integradores muertos de un regulador llamándolos simetrías.

### Corolario T1a — Reubicación numérica de la simetría

Para $\beta>0$,

$$A^{\sharp}=A-\beta UU^T$$

satisface

$$\det(sI-A^{\sharp})=(s+\beta)^{r_g}\det(sI-A_q). \tag{7}$$

Sólo se desplaza el subespacio de simetría conocido. No se está añadiendo una barra infinita ni amortiguamiento físico. Esta identidad puede facilitar resolventes en $s=0$ **si las actualizaciones preservan el mismo subespacio de simetría**.

### Ejemplo verificado V1

$$A(\varepsilon)=\begin{bmatrix}0&1&0\\0&\varepsilon&0\\0&0&-2\end{bmatrix},\qquad R=e_1.$$

El cociente es $\operatorname{diag}(\varepsilon,-2)$. Para $\varepsilon=10^{-7}$ el sistema relativo es inestable, aunque el filtro $|\lambda|>10^{-3}$ lo ocultaría. Para $\varepsilon=0$ queda un cero físico; el proyector basado en emparejar autovectores derechos/izquierdos del cero puede fallar porque la matriz es defectiva. (5) no necesita ese emparejamiento.

## 4.1 Definiciones de seguridad

Se usa

$$\alpha_q(S;\theta,\xi)=\max\Re\sigma(A_q(S;\theta,\xi)).$$

Un estado es estrictamente estable si $\alpha_q<0$. En una zona regular sin autovalores sobre el eje imaginario, el conteo

$$N_+(S)=\#\{\lambda\in\sigma(A_q(S)):\Re\lambda>0\} \tag{8}$$

certifica estabilidad cuando es cero. **En la frontera, $N_+=0$ no implica estabilidad asintótica.** Numéricamente se usan `STABLE`, `UNSTABLE`, `BOUNDARY_OR_UNRESOLVED`, `INFEASIBLE`.

Para diseño con tasa de decaimiento $\gamma>0$, se aplica lo mismo a $A_q+\gamma I$. Todos los conteos son con multiplicidad algebraica y por autovalor: un par conjugado aporta dos. Los pares y polos reales se pueden reportar por separado, sin mezclarlos en un único entero ambiguo.

# 5. Objeto combinatorio: mínimo, no «causa irreducible»

Sobre una familia de candidatos factible y completa, y fuera de la frontera,

$$\kappa(\theta)=\min\{|S|:N_+(S;\theta)>0\},\quad\min\varnothing=\infty. \tag{9}$$

$$\mathcal H(\theta)=\min_{\subseteq}\{S:N_+(S;\theta)>0\}. \tag{10}$$

Si hay portfolios inviables, se conservan como una categoría distinta. La familia de restricciones de factibilidad debe formar parte del problema; no se cuentan inviables como estables o inestables.

Un conteo mínimo no es una demostración de interacción de orden alto en una función escalar. La descomposición Möbius

$$\mu_f(S)=\sum_{R\subseteq S}(-1)^{|S|-|R|}f(R) \tag{11}$$

es exacta para cualquier función de conjuntos $f$. Si $f=\max\Re\lambda$, la maximización entre modos puede generar términos altos aunque las ramas subyacentes tengan dependencia sencilla. El cambio a envolventes evita falsos negativos de tracking, pero **no transforma automáticamente los coeficientes Möbius de esa envolvente en mecanismos físicos de cuarto orden**.

## Proposición T2 — No-go de orden fijo, con interacción aditiva

Para cualquier $m\ge2$, el sistema escalar

$$\dot x=\left[-(m-\tfrac12)+\sum_{i=1}^m\delta_i\right]x \tag{12}$$

es estable en todo subconjunto propio y es inestable con las $m$ acciones. Así $\kappa=m$, aunque todos los coeficientes Möbius de orden mayor que uno de la abscisa son cero.

**Demostración.** Si $|S|\le m-1$, el autovalor es como máximo $-1/2$. Si $|S|=m$, vale $+1/2$. La expresión es afín en los indicadores. $\square$

Esto es un contraejemplo elemental y no se presenta como un teorema fundamental nuevo. Su utilidad es impedir la confusión entre cardinalidad mínima y no aditividad genuina.

# 6. Representación de puertos y contabilidad de polos

En equilibrio común, cada opción de dispositivo tiene realización incremental

$$\dot x_i=A_i x_i+B_i v_i,\qquad i_i=C_i x_i+D_i v_i,$$

$$Y_i(s)=D_i+C_i(sI-A_i)^{-1}B_i. \tag{13}$$

Se reúnen los dispositivos seleccionados en $A_d,B_d,C_d,D_d$ y se toma

$$\mathcal P_S(s)=\begin{bmatrix}sI-A_{d,S}&-B_{d,S}\\C_{d,S}&Y_{\rm net}+D_{d,S}\end{bmatrix}.$$

Su complemento de Schur es

$$T_S(s)=Y_{\rm net}+D_{d,S}+C_{d,S}(sI-A_{d,S})^{-1}B_{d,S}. \tag{14}$$

Luego

$$\det\mathcal P_S(s)=h_S(s)\det T_S(s),\qquad h_S(s)=\det(sI-A_{d,S}). \tag{15}$$

Esta igualdad es meromorfa antes de cancelar factores y analítica para el determinante completo. Un polo de $T_S$ puede cancelarse con un cero de $h_S$; un modo interno puede no verse en la transferencia. Por ello, alrededor de un contorno $\Gamma$ admisible,

$$N_{\rm full}(S)=N_{h_S}+\operatorname{wind}_\Gamma\det T_S. \tag{16}$$

Y para comparar portfolios,

$$N_{\rm full}(S)-N_{\rm full}(0)=\operatorname{wind}_\Gamma\frac{\det T_S}{\det T_0}+N_{h_S}-N_{h_0}. \tag{17}$$

Se presupone una orientación antihoraria del contorno; las convenciones horarias de Nyquist requieren invertir el signo correspondiente. Se incluyen los factores de simetría retirados y las constantes de eliminación algebraica en la identidad completa.

**Punto crítico:** la razón de determinantes de puerto no siempre equivale por sí sola a la razón de polinomios característicos cuando cambia el tipo/número de estados de dispositivos. Los términos de (17) son obligatorios.

## 6.1 Acciones localizadas

Con matrices de inserción $E_i$ y políticas congeladas,

$$T_S=T_0+\sum_{i\in S}E_i\Delta Y_i(s)E_i^T.$$

Definiendo

$$U=[E_1\cdots E_m],\quad V^T=\operatorname{col}(E_i^T),\quad C(s)=\operatorname{blkdiag}(\Delta Y_i(s)),$$

$$K=V^TT_0^{-1}U,\qquad M=CK,$$

se obtiene

$$\frac{\det T_S}{\det T_0}=\det(I+M_{SS}). \tag{18}$$

Cada acción local tiene hasta dos canales de tensión/corriente en el modelo equilibrado dq/xy. Los bloques $C_i$ pueden ser racionales en $s$.

Con

$$D=\operatorname{blkdiag}(I+M_{ii}),\qquad Q=D^{-1}(I+M)-I, \tag{19}$$

$$\det(I+M_{SS})=\prod_{i\in S}\det(I+M_{ii})\det(I+Q_{SS}). \tag{20}$$

El corte de bloques físicos debe permanecer fijo. Una similitud por bloques produce $Q'=S^{-1}QS$, pero reagrupar o mezclar proyectos distintos cambia la definición de «efecto individual». Ambas cosas no deben confundirse.

## 6.2 Métricas e invariancia

El espectro y el determinante son invariantes bajo similitud. La norma euclídea y el menor valor singular no lo son bajo transformaciones no unitarias. Eso no vuelve inútiles a todas las normas: una norma física ponderada es coherente si su métrica se transforma por congruencia.

Para $x_i=S_i x'_i$, la energía $x_i^*W_i x_i$ exige $W'_i=S_i^*W_iS_i$. La ganancia inducida entre los espacios métricos así transformados se conserva.

El ejemplo

$$Q_L=\begin{bmatrix}0&L\\0&0\end{bmatrix}$$

tiene $\operatorname{dist}(-1,\sigma Q_L)=1$ para cualquier $L$, mientras que

$$\sigma_{\min}(I+Q_L)\sim L^{-1}.$$

Una gran distancia espectral a −1 no excluye sensibilidad no normal. El margen de cierre será un **descriptor local de frontera**, no un radio universal de robustez.

# 7. Extensión al cruce aperiódico: qué queda demostrado

T1 elimina la simetría conocida sin borrar un polo físico que cruza $s=0$. Esto cierra el problema del conteo en espacio de estados. **No garantiza todavía que el operador de puertos original se vuelva invertible en el origen.**

La extensión de puertos debe hacerse en este orden:

1. Construir el cociente físico del descriptor/DAE.
2. Derivar una nueva partición de Schur de ese cociente, sin proyectar arbitrariamente una fila y una columna del viejo $T$.
3. Verificar su identidad de determinantes y su contabilidad de polos.
4. Usar $Q^q(0)$ sólo si la referencia y los factores individuales son regulares allí.

Si la referencia tiene un modo físico nulo o el bloque oculto conserva un integrador, la singularidad no es una gauge a borrar. Puede ser necesario usar una factorización coprima, una realización estable alternativa o el descriptor completo. [L4] es un antecedente especialmente relevante para representaciones de puertos con polos inestables.

**Resultado E prioritario:** intentar certificar los 194 cruces aperiódicos Kundur reportados en G1, conservando el modo real peligroso y demostrando cuáles son accesibles a la nueva reducción. No prometer 194/194 antes de ejecutar.

# 8. Certificado positivo de composabilidad por pequeña ganancia

## Teorema T3 — Certificado para todas las subfamilias

Supóngase, en una realización común y una política fija:

1. El baseline y todas las acciones individuales son internamente estables después del cociente de simetría.
2. Existe un $Q(s)$ propio y estable, con diagonal de bloques nula, tal que la selección de un portfolio corresponde a $Q_{SS}$ y la factorización característica/polos está verificada.
3. Para normas físicas fijas se conocen cotas rigurosas

$$\bar r_{ij}\ge \|Q_{ij}\|_\infty,\quad \bar r_{ii}=0.$$

4. La matriz no negativa $\bar R=[\bar r_{ij}]$ cumple $\rho(\bar R)<1$.

Entonces todos los portfolios de la familia son internamente estables y $\mathcal H=\varnothing$.

**Demostración.** Por Perron–Frobenius/propiedades de matrices no negativas, existe $d>0$ y $q<1$ con $\bar R d\le qd$. Por ejemplo, $(I-\bar R)^{-1}{\bf1}>0$ proporciona pesos con desigualdad estricta. En la norma de señales

$$\|x\|_d=\max_i\|x_i\|_{L_2}/d_i,$$

cada $Q_{SS}$ tiene ganancia como máximo $q$. La serie de Neumann establece un inverso estable de $I+Q_{SS}$; las hipótesis de estabilidad interna y la contabilidad característica excluyen modos ocultos peligrosos. La desigualdad se conserva al eliminar bloques porque todos los $\bar r_{ij}$ son no negativos. $\square$

**Origen matemático:** adaptación de pequeña ganancia y comparación positiva, no un criterio nuevo de estabilidad por sí solo [L2, L3, L7].

### Lo que puede fallar

Conocer sólo ganancias en $j\omega$ sin verificar estabilidad de $Q$ es insuficiente. Por ejemplo $Q(s)=0.1/(s-1)$ cumple $\sup_\omega|Q(j\omega)|=0.1$, pero $1+Q$ se anula en $s=0.9$. Este contraejemplo está en el script.

Si el certificado falla, el resultado es **UNKNOWN**, no «inestable». Si los puertos no son estables se debe usar una transformación válida como las estudiadas en [L4], otra factorización, o el certificado de espacio de estados de la Sección 9.

## Corolario T3a — Cota inferior del orden mínimo

Para un tamaño $r$, defínase

$$c_r(d)=\max_i\ \operatorname{TopSum}_{r-1}\left\{\bar r_{ij}\frac{d_j}{d_i}:j\ne i\right\}. \tag{21}$$

Si $c_r(d)<1$, todos los portfolios con $|S|\le r$ son estables y

$$\boxed{\kappa\ge r+1.} \tag{22}$$

**Demostración.** Una fila activa puede tener como máximo $r-1$ bloques fuera de la diagonal. Su suma ponderada no supera (21). Aplicar la demostración de T3 a cada selección. $\square$

Este certificado obtiene una respuesta útil incluso cuando no puede certificar todo el cubo. El óptimo de pesos puede buscarse por bisección y programación geométrica/convexa cuando la formulación escogida lo permita; la factibilidad numérica debe comprobarse sobre las desigualdades originales.

### Ejemplo V2

Con seis canales, $\bar R_{ij}=0.3$ para $i\ne j$. El certificado global falla porque $\rho(\bar R)=1.5$. Sin embargo, $c_4=0.9<1$, de modo que ninguna selección de hasta cuatro canales falla. En el sistema $Q(s)=-\bar R/(s+1)$, la primera selección inestable tiene cinco canales: la cota $\kappa\ge5$ es exacta en este ejemplo.

### Frecuencia continua, no raster

Una rejilla finita da cotas inferiores del máximo de ganancia, no superiores. Un certificado debe emplear normas $H_\infty$ verificadas, bounded-real/KYP cuando proceda, aritmética de intervalos o una cota derivativa entre muestras y una cota de cola a alta frecuencia. Sin esto se etiqueta `SAMPLED_SCREEN`, no `CERTIFICATE`.

# 9. Certificado de un subcubo mediante Lyapunov

La pequeña ganancia puede resultar demasiado conservadora. Una segunda puerta opera en espacio de estados y permite podar **subcubos completos**.

Considérese una rama con decisiones forzadas $F$, decisiones libres $J$ y

$$A(\delta)=A_F+\sum_{i\in J}\delta_i A_i,\quad\delta_i\in\{0,1\}. \tag{23}$$

Debe existir una representación común exacta. No se presupone para SG–GFL sólo por compartir el flujo AC.

## Teorema T4 — Mayorante LMI para un subcubo

Si existen $P\succ0$, matrices simétricas $X_i\succeq0$ y $\epsilon>0$ tales que

$$X_i\succeq A_i^TP+PA_i,$$

$$A_F^TP+PA_F+2\gamma P+\sum_{i\in J}X_i\preceq-\epsilon I, \tag{24}$$

entonces todos los vértices de esa rama cumplen $\alpha(A)<-\gamma$. De hecho, la misma conclusión vale para el cubo continuo $0\le\delta_i\le1$.

**Demostración.** Como $0\le\delta_i\le1$,

$$\delta_i(A_i^TP+PA_i)\preceq\delta_iX_i\preceq X_i.$$

Sumar y aplicar el criterio de Lyapunov a $A(\delta)+\gamma I$. $\square$

Una versión barata fija $P$ resolviendo Lyapunov en $A_F+\gamma I$ y usa $X_i$ igual a la parte positiva de $A_i^TP+PA_i$. La versión SDP busca simultáneamente $P,X_i$. El test barato se implementó en `theory_checks.py`.

## Corolario T4a — Transferencia al modelo reequilibrado

Si $\widehat A$ tiene certificado

$$\widehat A^TP+P\widehat A+2\gamma P\preceq-\epsilon I$$

y el modelo exacto reequilibrado satisface $A=\widehat A+E$, entonces basta

$$2\|P\|\,\|E\|<\epsilon \tag{25}$$

para conservar la tasa de decaimiento certificada.

**Demostración.** $E^TP+PE\preceq2\|P\|\|E\|I$. $\square$

Esta es una manera honesta de intentar recuperar parte de la aceleración fuera del equilibrio común. Medir $E$ después de construir ambos modelos no ahorra ese primer cálculo, pero una cota analítica/intervalar válida sobre una celda puede permitir reutilización posterior. Sin cota no se reutiliza el operador.

# 10. Certificados que sí respetan la discreción binaria

T3 y T4 pueden certificar también interpolaciones continuas. Eso los hace seguros, pero puede volverlos inútiles cuando sólo se implementan vértices binarios estables y entre ellos hay una interpolación inestable.

## Ejemplo V3 — El cubo continuo exige algo que no pedimos

$$A_0=\begin{bmatrix}-1&4\\0&-1\end{bmatrix},\qquad
A_1=\begin{bmatrix}-1&0\\4&-1\end{bmatrix}.$$

Ambos son Hurwitz con autovalores $-1$. Su promedio tiene autovalores $1,-3$ y es inestable. Por tanto ningún certificado que garantice estabilidad de toda la envolvente convexa podrá certificar este par de opciones, aunque ambas opciones finales sean seguras.

Defínase

$$P_0=\begin{bmatrix}1/2&1\\1&9/2\end{bmatrix},\quad
P_1=\begin{bmatrix}9/2&1\\1&1/2\end{bmatrix},$$

$$A(d)=(1-d)A_0+dA_1,\quad P(d)=(1-d)P_0+dP_1.$$

Entonces

$$-(A(d)^TP(d)+P(d)A(d))-I
=32(d^2-d)\begin{bmatrix}0&1\\1&0\end{bmatrix}. \tag{26}$$

Sobre $d\in\{0,1\}$ el término derecho es cero, $P(d)\succ0$ y existe un certificado de estabilidad de las dos opciones. La identidad fue verificada simbólicamente.

## Teorema T5 — Completitud finita de Lyapunov sobre el cubo booleano

Sea $A(\delta)$ una familia real de dimensión común definida en $\{0,1\}^m$. Todos sus vértices son Hurwitz si y sólo si existe una matriz polinómica multilineal simétrica $P(\delta)$, de grado a lo sumo $m$, tal que en todos los vértices

$$P(\delta)\succ0,\qquad
A(\delta)^TP(\delta)+P(\delta)A(\delta)\prec0. \tag{27}$$

Además, estas desigualdades admiten certificados de suma de cuadrados matricial módulo el ideal booleano

$$\mathcal I_B=\langle\delta_i^2-\delta_i:i=1,\ldots,m\rangle$$

cuando se permite grado suficiente.

**Demostración constructiva.** Si cada vértice $v$ es estable, resolver

$$A(v)^TP_v+P_vA(v)=-I.$$

Definir los idempotentes de interpolación

$$\chi_v(\delta)=\prod_{i:v_i=1}\delta_i\prod_{i:v_i=0}(1-\delta_i).$$

Entonces $P(\delta)=\sum_v\chi_v(\delta)P_v$ es multilineal, tiene grado máximo $m$ y satisface $P(v)=P_v$. Como hay finitos vértices, existe $\epsilon>0$ tal que $P_v-\epsilon I\succ0$ para todos. Si $L_v^TL_v=P_v-\epsilon I$,

$$P(\delta)-\epsilon I\equiv\sum_v\chi_v(\delta)^2L_v^TL_v\pmod{\mathcal I_B}.$$

La derivada negativa vale $I$ en cada vértice y admite el mismo tipo de representación. Un polinomio que se anula en todo el cubo tiene residuo multilineal nulo y pertenece a $\mathcal I_B$. La implicación inversa es el criterio de Lyapunov en cada vértice. $\square$

**Alcance y novedad.** Es una especialización de Lyapunov dependiente de parámetros y positividad en conjuntos finitos, con herramientas conocidas [L7–L8]. No demuestra un algoritmo polinómico. La construcción explícita enumera $2^m$ vértices. La hipótesis de investigación es que **grados bajos y estructura local** basten en muchas familias SG–IBR y certifiquen opciones que una caja continua de incertidumbre no logra certificar.

## 10.1 Jerarquía computable

Buscar $P_d(\delta)=\sum_{|I|\le d}P_I\prod_{i\in I}\delta_i$ y matrices SOS $S_0,S_1$ con

$$P_d-\epsilon I=S_0+\sum_i(\delta_i^2-\delta_i)R_i,$$

$$-(A^TP_d+P_dA)-\epsilon I=S_1+\sum_i(\delta_i^2-\delta_i)T_i. \tag{28}$$

Los multiplicadores $R_i,T_i$ son matrices polinómicas simétricas libres; el ideal de igualdad no requiere positividad. Para certificar sólo $\sum_i\delta_i\le r$, se añade el localizador SOS $(r-\sum_i\delta_i)S_r$ y los demás generadores necesarios para la región elegida.

El modelo $A(\delta)$ puede ser multilineal en vez de afín. Obtenerlo ajustando todos los vértices puede destruir cualquier ventaja computacional. Por ello, el experimento debe distinguir: **prueba de principio completa en cuatro acciones** frente a **realización local simbólica sin enumeración** en casos grandes.

No aplicar (28) a controladores que conmutan arbitrariamente entre vértices: $P$ depende de la decisión fija. La Sección 15 trata las condiciones adicionales para cambios en el tiempo.

# 11. Poda exacta sin monotonía y cotas de cardinalidad

## 11.1 Certificado homotópico de una celda

Sea un pencil matricial analítico $F_0(s)$ no singular sobre un contorno cerrado $\Gamma$, y una familia $F(s)=F_0(s)+E(s)$ analítica en y dentro del contorno. Si

$$\sup_{s\in\Gamma}\|F_0(s)^{-1}E(s)\|<1, \tag{29}$$

entonces $F$ y $F_0$ tienen el mismo número de ceros de determinante dentro de $\Gamma$.

**Demostración.** $F_t=F_0(I+tF_0^{-1}E)$ no es singular en el contorno para $0\le t\le1$. El número entero de ceros de $\det F_t$ permanece constante por el principio del argumento. $\square$

Para matrices racionales se usa el pencil analítico completo o se prueba que la contabilidad de polos permanece fija. No se ignoran cambios de polos ocultos de dispositivos reemplazados.

Si el ancla es estable y la cota se cumple uniformemente sobre todas las decisiones libres y toda la celda de $\theta$, se certifica una región entera. Si el ancla es inestable, se certifica que la región conserva ese conteo, lo que también puede ser útil.

## 11.2 Algoritmo exacto

Una rama del árbol se describe por $(F_1,F_0,J)$: acciones forzadas a uno, a cero y aún libres.

1. Validar factibilidad y representación de la rama.
2. Intentar T3, T4, T5 o (29), en ese orden de costo según el modelo.
3. Si la rama se certifica estable, podarla.
4. Si no se certifica, dividir una decisión binaria o evaluar una hoja exactamente.
5. Al encontrar un conjunto inestable, registrar un testigo. Para el objetivo «todos los testigos mínimos», sus superconjuntos ya no pueden ser mínimos. Esa poda **no etiqueta su estabilidad**.
6. Probar minimalidad por eliminación exacta de todos los subconjuntos relevantes, por certificados o por un recorrido de cardinalidad creciente. No basta un greedy que quite una acción a la vez.

El ejemplo $U=\{\{1\},\{1,2,3\}\}$ tiene un triple inestable cuyas tres eliminaciones inmediatas son estables, aunque contiene un singleton inestable. Un algoritmo greedy daría un falso testigo mínimo.

## 11.3 Certificado de cardinalidad óptima

Mantener cotas

$$\underline\kappa\le\kappa\le\overline\kappa.$$

Un testigo inestable de tamaño $b$ establece $\overline\kappa\le b$. Certificar todos los conjuntos de tamaño menor que $b$ establece $\underline\kappa\ge b$. Sólo cuando ambas coinciden se declara cardinalidad mínima global.

## 11.4 Límite de complejidad

No existe aquí una promesa de búsqueda subexponencial. El sistema escalar

$$A(S)=|S|-(r-\tfrac12)$$

tiene como hiperaristas todos los subconjuntos de tamaño $r$: $\binom mr$ resultados. Para $r\approx m/2$, sólo escribir la salida requiere tamaño exponencial. Sin estructura adicional, un oráculo arbitrario de estabilidad puede ocultar el único conjunto inestable en cualquier vértice no evaluado.

La contribución algorítmica debe ser **reducción medida de evaluaciones con garantías**, no una afirmación imposible sobre todos los casos.

# 12. Del hipergrafo al diseño: el significado exacto de «conservador»

Sea $\mathcal F_{\rm safe}=\{S:\alpha_q(S)<0\}$. No tiene por qué ser hereditaria. Definamos su núcleo hereditario

$$\mathcal C=\{S:\text{todo }R\subseteq S\text{ es factible y estable}\}. \tag{30}$$

## Teorema T6 — Hipergrafo y seguridad independiente del orden

Si todos los portfolios considerados son factibles y el baseline es estable,

$$S\in\mathcal C\quad\Longleftrightarrow\quad
\text{ninguna }H\in\mathcal H\text{ está contenida en }S. \tag{31}$$

$\mathcal C$ es la mayor subfamilia hereditaria de $\mathcal F_{\rm safe}$. Además, para una política fija, (31) equivale a que **cada secuencia de incorporaciones de una acción por etapa hasta $S$ tiene todos sus estados estacionarios intermedios estables**.

**Demostración.** Si $S$ contiene un $H$ inestable, no todos sus subconjuntos son estables. Si existe un subconjunto inestable de $S$, uno de sus subconjuntos inestables de cardinalidad mínima pertenece a $\mathcal H$ y está contenido en $S$. La heredabilidad es inmediata. Cualquier subfamilia hereditaria segura sólo puede contener conjuntos cuyos subconjuntos son seguros. Finalmente, cada $R\subseteq S$ puede ser prefijo de alguna permutación de $S$. $\square$

Así, evitar hiperaristas no es sólo una aproximación conservadora sin interpretación. Es **exacto para el requisito operacional de seguridad ante cualquier orden de implementación parcial**, siempre a nivel estacionario small-signal.

## 12.1 Optimización de reemplazos

Con $p_i$ en MW correctamente definidos y $x_i\in\{0,1\}$,

$$\max\sum_i p_i x_i$$

sujeto a

$$\sum_{i\in H}x_i\le |H|-1\qquad\forall H\in\mathcal H. \tag{32}$$

Esto resuelve el máximo plan seguro para cualquier subconjunto/prefijo si el hipergrafo es completo y no cambia la política. No debe afirmarse que encuentra todo endpoint estable. Con un hipergrafo incompleto, (32) no certifica seguridad frente a hiperaristas aún desconocidas; hacen falta búsqueda/certificación de separación.

Restaurar soporte síncrono cambia el modelo y generalmente cambia $\mathcal H$. Un hitting set calculado sobre el viejo hipergrafo no prueba que un condensador elegido haya «borrado» esas amenazas: hay que evaluar o certificar el modelo modificado.

## 12.2 Tres nociones de planificación

$$\text{seguro para cualquier orden}\Rightarrow
\text{existe un orden seguro}\Rightarrow
\text{endpoint estable}. \tag{33}$$

Las implicaciones inversas fallan. El ejemplo diagonal

$$A_0=-I_2,\quad A_1=\operatorname{diag}(2,-3),\quad A_2=\operatorname{diag}(-3,2)$$

tiene baseline estable, ambos singles inestables y endpoint de dos acciones $-2I_2$ estable. No existe camino monotónico de una acción por etapa a ese endpoint sin atravesar una configuración inestable. El ejemplo está ejecutado.

Buscar un camino permitido en el lattice es útil, pero cada arista seguirá necesitando una validación dinámica de transición si se pretende operación real y no sólo planificación estacionaria.

# 13. Robustez: cuantificadores antes de porcentajes

Sea $\mathcal U$ un conjunto de incertidumbre física. Defínase

$$\kappa_{\rm wc}(\theta)=\min_{\xi\in\mathcal U}\kappa(\theta,\xi). \tag{34}$$

Para el diseño seguro frente a cualquier condición,

$$\mathcal H^{\exists}(\theta)=\min_{\subseteq}\{S:\exists\xi\in\mathcal U,\ \alpha_q(S;\theta,\xi)\ge0\}. \tag{35}$$

Cuando la definición de inseguridad incluye marginalidad, el mismo criterio se conserva en todas las ecuaciones. Para un conjunto finito de escenarios, el lado derecho coincide con tomar los mínimos de la unión de hiperaristas de todos los escenarios. Para incertidumbre continua, una muestra sólo aproxima la unión; no la certifica.

De forma diferente,

$$\mathcal H^{\forall}(\theta)=\min_{\subseteq}\{S:\forall\xi\in\mathcal U,\ \alpha_q(S;\theta,\xi)\ge0\}. \tag{36}$$

No es, en general, la intersección de las hiperaristas mínimas de cada escenario. La operación «tomar mínimos» y los cuantificadores no conmutan libremente.

Los certificados T3/T4/(29) pueden hacerse uniformes en $\xi$ mediante cotas verificadas. La jerarquía SOS puede incorporar restricciones polinómicas de $\mathcal U$, pero su costo y su conservadurismo son parte del experimento.

**Estadística.** Una tasa 140/140 no tiene intervalo [1,1]. Los intervalos binomiales exactos requieren ensayos Bernoulli iid o hipótesis explícitas. Un único Latin Hypercube introduce dependencia entre muestras; para inferencia de probabilidad usar campañas iid o replicaciones independientes del diseño y un estimador adecuado. La tasa de éxito en una malla determinista no es una probabilidad operacional.

# 14. Teorema local del margen de cierre

La coincidencia de frecuencia entre puertos y DAE no tiene que quedar como una correlación.

## Teorema T7 — Ley local del retorno cerca de un cruce simple

Sea $Q(s,t)$ analítico en $s$ y $C^2$ en un parámetro real $t$. Supóngase que una rama simple $\mu(s,t)$ cumple

$$f(s,t)=1+\mu(s,t),\quad f(j\omega_*,t_*)=0,\quad f_s(j\omega_*,t_*)\ne0,$$

y que el cero correspondiente $\lambda(t)=\alpha(t)+j\omega_\lambda(t)$ cruza transversalmente, $\alpha'(t_*)\ne0$. En una banda local que contiene sólo esa rama crítica,

$$m_{\rm loc}(t)=\min_{\omega\text{ local}}|1+\mu(j\omega,t)|$$

satisface

$$m_{\rm loc}(t)=|f_s(j\omega_*,t_*)|\,|\alpha(t)|+O((t-t_*)^2), \tag{37}$$

$$\omega_{\min}(t)=\omega_\lambda(t)+O((t-t_*)^2). \tag{38}$$

**Demostración.** Por función implícita existe la rama $\lambda(t)$. Localmente,

$$f(s,t)=a(s,t)(s-\lambda(t)),\quad a(j\omega_*,t_*)=f_s(j\omega_*,t_*).$$

En $s=j\omega$, $|s-\lambda|^2=\alpha^2+(\omega-\omega_\lambda)^2$. El factor $a$ es suave y no nulo. Minimizar desplaza la frecuencia óptima sólo a segundo orden en $\alpha$; la transversalidad hace $\alpha=O(t-t_*)$. Sustituir da (37)–(38). $\square$

La derivada es computable en el operador pequeño:

$$f_s=\frac{w^*Q_s v}{w^*v},\quad
\frac{d\lambda}{dt}=-\frac{f_t}{f_s}. \tag{39}$$

Así,

$$\frac{m_{\rm loc}}{|f_s|}\approx|\alpha|$$

produce una escala local en $\mathrm{s}^{-1}$. Es una aproximación de distancia modal local, **no el lado estable/inestable**, que requiere conteo/orientación. Tampoco es una cota global de robustez.

**Límites.** Repetición/defectividad de $\mu$, $f_s=0$, cruce tangencial, singularidad del baseline o normalización, otra frecuencia con un mínimo menor y modos no visibles. Un mínimo global puede cambiar de rama. El caso de lóbulo/Fold de F7 exige tratar por separado la falta de transversalidad de la trayectoria elegida.

El resultado es una aplicación local de perturbación analítica e implicit-function theory. Su valor de publicación dependerá de mostrar que proporciona una calibración o cota útil frente al margen sin normalizar, no de reclamar prioridad de la expansión de Taylor.

# 15. Extensión no lineal y transición híbrida: una primera garantía válida

El sistema (1), restringido a la variedad algebraica regular y expresado en coordenadas físicas relativas, admite localmente

$$\dot e=A_qe+r(e),\qquad \|r(e)\|\le L\|e\|^2 \quad(\|e\|\le r_0). \tag{40}$$

La cota $L$ debe provenir de derivadas/Hessianos y límites válidos en la región; un muestreo aleatorio de Hessianos no basta como cota certificada.

## Teorema T8 — Entorno local no lineal uniforme

Si para cada portfolio admitido existe $P_S\succ0$ con

$$A_{q,S}^TP_S+P_SA_{q,S}\preceq-2\beta P_S,$$

entonces dentro de

$$\|e\|<r_S:=\min\left\{r_0,\frac{\beta\lambda_{\min}(P_S)}{\|P_S\|L_S}\right\}, \tag{41}$$

$V_S=e^TP_Se$ decrece estrictamente. Una subnivelación $V_S\le c_S$ contenida en esa bola y en el dominio del mismo active set es positivamente invariante.

**Demostración.**

$$\dot V_S\le-2\beta\lambda_{\min}(P_S)\|e\|^2+2\|P_S\|L_S\|e\|^3<0.$$

Elegir $c_S<\lambda_{\min}(P_S)r_S^2$ garantiza que toda la elipsoide queda dentro de la bola. $\square$

Se intersecta la región con límites de corriente, tensión, Q disponible, ángulos permitidos y distancia a singularidad de $g_z$. No es estabilidad transitoria global ni cobertura de una falla severa; sí es un puente demostrable desde la certificación espectral a perturbaciones finitas pequeñas.

## 15.1 Cambio de portfolio

Una transición tiene un reset/mapa de inicialización

$$e^+=\mathcal R_{S\to S'}(e^-),$$

que incluye cambio de equilibrio y consistencia algebraica. Para garantizar una secuencia, exigir

$$\mathcal R_{S\to S'}(\{V_S\le c_S\})\subseteq\{V_{S'}\le c_{S'}\}. \tag{42}$$

Si sólo se dispone de $V_{S'}(e^+)\le\mu V_S(e^-)$ y $\dot V\le-2\beta V$ entre cambios, un tiempo de permanencia mayor que $\log\mu/(2\beta)$ es suficiente en el caso homogéneo sin término constante por desplazamiento de equilibrio. Con desplazamiento, (42) o una cota afín explícita es necesaria.

El script incluye dos matrices individualmente estables cuya conmutación periódica es inestable. Esto impide afirmar que un camino de configuraciones Hurwitz equivale automáticamente a una transición híbrida segura.

# 16. Qué se comprobó en esta sesión

`theory_checks.py` ejecutó **404 comprobaciones/assertions** deterministas. El archivo JSON contiene parámetros y resultados. No son 404 ensayos de redes eléctricas ni una validación estadística.

| Comprobación | Resultado |
|---|---|
| Cociente de simetría | Conserva polos físicos $10^{-7}$, cero físico y negativos; rechaza retirar una cadena de Jordan como dos simetrías. |
| Identidad de puertos y similitud por bloques | Coinciden dentro de tolerancias de doble precisión en los puntos sintéticos elegidos. |
| Margen espectral frente a no normalidad | $d_{-1}=1$ constante aunque el menor singular cae aproximadamente como $1/L$. |
| Certificado de cardinalidad | Certifica todas las selecciones de hasta cuatro; el mínimo inestable real es cinco. |
| Certificado booleano | Identidad (26) exacta simbólicamente; dos vértices estables con promedio inestable. |
| Búsqueda con poda segura | Recupera exactamente el hipergrafo en cinco familias diagonales afines de 512 portfolios con restabilización. Usó 66–319 consultas espectrales frente a 512, pero añadió 59–240 soluciones de Lyapunov. No se reclama speedup por ese conteo. |
| Núcleo hereditario | Coincide con evitar hiperaristas en ocho funciones booleanas no monótonas aleatorias. |
| Seguridad de endpoint y secuencia | Contraejemplo explícito con endpoint estable y ambos singles inestables. |
| Margen local de cierre | Coincide con la ley lineal cerca de un cruce oscilatorio sintético. |
| Conmutación | Dos modos estables producen multiplicador de Floquet mayor que uno bajo conmutación. |
| Certificado no lineal | Comprobación directa sobre $\dot x=-x+x^2$. |

El banco sintético es un punto de partida de regresión, no un benchmark favorable elegido para demostrar superioridad computacional. El certificado puede no aceptar ningún subcubo difícil de IEEE-39; eso sería un resultado a reportar.

# 17. Experimentos para cerrar el artículo

La especificación ejecutable en lenguaje natural está en `PROMPT_CLAUDE_CODE.md`. El siguiente resumen define el orden lógico y la evidencia que falta.

| ID | Pregunta | Experimento obligatorio | Resultado interpretable |
|---|---|---|---|
| C00 | ¿Hay errores de unidades o filtros de modos? | Auditoría MW/MVA, polaridad, bases, simetrías, conteo sin disco excluido | Descartar errores de reporte y descubrir si cambian $\mathcal H,\kappa$. |
| C01 | ¿El cociente cierra los cruces aperiódicos? | 194 cruces Kundur; full descriptor, cociente, nueva reducción a cero | Cobertura real por puertos, excepciones e identidad de polos. |
| C02 | ¿Es exacta la familia binaria común? | Comparar realizaciones y actualizaciones en los 16/512/4096 vértices disponibles | Aplicabilidad o rechazo de afinidad; no confundir con F7. |
| C03 | ¿T3/T3a certifican algo útil? | Q estable, cotas de ganancia; todos los subsets y cardinalidades | Cotas $\underline\kappa$, fracción certificada, conservadurismo. |
| C04 | ¿T4 y T5 superan la relajación continua? | Lyapunov común, subcubo, SOS booleano grado 1–3 en cuatro acciones | Casos binariamente seguros con caja no certificable, sin falsos seguros. |
| C05 | ¿La poda ahorra realmente? | Búsqueda certificada vs enumeración optimizada con caché de factores | Hipergrafo idéntico; consultas, setup, CPU, memoria y gap. |
| C06 | ¿La incertidumbre cambia la cardinalidad mínima? | Celdas físicas y cotas uniformes; escenarios sólo como contraste | $\kappa_{wc}$ certificado o intervalo con cobertura explícita. |
| C07 | ¿La estructura mejora una decisión? | MILP order-agnostic, camino order-aware y óptimo exhaustivo pequeño | MW reales, soporte, costo normalizado y garantía distinguible. |
| C08 | ¿El margen local se calibra? | Fronteras simples en tres redes, sin usar frecuencias DAE para buscarlas | Error de $m/|f_s|$ y frecuencia; fallos cerca de folds. |
| C09 | ¿Se puede reutilizar fuera del equilibrio común? | Bound (25)/(29) vs reequilibración exacta | Celdas autorizadas y rechazadas; cero errores de certificación. |
| C10 | ¿Existe entorno no lineal certificado? | Hessianos/guards y elipsoides de T8; TDS interiores/exteriores | Radio local, límites activos; no confundir con EMT. |
| C11 | ¿Es segura la implementación secuencial? | Dos rutas espectralmente estables, resets y rampas; verificar (42) | Distinguir endpoint, todos los prefijos y transición dinámica. |
| C12 | ¿La contribución supera el estado del arte? | Baselines de [L2–L8] adaptados bajo hipótesis compatibles | Beneficio y límites frente a Nyquist exhaustivo y robust control. |

# 18. Prioridad de publicación

El artículo aplicado existente puede mantener su aporte empírico si se conservan sus condiciones y resultados negativos. Para una contribución más fuerte, priorizar:

**Pieza A: conteo físico completo, incluido el origen.** Sustituye una exclusión numérica que puede ocultar inestabilidades por una reducción de simetría demostrada.

**Pieza B: certificados de familias discretas y búsqueda con prueba.** El resultado operacional no es «más puntos en un mapa», sino un límite inferior garantizado de cardinalidad y una lista de testigos mínimos obtenida sin pruebas inválidas de monotonía.

**Pieza C: diseño con una especificación exacta.** La seguridad independiente del orden se conecta de forma precisa con evitar hiperaristas. Se compara con planes que necesitan un orden particular y con endpoints seguros pero inalcanzables mediante pasos seguros.

**Pieza D: puente local no lineal.** Se añade sólo cuando hay cotas verificadas y no compromete la claridad del paper principal.

No hay garantía de aceptación en TPWRS ni prueba de prioridad de la combinación. La viabilidad editorial dependerá de si C03–C07 producen ventaja sobre baselines fuertes y si C00–C02 preservan los resultados físicos. No anunciar tres artículos por adelantado: cerrar primero uno con tres contribuciones comprobables.

# Referencias primarias y documentos

## Literatura consultada

[L1] V. Katewa y F. Pasqualetti, *On the Real Stability Radius of Sparse Systems*, arXiv:1810.10578. <https://arxiv.org/abs/1810.10578>

[L2] L. Huang et al., *Gain and Phase: Decentralized Stability Conditions for Power Electronics-Dominated Power Systems*, arXiv:2309.08037. <https://arxiv.org/abs/2309.08037>

[L3] L. Hallinan e I. Lestas, *Decentralised Plug-and-Play Stability Conditions for AC Grids—Part I*, arXiv:2609.07315, 7 septiembre 2026. <https://arxiv.org/abs/2609.07315>

[L4] L. Hallinan e I. Lestas, *Decentralised Plug-and-Play Stability Conditions for AC Grids—Part II: Unstable Subsystems*, arXiv:2609.07331, 7 septiembre 2026. <https://arxiv.org/abs/2609.07331>

[L5] Y. Zhu, Y. Gu, Y. Li y T. C. Green, *Participation Analysis in Impedance Models: The Grey-Box Approach for Power System Stability*, arXiv:2102.04226. <https://arxiv.org/abs/2102.04226>

[L6] P. Rodríguez-Ortega, D. Moutevelis, J. Roldán-Pérez y M. Prodanović, *Equivalent modelling for the fundamental frequency dynamic variation: State-space, impedance, and power-frequency representations*, arXiv:2305.19655v2, 11 de junio de 2026; Electric Power Systems Research 258, 113098. <https://arxiv.org/abs/2305.19655v2>

[L7] S. Boyd, L. El Ghaoui, E. Feron y V. Balakrishnan, *Linear Matrix Inequalities in System and Control Theory*, SIAM, 1994. <https://web.stanford.edu/~boyd/lmibook/>

[L8] P. A. Parrilo, *Semidefinite programming relaxations for semialgebraic problems*, Mathematical Programming 96, 293–320, 2003. Copia del autor: <https://www.mit.edu/~parrilo/pubs/files/SDPrelaxations.pdf>

[L9] H. Yazdani y S. Lotfifard, *Generalized Nyquist Criterion Limitations and Misconceptions for Frequency Domain Stability Analysis of Inverter-based Resources Integrated Power Grids*, arXiv:2608.07785. <https://arxiv.org/abs/2608.07785>

[L10] ANDES, documentación oficial de modelos de excitación. La página consultada se identifica como documentación 1.9.3; no se supone identidad con la instalación 2.0.0 del repositorio sin verificar el código fijado. <https://docs.andes.app/en/stable/groupdoc/Exciter.html>

## Documentos del proyecto leídos

[D1] `TPWRS_FINAL_EVIDENCE_TABLE.md`, versión suministrada el 10 septiembre 2026.

[D2] `G5_NOVELTY_STATEMENT.md`.

[D3] `G1_WHOLE_RHP_COMPOSABILITY.md`.

[D4] `G3_IEEE68.md`.

[D5] `F10_F11_BASELINES_AND_SCALING.md`.

[D6] `F1_ANDES_MODEL_RECONCILIATION.md`.

[D7] `PHASE_E_REPORT.md`, usado sólo para la alerta de unidades, no como fuente de claims vigentes.

Los enlaces de literatura identifican fuentes; los ejemplos, demostraciones y protocolos desarrollados en esta nota no se atribuyen como resultados empíricos de esas fuentes.
