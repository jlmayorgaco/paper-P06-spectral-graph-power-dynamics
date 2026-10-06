# Co-diseño SG → GFL: ecuaciones, fronteras algebraicas, robustez y validación incremental

**Base experimental:** ExpN, PowerDynamics 5.0.0, NetworkDynamics 1.3.0, Julia 1.11.9.
**Estado de este documento:** especificación matemática y de implementación. No contiene una nueva ejecución IEEE-39 ni un nuevo óptimo robusto/global.

## 0. Qué se conserva y qué se intenta demostrar

ExpN reporta un candidato nominal local con soporte SG {38}, 1.1243444776535036 MW síncronos retenidos, Kp=7.853981633974483 y Ki=986.9604401089358 en los diez GFL. Su reparto P/Q y sus matrices coinciden con una construcción independiente de PowerDynamics. Su certificado es local dentro de ese soporte y del rectángulo de ganancias, evaluado numéricamente. No hay certificado global ni de incertidumbre robusta. El margen adicional respecto a −0.05 s⁻¹ es diminuto.

Se conserva el modelo de ExpN, no los ensamblajes mixtos anteriores a la corrección P/Q. Las nuevas identidades algebraicas, la reducción GSP, el diseño robusto y la certificación global se validan por separado. Una etapa opcional que no funcione no invalida retrospectivamente la reproducción de ExpN.

La jerarquía de resultados es:

1. identidad de ecuaciones/puertos/realizaciones;
2. fronteras algebraicas exactas condicionadas a parámetros y soporte;
3. candidato continuo factible y certificado local;
4. candidato robusto para una incertidumbre definida;
5. cotas globales, solamente cuando sean matemáticamente válidas;
6. simulación no lineal independiente, que no sustituye ningún certificado anterior.

## 1. Variables, objetivo y unidades

Sean N=39 barras y n=10 emplazamientos de generación G={30,…,39}. Para cada i:

\[
\epsilon_i=1-\rho_i,\quad 0\le\rho_i\le1,\qquad
Z=(\epsilon,k_p,k_i).
\]

La minimización es

\[
J(Z)=w^T\epsilon,\qquad P_{\mathrm{reemplazado}}=\sum_iw_i-J.
\tag{1}
\]

En la reproducción de ExpN, w_i=P_{G,i}^0 es despacho inicial en MW. Para optimizar capacidad instalada activa se requieren capacidades C_i en MW verificadas y se cambia w; los MVA nominales Sni no son automáticamente MW instalados. Cambiar w define otro problema. No reutilizar el porcentaje de despacho como porcentaje de capacidad instalada.

Normalizar para optimización:

\[
u_{p,i}=K_{p,i}/K_{p0,i},\qquad u_{I,i}=K_{i,i}/K_{i0,i},
\quad u_{p,i},u_{I,i}\in[0.25,4].
\tag{2}
\]

El objetivo histórico usa base de sistema S_b=100 MVA, f_b=60 Hz y ω_b=2πf_b. Cada tensión usa su base física por barra. Las corrientes, potencias, parámetros y fuentes de perturbación deben documentar conversión a pu. Hz, rad/s y pu de velocidad no son intercambiables.

Cada barra puede tener SG sola, SG+GFL, o GFL solo. Hay hasta 3^10 combinaciones de presencia de ambos tipos. Cerca de un candidato con J<min_i w_i, no puede haber ε_i=1; allí todos los GFL están presentes y bastan las 2^10 máscaras de SG. Fuera de esa región no confundir ambos conteos.

## 2. Red, cargas y contrato P/Q

Usar inyección positiva hacia la red. V_j e I_j son fasores en la misma base de sistema:

\[
S_j=V_j\overline I_j.
\]

Para una rama π con admitancia serie y_l, admitancia shunt total y_sh y tap complejo a en el extremo from:

\[
\begin{bmatrix}I_f\\I_t\end{bmatrix}=
\begin{bmatrix}
(y_l+y_{sh}/2)/|a|^2&-y_l/\bar a\\
-y_l/a&y_l+y_{sh}/2
\end{bmatrix}
\begin{bmatrix}V_f\\V_t\end{bmatrix}.
\tag{3}
\]

Es una representación estándar; el ensamblaje debe verificar orientación/taps/shunts contra la implementación congelada, no sustituirla sin comparación. En rectangular, para Y=G+jB:

\[
\mathcal Y=\begin{bmatrix}G&-B\\B&G\end{bmatrix}
\]

si se apilan primero partes reales y luego imaginarias. Un apilado intercalado necesita permutación explícita.

Para consumo ZIP positivo, con r_j=|V_j|/V_{ref,j}:

\[
P_{L,j}=P_{L,ref,j}(z_{p,j}r_j^2+i_{p,j}r_j+p_{p,j}),
\]
\[
Q_{L,j}=Q_{L,ref,j}(z_{q,j}r_j^2+i_{q,j}r_j+p_{q,j}),\quad
I_{L,j}=\overline{(P_{L,j}+jQ_{L,j})/V_j}.
\tag{4}
\]

Los coeficientes y referencias son los inicializados de ExpN. En el caso de impedancia constante, la carga puede incorporarse a Y. No convertir una carga P o I en Z sin declararlo. En barras 31 y 39 conservar por separado carga y generador; el neto no determina por sí solo ambos componentes.

El contrato de sustitución es:

\[
S_{SG,i}^0=\epsilon_iS_{G,i}^0,\qquad
S_{F,i}^0=(1-\epsilon_i)S_{G,i}^0,\qquad V_i^0\ \text{fijo}.
\tag{5}
\]

Por tanto

\[
I_{SG,i}^0=\epsilon_i\overline{S_{G,i}^0/V_i^0},\quad
I_{F,i}^0=(1-\epsilon_i)\overline{S_{G,i}^0/V_i^0}.
\tag{6}
\]

KCL conserva el mismo equilibrio externo:

\[
Y_{net}V+I_L(V)-\sum_i\Pi_i(I_{SG,i}+I_{F,i})=0,
\tag{7}
\]

con Π_i inserción nodal. Estas identidades no garantizan que las ecuaciones internas admitan cualquier trim: hay que resolverlas y comprobar límites.

## 3. SG detallada: no sustituir Sauer–Pai por un swing de dos estados

La fuente consultada es `PowerDynamics.jl/v5.0.0/src/Library/Machines/SauerPaiMachine.jl`. Las fórmulas siguientes usan su convención rotor dq, distinta de la d-alineada del GFL.

Defina

\[
a_{d1}={X_d''-X_{ls}\over X_d'-X_{ls}},\quad
a_{d2}={X_d'-X_d''\over(X_d'-X_{ls})^2},
\]

con definiciones análogas para q. Sean E_d',E_q',ψ_d'',ψ_q'',δ,ω los seis estados cuando no se habilitan dinámicas de estator. Las ecuaciones son:

\[
\dot\delta=\omega_b(\omega-\omega_{frame}),\qquad
2H\dot\omega=\tau_m-\tau_e-D(\omega-1),
\tag{8}
\]

\[
\tau_e=\psi_d I_q-\psi_q I_d,
\]

\[
T_{d0}'\dot E_q'=-E_q'-(X_d-X_d')
[I_d-a_{d2}\psi_d''-(1-a_{d1})I_d+a_{d2}E_q']+v_f,
\]

\[
T_{q0}'\dot E_d'=-E_d'+(X_q-X_q')
[I_q-a_{q2}\psi_q''-(1-a_{q1})I_q-a_{q2}E_d'],
\]

\[
T_{d0}''\dot\psi_d''=-\psi_d''+E_q'-(X_d'-X_{ls})I_d,
\]
\[
T_{q0}''\dot\psi_q''=-\psi_q''-E_d'-(X_q'-X_{ls})I_q.
\tag{9}
\]

Las relaciones algebraicas son

\[
0=R_s I_d+\omega\psi_q+V_d,\qquad
0=R_s I_q-\omega\psi_d+V_q,
\]

\[
\psi_d=-X_d''I_d+a_{d1}E_q'+(1-a_{d1})\psi_d'',
\]
\[
\psi_q=-X_q''I_q-a_{q1}E_d'+(1-a_{q1})\psi_q''.
\tag{10}
\]

Con

\[
T_{SG}(\delta)=\begin{bmatrix}\sin\delta&-\cos\delta\\\cos\delta&\sin\delta\end{bmatrix},
\]

\[
v={V_n\over V_b}T_{SG}^T[V_d,V_q]^T,
\quad [I_d,I_q]^T={I_b\over I_n}T_{SG}i_{SG}.
\tag{11}
\]

El cociente de bases de corriente debe calcularse con los S,V correspondientes. Sn=ε Sn0, con H fijo, escala HSn una sola vez. No multiplicar además H por ε.

### AVRTypeI y TGOV1

Para las máquinas controladas, conservar los bloques originales. En su rama interior suave:

\[
T_r\dot v_m=|v_{dq}|-v_m,
\quad T_f\dot v_{fb}=K_f\dot v_f-v_{fb},
\]
\[
T_a\dot v_r=\mathcal A(v_r,K_a(v_{ref}-v_m-v_{fb});v_{r,min},v_{r,max}),
\]
\[
T_e\dot v_f=v_r-\Phi(v_f)-K_e v_f.
\tag{12}
\]

Aquí A(z,r;l,u)=0 si (z>u y r>z) o (z<l y r<z); en otro caso A=r-z, reproduciendo la lógica anti-windup consultada. En TDS usar exactamente la lógica congelada, no una suavización nueva. Si el equilibrio cae en una frontera no suave, las derivadas clásicas requieren tratar la rama activa.

Para la saturación cuadrática del ejemplo IEEE-39, Φ(vf)=|vf| SE(|vf|). Si SE1≠SE2:

\[
a=\sqrt{SE_1E_1/(SE_2E_2)},\quad
A_s=E_2-(E_1-E_2)/(a-1),
\quad B_s=SE_2E_2(a-1)^2/(E_1-E_2)^2,
\]

\[
SE(u)=\begin{cases}0,&u\le A_s,\\B_s(u-A_s)^2/u,&u>A_s.\end{cases}
\tag{13}
\]

Si SE1=SE2, la función congelada devuelve esa constante. No evaluar puntos singulares artificiales fuera del dominio físico.

TGOV1 usa

\[
\Delta\omega=\omega-\omega_{ref},\quad r_g=(p_{ref}-\Delta\omega)/R,
\]
\[
T_1\dot x_{g1}=\mathcal A(x_{g1},r_g;V_{min},V_{max}),
\quad T_3\dot x_{g2}=x_{g1}+T_2\dot x_{g1}-x_{g2},
\]
\[
\tau_m\omega=x_{g2}-D_T\Delta\omega.
\tag{14}
\]

No sustituir esta relación de potencia/torque por τm=xg2 fuera de ω=1. En bus 39, el modelo congelado no tiene estos controles: v_f y τ_m son sus referencias fijas inicializadas.

## 4. GFL de nueve estados y su trim P/Q explícito

Estas fórmulas transcriben en notación vectorial el bloque `SimpleGFLDC` de v5.0.0; los parámetros numéricos son los de ExpN, no los defaults del tutorial. Se revisaron los bloques LFilter, CC1 y PLL_LPF. La documentación textual de algún eje no sustituye el código de conexión: la referencia activa de corriente entra al eje d en el código consultado.

Defina

\[
R(\theta)=\begin{bmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{bmatrix},\qquad
W=\begin{bmatrix}0&1\\-1&0\end{bmatrix},
\]
\[
v_{dq}=R^Tv,\quad i_{dq}=R^Ti_f,\quad e=e_2^Tv_{dq}.
\]

Los estados son θ,ν,ξ (PLL), i_f∈R² (filtro), γ∈R² (CC), vdc,ηdc (DC): nueve en total.

\[
\dot\theta=\nu,\qquad
\dot\xi=K_i e,\qquad
\tau\dot\nu=\xi+K_p e-\nu.
\tag{15}
\]

\[
i_{ref}=\begin{bmatrix}k_{pv}(V_{dc}-v_{dc})+\eta_{dc}\\ i_{q,set}\end{bmatrix},\quad
\dot\gamma=i_{ref}-i_{dq},
\]

\[
v_{I,dq}=-F_cX_fWi_{dq}+K_{PC}(i_{ref}-i_{dq})+K_{IC}\gamma+F_vv_{dq},
\quad v_I=Rv_{I,dq},
\]

\[
{X_f\over\omega_b}\dot i_f=v_I-v-R_fi_f+\omega_{frame}X_fWi_f,
\tag{16}
\]

\[
\dot\eta_{dc}=k_{iv}(V_{dc}-v_{dc}),\qquad
C_{dc}\dot v_{dc}={v_{I,dq}^Ti_{dq}-P_{dc}\over v_{dc}}.
\tag{17}
\]

El signo de (17) es la convención del modelo consultado. No cambiarlo basándose en otra convención de potencia. Las ganancias CC y DC no son las ganancias PLL a optimizar.

Para el agregado, i_F=ρ i_f cuando i_f representa el módulo equivalente de despacho completo en la convención validada de ExpN. Es obligatorio conservar la interpretación de energía/potencia total del agregado y no contar las pérdidas dos veces.

### Trim cerrado

En el módulo equivalente, sea i_f^0=conj(S_G^0/V^0) en las bases correctas. Entonces:

\[
\theta^0=\arg V^0,\quad \nu^0=\xi^0=0,\quad
v_{dc}^0=V_{dc},\quad
\eta_{dc}^0=i_{d}^0,\quad i_{q,set}=i_q^0,
\]

\[
v_I^0=v^0+R_fi_f^0-\omega_{frame}X_fWi_f^0,
\]

\[
\gamma^0={1\over K_{IC}}
\left(R^Tv_I^0+F_cX_fWi_{dq}^0-F_vv_{dq}^0\right),
\]

\[
P_{dc}=(v_I^0)^Ti_f^0.
\tag{18}
\]

El error PLL es cero; por eso Kp,Ki no aparecen en esos valores de equilibrio. La invariancia bajo ρ depende del contrato de agregado. ExpN la verificó; no extrapolarla a otros tipos de inversor, otro reparto Q o puntos operativos sin nueva derivación. En unidades físicas convertir corriente/potencia mediante bases antes de usar (18).

Cada PLL opera con medidas locales. La red completa se usa offline para calcular sus ganancias; no se introduce comunicación online al parametrizar ganancias mediante GSP.

## 5. Linealización, reducción y polos internos

Agrupe estados diferenciales x y algebraicos v. En un equilibrio regular:

\[
E_x\delta\dot x=f_x\delta x+f_v\delta v+f_d d,
\quad 0=g_x\delta x+g_v\delta v+g_d d.
\tag{19}
\]

El pencil completo es

\[
\mathscr P(s,Z)=
\begin{bmatrix}sE_x-f_x&-f_v\\-g_x&-g_v\end{bmatrix}.
\tag{20}
\]

Para gv invertible, sin formar inversas explícitas:

\[
A_r=E_x^{-1}(f_x-f_vg_v^{-1}g_x),\quad
B_r=E_x^{-1}(f_d-f_vg_v^{-1}g_d).
\tag{21}
\]

Para salida y=h_x dx+h_v dv+h_d d:

\[
C_r=h_x-h_vg_v^{-1}g_x,\quad D_r=h_d-h_vg_v^{-1}g_d.
\tag{22}
\]

Eliminar solo el gauge probado mediante una transformación explícita. No borrar polos pequeños por umbral. Para equivalencias descriptor generales pueden intervenir transformación de estados y transformación de ecuaciones distintas: E_an=L E_PD T, A_an=L A_PD T. No exigir sin justificación una semejanza única para ambos.

El puerto del dispositivo, con corriente saliente hacia red, es

\[
Y_a(s)=D_a+C_a(sM_a-A_a)^{-1}B_a.
\tag{23}
\]

Bajo el contrato de ExpN:

\[
Y_i^{mix}=\epsilon_iY_{SG,i}^{(0)}+(1-\epsilon_i)Y_{F,i}^{(0)}(s,k_i).
\tag{24}
\]

Con YL incluido en la red si es carga Z, el operador nodal es

\[
T(s,Z)=\mathcal Y_{net}+\mathcal Y_L-
\sum_i\Pi_iY_i^{mix}\Pi_i^T.
\tag{25}
\]

La igualdad por Schur entre det(P) y det(T) incluye determinantes de bloques eliminados. T es meromorfo. Por ello NO se iguala sin más “todos los polos” con “ceros de det(T)”: reconstruir multiplicidades y cancelaciones o usar el pencil completo como árbitro. En particular, C=I+ΔG es I cuando Δ=0 y no exhibe por sí solo los polos del sistema de referencia.

## 6. Dos identidades algebraicas que pueden acelerar el cálculo

### 6.1. Retención SG de un puerto: ecuación cuadrática condicionada

Con K fijo, defina T_F=operador all-GFL y D_i=Y_F,i−Y_SG,i. Entonces:

\[
T=T_F+\epsilon_i\Pi_iD_i\Pi_i^T.
\]

Cuando TF(s) y los bloques eliminados requeridos sean regulares:

\[
\det T=\det T_F\det(I_2+\epsilon_iN_i),\quad
N_i=D_i\Pi_i^TT_F^{-1}\Pi_i.
\tag{26}
\]

Para una frontera real sb=−σ:

\[
1+b_i\epsilon_i+c_i\epsilon_i^2=0,
\quad b_i=\operatorname{tr}N_i(sb),\quad c_i=\det N_i(sb).
\tag{27}
\]

Resolver todas las raíces reales en [0,1], incluyendo degeneración c≈0 con tratamiento bien condicionado. También pueden calcularse como −1/λ(Ni), verificando realidad/residuo. Evitar cancelación en la fórmula cuadrática. Después: verificar todos los polos, clasificación de soporte y restricciones. Esto es una raíz exacta de una frontera, NO una prueba de óptimo global.

Para dos puertos, formar Π=[Πi Πj] y Dε=diag(εiDi,εjDj):

\[
\det T=\det T_F\det(I_4+D_\epsilon\Pi^TT_F^{-1}\Pi).
\tag{28}
\]

No sumar dos raíces de un puerto: hay términos cruzados. Para frontera compleja, resolver Re F=Im F=0 con s=−σ+jω, además de las condiciones necesarias de diseño.

### 6.2. Un par Kp/Ki: dependencia afín del determinante

Esta identidad requiere demostrar sobre el pencil completo (o quotient válido de dimensión fija):

\[
\mathscr P(s,K)=\mathscr P_r(s)
-(\Delta K_{p,i}h_{p,i}+\Delta K_{i,i}h_{I,i})q_i^T.
\tag{29}
\]

Se usa una referencia Kr admisible; las diferencias no obligan a usar K=0. La estructura de (15) y el error común motivan (29), pero debe auditarse tras las conversiones de coordenadas y ensamblaje. Si E también cambia o la actualización no es rango uno, no afirmar la forma afín sin una nueva derivación.

La identidad de adjugadas da, incluso donde Pr sea singular:

\[
p(s,K)=p_r(s)+\Delta K_{p,i}p_{p,i}(s)+\Delta K_{i,i}p_{I,i}(s),
\]
\[
p_r=\det\mathscr P_r,\quad
p_{p,i}=-q_i^T\operatorname{adj}(\mathscr P_r)h_{p,i},\quad
p_{I,i}=-q_i^T\operatorname{adj}(\mathscr P_r)h_{I,i}.
\tag{30}
\]

No computar adjugadas densas en producción. Usar sistemas bordeados, LU con escalas comunes o representación polinómica estable. No expandir innecesariamente el polinomio de ~100+ estados en monomios Float64.

En sb=−σ+jω, ω>0, cuando Q sea invertible:

\[
Q_i(\omega)\begin{bmatrix}\Delta K_{p,i}\\\Delta K_{i,i}\end{bmatrix}=
-\begin{bmatrix}\Re p_r(sb)\\\Im p_r(sb)\end{bmatrix},
\quad
Q_i=\begin{bmatrix}\Re p_{p,i}&\Re p_{I,i}\\\Im p_{p,i}&\Im p_{I,i}\end{bmatrix}.
\tag{31}
\]

Para ω=0, todos los coeficientes son reales: Q tiene rango como máximo uno. Usar la recta

\[
p_r(-\sigma)+\Delta K_p p_p(-\sigma)+\Delta K_i p_I(-\sigma)=0.
\tag{32}
\]

No invertir una matriz 2×2 singular. Tampoco omitir líneas singulares, límites del rectángulo, pérdida de regularidad ni cruces de polos de bloques eliminados. El lado estable se identifica con el espectro completo/conteo, no con el dibujo de la curva.

Para los m PLL simultáneamente:

\[
\mathscr P=\mathscr P_r-H(\Delta k)Q^T,
\quad
p=p_r\det(I_m-Q^T\mathscr P_r^{-1}H).
\tag{33}
\]

Cada columna de H depende linealmente del par de ganancias de su dispositivo. El determinante es afín en cada par con los demás pares fijos, pero aparecen productos entre dispositivos. Diez PLL con dos ganancias comunes NO constituyen una única actualización de rango uno. La reducción a un determinante m×m acelera evaluación; no elimina automáticamente las veinte variables de diseño.

## 7. Sensibilidades, KKT y cómo calcular el candidato

Sea y*(z) el conjunto de estados y referencias de trim fijado por h(y*,z)=0. Cuando hy es invertible:

\[
y_z^*=-h_y^{-1}h_z,\qquad
{dT\over dz}=T_z+T_y y_z^*.
\tag{34}
\]

En ExpN, las coordenadas normalizadas nominales hacen yz*=0 para los parámetros de diseño considerados. Verificarlo, no añadir artificialmente un término no existente ni omitirlo al cambiar el punto operativo.

Para un polo simple:

\[
\lambda_z=-{l^H(dT/dz)r\over l^HT_sr}.
\tag{35}
\]

En descriptor: λz=lᴴ(Az−λEz)r/(lᴴEr). En ODE: λz=lᴴAzr/(lᴴr). No confundir autovector de T(λ) con autovector matricial al valor λ. No usar complex-step estándar sobre max, abs, conjugado, Re o eigenvalores reordenados; usar AD real y diferencias centradas con rastreo de rama, escalamiento y precisión apropiada.

Defina g_j=Re λj+σ. Para ramas simples distintas que empatan en parte real, mantener todos los g_j activos. Si hay un polo múltiple/defectivo físico, usar subespacios/funciones características y un procedimiento no suave; no dividir por un denominador casi cero.

Para z normalizado y restricciones g≤0 (incluye bounds), el problema local es:

\[
\mathcal L=J+\mu^Tg,\quad
\nabla J+J_g^T\mu=0,\quad g\le0,\quad\mu\ge0,\quad\mu\odot g=0.
\tag{36}
\]

Una implementación efectiva usa un paso SQP/active-set propio con curvatura regularizada:

\[
\min_d\ c^Td+\tfrac12d^TB_kd,
\quad g+J_gd\le0,\quad \|d\|_\infty\le\Delta_k.
\tag{37}
\]

Para un conjunto activo dado, resolver su sistema KKT; retirar multiplicadores de signo incorrecto, añadir restricciones violadas, aplicar búsqueda de paso de mérito y actualizar curvatura. Cerca del óptimo se puede usar Newton sobre las ecuaciones KKT:

\[
\begin{bmatrix}H_\mathcal L&J_A^T\\J_A&0\end{bmatrix}
\begin{bmatrix}\Delta z\\\Delta\mu\end{bmatrix}
=-\begin{bmatrix}\nabla\mathcal L\\g_A\end{bmatrix}.
\tag{38}
\]

Las fronteras (27), (31), (32) y (33) proporcionan candidatos, correcciones y pasos admisibles. El corrector (37)–(38) modifica ρ y las ganancias conjuntamente. Un pseudoinverso de mínima norma que arregla un polo no sustituye la minimización de J.

Reutilizar abiertamente ExpN como incumbent. No reintroducir una ficción de búsqueda ciega. Comprobar primal, dual, complementariedad, LICQ y condiciones de segundo orden en coordenadas escaladas y originales. Distinguir certificado numérico Float64 de prueba con inclusión intervalar. No exigir a priori que el nuevo punto sea mejor: reproducirlo por una ecuación reducida con menos coste es otro resultado válido.

## 8. Globalidad: un módulo separado

Si U es el coste de un candidato nominal factible, cualquier mejora está en

\[
\mathcal B(U)=\{\epsilon\ge0:\ w^T\epsilon<U\},
\quad \epsilon_i<U/w_i.
\tag{39}
\]

Con U≈1.124344 MW y los pesos nominales, todos los GFL están presentes. Persisten máscaras de SG y ramas desconectadas. No asumir dominancia de barra 38 ni monotonía global de Kp/Ki a partir de signos locales.

Para una certificación global se requiere una cota inferior válida L y una superior U:

\[
L\le J_{global}^*\le U,\qquad U-L\le\varepsilon_{MW}.
\tag{40}
\]

L=0 es una cota trivial honesta, no evidencia de una técnica global nueva. La factibilidad de U necesita la misma definición de incertidumbre y restricciones que el problema certificado. El U nominal no sirve para restringir el dominio robusto si no satisface robustez.

Un branch-and-bound certificable puede subdividir cajas de parámetros por soporte, usar inclusión espectral uniforme/conteo para descartar cajas incompatibles con estabilidad, y usar cotas del objetivo. Las cajas con ceros de denominadores, cambio de arquitectura o colisiones no se eliminan por fracaso de Newton. Si no se resuelven, se reportan como abiertas. Las relajaciones conservadoras internas de estabilidad no son, por sí solas, cotas inferiores de esta minimización.

## 9. GSP y Beyond Nodal Damping

### Control local con ganancias diseñadas mediante grafo

Comenzar con k_p=kp_bar·1, k_i=ki_bar·1: esa familia contiene exactamente el candidato de ExpN. Posteriormente:

\[
k_p=k_{p,r}+\Phi_da_p,\qquad k_i=k_{i,r}+\Phi_da_i.
\tag{41}
\]

Construir Φd con palabras de operadores SG,SB normalizados sobre el grafo físico y características exógenas congeladas. No llamarlos Laplacianos PSD si no cumplen esas propiedades. No usar eigenvalues negativos “corregidos” por valor absoluto. No equiparar un salto en la red Kron densa con un salto físico.

Conservar rango real y número de parámetros; una base de rango 10 sobre diez nodos no comprime. Si se usa autoridad calculada del modelo como feature, denominarlo physics-informed, no evidencia independiente de que el grafo solo predijo esa autoridad. Comparar con bases de igual dimensión y evaluar fuera de las condiciones utilizadas para construirlas.

La parametrización de ganancias es offline. K=g(S) aplicado online a errores de otros PLL sería otro controlador, con comunicación, retardos y pérdidas; no introducirlo sin una rama experimental independiente.

### Self-energy y cota de truncación

Para una partición fija:

\[
F=T_{cc}+\Gamma,\quad \Gamma=-T_{cr}T_{rr}^{-1}T_{rc}.
\tag{42}
\]

La derivada total es

\[
\Gamma_z=-T_{cr,z}RT_{rc}
+T_{cr}RT_{rr,z}RT_{rc}
-T_{cr}RT_{rc,z},\quad R=T_{rr}^{-1}.
\tag{43}
\]

Las contribuciones directa/colectiva dependen del pivote y de las coordenadas. Reportar la sensibilidad total y reconstrucción; nunca vender “100% self-energy” como porcentaje físico invariante.

Con Trr=D−W, Q=D⁻¹W, y ||Q||≤q<1:

\[
\Gamma^{(d)}=-\sum_{m=0}^{d}T_{cr}Q^mD^{-1}T_{rc},
\]
\[
\|\Gamma-\Gamma^{(d)}\|\le
\|T_{cr}\|\|D^{-1}\|\|T_{rc}\|{q^{d+1}\over1-q}.
\tag{44}
\]

Si q≥1, no aplicar esta cota. Mantener inversa exacta o ampliar bloques que retengan ciclos relevantes. Path-sums motiva resumas exactas, pero no prometer complejidad pequeña para todo grafo.

Para T=That+E, analíticos dentro de un contorno cerrado C, sin singularidades de That en C, una cota uniforme

\[
\sup_{s\in C}\|\widehat T(s)^{-1}E(s)\|<1
\tag{45}
\]

conserva el número de ceros interiores por homotopía. Para cierres meromorfos contabilizar polos o trabajar sobre el pencil completo. Un chequeo finito en C es evidencia, no una cota uniforme. Para concluir estabilidad se necesita además cubrir todos los posibles polos del semiplano relevante, incluyendo una cota de alta frecuencia/infinito.

## 10. Robustez: certificado y optimización no son lo mismo

Declare incertidumbre estática (o instantánea con norma acotada) en coordenadas fijadas:

\[
A_\Delta=A_r+B_w\Delta C_w,\qquad\|\Delta\|_2\le\beta.
\tag{46}
\]

Bw,Cw y las unidades se congelan; no normalizar cada candidato para mejorar artificialmente β. Si la perturbación física requiere una LFT más general, construirla y declarar sus hipótesis. No llamar porcentaje a β normalizado. No aplicar automáticamente el desplazamiento espectral a una incertidumbre dinámica cuya norma exponencialmente ponderada no ha sido especificada.

Sea Aσ=Ar+σI. Para D=0:

\[
G_\sigma(s)=C_w(sI-A_\sigma)^{-1}B_w,
\quad A_\sigma\ \text{Hurwitz},\quad\beta\|G_\sigma\|_\infty<1.
\tag{47}
\]

Un certificado suficiente es X=Xᵀ>0 tal que

\[
\begin{bmatrix}
A_\sigma^TX+XA_\sigma+C_w^TC_w&XB_w\\
B_w^TX&-\gamma^2I
\end{bmatrix}\prec0,\qquad\beta\gamma<1.
\tag{48}
\]

Para D≠0 usar la versión completa bounded-real y comprobar well-posedness. Un solver SDP puede usarse para CERTIFICAR un Z fijo: no se presenta como el método primario de diseño. Validar residuos, escalas y holguras; una prueba formal requiere enclosures, no únicamente éxito del solver.

El máximo sobre muestras de ω es una cota INFERIOR de ||G||∞. Su inverso es una cota SUPERIOR de la tolerancia, no un radio certificado. Obtener γL≤||G||∞≤γU con un test bounded-real/Hamiltoniano/inclusión validada. Reportar β_cert=1/γU.

Para pico simple aislado:

\[
\partial_z\bar\sigma(G(j\omega_*))=
\Re[u^H G_z(j\omega_*)v].
\tag{49}
\]

Si empatan frecuencias o singular values, mantener restricciones/subgradientes múltiples. Añadir gR=βγ−1 al co-diseño; no asumir que robustez sale gratuitamente del óptimo nominal.

## 11. Transitorios, RoCoF, frecuencia, corrientes y límites de la afirmación

Fijar salida: COI SG, frecuencia de tensión de barras y/o PLL, con sus filtros y unidades. No asumir que son iguales. Para Ar estable con gauge eliminado y salida sin feedthrough:

\[
\dot x=A_rx+B_\ell\Delta P\,u(t),\quad f=C_fx,
\]
\[
R_\ell(t)=\Delta P C_fe^{A_rt}B_\ell,\quad
f_\ell(t)=\Delta P C_f\int_0^t e^{A_r\tau}B_\ell d\tau.
\tag{50}
\]

Un pulso rectangular de duración Tp es la diferencia de dos respuestas al escalón. No extrapolar capacidad de pulso como capacidad de escalón sostenido. Evaluar la integral con exponencial aumentada, evitando Ar⁻¹ cuando sea mal condicionada.

Si Df≠0, un escalón puede producir un salto; su derivada ideal incluye una distribución. Modelar el filtro/ventana de medición de RoCoF antes de atribuir un pico finito.

Defina para el conjunto de ubicaciones y perfiles declarado:

\[
G_R=\sup_{\ell,t\ge0}|R_\ell(t;1\mathrm{MW})|,\qquad
G_f=\sup_{\ell,t\ge0}|f_\ell(t;1\mathrm{MW})|.
\]

\[
\Delta P_{\max}^{lin}=\min(R_{max}/G_R,F_{max}/G_f).
\tag{51}
\]

Esta igualdad es del modelo lineal y perfiles escogidos; una cota superior de GR,Gf produce una capacidad garantizada conservadora para ese modelo. No constituye una garantía no lineal/EMT/current-limit.

Con polos simples y expansión modal completa:

\[
q_{j\ell}={(C_fr_j)(l_j^HB_\ell)\over l_j^Hr_j},\quad
R_\ell(t)=\Delta P\sum_jq_{j\ell}e^{\lambda_jt}.
\tag{52}
\]

Para polos repetidos usar bloques Schur/Jordan/exponenciales completos; no dividir por un condicionamiento degenerado. El polo más a la derecha no es necesariamente el que domina RoCoF. Para pares complejos ζj=−Reλj/|λj|; el polo real activo de ExpN no es una oscilación.

Una cota transitoria alternativa: si AᵀX+XA≤−2aX, X>0,

\[
|C_fe^{At}B_\ell|\le c_\ell e^{-at},\quad
c_\ell=\|C_fX^{-1/2}\|\|X^{1/2}B_\ell\|.
\tag{53}
\]

Da Rpeak≤|ΔP|cℓ y fpeak≤|ΔP|cℓ/a. Es conservadora; si X es común para todo el conjunto de incertidumbre, la cota se extiende a esa familia bajo las hipótesis declaradas. Para cerrar un horizonte finito T, usar también el estado x(T) y la contracción/tail para excluir un pico posterior.

No imponer automáticamente f0ΔP/(2Σ εHSn) como RoCoF de este modelo. Esa fórmula necesita hipótesis adicionales sobre qué almacenes de energía responden inicialmente y sobre carga/voltaje. Comparar primero con el límite exacto Cf Bℓ y con los estados del GFL. Las cargas Z y cambios algebraicos pueden invalidar una identificación ingenua del desequilibrio con potencia instantánea absorbida por SG.

El modelo instalado no tiene un limitador duro GFL validado. Se puede medir corriente y headroom con rating conocido, pero no afirmar seguridad ante saturaciones no incluidas. Añadir un limitador cambia el modelo y exige otra validación.

## 12. Secuencia incremental propuesta

P0: reproducir ExpN sin cambios. P1: validar rango uno y fronteras condicionadas. P2: recuperar el candidato de un SG mediante la raíz de rango dos y contar todo el espectro. P3: co-diseño nominal continuo, uniforme/graph/libre, con KKT real y presupuesto computacional. P4: certificado de robustez y co-diseño para β declarado; añadir transitorios obligatorios como restricciones. P5: validación PD y TDS de candidatos congelados. P6: mejora global por cotas dentro del presupuesto del incumbent, y reducción GSP certificada, como módulos separados que no bloquean resultados ya válidos.

En cada fase guardar entradas/hashes, ecuaciones usadas, pruebas que pueden fallar, coste real, resultados y alcance. No nuevos barridos de 1024 soportes como sustituto del optimizador. No imponer bus 38, la esquina de ganancias ni mejoras numéricas. No bajar tolerancias después de ver resultados sin crear otra versión de la especificación.

## 13. Referencias y procedencia

**Fuente experimental suministrada:** REPORT_EXP_N.md, archivo de conversación “Texto pegado(10).txt”. Los resultados de ExpN son reportados, no reejecutados aquí.

**Fuentes de ecuaciones consultadas (tag explícito):**
- https://github.com/JuliaEnergy/PowerDynamics.jl/blob/v5.0.0/src/Library/Machines/SauerPaiMachine.jl
- https://github.com/JuliaEnergy/PowerDynamics.jl/blob/v5.0.0/src/Library/Controls/AVRs.jl
- https://github.com/JuliaEnergy/PowerDynamics.jl/blob/v5.0.0/src/Library/Controls/Govs.jl
- https://github.com/JuliaEnergy/PowerDynamics.jl/blob/v5.0.0/src/Library/Renewables/ComposableInverter.jl
- https://github.com/JuliaEnergy/PowerDynamics.jl/blob/v5.0.0/src/Library/building_blocks.jl
- https://juliadynamics.github.io/NetworkDynamics.jl/dev/API/ (solo referencia; verificar versión instalada).

**Herramientas matemáticas externas, no prueba de novedad del proyecto:**
- Henrion y Šebek, Plane geometry and convexity of polynomial stability regions, arXiv:0801.2499.
- Bindel y Hood, Localization theorems for nonlinear eigenvalue problems, arXiv:1303.4668.
- Giscard, Thwaite y Jaksch, Evaluating Matrix Functions by Resummations on Graphs, arXiv:1112.1588.
- Butler, Parada-Mayorga y Ribeiro, Convolutional Learning on Multigraphs, arXiv:2209.11354.
- Boyd et al., Linear Matrix Inequalities in System and Control Theory, https://web.stanford.edu/~boyd/lmibook/.

Las aplicaciones de fronteras condicionadas, la combinación con el reparto P/Q, los tests GSP y la secuencia de co-diseño aquí especificados son propuestas a validar. No se atribuye a esos artículos haber resuelto nuestro óptimo SG→GFL.
