# Beyond Nodal Damping — Ajuste mínimo de PLL con márgenes colectivos

Fecha: 5 de octubre de 2026.
Estado: formulación y demostraciones en papel; NO es una nueva campaña experimental.

## 0. Objeto científico

A un reemplazo SG→GFL prescrito, encontrar el menor número de PLL cuyos pares (Kp, Ki) deben modificarse para recuperar un contrato de margen espectral y comportamiento temporal. Mantener los PLL implementados localmente; usar el modelo completo durante el diseño offline. La sparsity se aplica a los CAMBIOS de ganancias, no a apagar PLL ni a añadir enlaces de comunicación.

No se afirma que esta formulación sea un teorema de optimalidad global de reemplazo. Tampoco se atribuye prioridad bibliográfica a Schur, sensibilidad modal, funciones soporte, sparsity por grupos o asignación de polos.

### Procedencia

Base del proyecto: `Beyond_Nodal_Damping_Network_Moment_Laws_20261004.pdf`, secciones 2, 3.5, 8, 9 y apéndices B, G, I. El apéndice B contiene el resultado de delay finito; G el retorno de un PLL y la asignación de un polo; I el retorno all-PLL. Las proposiciones sobre cambios de ganancias y soporte ordenado que siguen se derivan aquí bajo hipótesis explícitas; su utilidad en IEEE-39 no está demostrada todavía.

Antecedentes primarios a comparar:
- Dörfler, Jovanović, Chertkov y Bullo, *Sparsity-Promoting Optimal Wide-Area Control of Power Networks*, arXiv:1307.4342; IEEE TPWRS 2014.
- Lin, Fardad y Jovanović, *Design of Optimal Sparse Feedback Gains via the Alternating Direction Method of Multipliers*, arXiv:1111.6188; IEEE TAC 2013.
- Huang, Xin, Dong y Dörfler, *Impacts of Grid Structure on PLL-Synchronization Stability of Converter-Integrated Power Systems*, arXiv:1903.05489.
- Michiels y Gumussoy, *Eigenvalue based algorithms and software for the design of fixed-order stabilizing controllers for interconnected systems with time-delays*, arXiv:2003.05496; capítulo Springer 2014.
- Gumussoy y Michiels, *Fixed-order strong H-infinity control of interconnected systems with time-delays*, arXiv:2003.05728.

## 1. Contrato físico

En un soporte de estados fijo:

\[
E\dot x=f_0(x,z;\rho)+\sum_i b_i(k_i)e_i(t-\tau_i),\qquad 0=g(x,z;\rho,d),
\]

\[
\dot\theta_i=\omega_{p,i},\quad t_{f,i}\dot\omega_{p,i}=-\omega_{p,i}+\xi_i+k_{p,i}e_i(t-\tau_i),\quad \dot\xi_i=k_{I,i}e_i(t-\tau_i).
\]

Aquí \(k_i=(k_{p,i},k_{I,i})^T\). Retener red AC con pérdidas, corrientes, DC, P/Q y SG de la implementación canónica. No sustituir el Jacobiano por un Laplaciano lossless.

Primera campaña: rho, despacho, planta y punto de equilibrio prescritos. Ganancias únicamente en las ecuaciones PLL declaradas. Latencias exógenas: no se optimizan para fabricar una mejora. No afirmar que 40 ms sea el delay digital típico del PLL. Identificar dónde se aplica el delay; un retardo de detector no es intercambiable con PWM, filtro o retardo de una señal remota.

Con reducción algebraica regular:

\[
\Delta(s;k,\vartheta)=sE-A_0-\sum_i b_i(k_i)c_i^T e^{-s\tau_i},\qquad b_i=b_{p,i}k_{p,i}+b_{I,i}k_{I,i}.
\]

En el modelo exportado E=I. Si A0, E o c dependen de las ganancias, utilizar sus derivadas completas: las identidades de rango uno siguientes dejan de aplicarse sin modificación. Cambiar rho o topología puede mover el equilibrio y requiere reensamblaje y derivadas totales.

## 2. Problema exacto de diseño

Elegir escalas positivas fijas \(S_i=\operatorname{diag}(s_{p,i},s_{I,i})\). Definir

\[
k_i=k_i^0+S_i u_i,\qquad \operatorname{supp}_B(u)=\{i:u_i\ne0\}.
\]

Las escalas se registran ANTES de comparar; coste de ganancia no equivale a energía de actuación. Definir \(\mathcal U_i\) como límites físicos normalizados y presupuesto de ajuste. Todos los conjuntos contienen cero.

\[
\min_u^{\mathrm{lex}}\left(\#\operatorname{supp}_B(u),\frac12\sum_i\|u_i\|_2^2\right)
\]

sujeto a ganancias admisibles y, para escenarios \(\vartheta\in\mathcal V\),

\[
\alpha_\perp(k^0+Su;\vartheta)\le-\sigma_*.
\]

Imponer también los límites temporales prescritos para cada evento. Si el objetivo incluye desempeño entrada-salida, añadir

\[
\sup_{\vartheta\in\mathcal V}\|\mathcal T_{zd}(s;k,\vartheta)\|_\infty\le\gamma_*,
\quad
\mathcal T_{zd}=W_z(C_o\Delta^{-1}B_d+D_{od})W_d.
\]

Definir puertos, pesos, unidades y filtrado; tratar apropiadamente el modo rotacional. Una malla de frecuencias no certifica la norma H-infinity. Una colección finita de escenarios no certifica un continuo de incertidumbre. Para el póster, el núcleo mínimo es margen de polos + eventos; H-infinity es una extensión registrada separadamente.

## 3. El puente físico de BND

Retener velocidades SG nu y eliminar estados h donde el bloque oculto es invertible. En coordenadas exportadas E=I:

\[
Z(s)=M\left(\Delta_{\nu\nu}-\Delta_{\nu h}\Delta_{hh}^{-1}\Delta_{h\nu}\right),\quad d=Z(s)\nu.
\]

\[
D_H(\omega)=\frac{Z(j\omega)+Z(j\omega)^H}{2},\qquad
\langle d^T\nu\rangle=\frac12\hat\nu^HD_H\hat\nu.
\]

Es suministro mecánico incremental en los puertos declarados, no una función de Lyapunov no lineal ni la potencia AC total. No tomar la parte Hermitiana de cualquier matriz característica y llamarla damping.

### Proposición A: un retuning local también admite actualización de rango uno

A rho, equilibrio y tau fijos, modificar solo el PLL i. Definir

\[
p_i(s)=-e^{-s\tau_i}(b_{p,i}\Delta k_{p,i}+b_{I,i}\Delta k_{I,i}),\qquad
\Delta_{new}=\Delta_0+p_i c_i^T.
\]

Sea \(D=\Delta_{hh,0}\). Particionar p y c en nu,h. Si D y \(D+p_hc_h^T\) son invertibles:

\[
a_i=M\frac{p_\nu-\Delta_{\nu h,0}D^{-1}p_h}{1+c_h^TD^{-1}p_h},\qquad
v_i^H=c_\nu^T-c_h^TD^{-1}\Delta_{h\nu,0}.
\]

Entonces

\[
\delta Z=a_i v_i^H,\qquad \delta D_H=\tfrac12(a_i v_i^H+v_i a_i^H).
\]

Demostración: sustituir Sherman–Morrison para \((D+p_hc_h^T)^{-1}\) en ambos complementos de Schur y agrupar. No usar esta identidad si el cambio de ganancia modifica canales fuera del detector declarado.

A una frecuencia fija, si a y v son linealmente independientes,

\[
\ell_\pm=\frac{\Re(v^Ha)\pm\sqrt{\|a\|^2\|v\|^2-[\Im(v^Ha)]^2}}{2},
\quad
\ell_+\ell_-=-\tfrac14(\|a\|^2\|v\|^2-|v^Ha|^2)<0.
\]

Luego el incremento tiene firma (+,-,0,...,0). Si son colineales puede tener rango uno, cero o signo diferente. El rango dos es del INCREMENTO, no del operador total ni del número de modos dinámicos. No garantiza que algún polo mejore o empeore.

La misma construcción cubre un cambio conjunto de ganancias y delay con c fijo sustituyendo
\(p_i=e^{-s\tau_i^0}b_i^0-e^{-s\tau_i^1}b_i^1\).

## 4. Reducción exacta por soporte: m PLL, determinante m por m

Para soporte S de cardinal m, apilar \(P_S=[p_i]_{i\in S}\) y \(C_S^T=[c_i^T]_{i\in S}\). Entonces

\[
\Delta_u=\Delta_0+P_SC_S^T,
\qquad
\frac{\det\Delta_u}{\det\Delta_0}=\det\left(I_m+C_S^T\Delta_0^{-1}P_S\right).
\]

Es exacta donde \(\Delta_0\) sea invertible. Con un contorno común regular, el winding del cociente da el cambio de conteo de raíces. Los polos de la referencia NO se ignoran. Una reducción m por m no reduce el número de estados físicos a m ni convierte un DDE en un sistema con m polos.

Cada evaluación necesita solves de la planta completa; éstos pueden reutilizarse por frecuencia y configuración. Medir aceleración, no inferirla del tamaño del determinante. Con un solo PLL, dos ganancias reales comparten UN detector, por eso la actualización tiene rango a lo sumo uno.

## 5. Autoridad modal exacta

Para una raíz simple \(\lambda_r\) y vectores izquierdo/derecho \(\ell_r,v_r\):

\[
\lambda_{r,p}=-\frac{\ell_r^H\Delta_p(\lambda_r)v_r}{\ell_r^H\Delta_s(\lambda_r)v_r}.
\]

En el contrato anterior:

\[
\Delta_s=E+\sum_i\tau_i b_ic_i^T e^{-s\tau_i}.
\]

No omitir el término de delay del denominador. Definir

\[
g_{ri}=\begin{bmatrix}\Re\lambda_{r,k_{p,i}}\\\Re\lambda_{r,k_{I,i}}\end{bmatrix},\qquad h_{ri}=S_i^Tg_{ri}.
\]

Con el retorno PLL del modelo canónico:

\[
\lambda_{r,k_{p,i}}=\lambda_r\lambda_{r,k_{I,i}},
\qquad
g_{ri}=\begin{bmatrix}\Re(\lambda_r d_{ri})\\\Re d_{ri}\end{bmatrix},\quad d_{ri}=\lambda_{r,k_{I,i}}.
\]

El residuo colectivo \(d_{ri}\) incluye red, otros controladores y fases. Schur exacto debe dar el MISMO gradiente que la formulación completa. Una diferencia persistente es un bug/model mismatch, no una ventaja BND.

## 6. Proposición B: predictor sparse analítico para una restricción

Una raíz viola el objetivo por \(d=\Re\lambda_0+\sigma_*>0\). A primer orden:

\[
\sum_i h_i^Tu_i\le-d.
\]

Si se permite actuar solo en i, el ajuste de norma mínima sin límites activos es

\[
u_i^*=-\frac{d}{\|h_i\|_2^2}h_i,\qquad
\|u_i^*\|_2=\frac{d}{\|h_i\|_2}.
\]

Prueba: Cauchy–Schwarz da \(d\le\|h_i\|\|u_i\|\); la igualdad se obtiene en la dirección opuesta a h. Si h_i=0 y d>0, este predictor no ofrece autoridad. La fórmula no ignora bounds: debe pertenecer a \(\mathcal U_i\).

Para bolas independientes \(\|u_i\|\le r_i\), definir \(A_i=r_i\|h_i\|\). Un soporte S satisface esta única restricción lineal si y solo si

\[
\sum_{i\in S}A_i\ge d.
\]

Ordenando \(A_{(1)}\ge\cdots\ge A_{(n)}\):

\[
m_{lin}=\min\{m:\sum_{j=1}^{m}A_{(j)}\ge d\}.
\]

Es mínimo soporte del MODELO LINEAL DE UNA RESTRICCIÓN, no del DDE no lineal. Con otras formas de bounds usar la función soporte exacta \(A_i=h_{\mathcal U_i}(-h_i)\). Un presupuesto global adicional puede impedir suficiencia; la suma continúa siendo cota superior de autoridad y sirve para exclusión.

## 7. Corrector analítico de un PLL, sin congelar el patrón global

Abrir únicamente el PLL candidato i y conservar toda la red restante. El retorno es

\[
e_i=G_i(s;\rho,k_{-i})\theta_i,
\quad
F_i=s^2(1+t_{f,i}s)-(sk_{p,i}+k_{I,i})e^{-s\tau_i}G_i.
\]

Para \(\lambda_*=-\sigma_*+j\omega_*\), \(\omega_*>0\), y retorno regular no nulo:

\[
W_i=\frac{\lambda_*^2(1+t_{f,i}\lambda_*)e^{\lambda_*\tau_i}}{G_i(\lambda_*;\rho,k_{-i})},
\quad
k_{p,i}=\frac{\Im W_i}{\omega_*},\quad
k_{I,i}=\Re W_i+\sigma_*k_{p,i}.
\]

Es exacto para imponer ese polo; no para asegurar dominancia, todos los polos o eventos. Es un resultado existente del apéndice G, no una nueva prueba experimental. Proponer una frecuencia mediante el predictor; si se busca optimizarla, registrar una búsqueda unidimensional y no llamarla solución cerrada global. Fallar para una frecuencia no excluye otras frecuencias, otros polos o componentes factibles.

No confundir esta construcción single-PLL con el mapa all-PLL que fija también un patrón q y puede cambiar veinte ganancias. La antigua campaña 30–37 preservaba solo la parte real y prescribía Kp; es otro contrato de diseño.

## 8. Varios modos y escenarios: subproblema convexo

Apilar las restricciones \(g_r=\Re\lambda_r+\sigma_*\), o medias de clusters como restricciones NECESARIAS con conteos conservados. Incluir modos sanos cercanos al guard, no únicamente el objetivo.

\[
\min_u\ \frac12\|u\|_2^2+\gamma\sum_i w_i\|u_i\|_2
\quad\text{s.a.}\quad g^0+Hu\le0,\quad u_i\in\mathcal U_i,\quad\|u\|\le r.
\]

Es un programa convexo con normas de segundo orden (SOCP con epígrafes), no un QP puro salvo reformulación/caso especial. w_i son costes prescritos o unos. No multiplicar la sensibilidad por un índice arbitrario de ciclos y atribuirle la optimalidad de esta derivación.

Fijar el soporte elegido y refinar con el modelo exacto. La penalización de grupos propone sparsity; NO demuestra cardinalidad mínima global. Gestionar infeasibilidad local mediante restauración, sin aceptar violaciones del problema físico.

## 9. Proposición C: certificado finito de imposibilidad con m PLL

Esta es la extensión de interés al certificado regional del manuscrito.

Supóngase que en TODO el dominio declarado se conoce

\[
g(u)=g^0+Hu+r(u),\qquad |r_j(u)|\le\beta_j.
\]

Las restricciones g usadas deben ser necesarias para el contrato real. Pueden ser raíces simples continuables sin colisiones; alternativamente medias de clusters con contornos/conteos verificados. Una media negativa NO garantiza que todas las raíces sean seguras.

Para \(\eta\ge0\), definir la autoridad de cada PLL sobre la combinación de restricciones:

\[
A_i(\eta)=\max_{u_i\in\mathcal U_i}[-\eta^TH_i u_i].
\]

Para bolas, \(A_i(\eta)=r_i\|H_i^T\eta\|_2\). Ordenar los A de mayor a menor. Si

\[
\eta^T(g^0-\beta)>\sum_{j=1}^{m}A_{(j)}(\eta),
\]

NINGÚN ajuste de a lo sumo m PLL dentro de ese dominio puede satisfacer el contrato.

Prueba: una solución física cumple \(g^0+Hu\le\beta\). Multiplicar por eta da \(\eta^T(g^0-\beta)\le-\eta^THu\). Un soporte S de cardinal a lo sumo m puede ofrecer a lo sumo \(\sum_{i\in S}A_i\), que no supera la suma de los m mayores. Contradicción.

La desigualdad es homogénea: puede normalizarse \(\mathbf1^T\eta=1\). Un certificado hallado con optimización flotante necesita verificación con errores. La ausencia de certificado NO establece factibilidad.

Las cotas beta deben ser uniformes y rigurosas. Por ejemplo, si \(\|\nabla^2g_j(u)\|_2\le L_j\) en una bola \(\|u\|\le r\), se permite \(\beta_j=L_jr^2/2\), incorporando errores de referencia/gradiente. Hessianas muestreadas no certifican ese supuesto. Cerca de raíces múltiples abandonar el modelo de raíz simple o subdividir/cambiar a clusters.

Si se excluyen todos los soportes de cardinal <=m-1 y se valida uno de cardinal m, se obtiene cardinalidad mínima dentro del MISMO dominio y contrato, al nivel de certificación alcanzado. Una solución nominal no basta para cerrar un upper bound robusto continuo.

Caso fácil: si el diseño sin retuning incumple y un único PLL satisface el contrato completo, la cardinalidad mínima es uno; no hace falta excluir todos los otros PLL para probarlo. Ello no prueba mínimo esfuerzo ni unicidad.

## 10. Robustez de delay: distinguir certificar de muestrear

Los escenarios de delay finitos son evidencia. Para un certificado continuo, usar una incertidumbre explícita que contenga la referencia y todos los segmentos de homotopía.

Con \(\Delta_*(s)\) del diseño candidato y \(\mathcal E(s,\vartheta)=\Delta(s;\vartheta)-\Delta_*(s)\), una condición suficiente sobre el contorno que encierra la región peligrosa es

\[
\sup_{s\in\Gamma,\vartheta\in\mathcal V}\|\Delta_*^{-1}(s)\mathcal E(s,\vartheta)\|_2<1.
\]

Debe existir conteo inicial cero, regularidad del modelo/gauge y cota de cola. Para un retarded DDE reducido E=I, con tau no negativas y objetivo Re(s)>=-sigma:

\[
R>\|A_0\|_2+\sum_i\|b_i\|_2\|c_i\|_2 e^{\sigma\tau_{i,max}}
\]

excluye raíces peligrosas con |s|>R. No aplicar esta cola automáticamente a un descriptor con restricciones retardadas/estructura neutral. Una malla no demuestra el supremo; emplear intervalos o cotas de variación con redondeo controlado. El fallo de la condición suficiente no prueba inestabilidad.

No tratar sigma_min(I-Rc) en puertos arbitrariamente escalados como un certificado de robustez físico. Se requieren incertidumbre/normalización, conteos y frecuencia/cola completos.

## 11. Papel de GSP

Con Laplaciano físico simétrico verificado y M positivo sobre las velocidades retenidas:

\[
L_M=M^{-1/2}L_PM^{-1/2}=U\Lambda U^T,\quad
\widehat D=U^TM^{-1/2}D_HM^{-1/2}U.
\]

Usar esta base para explicar el alcance espacial de acciones y para probar aproximaciones. Retener toda la base no cambia las decisiones exactas ni las sensibilidades. No identificar los autovectores de D_H con los modos completos del DDE.

Un conmutador no nulo es una señal de falta de invariancia de subespacios, NO un criterio de inestabilidad. Un damping diagonal nodal heterogéneo puede mezclar modos del grafo; un damping no nodal puede ser beneficioso. Minimizar offdiagonal o conmutador sin objetivo de desempeño no es una ley de robustez.

Si se trunca la base, el resultado interesante es preservar una DECISIÓN (soporte suficiente/insuficiente, margen, ranking) con error acotado; no únicamente producir una figura de Fourier.

## 12. Qué pertenece ya a los archivos y qué falta

El caso registrado 30–37 tiene raíz conjunta -0.04139961 /s y guard -0.05 /s. Una corrección de dos Ki a Kp fijo alcanza -0.06000000 /s manteniendo 7.9 MW adicionales. Es un resultado reportado en el manuscrito, no reejecutado aquí. La raíz conjunta viola margen pero no es inestable por su signo. La campaña anterior también contiene un comparador complex-pole que satisface margen; no ocultarlo.

De esos números, el desplazamiento mínimo de parte real hasta el guard es 0.00860039 /s; hasta el objetivo -0.06 es 0.01860039 /s. Son aritmética de los resultados existentes, no predicciones nuevas.

Todavía NO están establecidos: reparación single-PLL del caso, soporte mínimo dos, ventaja de esfuerzo frente a un optimizador sparse del modelo completo, robustez continua de la nueva síntesis, validación EMT/fabricante ni óptimo de rho.

## 13. Resultado de póster que buscar

Afirmación candidata fuerte: “El cierre dinámico de la red permite seleccionar y verificar el mínimo conjunto de PLL que debe retunearse para recuperar el margen, sin reducir el reemplazo y sin añadir comunicación”.

El resultado ideal consta de m_L <= m* <= m_U y, cuando cierre, m_L=m_U. Complementarlo con coste normalizado, margen espectral, incertidumbre de delay declarada y eventos. Si solo hay m_U, decir “soporte factible de m_U PLL”. Si solo hay una reducción de coste, reportar ese resultado y no inventar minimalidad.

La comparación contra sensibilidad y optimización sparse del MODELO COMPLETO es obligatoria: Schur es una representación equivalente, no información física adicional. Puede ofrecer explicación, reutilización computacional, restricciones exactas por soporte o certificados; no superioridad mágica.
