from pathlib import Path
import hashlib
import json
import tomllib

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports/nonlinear_codesign_20261001/full_gfl_sector_contract"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda name:json.loads((OUT/name).read_text(encoding="utf-8"))
s=read("synthesis_final.json");v=read("verification_final.json")
n=read("nonlinear_validation.json");w=read("network_closure_witness.json")
j=tomllib.loads((OUT/"julia_mapping_audit.toml").read_text(encoding="utf-8"))
jn=tomllib.loads((OUT/"julia_network_witness_audit.toml").read_text(encoding="utf-8"))
models=tomllib.loads((OUT/"models.toml").read_text(encoding="utf-8"))
assert v["status"]=="ALL_LOCAL_CERTIFICATES_VERIFIED"
assert s["models_sha256"]==sha(OUT/"models.toml")
assert v["synthesis_sha256"]==n["synthesis_sha256"]==w["synthesis_sha256"]==sha(OUT/"synthesis_final.json")
assert n["verification_sha256"]==w["verification_sha256"]==sha(OUT/"verification_final.json")
assert jn["witness_csv_sha256"]==sha(OUT/"closure_witnesses.csv")
assert jn["joint_witness_csv_sha256"]==sha(OUT/"joint_closure_witness.csv")
assert abs(jn["joint_input_usage"]-w["joint_witness"]["joint_input_usage"])<1e-9
for rel,digest in models["source_hashes"].items():assert sha(ROOT/rel)==digest
for name,digest in [("full_gfl_sector_contract.py",s["source_sha256"]),
                    ("verify_full_gfl_sector_contract.py",v["script_sha256"]),
                    ("validate_full_gfl_sector_contract.py",n["script_sha256"]),
                    ("audit_full_gfl_mapping.jl",j["script_sha256"]),
                    ("gfl_network_closure_witness.py",w["script_sha256"]),
                    ("audit_gfl_network_witness.jl",jn["script_sha256"])]:
    assert sha(ROOT/"experiments/nonlinear_codesign_20261001"/name)==digest,name
cert=v["certificates"];joint=w["joint_witness"];p=s["specification"]
rows="\n".join(f"  {c['bus']:2d}  {c['Kp']:12.7f}  {c['Ki']:12.7f}  {c['epsilon_actual_model_verified']:.10f}  {100*c['ideal_fixed_gain_relative_gap']:7.2f}%" for c in cert)
report=f"""GFL COMPLETO: CERTIFICADOS LOCALES Y OBSTRUCCIÓN DE SU COMPOSICIÓN
2026-10-01. Avance hacia el co-diseño no lineal SG -> GFL.

Este informe documenta la etapa de contratos locales. El avance posterior sobre
un certificado conjunto microscópico se encuentra en
../coupled_network_contract/CERTIFICADO_ACOPLADO_Y_LIMITES_ES.txt.

RESULTADO ACTUAL

Se verificaron contratos de invariancia para los nueve estados físicos/de control
de cada uno de los diez GFL, más un estado de medida de RoCoF. Cada contrato cubre
toda señal medible dentro de una envolvente conjunta de tensión, frecuencia de
referencia y potencia DC. Se conservaron el seno, los productos de corriente,
la dinámica DC y el efecto de la frecuencia sobre el marco de corriente.

La representación usa cinco parámetros sectoriales y 32 vértices por dispositivo.
La verificación racional de {32*len(cert)} desigualdades matriciales y de sus cotas
de salida fue exitosa. Un margen adicional cubre el residuo de equilibrio producido
por los valores numéricos congelados del modelo, sin suponer que ese residuo es cero.

Estas son garantías LOCALES CONDICIONADAS A SUS ENTRADAS. Un cálculo independiente
demuestra que el producto de estas regiones no conserva sus hipótesis de entrada
al conectarlas a la red del candidato actual. Por tanto, aún no hay certificado
de seguridad de toda la red ni capacidad óptima de sustitución demostrada.

1. COORDENADAS QUE CONSERVAN EL CONVERTIDOR COMPLETO

Sea delta=theta_PLL-theta_ref. Conservamos la frecuencia física omega_PLL.
Los estados incrementales, en el orden usado por el código, son
  z=[delta,omega,xi,chi,Delta gamma_q,Delta gamma_d,
     Delta i_d,Delta i_q,Delta v_dc_i,Delta v_dc]^T.
Los dos últimos son integrador del regulador DC y tensión DC; no se suprime
ninguno. Las corrientes d/q giran con el PLL y theta_ref sólo es una coordenada
de análisis, no una señal remota añadida al controlador.

Si u_ref=[V0,0]+Delta u_ref, el cambio
  Delta u_PLL=R(-delta) Delta u_ref
preserva exactamente la norma de la entrada. Así,
  u_d=V0 cos(delta)+Delta u_PLL,d,
  u_q=-V0 sin(delta)+Delta u_PLL,q.
La envolvente general es
  ||Delta u_ref||^2/epsilon^2
  +nu^2/(kappa_f epsilon)^2
  +Delta Pdc^2/(kappa_P epsilon)^2 <= 1,
donde nu=theta_ref_dot. Las ponderaciones y límites son parámetros; no se fija
una forma de evento, soporte nodal o duración. La fuente DC sigue siendo parte
del modelo declarado, no una reserva energética inferida de los Kp,Ki del PLL.

Los controles de corriente son
  e_d=k_dc,p Delta vdc+Delta v_dc_i-Delta i_d,
  e_q=-Delta i_q,
  Delta v_inv,d=k_cc,p e_d+k_cc,i Delta gamma_d,
  Delta v_inv,q=k_cc,p e_q+k_cc,i Delta gamma_q.
Estas expresiones suponen trim ideal; los pequeños términos constantes reales
se incluyen en la cota de residuo descrita en §3.

La transformación exacta de corriente da
  i_d_dot=(omega_base/Xf)(v_inv,d-u_d-Rf i_d)
                              +(omega_base omega_frame+omega) i_q,
  i_q_dot=(omega_base/Xf)(v_inv,q-u_q-Rf i_q)
                              -(omega_base omega_frame+omega) i_d.
La ecuación DC conserva
  Cdc vdc vdc_dot=Pdc+Delta Pdc-v_inv,d i_d-v_inv,q i_q.
En particular, para E=Cdc vdc^2/2+Xf(i_d^2+i_q^2)/(2 omega_base),
  E_dot=Pdc+Delta Pdc-u_d i_d-u_q i_q-Rf(i_d^2+i_q^2).
Los términos de rotación se cancelan en la energía. No se cambiaron signos de
potencia ni se retiraron estados para obtener el certificado.

2. REPRESENTACIÓN SECTORIAL EXACTA, SIN TRUNCAR LA PLANTA

En |delta|<=delta_bar y |Delta vdc|<=v_bar<Vdc0, definimos
  q1=sin(delta)/delta,
  q2=(cos(delta)-1)/delta,
  q3=Delta i_d, q4=Delta i_q,
  q5=1/(Vdc0+Delta vdc),
con los valores continuos en delta=0. Son válidas las cotas
  1-delta_bar^2/6 <= q1 <= 1,
  -delta_bar/2 <= q2 <= delta_bar/2,
  1/(Vdc0+v_bar) <= q5 <= 1/(Vdc0-v_bar),
y las cotas declaradas de q3,q4.

Las ecuaciones incrementales tienen la forma exacta
  zdot=A(q;Kp,Ki)z+epsilon B(q)eta+r(z),  ||eta||<=1.
Por ejemplo, el producto omega Delta i_q se conserva como un coeficiente q4
de la columna omega, y la fila DC contiene
  -q5/Cdc [h_linear+q3 c_d+q4 c_q],
donde Delta v_inv,d=c_d z y Delta v_inv,q=c_q z.
La potencia DC perturbada entra mediante q5 kappa_P/Cdc.

A y B son multiafines en los cinco q: afines en cada uno cuando los demás
se fijan. Por interpolación multilineal pertenecen EXACTAMENTE a la envolvente
convexa de los 2^5 vértices. Se ignoran correlaciones entre q para obtener una
condición suficiente conservadora, pero no se aproxima sin(delta) por delta,
no se omiten términos cuadráticos y no se reemplaza el convertidor por su Jacobiano.

3. INVARIANCIA Y TRATAMIENTO DEL TRIM NUMÉRICO

Con y=S^-1 z, alpha>0 y V=y^T X^-1 y, el SDP a ganancias fijas busca X>0,
epsilon y las mismas desigualdades de disipatividad usadas en el contrato PLL:
  [[Atilde_v X+X Atilde_v^T+2 alphatilde X, epsilon Btilde_v],
   [epsilon Btilde_v^T, -2 alphatilde I]] <= -g I
para TODOS los 32 vértices, más contención de salidas mediante C_j X C_j^T<=1-g.
Atilde=t0 S^-1 A S, Btilde=t0 S^-1 B, alphatilde=t0 alpha.
Las salidas aseguran el dominio angular, los incrementos d/q, vdc positivo,
frecuencia del PLL y su RoCoF filtrado. También se imponen X>=gI, trace(X)<=T
y epsilon<=E para definir la familia de búsqueda y acotar residuos duales.

La congruencia con diag(X^-1,I) da, en el modelo incremental ideal,
  Vdot<=-2 alpha V+2 alpha||eta||^2.
La contención en el dominio que define q cierra la prueba sectorial por invariancia.

En el modelo congelado el trim no es matemáticamente cero. El verificador
reconstruye r_gamma, r_current y r_DC con fracciones de los parámetros guardados.
Por ejemplo, si los voltajes reales del controlador en el origen son vhat_d0,
vhat_q0, el residuo de la fila DC se acota uniformemente usando
  (|Pdc-vhat_d0 id0-vhat_q0 iq0|
   +|vhat_d0-v_d0| i_d_bar+|vhat_q0-v_q0| i_q_bar)
     /[Cdc(Vdc0-v_bar)].
El efecto sobre almacenamiento está acotado por
  M <= ||S^-1 r||/sqrt(g),
pues X>=gI. Se usa un radio físico epsilon_actual=beta epsilon_SDP, beta=0.999.
Entonces, en V=1,
  Vdot <= -2 alpha(1-beta^2)+2M < 0.
La condición se verificó con fracciones, comparando cuadrados de cantidades
positivas. El mayor M fue {max(c['storage_trim_disturbance_bound'] for c in cert):.9g},
frente a un margen reservado alpha(1-beta^2)={cert[0]['inward_margin_reserved']:.9g}.
Así se certifica invariancia del modelo con sus parámetros numéricos reales,
no sólo de una versión a la que se le restó silenciosamente el error de trim.
No se afirma convergencia exactamente al origen exportado si su residuo no es cero.

4. RESULTADOS LOCALES Y LÍMITES DE OPTIMALIDAD

La instancia usa alpha={p['decay_rate']} s^-1, tau_m={p['measurement_tau']} s,
fmax={p['frequency_max_hz']} Hz y rmax={p['rocof_max_hz_s']} Hz/s.
El dominio sectorial usa delta_bar={p['angle_box']} rad,
|Delta i_d|,|Delta i_q|<={p['id_box']} pu y |Delta vdc|<={p['vdc_box']}.
Estas cotas incrementales NO son ratings de hardware verificados. Las restricciones
de placa y disponibilidad de energía/potencia deben agregarse antes de un diseño físico.
Las ganancias son las del candidato de red previamente encontrado, sin retuning nuevo.

  bus       Kp            Ki          epsilon_actual      brecha dual
{rows}

Los radios de tensión van de {min(c['epsilon_actual_model_verified'] for c in cert):.9g}
a {max(c['epsilon_actual_model_verified'] for c in cert):.9g} pu. Son pequeños.
Cada radio pertenece a una elipsoide conjunta de entradas; no autoriza alcanzar
simultáneamente los máximos de tensión, frecuencia de referencia y potencia DC.

Todos los contratos anteriores pasan verificación racional de factibilidad.
El solver devolvió estados optimal_inaccurate; se aceptaron únicamente tras esa
verificación independiente. Las cotas superiores duales rigurosas del SDP dejan
brechas de {100*min(c['ideal_fixed_gain_relative_gap'] for c in cert):.2f}% a
{100*max(c['ideal_fixed_gain_relative_gap'] for c in cert):.2f}%: NO hay prueba de
optimalidad estrecha para estos SDP completos. El resultado del PLL aislado no
se extrapola. Una reformulación experimental por complemento de Schur no logró
resolver sus dos casos de prueba y se conserva como diagnóstico, fuera del resultado.

5. CONTRASTES DE MODELO E IMPLEMENTACIÓN

Julia comprobó {j['cases']} puntos contra los nueve estados del GFL instalado:
  error máximo de RHS = {j['max_absolute_rhs_error']:.9g};
  error de identidad energética = {j['max_energy_identity_error']:.9g}.
Se verificaron {sum(c['samples'] for c in n['point_checks'])} puntos no lineales y su
interpolación multiafín, así como {len(n['trajectories'])} trayectorias, incluidas
condiciones iniciales cerca de la frontera de corriente. Todas pasaron.
Estos ensayos comprueban implementación; la garantía para toda señal admisible
procede de la prueba sectorial y las matrices verificadas, no de las trayectorias.

6. POR QUÉ ESTOS CERTIFICADOS LOCALES NO CIERRAN LA RED

La red sólo se usa aquí para una prueba algebraica independiente de cualquier
forma temporal de perturbación. En la rebanada de estados SG nominales y
delta_GFL=0, la respuesta de tensión a corrientes GFL es EXACTAMENTE afín:
  G(rho)=Y_Kron+sum_i(1-rho_i) S_i D_SG,i S_i^T,
  Z(rho)=-G(rho)^-1,
  Delta v_i=sum_j rho_j Z_ij R(theta_j0) Delta i_dq,j.
La fórmula conserva la admitancia de red y el efecto de cada SG retenido.
No es una simulación temporal ni una linealización respecto a esas corrientes.

Sea X_j la elipsoide local y e_delta el selector del ángulo. Al fijar delta=0,
  Xbar_j=X_j-X_j e_delta e_delta^T X_j/(e_delta^T X_j e_delta).
Sea C_j el mapa de y_j a sus incrementos de corriente. Definimos
  T_ij=rho_j Z_ij R(theta_j0) C_j,
  H_ij=T_ij Xbar_j T_ij^T,
  gamma_ij=sqrt(lambda_max(H_ij))/epsilon_i.
gamma_ij es el máximo exacto de uso de la cota de tensión receptora en esa
rebanada, variando sólo el estado del GFL j. Depende de las interacciones de red,
rho, ganancias y forma del certificado; no sólo de damping o grado del nodo.
Los rótulos de buses no intervienen en la fórmula general.

En el candidato de red actual, {w['violating_donor_receiver_pairs']} de {w['total_pairs']}
pares producen un testigo que excede una hipótesis de entrada. El mayor uso
construido fue {w['worst']['witness_input_usage']:.6f} veces la cota.

TESTIGO COLECTIVO MÁS EXIGENTE. Se contrajeron todos los estados por un factor
común y se escogieron sus direcciones mediante la función soporte
  max_(||d||=1) sum_j sqrt(d^T H_ij d).
La optimización angular sólo construye un testigo; no se afirma que haya hallado
el máximo global de esa función. Los estados concretos guardados cumplen:
  mayor V_j local = {joint['max_local_storage']:.12g} < 1;
  al aplicar cada variación por separado, el mayor uso en toda la red es
     {joint['largest_single_input_usage']:.12g} < 1;
  al aplicarlas conjuntamente, un receptor usa
     {joint['joint_input_usage']:.12g} veces su cota de entrada de tensión.
No se eligieron eventos temporales para construir este resultado.

La DAE independiente de Julia reprodujo los {jn['cases']} testigos simples y el
colectivo, con error máximo de tensión {jn['max_absolute_voltage_prediction_error']:.9g}.
La figura full_gfl_composition.png/pdf/svg presenta las regiones locales y
este testigo colectivo; figure_data.json conserva los valores y hashes de origen.
La coalición respeta los conjuntos locales de estado pero rompe una suposición
necesaria de sus certificados al interconectarlos.

INTERPRETACIÓN PRECISA. Esto refuta la suficiencia de ESA composición por producto
de elipsoides y cotas de entrada. No prueba inestabilidad, no prueba que el estado
sea alcanzable desde el equilibrio y no limita la capacidad física de sustitución.
Un certificado acoplado puede excluir esos estados o conservar las correlaciones
que el producto independiente pierde. El resultado justifica investigar ese
certificado; no constituye por sí mismo una nueva teoría de estabilidad de redes.

7. ECUACIÓN PARA LA SIGUIENTE SÍNTESIS ACOPLADA

El paso necesario es mantener la DAE y sus términos cruzados. Para una inclusión
exacta levantada del conjunto SG/GFL, escribimos
  zdot=A(q,p)z+B(q,p)v+B_w(q,p)w,
  0=C(q,p)z+Y(q,p)v+G_w(q,p)w.
Con un almacenamiento cuadrático V=z^T Pz, una condición suficiente es
  [[A^T P+PA+2 alpha P, PB, PB_w],
   [B^T P,                 0, 0],
   [B_w^T P,               0,-mu R_w]]
    +He(N^T [C,Y,G_w]) <=0,
para todos los parámetros q admisibles, donde He(M)=M+M^T.
El término agregado vale cero sobre la variedad algebraica. El almacenamiento
puede incluir bloques cruzados determinados por el grafo. P,N son variables;
rho entra en los puertos SG/GFL y Kp,Ki en los bloques de control.

Para p,alpha y la envolvente de entrada fijos, estas condiciones son afines en
P,N. Si la matriz se descompone como H0+sum_i H_i(q_i), una relajación suficiente
evita enumerar el producto exponencial de vértices: buscar Q_i con
  H_i(q_i_vertex)<=Q_i para cada vértice LOCAL,
  H0+sum_i Q_i<=0.
La estructura y sparsidad elegidas afectan el conservadurismo. Este esquema
requiere todavía levantar los SG, fijar la referencia angular, imponer contención
de dominio/condiciones iniciales y verificar las salidas. No está implementado
ni resuelto para la red completa en esta entrega.

El co-diseño exterior actualizará rho y ganancias manteniendo un certificado
verificado. En la rebanada algebraica anterior, por ejemplo,
  dZ/d rho_k=-Z S_k D_SG,k S_k^T Z,
lo que da una sensibilidad analítica de las interacciones. Ese gradiente puede
guiar una búsqueda; no convierte gamma<=1 en prueba suficiente de toda la dinámica.

8. MEDICIÓN DE FRECUENCIA DE BUS SIN ELIMINAR LA FRECUENCIA COMÚN

Los certificados actuales usan omega_PLL; el objetivo de red puede exigir la
frecuencia medida de la tensión de cada bus, que es otra salida. Una realización
exacta y finita, para una medición filtrada declarada, es la siguiente.
Escribamos la fase de bus como alpha(t)+phi_i(v), con nu=alpha_dot y una rama
angular regional consistente. Si q_abs sigue esa fase mediante un filtro tau_f,
definimos eta_i=q_abs-alpha. Entonces
  eta_i_dot=(phi_i(v)-eta_i)/tau_f-nu,
  fhat_i=(phi_i(v)-eta_i)/(2 pi tau_f),
  chi_i_dot=(fhat_i-chi_i)/tau_r,
  rhat_i=(fhat_i-chi_i)/tau_r.
Estas ecuaciones conservan la frecuencia común y permiten saltos finitos de la
fase algebraica sin diferenciar la entrada idealmente. Sus estados se añaden al
certificado; tau_f,tau_r son especificaciones, no constantes de la teoría.
phi_i requiere tensión no nula y una rama regional definida. Una ventana móvil
exacta necesita un funcional con memoria, como ya señala la teoría principal.

9. ESTADO DEL OBJETIVO COMPLETO

Hecho: inclusión no lineal exacta del GFL completo, certificados locales para los
diez dispositivos, verificación racional, auditoría de energía/modelo, prueba
constructiva de que las cotas locales independientes no cierran la red actual.
Pendiente: certificado acoplado SG/GFL/red, límites reales de hardware/energía,
medición de frecuencia/RoCoF elegida para todos los buses y co-diseño de rho,Kp,Ki
con evidencia de optimalidad adecuada. No se actualizó el póster como si esos
resultados pendientes estuvieran demostrados.

REPRODUCCIÓN PRINCIPAL
  julia --project=. experiments/nonlinear_codesign_20261001/export_gfl_sector_models.jl
  python experiments/nonlinear_codesign_20261001/full_gfl_sector_contract.py --label final
  python experiments/nonlinear_codesign_20261001/verify_full_gfl_sector_contract.py --label final
  python experiments/nonlinear_codesign_20261001/validate_full_gfl_sector_contract.py
  julia --project=. experiments/nonlinear_codesign_20261001/audit_full_gfl_mapping.jl
  julia --project=. experiments/nonlinear_codesign_20261001/export_gfl_network_closure.jl
  python experiments/nonlinear_codesign_20261001/gfl_network_closure_witness.py
  julia --project=. experiments/nonlinear_codesign_20261001/audit_gfl_network_witness.jl
  python experiments/nonlinear_codesign_20261001/render_full_gfl_composition.py
  python experiments/nonlinear_codesign_20261001/write_full_gfl_contract_report.py

Los antecedentes de disipatividad DAE/composición están en THEORY_GRAPH_ROBUST_CONTROL.txt.
Las pruebas sectoriales, de composición y los complementos de Schur son ingredientes
clásicos. La contribución aplicada deberá sustentarse en la síntesis acoplada útil,
no en reclamar novedad por escribir una desigualdad de Lyapunov o usar un SDP.
"""
(OUT/"GFL_COMPLETO_Y_COMPOSICION_ES.txt").write_text(report,encoding="utf-8")
manifest={path.name:sha(path) for path in OUT.iterdir() if path.is_file() and path.name!="manifest.json"}
(OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
print(json.dumps(dict(report=str(OUT/"GFL_COMPLETO_Y_COMPOSICION_ES.txt"),local_certificates=len(cert),
                      collective_input_usage=joint["joint_input_usage"],status="USEFUL_FULL_NETWORK_CERTIFICATE_STILL_PENDING"),indent=2))
