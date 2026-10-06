"""Generate the report and merge the theory extension into the open source."""
from pathlib import Path
import json
import re
from checks import OUT,THEORY,ROOT,sha,save,freeze

freeze()
s=json.loads((OUT/"SUMMARY.json").read_text())
d=json.loads((OUT/"DESCRIPTOR_SUMMARY.json").read_text())
snapshot=OUT/"source_snapshot/THEORY_before_regional.tex"
provenance=OUT/"DOCUMENT_PROVENANCE.json"
if provenance.exists():
    assert sha(THEORY)==json.loads(provenance.read_text())["theory_after_sha256"],"Manual edit since integration; inspect before merging."
else:
    assert sha(THEORY)==sha(snapshot),"Open theory changed; merge manually."
section=(OUT/"SECTION.tex").read_text(encoding="utf-8").replace(r"\mathscr",r"\mathcal")
(OUT/"SECTION.tex").write_text(section,encoding="utf-8")
def sci(x):
    a,b=f"{x:.3e}".split("e")
    return a+r"\times10^{"+str(int(b))+"}"
checks=r"""
\paragraph{Executed small checks and a negative bound result.}
Exact symbolic checks passed @NS@ identities; the two-state pure-DDE
fixture passed @NT@ contour/root checks. Its uniform remainder bound
is proved analytically over $p\in[-1/100,1/100]$ and the entire circle,
not inferred from the sampled roots. The exact bound is $36/28525$.
The IEEE--39 parameter-factorization checks passed at the three
preregistered frequencies. The largest lifted determinant/equivalence
discrepancy was $@ERR@$.

The first scalar product-of-norms bound was inconclusive at every
point: its values ranged from @LOOSEMIN@ to @LOOSEMAX@.
Keeping those outputs, the affine descriptor and positive-vector
envelope were evaluated on exactly the same parameter box, gains,
delay vector and frequency points. Their pointwise upper expressions
ranged from @TIGHTMIN@ to @TIGHTMAX@, all below one.
This is a numerical comparison of two analytical envelopes; it is
not a controller-performance or certified stability comparison.
The heterogeneous delay vector used for these matrix checks was
fixed beforehand. No new equilibrium, root-count or nonlinear
event validation was run for that vector.

\status{BLOCKED}
The continuous contour enclosures and baseline cluster integrals
required for an IEEE--39 regional MW bound remain unexecuted.
There is no new maximum replacement or new feasible design.
The initial loose envelope is retained as a negative result.

"""
for key,value in {
    "NS":str(s["symbolic_passed"]),"NT":str(s["toy_DDE_passed"]),
    "ERR":sci(d["maximum_algebra_error"]),
    "LOOSEMIN":f"{s['IEEE39_scalar_envelope_min']:.3f}",
    "LOOSEMAX":f"{s['IEEE39_scalar_envelope_max']:.3f}",
    "TIGHTMIN":f"{d['weighted_box_bound_min']:.6f}",
    "TIGHTMAX":f"{d['weighted_box_bound_max']:.6f}"}.items():
    checks=checks.replace("@"+key+"@",value)
assert "@" not in checks
section += checks
source=snapshot.read_text(encoding="utf-8")
marker=r"\section{Claim ledger, novelty and missing gates}"
assert source.count(marker)==1
source=source.replace(marker,"% BEGIN REGIONAL CERTIFICATE EXTENSION\n"+section+
    "% END REGIONAL CERTIFICATE EXTENSION\n\n"+marker)
abstract_anchor="demonstrate a replacement maximum or superiority over fixed gains.\n"
assert source.count(abstract_anchor)==1
source=source.replace(abstract_anchor,abstract_anchor+
    "A further regional construction retains algebraic variables to obtain an\n"
    "affine descriptor perturbation, bounds finite root-cluster changes by a\n"
    "matrix-logarithm remainder, and yields a necessary replacement LP with\n"
    "free modal patterns. Its IEEE--39 continuous-domain certificate remains\n"
    "unevaluated.\n")
ledger_anchor="First-ever result or superiority over literature & BLOCKED: not established\\\\"
assert source.count(ledger_anchor)==1
ledger_end=source.index(r"\end{center}",source.index(marker))+len(r"\end{center}")
source=source[:ledger_end]+r"""

The regional extension has separate statuses:
affine descriptor and finite cluster identity are EXACT\_IDENTITY;
the necessary LP and remainder bound are PROVED\_REDUCED\_MODEL
under the stated uniform hypotheses; pointwise matrix checks are
NUMERICALLY\_VALIDATED. A numerical IEEE--39 regional MW certificate
is BLOCKED until the continuous enclosures are evaluated.
"""+source[ledger_end:]
bib=r"""
\bibitem{beyn} W.-J. Beyn,
\emph{An integral method for solving nonlinear eigenvalue problems},
arXiv:1003.1580, 2010.
\url{https://arxiv.org/abs/1003.1580}.
\bibitem{lessard} J.-P. Lessard and J. D. Mireles James,
\emph{A functional analytic approach to validated numerics for
eigenvalues of delay equations}, Journal of Computational Dynamics,
7(1):123--158, 2020.
\url{https://doi.org/10.3934/jcd.2020005}.
\bibitem{michielsdesign} W. Michiels and S. Gumussoy,
\emph{Eigenvalue based algorithms and software for the design of
fixed-order stabilizing controllers for interconnected systems with
time-delays}, in Delay Systems: Advances in Delays and Dynamics,
Springer, 2014; author version arXiv:2003.05496.
\url{https://doi.org/10.1007/978-3-319-01695-5_18}.
"""
assert source.count(r"\end{thebibliography}")==1
source=source.replace(r"\end{thebibliography}",bib+r"\end{thebibliography}")
# Structural checks supplement, and do not replace, native compilation.
labels=re.findall(r"\\label\{([^}]+)\}",source)
assert len(labels)==len(set(labels)),"Duplicate LaTeX label"
refs=re.findall(r"\\(?:eqref|ref)\{([^}]+)\}",source)
assert set(refs)<=set(labels),"Unresolved cross-reference"
keys=set(re.findall(r"\\bibitem\{([^}]+)\}",source))
cites={k for group in re.findall(r"\\cite\{([^}]+)\}",source) for k in group.split(",")}
assert cites<=keys,"Unresolved bibliography"
assert r"\mathscr" not in source
for env in ["equation","align","align*","theorem","proof","bmatrix","enumerate"]:
    assert source.count("\\begin{"+env+"}")==source.count("\\end{"+env+"}"),env
THEORY.write_text(source,encoding="utf-8")
save("DOCUMENT_PROVENANCE.json",{
    "theory_before_sha256":sha(snapshot),"theory_after_sha256":sha(THEORY),
    "snapshot":str(snapshot),"theory":str(THEORY),
    "summary_sha256":sha(OUT/"SUMMARY.json"),
    "descriptor_summary_sha256":sha(OUT/"DESCRIPTOR_SUMMARY.json"),
    "section_sha256":sha(OUT/"SECTION.tex"),
    "scope":"append regional finite-remainder theory and generated check results; preserve previous experiments",
    "latex_structure_checks_pass":True})
report=f"""CIERRE TEÓRICO REGIONAL — AVANCE Y LÍMITES

Se añadió la sección 'A finite regional replacement bound with moving modal
patterns' al THEORY.tex abierto. No se ejecutaron simulaciones de red ni
optimización. Los resultados congelados y su procedencia se preservaron.

PROBLEMA QUE SE RESUELVE MATEMÁTICAMENTE
La aproximación anterior g=g0+J*d+r dejaba sin construir la cota de r.
Ahora se deriva una identidad de la media de raíces dentro de un contorno:

  z(p)-z(p0) = -1/(2*pi*i*m) integral_Gamma tr(log(I+H(s,d))) ds.

Se permiten movimientos y colisiones internas de los modos; no se fijan
autovectores. Se requiere conteo constante y regularidad en todo el contorno.
Una media estable es NECESARIA, no suficiente, para estabilidad de cada raíz.

MEJORA ESPECÍFICA DEL CONTRATO IEEE-39 INSPECCIONADO
Retener voltajes y detectores algebraicos produce una matriz descriptor F
de 234 coordenadas equivalente al modelo de 204 estados:
det(F)=det(G)*det(Delta), con G regular.
Su cambio es exactamente U*Theta(d)*V(s), afín en rho,Kp,Ki ordinarios.
La matriz T=V*F0^(-1)*U tiene 40 canales y H=Theta*T.
No hay resto de parametrización matricial: H_rem=0. El resto finito
proviene de la serie del logaritmo, que se acota explícitamente:

  beta = longitud(Gamma)*40/(2*pi*m)*[-log(1-kappa)-kappa].

La fórmula conserva retrasos puros fijos, incluso heterogéneos. No convierte
el modelo físico no lineal en lineal ni demuestra estabilidad de gran señal.

QUÉ PUEDE CERTIFICAR
Con enclosures rigurosos de kappa<1 y de los coeficientes, el LP necesario
maximiza P0'*delta_rho sobre TODAS las direcciones de la caja regional y
las ganancias admisibles. Su dual elimina las ganancias mediante soportes
de sus intervalos. Un testigo positivo puede excluir un reemplazo concreto
en la región, incluso permitiendo que cambien los patrones modales.
Esto daría P_F <= P*_D <= P_U,D. P_F debe pertenecer al mismo dominio,
tener el mismo vector de retrasos y pasar espectro y eventos completos.
La factibilidad numérica no se convierte en un certificado por combinarla
con una cota superior. Un límite regional no es un máximo global.

COMPROBACIONES EJECUTADAS
Identidades simbólicas: {s['symbolic_passed']}/{s['symbolic_checks']}.
Fixture DDE de dos estados: {s['toy_DDE_passed']}/{s['toy_DDE_checks']}.
Discrepancia máxima raíz-integral: {s['toy_max_moment_error']:.6e}.
El ejemplo pequeño tiene cota uniforme demostrada 36/28525,
sin estimarla a partir de muestras. No es una red de potencia.
Factorizaciones IEEE-39: {s['IEEE39_algebra_passed']}/{s['IEEE39_algebra_points']}.
Descriptor afín y equivalencia de determinantes: {d['points_passed']}/3.
Discrepancia máxima del descriptor: {d['maximum_algebra_error']:.6e}.

RESULTADO NEGATIVO CONSERVADO
La primera cota escalar por productos de normas estuvo entre
{s['IEEE39_scalar_envelope_min']:.6f} y {s['IEEE39_scalar_envelope_max']:.6f},
por encima del requisito kappa<1 en los tres puntos.
La reformulación afín y la norma ponderada dieron expresiones superiores
puntuales entre {d['weighted_box_bound_min']:.6f} y {d['weighted_box_bound_max']:.6f},
en la misma caja y los mismos puntos. El refinamiento cambió el análisis,
no el experimento ni el controlador. Esto es NUMERICALLY_VALIDATED.
La primera ejecución del script del descriptor pasó el álgebra y falló al
serializar un entero NumPy; se corrigió la conversión a int y se repitieron
los mismos puntos. No se cambió ninguna tolerancia o definición científica.

QUÉ FALTA, EXACTAMENTE
1. Elegir contornos de grupos físicos relevantes bajo un vector de tau fijo.
2. Verificar inversas, exponenciales, conteos e integrales sobre todos sus
   paneles continuos con redondeo exterior. Hay una fórmula explícita para
   acotar entre puntos: no basta incrementar la densidad de una malla.
3. Construir y resolver el LP superior con los errores certificados incluidos.
4. Comprobar que limita la región más que sus propios bordes; de lo contrario
   es una cota válida pero inútil para explicar soporte síncrono.
5. Identificar si el límite es colectivo/modal, una ganancia saturada o un
   evento no lineal. No presuponer que los PLL son el cuello de botella.
6. Contrastar con métodos existentes bajo el mismo contrato para una novedad
   de aplicación demostrable. La dualidad, los contornos, las expansiones
   log-det y las cotas de matrices positivas son herramientas conocidas.

CONEXIÓN CON BEYOND NODAL DAMPING
tr(H^2) contiene productos H_ij*H_ji; órdenes superiores contienen recorridos
cerrados entre puertos dinámicos. Su signo es complejo/frecuencial y depende
de la red completa. No se denominan ciclos físicos de líneas ni Laplaciano.
La expansión permite acotar cuánto se pierde al omitir acoplamientos;
todavía no demuestra que esos términos limiten reemplazo en IEEE-39.

ESTADO FINAL
Teoremas/identidades condicionales: derivados y documentados.
Cota uniforme analítica del ejemplo pequeño: demostrada.
Identidades del modelo IEEE-39: verificadas numéricamente.
Cota superior IEEE-39 en MW: NO OBTENIDA.
Óptimo global/nuevo diseño/ventaja de reemplazo: NO OBTENIDOS.
Novedad bibliográfica: NO ESTABLECIDA.
No se necesita otra simulación pesada para el siguiente paso: falta cerrar
la verificación de contorno y evaluar si el LP es informativo.
"""
(OUT/"REPORT_ES.txt").write_text(report,encoding="utf-8")
save("STATUS.json",{
    "theory":"CONDITIONAL_REGIONAL_BOUND_DERIVED",
    "parameter_lift":"EXACT_IDENTITY",
    "toy_bound":"PROVED_REDUCED_MODEL",
    "IEEE39_algebra":"NUMERICALLY_VALIDATED",
    "naive_scalar_envelope":"NEGATIVE_RESULT",
    "IEEE39_regional_MW_bound":"BLOCKED_CONTINUOUS_CONTOUR_ENCLOSURES",
    "poster_claim_maximum_replacement":False,"novelty_established":False})
save("FINAL_MANIFEST.json",{"source":str(THEORY),"source_sha256":sha(THEORY),
    "artifacts":{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob("*"))
       if p.is_file() and p.name!="FINAL_MANIFEST.json" and "__pycache__" not in p.parts}})
print("Merged the regional proof and generated evidence into the existing THEORY.tex.")
