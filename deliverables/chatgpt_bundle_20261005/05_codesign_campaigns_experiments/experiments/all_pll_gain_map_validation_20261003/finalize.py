"""Aggregate the frozen one-step trial. No model evaluation or retuning."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import tomllib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
THEORY = ROOT / "experiments/theory_collective_damping_20261003/THEORY.tex"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

def read(name):
    return pd.read_csv(OUT / name)

def yes(series):
    return series.astype(str).str.lower().eq("true")

def verify_inputs():
    manifest = json.loads((OUT / "INPUT_MANIFEST.json").read_text())
    assert sha(OUT / "PROTOCOL.txt") == manifest["protocol_sha256"]
    result = {}
    for relative, expected in manifest["inputs"].items():
        path = ROOT / relative
        if path == THEORY and (OUT / "source_snapshot/THEORY_before_validation.tex").exists():
            path = OUT / "source_snapshot/THEORY_before_validation.tex"
        actual = sha(path)
        assert actual == expected, f"Frozen input changed: {relative}"
        result[relative] = {"sha256": actual, "verified_path": str(path)}
    lock = json.loads((OUT / "PREDICTION_LOCK.json").read_text())
    for name, expected in lock["trial_hashes"].items():
        assert sha(OUT / name) == expected, f"Prediction changed: {name}"
    dump("PRESERVATION_CHECK.json", {"inputs": result, "prediction_lock_pass": True,
         "protocol_pass": True})
    return lock

def main():
    lock = verify_inputs()
    gain = read("TABLE_01_GAIN_RECONSTRUCTION.csv")
    der = read("TABLE_02_GAIN_DERIVATIVES.csv")
    design = read("TABLE_03_GAIN_PREDICTION.csv")
    compat = read("TABLE_05_FIXED_PATTERN_COMPATIBILITY.csv")
    independent = read("TABLE_06_INDEPENDENT_MODEL.csv")
    contour = read("TABLE_07_COMPLETE_SPECTRAL_REGION.csv")
    events = read("TABLE_08_NONLINEAR_EVENTS.csv")
    poles = read("TABLE_09_INDEPENDENT_ROOTS.csv")
    baseline = pd.read_csv(ROOT / "experiments/graph_gsp_codesign_20261003/eval/baseline/events.csv")
    candidate = tomllib.loads((OUT / "designs/analytic.toml").read_text())
    expected_events = {(8, -100.), (16, 100.), (16, -100.), (29, 100.), (29, -100.)}
    assert len(events) == 5 and set(zip(events.bus, events.delta)) == expected_events
    assert set(independent.design) == {"baseline", "analytic", "fixed"}
    assert set(contour.design) == {"baseline", "analytic", "fixed"}
    assert len(gain) == 11 and len(der) == 9
    target = poles[(poles.design == "analytic") & yes(poles.target)].iloc[0]
    other = poles[(poles.design == "analytic") & ~yes(poles.target)]
    margin_pass = contour.status.eq("PASS_NUMERICAL") & contour["count"].eq(0)
    event_pass = ((events.F <= .5) & (events.R <= .5) & (events.Vmin >= .9) &
                  (events.Vmax <= 1.1) & (events.slack >= .002) & yes(events["pass"]))
    gates = {
        "reconstruction": bool(yes(gain.gain_pass).all() and (gain.gain_max_relative_error <= 1e-5).all()),
        "derivatives": bool(yes(der.derivative_pass).all()),
        "gain_bounds": bool(lock["gain_bounds_pass"]),
        "independent_equilibrium_and_jacobian": bool(yes(independent["pass"]).all()),
        "target_eigenpair": bool(target.target_pole_error <= 1e-5 and target.tracked_pattern_MAC >= .999999),
        "tracked_root_residuals": bool((poles.full_residual <= 1e-9).all()),
        "full_margin_numerical": bool(margin_pass.all()),
        "five_frozen_events": bool(event_pass.all()),
    }
    trial = independent[independent.design == "analytic"].iloc[0]
    before = independent[independent.design == "baseline"].iloc[0]
    fworst = events.loc[events.F.idxmax()]
    rworst = events.loc[events.R.idxmax()]
    sworst = events.loc[events.slack.idxmin()]
    summary = {
        "status": "TARGETED_VALIDATION_PASSED" if all(gates.values()) else "TARGETED_VALIDATION_FAILED",
        "claim_status": "NUMERICALLY_VALIDATED", "gates": gates,
        "scope": "one preregistered trial; all delays fixed and uniform; no optimization",
        "formula": "y=G_PLL(lambda,rho)q; W_i=lambda^2(1+t_fi lambda)exp(lambda tau_i)q_i/y_i; Kp_i=Im(W_i)/Im(lambda); Ki_i=Re(W_i)-Re(lambda)Kp_i",
        "rho": candidate["rho"], "Kp": candidate["Kp"], "Ki": candidate["Ki"],
        "tau_ms": 1000*candidate["tau"],
        "baseline_GFL_MW": float(before.GFL_MW), "GFL_MW": float(trial.GFL_MW),
        "retained_SG_MW": float(trial.SG_MW), "GFL_percent": 100*float(trial.GFL_MW)/lock["P_total_MW"],
        "added_GFL_MW": float(trial.GFL_MW-before.GFL_MW),
        "gain_reconstruction_max_relative": float(gain.gain_max_relative_error.max()),
        "derivative_max_relative": float(der.relative_inf_error.max()),
        "derivative_max_absolute_inf": float(der.absolute_inf_error.max()),
        "hidden_block_max_condition": float(gain.hidden_condition.max()),
        "equilibrium_max_inf": float(independent.equilibrium_inf.max()),
        "independent_jacobian_max_relative": float(independent.jacobian_relative.max()),
        "target_lambda_real": lock["lambda_real"], "target_lambda_imag": lock["lambda_imag"],
        "target_frequency_Hz": lock["lambda_imag"]/(2*np.pi),
        "target_pole_error": float(target.target_pole_error),
        "target_pattern_MAC": float(target.tracked_pattern_MAC),
        "tracked_pattern_min_MAC": float(poles.tracked_pattern_MAC.min()),
        "other_modes_max_relative_prediction_error": float(other.relative_prediction_error.max()),
        "other_modes_max_absolute_prediction_error": float(other.absolute_prediction_error.max()),
        "full_margin_per_second": .05, "interval_spectral_certificate": False,
        "max_frequency_Hz": float(fworst.F), "worst_frequency_event": [int(fworst.bus),float(fworst.delta)],
        "max_RoCoF_Hz_s": float(rworst.R), "worst_RoCoF_event": [int(rworst.bus),float(rworst.delta)],
        "min_SG_actuator_slack": float(sworst.slack), "worst_actuator_event": [int(sworst.bus),float(sworst.delta)],
        "min_voltage_pu": float(events.Vmin.min()), "max_voltage_pu": float(events.Vmax.max()),
        "nonlinear_events_passed": int(event_pass.sum()), "nonlinear_events_total": len(events),
        "nonlinear_comparator": "NOT_RUN; baseline is historical 87.5%, not same-rho fixed-gain comparator",
        "fixed_gain_comparator_margin_pass": bool(margin_pass[contour.design=="fixed"].iloc[0]),
        "fixed_pattern_conflict_min_relative": float(compat[compat.design=="replacement"].relative_residual.min()),
        "fixed_pattern_conflict_max_relative": float(compat[compat.design=="replacement"].relative_residual.max()),
        "optimality": "NOT_OPTIMIZED", "superiority": "NOT_ESTABLISHED",
        "heterogeneous_delay_validation": "NOT_RUN", "neighbor_control_validation": "NOT_RUN",
        "poster_ready_for_maximum_or_superiority_claim": False,
    }
    dump("SUMMARY.json", summary)
    dump("STATUS.json", {k: summary[k] for k in [
        "status","claim_status","gates","optimality","superiority",
        "heterogeneous_delay_validation","neighbor_control_validation",
        "poster_ready_for_maximum_or_superiority_claim"]})
    if not all(gates.values()):
        raise RuntimeError("Validation failed; positive report and manuscript generation disabled.")
    s = summary
    report = f"""VALIDACIÓN DEL MAPA ANALÍTICO DE GANANCIAS PLL — IEEE-39

Resultado: {s['status']}.
Se ejecutó UN paso preregistrado: rho_i=0.875 -> 0.876 en los diez sitios.
Tau_i={s['tau_ms']:.0f} ms, fijos e iguales. No hubo optimización ni búsqueda
del tamaño de paso. La predicción y los diseños se sellaron antes de la
validación independiente. Todos los hashes de entradas congeladas coinciden.

ECUACIÓN PUESTA A PRUEBA
y = G_PLL(lambda,rho) q
W_i = lambda^2 (1+t_fi lambda) exp(lambda tau_i) q_i/y_i
Kp_i = Im(W_i)/Im(lambda)
Ki_i = Re(W_i) - Re(lambda) Kp_i
G_PLL es el retorno de la red completa con todos los PLL abiertos; depende
de rho y de la frecuencia, pero no de las ganancias que se están calculando.
El modelo conserva los estados ocultos de red/dispositivos. Se requieren un
bloque oculto regular, un polo no real y denominadores y_i distintos de cero.
Es una identidad del modelo linealizado con retraso puro. No es una solución
analítica de toda la dinámica no lineal ni de la maximización de reemplazo.

VERIFICACIÓN
11 pares reconstruyen las ganancias conocidas: máximo error relativo
{s['gain_reconstruction_max_relative']:.6e}, límite 1e-5.
9 derivadas vectoriales (3 modos x buses 30,34,38) frente a diferencias
centradas: máximo error relativo {s['derivative_max_relative']:.6e}.
Máximo error absoluto infinito {s['derivative_max_absolute_inf']:.6e};
se conserva este dato aunque la magnitud de las derivadas sea grande.
Condición máxima del bloque oculto: {s['hidden_block_max_condition']:.6e}.
Por ello, un residuo pequeño no debe confundirse con precisión física de
trece cifras: son comprobaciones de consistencia en coma flotante.
Residuo de equilibrio máximo: {s['equilibrium_max_inf']:.6e}.
Jacobiano paramétrico frente a ForwardDiff de las ecuaciones no lineales:
error relativo máximo {s['independent_jacobian_max_relative']:.6e}.
Objetivo: lambda={s['target_lambda_real']:.10f}+j{s['target_lambda_imag']:.10f} /s
({s['target_frequency_Hz']:.6f} Hz). Es un modo PLL, no el polo dominante global.
Error de conservación del polo: {s['target_pole_error']:.6e} /s;
MAC del patrón PLL: {s['target_pattern_MAC']:.15f}.
En los otros modos, la predicción de primer orden tiene error relativo
máximo {100*s['other_modes_max_relative_prediction_error']:.3f}% y absoluto
{s['other_modes_max_absolute_prediction_error']:.6e} /s. El límite del 1%
preregistrado se refiere a DERIVADAS, no al error de un paso finito.

ESPECTRO COMPLETO EN LA REGIÓN DE SEGURIDAD
El conteo numérico por argumento/tra za y fase, con exponenciales exactas,
rotación eliminada y cota de cola, devuelve cero raíces a la derecha de
Re(s)=-0.05 /s para baseline, diseño analítico y comparador de ganancias fijas.
No se usa Padé como verdad espectral. El conteo es numérico, no un certificado
con aritmética de intervalos. No enumera el espectro infinito completo.
TABLE_07 conserva errores de cuadratura, pasos de fase y radios utilizados.

VALIDACIÓN NO LINEAL
Método de pasos, Rodas5P, tol=1e-9, paso máximo 0.01 s, horizonte 60 s.
Eventos Z-load congelados: (8,-100), (16,+100), (16,-100), (29,+100), (29,-100) MW.
Métricas de fase de tensión en 39 buses, ventana 0.5 s.
Eventos que pasan: {s['nonlinear_events_passed']}/{s['nonlinear_events_total']}.
Máximo |Delta f| = {s['max_frequency_Hz']:.9f} Hz; límite 0.5 Hz.
Máximo RoCoF = {s['max_RoCoF_Hz_s']:.9f} Hz/s; límite 0.5 Hz/s.
Tensión = [{s['min_voltage_pu']:.9f},{s['max_voltage_pu']:.9f}] pu; límites [0.9,1.1].
Holgura SG mínima = {s['min_SG_actuator_slack']:.9f}; mínimo exigido 0.002.
Cada restricción tiene su peor evento: frecuencia {s['worst_frequency_event']},
RoCoF {s['worst_RoCoF_event']}, actuador {s['worst_actuator_event']}.
No se afirma seguridad de un limitador de corriente no modelado ni robustez
para perturbaciones externas al conjunto ensayado.

DISEÑO FACTIBLE EN ESTE ENSAYO
GFL = {s['GFL_MW']:.9f} MW ({s['GFL_percent']:.6f}%).
SG retenida = {s['retained_SG_MW']:.9f} MW.
Incremento = {s['added_GFL_MW']:.9f} MW, o 0.1 punto porcentual.
Los vectores completos están en designs/analytic.toml y TABLE_03.
Este porcentaje fue prescrito, NO maximizado ni presentado como nuevo récord.

RESULTADOS QUE LIMITAN LA INTERPRETACIÓN
El comparador con ganancias fijas también pasa el margen espectral.
Este ensayo NO demuestra que retunar fuera necesario para aceptar el paso.
No se ejecutaron sus eventos no lineales: no se puede afirmar equivalencia
ni superioridad de seguridad entre ambos diseños a rho=0.876.
La gráfica de eventos compara el baseline histórico rho=0.875 con el nuevo
rho=0.876; no aísla el efecto causal del retuning.
La incompatibilidad de DOS polos Y patrones congelados produce residuos
relativos entre {s['fixed_pattern_conflict_min_relative']:.6e} y
{s['fixed_pattern_conflict_max_relative']:.6e}. Es una restricción de esa
asignación; no prueba que los PLL locales sean incapaces de estabilizar.
No se probaron retrasos heterogéneos, enlaces vecinos, compresión GSP,
convergencia de una iteración, máximo global ni una cota de optimalidad.

QUÉ CIERRA Y QUÉ SIGUE
Se cierra la trazabilidad fórmula -> ganancias -> polo/patrón -> validación
espectral y cinco eventos de un candidato. La ley utiliza el retorno
colectivo de la red, aunque su implementación siga siendo PLL local.
La ventaja práctica de diseño sigue abierta. Antes de un barrido grande,
formular un problema de comparación que permita mover los patrones modales
y use restricciones y presupuestos idénticos. Un conflicto con patrones
arbitrariamente congelados no justifica cambiar la arquitectura.
No está listo un claim de máximo reemplazo o superioridad para el póster.
"""
    (OUT/"REPORT_ES.txt").write_text(report.replace("argumento/tra za","argumento/traza"), encoding="utf-8")
    claims = f"""C1 NUMERICALLY_VALIDATED — La identidad del mapa all-PLL reconstruye 11 pares
del catálogo con error relativo máximo {s['gain_reconstruction_max_relative']:.6e}.
C2 NUMERICALLY_VALIDATED — Un paso prescrito de rho=0.875 a 0.876 conserva el
polo PLL elegido y su patrón con ganancias locales calculadas analíticamente;
error del polo {s['target_pole_error']:.6e} /s, MAC {s['target_pattern_MAC']:.15f}.
C3 NUMERICALLY_VALIDATED — Ese candidato ({s['GFL_MW']:.6f} MW de GFL) pasa el
conteo numérico de seguridad y los cinco eventos declarados, con tau=40 ms.
C4 NEGATIVE_RESULT — La prueba no discrimina una ventaja espectral necesaria:
mantener las ganancias también pasa el margen de -0.05 /s para el mismo paso.
No claim de máximo, optimalidad, superioridad no lineal, novedad bibliográfica,
efecto heterogéneo o necesidad de comunicación entre vecinos está sustentado.
"""
    (OUT/"POSTER_CLAIMS.txt").write_text(claims, encoding="utf-8")
    plots(gain, design, other, baseline, events)
    tex_section(s, events)
    tracked = subprocess.run(["git","diff","--name-only"],cwd=ROOT,capture_output=True,text=True,check=True).stdout.splitlines()
    head = subprocess.run(["git","rev-parse","HEAD"],cwd=ROOT,capture_output=True,text=True,check=True).stdout.strip()
    source_paths = [
        ROOT/"experiments/physical_collective_damping_20261003/source_snapshot/Project.toml",
        ROOT/"experiments/physical_collective_damping_20261003/source_snapshot/Manifest.toml",
    ]
    dump("FINAL_MANIFEST.json", {"git_HEAD":head,"preexisting_tracked_dirty_paths":tracked,
        "python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,
        "matplotlib":matplotlib.__version__,
        "project_files":{str(p.relative_to(ROOT)):sha(p) for p in source_paths},
        "artifacts":{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob("*"))
          if p.is_file() and p.name not in {"FINAL_MANIFEST.json","THEORY_UPDATE_PROVENANCE.json"}
          and "__pycache__" not in p.parts}})
    print(json.dumps(summary, indent=2, ensure_ascii=False))

def plots(gain, design, other, baseline, events):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,
        "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":220})
    fig, ax = plt.subplots(1,3,figsize=(15,4.5),layout="constrained")
    ax[0].semilogy(gain.root,gain.gain_max_relative_error,"o",color="#087f8c")
    ax[0].axhline(1e-5,ls="--",color="#9c3b32",label="Declared tolerance")
    ax[0].set(xlabel="Catalog mode",ylabel="Maximum relative gain error",
              title="A  Known gains reconstructed")
    ax[0].legend(fontsize=9)
    x,y=other.predicted_delta_real.to_numpy(),other.actual_delta_real.to_numpy()
    lo,hi=min(x.min(),y.min()),max(x.max(),y.max())
    pad=.06*(hi-lo)
    ax[1].plot([lo-pad,hi+pad],[lo-pad,hi+pad],"--",color="0.6")
    ax[1].scatter(x,y,color="#087f8c",s=38)
    ax[1].set(xlabel="Predicted change in Re(pole) [1/s]",
              ylabel="Independent full-model change [1/s]",
              title="B  Other modes: finite-step prediction")
    ax[2].plot(design.bus,100*(design.Kp_exact/design.Kp_before-1),"o-",label="Kp")
    ax[2].plot(design.bus,100*(design.Ki_exact/design.Ki_before-1),"s-",label="Ki")
    ax[2].axhline(0,color="0.7",lw=.7)
    ax[2].set(xlabel="Generator bus",ylabel="Gain change [%]",title="C  Analytically calculated retuning")
    ax[2].legend()
    fig.suptitle("One prescribed replacement step: 87.5% to 87.6% GFL, fixed 40 ms delay",fontsize=14)
    fig.savefig(OUT/"FIG_01_ANALYTICAL_PREDICTION.png")
    plt.close(fig)
    baseline=baseline.set_index(["bus","delta"]).loc[list(zip(events.bus,events.delta))].reset_index()
    labels=[f"{int(b)} / {int(d):+d}" for b,d in zip(events.bus,events.delta)]
    fig,ax=plt.subplots(1,3,figsize=(15,4.8),layout="constrained")
    x=np.arange(len(events))
    for a,col,title,limit in zip(ax,["F","R","slack"],
               ["Frequency excursion [Hz]","Windowed RoCoF [Hz/s]","Normalized SG actuator slack"],
               [.5,.5,.002]):
        a.bar(x-.19,baseline[col],.36,color="#8e9ba8",label="87.5%: frozen baseline")
        a.bar(x+.19,events[col],.36,color="#087f8c",label="87.6%: analytical design")
        a.axhline(limit,ls="--",color="#9c3b32",label="Declared limit")
        a.set(xticks=x,xticklabels=labels,xlabel="Event: bus / load change [MW]",title=title)
    ax[0].legend(fontsize=8,loc="lower left")
    fig.suptitle("Five frozen nonlinear events — comparison changes both replacement and gains",fontsize=14)
    fig.savefig(OUT/"FIG_02_FROZEN_EVENT_VALIDATION.png")
    plt.close(fig)

def tex_section(s, events):
    def sci(x):
        a,b=f"{x:.3e}".split("e")
        return a+r"\times10^{"+str(int(b))+"}"
    rows="\n".join(f"{int(r.bus)} & {int(r.delta):+d} & {r.F:.6f} & {r.R:.6f} & {r.slack:.6f}"+r"\\" for r in events.itertuples())
    tex = r"""
% BEGIN GENERATED ALL-PLL VALIDATION
\section{Targeted validation of the all-PLL construction}\label{sec:mapvalidation}
\status{NUMERICALLY\_VALIDATED}
This subsequent experiment tests Eq.~\eqref{eq:allgain} with one
preregistered replacement step. It is not an optimization campaign.
All ten fractions change from $\rho_i=0.875$ to $0.876$ with
$\tau_i=40$ ms fixed. The target is the leading PLL-family pole
in the frozen 2--12 Hz catalog, not the globally rightmost pole:
\[
 \lambda_\star=@LRE@+\ii\,@LIM@\ \mathrm{s}^{-1}.
\]
Its PLL-angle pattern is held fixed and all twenty gains are calculated
from the full-grid return. The trial design and modal predictions
were hashed before independent full-model validation.

Known gains were reconstructed for all eleven nonreal catalog roots,
with maximum relative error $@GRE@$. Nine vector derivative tests,
using three modes and buses 30, 34 and 38, had maximum relative error
$@DRE@$ and maximum absolute infinity error $@DAE@$.
The largest hidden-block condition number was $@COND@$;
these residuals establish numerical consistency, not physical accuracy
to all the reported digits.
Direct automatic differentiation of the nonlinear equations gave
maximum relative Jacobian discrepancy $@JAC@$ and maximum equilibrium
infinity residual $@EQ@$.

The independent delayed full-model calculation preserved the target
pole within $@PE@\ \mathrm{s}^{-1}$ with PLL-pattern MAC $@MAC@$.
For the other modes, the largest relative finite-step prediction error
was @OTHER@\%; the preregistered 1\% criterion concerns derivative
validation, not a finite-step Taylor remainder bound.
A numerical argument-principle count, cross-checked against determinant
phase and closed by a tail bound after rotational deflation, found no
roots with $\Re s>-0.05\ \mathrm{s}^{-1}$ for the baseline, retuned
design, or fixed-gain comparator. Exact exponentials were used;
floating-point contour counts are not interval certificates.

The retuned design passed all five frozen nonlinear events, simulated
with the canonical method of steps (Rodas5P, tolerance $10^{-9}$,
maximum step $0.01$ s, horizon $60$ s). The 39-bus phase metrics use
$0.5$ s windows. Frequency and RoCoF limits are $0.5$ Hz and
$0.5$ Hz/s; normalized SG actuator slack must exceed $0.002$.
\begin{center}\small
\begin{tabular}{rrrrr}\toprule
Bus & Load change [MW] & $|\Delta f|$ [Hz] & RoCoF [Hz/s] & SG slack\\
\midrule
@EVENTROWS@
\bottomrule\end{tabular}
\end{center}
All voltage samples remained within $[@VMIN@,@VMAX@]$ pu
against the declared $[0.9,1.1]$ limits.
The candidate contains @GFL@ MW of GFL (@PERCENT@\%) and
@SG@ MW of retained SG: a prescribed increase of @ADDED@ MW.
The full vectors are in
\path{experiments/all_pll_gain_map_validation_20261003/designs/analytic.toml}.

\paragraph{Interpretation and negative evidence.}
This closes the numerical chain from an analytical gain construction
to a feasible full-model candidate for the declared events.
It does not establish a maximum or demonstrate that retuning is needed
for this step: the same replacement with unchanged gains also passes
the spectral margin. Its nonlinear events were not run, so nonlinear
superiority or equivalence is undetermined.
Two frozen eigenpairs/patterns produce nonzero assignment residuals
after replacement, but this does not prove that local PLLs cannot
stabilize the network when their modal patterns are allowed to change.
No heterogeneous-delay, neighboring-feedback, GSP-compression,
convergence, or global-optimality claim follows from this trial.
Data, prediction hashes, code, trajectories and claim limits are in
\path{experiments/all_pll_gain_map_validation_20261003/}.
% END GENERATED ALL-PLL VALIDATION

"""
    replacements={
        "LRE":f"{s['target_lambda_real']:.8f}","LIM":f"{s['target_lambda_imag']:.8f}",
        "GRE":sci(s["gain_reconstruction_max_relative"]),"DRE":sci(s["derivative_max_relative"]),
        "DAE":sci(s["derivative_max_absolute_inf"]),"COND":sci(s["hidden_block_max_condition"]),
        "JAC":sci(s["independent_jacobian_max_relative"]),"EQ":sci(s["equilibrium_max_inf"]),
        "PE":sci(s["target_pole_error"]),"MAC":f"{s['target_pattern_MAC']:.12f}",
        "OTHER":f"{100*s['other_modes_max_relative_prediction_error']:.3f}",
        "EVENTROWS":rows,"VMIN":f"{s['min_voltage_pu']:.6f}","VMAX":f"{s['max_voltage_pu']:.6f}",
        "GFL":f"{s['GFL_MW']:.6f}","PERCENT":f"{s['GFL_percent']:.3f}",
        "SG":f"{s['retained_SG_MW']:.6f}","ADDED":f"{s['added_GFL_MW']:.6f}"}
    for key,value in replacements.items():
        tex=tex.replace("@"+key+"@",value)
    assert "@" not in tex
    (OUT/"GENERATED_VALIDATION_SECTION.tex").write_text(tex,encoding="utf-8")

if __name__=="__main__":
    main()
