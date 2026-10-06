from __future__ import annotations

import csv
import json
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "poster" / "ias2026" / "compensability_ieee39_20261002"
EVENTS = OUT / "events.csv"
LINEAR = OUT / "linear_envelope.toml"
CANDIDATE = ROOT / "reports" / "nonlinear_codesign_20261001" / "candidate_final_physical.toml"

with EVENTS.open(newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
with LINEAR.open("rb") as f:
    envelope = tomllib.load(f)
with CANDIDATE.open("rb") as f:
    candidate = tomllib.load(f)

lookup = {(r["arm"], int(r["bus"]), float(r["delta_MW"])): r for r in rows}
events = []
for bus in (8, 16, 29):
    for delta in (100.0, -100.0):
        tuned = lookup[("tuned_PLL", bus, delta)]
        nominal = lookup[("nominal_PLL_same_rho", bus, delta)]
        events.append(
            {
                "bus": bus,
                "delta_MW": int(delta),
                "F_nominal_Hz": float(nominal["Fpeak_Hz"]),
                "F_tuned_Hz": float(tuned["Fpeak_Hz"]),
                "F_gain_from_tuning_Hz": float(nominal["Fpeak_Hz"])
                - float(tuned["Fpeak_Hz"]),
                "R_nominal_Hz_s": float(nominal["Rpeak_Hz_s"]),
                "R_tuned_Hz_s": float(tuned["Rpeak_Hz_s"]),
                "R_cost_of_tuning_Hz_s": float(tuned["Rpeak_Hz_s"])
                - float(nominal["Rpeak_Hz_s"]),
                "nominal_frequency_pass": nominal["F_pass"].lower() == "true",
                "tuned_frequency_pass": tuned["F_pass"].lower() == "true",
                "nominal_voltage_pass": nominal["V_pass"].lower() == "true",
                "tuned_voltage_pass": tuned["V_pass"].lower() == "true",
            }
        )

with (OUT / "comparison.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(events[0]))
    writer.writeheader()
    writer.writerows(events)

nominal = [r for r in rows if r["arm"] == "nominal_PLL_same_rho"]
tuned = [r for r in rows if r["arm"] == "tuned_PLL"]
F0 = max(float(r["Fpeak_Hz"]) for r in nominal)
F1 = max(float(r["Fpeak_Hz"]) for r in tuned)
R0 = max(float(r["Rpeak_Hz_s"]) for r in nominal)
R1 = max(float(r["Rpeak_Hz_s"]) for r in tuned)
nominal_F_fails = sum(r["F_pass"].lower() != "true" for r in nominal)
tuned_F_fails = sum(r["F_pass"].lower() != "true" for r in tuned)
nominal_R_fails = sum(r["R_pass"].lower() != "true" for r in nominal)
tuned_R_fails = sum(r["R_pass"].lower() != "true" for r in tuned)
nominal_V_fails = sum(r["V_pass"].lower() != "true" for r in nominal)
tuned_V_fails = sum(r["V_pass"].lower() != "true" for r in tuned)

# The two panels show the metric tradeoff at fixed rho and common event set.
labels = [f"{int(e['delta_MW']):+d} MW @ {e['bus']}" for e in events]
y = np.arange(len(events))
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.7), sharey=True)
colors = {"nominal": "#D19A2A", "tuned": "#007C70"}
for ax, nominal_key, tuned_key, limit, title, xlabel in [
    (
        axes[0],
        "F_nominal_Hz",
        "F_tuned_Hz",
        0.5,
        "Frequency peak",
        "Maximum |Δf| (Hz), all 39 buses",
    ),
    (
        axes[1],
        "R_nominal_Hz_s",
        "R_tuned_Hz_s",
        0.5,
        "RoCoF peak",
        "Maximum |RoCoF| (Hz/s), all 39 buses",
    ),
]:
    ax.axvline(limit, color="#B33C32", linestyle="--", linewidth=1.5, label="Limit 0.5")
    for i, event in enumerate(events):
        a = event[nominal_key]
        b = event[tuned_key]
        ax.plot([a, b], [i, i], color="#9BA5A3", linewidth=1.6, zorder=1)
        ax.scatter(a, i, color=colors["nominal"], marker="x", s=62, linewidth=2.0, zorder=3)
        ax.scatter(b, i, color=colors["tuned"], marker="o", s=42, zorder=4)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#DCE3E1", linewidth=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
axes[0].set_ylabel("Load-setpoint event")
axes[0].set_xlim(0.40, 0.535)
axes[1].set_xlim(0.08, 0.38)
fig.suptitle(
    "Same IEEE-39 replacement mix (90.785% GFL): PLL tuning lowers frequency peaks\n"
    "and raises RoCoF in every one of the six tested events",
    x=0.06,
    ha="left",
    fontsize=13,
    fontweight="bold",
)
handles = [
    plt.Line2D([0], [0], color=colors["nominal"], marker="x", linestyle="None",
               markersize=8, markeredgewidth=2, label="Nominal PLL"),
    plt.Line2D([0], [0], color=colors["tuned"], marker="o", linestyle="None",
               markersize=6, label="Co-designed PLL"),
    plt.Line2D([0], [0], color="#B33C32", linestyle="--", label="Limit"),
]
fig.legend(handles=handles, ncol=3, frameon=False, loc="lower center",
           bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.06, 1, 0.88))
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT / f"retuning_tradeoff.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)

report = f"""COMPENSABILIDAD PLL A MEZCLA SG/GFL FIJA — IEEE-39
Prueba de concepto, 2026-10-02

PREGUNTA Y MODELO

Se mantuvo fijo el vector rho de un candidato de 90.78544633% GFL
(497.84032003 MW SG retenidos). Se compararon exactamente el mismo punto
operativo, soporte de dispositivos, red, cargas y eventos con:
  (1) las ganancias Kp, Ki del candidato co-diseñado;
  (2) Kp=2*pi*5, Ki=(2*pi*5)^2/4 en todos los GFL.

El DAE de PowerDynamics se inicializó de nuevo para cada brazo. Se integraron
61 s incluyendo 1 s de prefalla, paso de salida 0.005 s, tolerancias 1e-9,
y se midieron frecuencia, RoCoF de ventana 0.5 s y tensión en los 39 buses.

RESULTADO DEL CONTRAFACTUAL NO LINEAL

Ambos brazos completaron los seis eventos. El tuning redujo el pico de
frecuencia en LOS SEIS:
"""
for e in events:
    report += (
        f"  bus {e['bus']}, carga {e['delta_MW']:+d} MW: "
        f"{e['F_nominal_Hz']:.6f} -> {e['F_tuned_Hz']:.6f} Hz "
        f"(mejora {e['F_gain_from_tuning_Hz']:.6f} Hz); "
                f"RoCoF {e['R_nominal_Hz_s']:.6f} -> {e['R_tuned_Hz_s']:.6f} Hz/s.\n"
    )
report += f"""
El peor pico de frecuencia pasó de {F0:.6f} a {F1:.6f} Hz. Con el mismo
límite exploratorio de 0.5 Hz, los PLL nominales fallaron {nominal_F_fails}/6
eventos y el candidato ajustado  {tuned_F_fails}/6.

El tuning subió el RoCoF máximo de {R0:.6f} a {R1:.6f} Hz/s. Ningún brazo
excedió el límite de 0.5 Hz/s ({nominal_R_fails}/6 frente a {tuned_R_fails}/6).
Todas las tensiones quedaron entre 0.9 y 1.1 pu ({nominal_V_fails}/6 y
{tuned_V_fails}/6 fallos, respectivamente).

Esto demuestra una COMPENSACIÓN condicionada en esta instancia: a igual
porcentaje de reemplazo, los PLL ajustados recuperan margen de frecuencia
frente a dos perturbaciones de carga que hacen fallar la referencia nominal.
El efecto no es gratuito: en este conjunto el pico de RoCoF aumentó en los
seis casos. La figura retuning_tradeoff muestra la relación de un vistazo.

ENVOLVENTE MATEMÁTICA LOCAL

Para la DAE index-1 xdot=f(x,z,p,w), 0=g(x,z,p,w), en una rama algebraica
regular, la eliminación local de z da:
  A = f_x - f_z g_z^(-1) g_x,
  B_eta = f_eta - f_z g_z^(-1) g_eta,
  Sdot_eta = A S_eta + B_eta,
para eta=rho o kappa=(log Kp,log Ki). Apilando restricciones normalizadas
en g(p)<=0, su mapa de compensación de primer orden es:
  maximize tau
  sujeto a g + J_rho * 1 * tau + J_kappa * delta_kappa <= 0,
           |delta_kappa_l| <= log(1.1), 0<=tau<=0.01.

Con 36 filas (6 modos y 5 métricas por cada uno de 6 eventos), el QP
regularizado sobre la tangente predijo tau={envelope['max_uniform_delta_rho']:.8g},
equivalente a {envelope['capacity_gain_percentage_points']:.6f} puntos
porcentuales adicionales de reemplazo bajo esa dirección uniforme. La
restricción activa es el polo de decaimiento. La sensibilidad de ese polo
analítica y la diferencia central fresca difieren {envelope['active_modal_gradient_finite_difference_relative_error']:.3e}
relativamente; la fila modal guardada coincide con el gradiente analítico
actualizado (error {envelope['stored_modal_row_vs_fresh_analytic_relative_error']:.3e}).

Ese número es una predicción local de la tangente; NO es una cota física.
Las filas de trayectoria proceden de una tangente SDIRK de paso 0.025 s y
su resto de segundo orden no está acotado. Por eso aún no certificamos que
ninguna ganancia del recuadro pueda superar la frontera no lineal.

FORMA CONDICIONAL DEL TEOREMA QUE FALTA VALIDAR

Sea h(p)<=0 una familia de restricciones de trayectoria diferenciables en
un recuadro B, y lambda>=0. Si se prueba en todo B la cota
  ||nabla^2(lambda' h)(p)||_2 <= L_lambda,
entonces, para cada delta en B,
  lambda'h(p+delta) >= lambda'h(p)
      + lambda'J delta - (L_lambda/2)||delta||_2^2.
Si el mínimo global de la parte derecha sobre B es estrictamente positivo,
no existe diseño factible en B: cualquier punto factible tendría
lambda'h<=0, contradicción. Para una cota rectangular en las ganancias,
el mínimo lineal es una función soporte que puede calcularse por coordenada.

La prueba es una aplicación de Taylor con resto y separación dual; por sí
sola no sería una novedad matemática. El hueco investigable es obtener una
cota de curvatura estrecha y reproducible a través del DAE completo y de sus
interacciones de red, y demostrar que esa exclusión finita anticipa la
frontera observada por PowerDynamics. Si se activan limitadores, hay que
certificar cada modo híbrido o tratar los eventos; la suavidad de un Taylor
único ya no se puede asumir.

ALCANCE Y SIGUIENTE PRUEBA

La comparación mantiene fija la rho del candidato, pero su contraparte de
ganancias nominales no es una búsqueda óptima sobre todo el recuadro. Los
seis escalones sólo falsifican/ilustran casos; no representan toda la clase
general de perturbaciones. Los límites de corriente del convertidor y de
energía DC no están certificados. 90.785% no es máximo global.

El siguiente paso para que esto sea un resultado de investigación es:
1. verificar cotas Hessianas/intervalares para las restricciones activas y
   mantener separadas las regiones de limitadores;
2. recorrer radios de recuadro y de sustitución, cotejando la predicción
   del no-go finito con optimización y trayectoria PowerDynamics;
3. buscar un certificado no local o un gap interior/exterior en una instancia
   reducida; después repetir la ablación de grafo y DAE completa en IEEE-39.

El experimento prueba que el tuning importa a porcentaje SG/GFL fijo. Aún
no prueba una teoría nueva, no demuestra imposibilidad global y no valida
robustez frente a cualquier perturbación.

ARCHIVOS Y REPRODUCCIÓN

Experimento PowerDynamics:
  julia --startup-file=no --project=. experiments/compensability_ieee39_20261002/pd_compare.jl
Sobre local y auditoría modal:
  julia --startup-file=no --project=. experiments/compensability_ieee39_20261002/linear_envelope.jl
Resumen y figura:
  C:/Users/walla/anaconda3/python.exe experiments/compensability_ieee39_20261002/summarize_results.py

Resultados: events.csv, comparison.csv, linear_envelope.toml y
retuning_tradeoff.png/.svg/.pdf en reports/poster/ias2026/compensability_ieee39_20261002/.
"""
(OUT / "RESULTADOS_ES.txt").write_text(report, encoding="utf-8")

print(json.dumps({
    "events": events,
    "frequency_peak_nominal": F0,
    "frequency_peak_tuned": F1,
    "rocof_peak_nominal": R0,
    "rocof_peak_tuned": R1,
    "frequency_failures_nominal": nominal_F_fails,
    "frequency_failures_tuned": tuned_F_fails,
    "local_capacity_gain_percentage_points": envelope["capacity_gain_percentage_points"],
}, indent=2, ensure_ascii=False))
