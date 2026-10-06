"""Record the scope and model compatibility of the read-only Q2B transfer.

This reads published artifacts and primary source, without rerunning or claiming
independent verification of Q2B's optimization/KKT/Riccati calculations.
"""
from pathlib import Path
import csv
import hashlib
import json
import tomllib

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "reports/experiment_Q2B/CERTIFIED_SEARCH"
OUT = ROOT / "reports/nonlinear_codesign_20261001"
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
read = lambda p: p.read_text(encoding="utf-8")
summary = json.loads(read(SRC / "TERMINAL_SUMMARY.json"))
candidate = tomllib.loads(read(SRC / "Z_LOCAL_SECURE_FINAL.toml"))
with (SRC / "PD/EVENTS.csv").open(encoding="utf-8", newline="") as f:
    events = {r["event"]: r for r in csv.DictReader(f)}

candidate_sha = sha(SRC / "Z_LOCAL_SECURE_FINAL.toml")
assert candidate_sha == "66513a8d3a1b4cc37c9d6f2b1a5824371fed2d0c5621d2ccdb7540490fe59124"
assert summary["status"] == "PASS_LOCAL_SECURE"
assert abs(summary["GFL_fraction"] - candidate["GFL_fraction"]) < 1e-14
assert summary["global_certified"] is False
assert float(events["bus8_100MW"]["Rpeak"]) > candidate["RoCoF_limit_Hz_s"]
assert float(events["bus29_100MW"]["Fpeak"]) > candidate["frequency_limit_Hz"]
assert events["bus16_100MW"]["frequency_pass"] == "true"
assert events["bus16_100MW"]["rocof_pass"] == "true"

model = ROOT / "src/pd39/model.jl"
reference = ROOT / "src/bnd_model_expN/PDReferenceN.jl"
validator = ROOT / "experiments/bnd_expQ2B/certified_search/validate_independent_pd.jl"
physical = ROOT / "experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl"
assert '"PDReferenceN.jl"' in read(validator)
assert "PD39.PD39Model.WeightedSimpleGFLDC(" in read(reference)
assert "(p_ac-P_dc)/v_dc_state" in read(model)
assert "iset_d_dc=(V_dc-v_dc_state)*kp_v_dc+v_dc_i" in read(model)
assert "(P_dc-p_ac)/v_dc_state" in read(physical)
assert "iset_d_dc=(v_dc_state-V_dc)*kp_v_dc+v_dc_i" in read(physical)

paths = [SRC / f for f in ("TERMINAL_SUMMARY.json", "Z_LOCAL_SECURE_FINAL.toml",
         "PD/EVENTS.csv", "RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md")]
paths += [model, reference, validator, physical,
          ROOT / "experiments/bnd_expQ2B/certified_search/LocalOracle.jl"]
data = dict(
    status="TRANSFER_SCOPE_AND_MODEL_DIFFERENCE_RECORDED",
    q2b_reported_status=summary["status"],
    q2b_reported_GFL_fraction=summary["GFL_fraction"],
    independently_reran_q2b_optimization=False,
    nonlinear_DC_model_identical=False,
    robust_all_input_capacity_bound=False,
    reported_primary_event_F_Hz=float(events["bus16_100MW"]["Fpeak"]),
    reported_primary_event_R_Hz_s=float(events["bus16_100MW"]["Rpeak"]),
    reported_external_bus8_R_Hz_s=float(events["bus8_100MW"]["Rpeak"]),
    reported_external_bus29_F_Hz=float(events["bus29_100MW"]["Fpeak"]),
    source_hashes={str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in paths},
    script_sha256=sha(Path(__file__)),
)
(OUT / "q2b_transfer_audit.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
report = f"""TRANSFERENCIA Q2B: ALCANCE Y COMPATIBILIDAD DE MODELO
2026-10-01. Lectura de artefactos; no se repitió el optimizador del otro chat.

1. RESULTADO NUEVO QUE DEBE RECONOCERSE
CERTIFIED_SEARCH declara PASS_LOCAL_SECURE, certificación NUMÉRICA local del
problema analítico congelado de un evento. Supera la tentativa fallida histórica
en LOCAL_CERTIFICATE. El nuevo candidato no debe etiquetarse con ese fallo previo.
Sustitución reportada: {100*summary['GFL_fraction']:.9f}%.
SG retenido: {summary['retained_MW']:.9f} MW.
Estacionariedad reportada: {summary['stationarity']:.9g}; primal cero.
El informe incluye LICQ y SOSC numéricos, no sólo un residuo KKT.
SHA candidato: {candidate_sha}.

Su optimalidad es local para restricciones analíticas del modelo nominal y el
evento especificados. La trayectoria no lineal se comprobó aparte; no demuestra
optimalidad del problema de trayectorias no lineales ni de nuestro problema robusto.
La auditoría presente confirma consistencia de los archivos leídos, no reproduce
independientemente sus certificados KKT/SOSC/Riccati ni sus simulaciones PD.

2. ALCANCE DE ENTRADA Y MEDICIÓN
Evento primario: cambio sostenido de Pset del ZIP de 100 MW en bus 16.
Medida: diferencias causales de fase, ventana 0.5 s, diez buses generadores.
PD reporta F={data['reported_primary_event_F_Hz']:.9f} Hz,
R={data['reported_primary_event_R_Hz_s']:.9f} Hz/s: pasa sus dos límites de 0.5.
Fuera del caso de diseño, la misma tabla conserva dos fallos:
  bus8/100 MW: R={data['reported_external_bus8_R_Hz_s']:.9f} Hz/s;
  bus29/100 MW: F={data['reported_external_bus29_F_Hz']:.9f} Hz.
No se usa el porcentaje como capacidad para cualquier señal o soporte de entrada.
El campo PowerDynamics_validation_at_freeze=NOT_RUN en el TOML describe el momento
de congelación; la tabla PD posterior y el resumen sí reportan la validación.

3. DIFERENCIA DE MODELO COMPROBADA EN EL CÓDIGO
validate_independent_pd.jl -> PDReferenceN.jl -> WeightedSimpleGFLDC en
src/pd39/model.jl. La realización transferida usa
  Cdc*vdc*vdcdot=Pac-Pdc,
  iref_d=kpdc*(Vref-vdc)+vdi,
  vdi_dot=kidc*(Vref-vdc).
La realización física aislada de esta investigación, PDPhysicalReference.jl, usa
  Cdc*vdc*vdcdot=Pdc-Pac,
  iref_d=kpdc*(vdc-Vref)+vdi,
  vdi_dot=kidc*(vdc-Vref).
Con Pac definido como potencia del inversor hacia su filtro, Pdc suministro positivo
y E=Cdc*vdc^2/2, la segunda ecuación DC conserva Edot=Pdc-Pac.
La auditoría energética previa está en dc_physics_audit.toml. Los Jacobianos
nominales de estas dos convenciones son similares mediante un cambio de signo
en la perturbación DC; sus polos iguales no validan identidad no lineal.
No se modificaron los archivos transferidos ni el modelo compartido.

4. QUÉ SE PUEDE REUTILIZAR
  - SQP con múltiples picos, corrección de segundo orden y región de confianza.
  - Disciplina de KKT, LICQ, SOSC y clasificación numérico/formal/global.
  - Tratamiento explícito de soporte y límites unilaterales al insertar SG.
  - Reconstrucción algebraica del límite derecho de un evento al exportar sensores.
La maquinaria temporal puede servir de contraste/falsificación del control robusto;
sus eventos no pasan a definir la teoría. Las ganancias y rho transferidas son
semillas candidatas que necesitan reequilibrio, reevaluación física y nueva
verificación de restricciones antes de entrar en una comparación.

5. QUÉ SIGUE SIN TRANSFERIRSE
No se transfiere una garantía para entradas generales, límite de corriente,
reserva energética, estabilidad no lineal regional ni máximo global.
Las cotas L=0, U=312.523288435 MW son del problema de minimización SG de Q2B;
no son cotas para nuestro problema de seguridad no lineal para toda señal.
Los diez contratos GFL aquí verificados, y su fallo de cierre independiente,
pertenecen al candidato físico anterior y no han sido recalculados para Q2B.

Trazabilidad: q2b_transfer_audit.json. Referencias primarias leídas:
  reports/experiment_Q2B/CERTIFIED_SEARCH/RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md
  reports/experiment_Q2B/CERTIFIED_SEARCH/TERMINAL_SUMMARY.json
  reports/experiment_Q2B/CERTIFIED_SEARCH/Z_LOCAL_SECURE_FINAL.toml
  reports/experiment_Q2B/CERTIFIED_SEARCH/PD/EVENTS.csv
"""
(OUT / "Q2B_CERTIFIED_SEARCH_TRANSFER_AUDIT.txt").write_text(report, encoding="utf-8")
print(json.dumps({k: data[k] for k in ("status", "nonlinear_DC_model_identical",
                                      "robust_all_input_capacity_bound")}, indent=2))
