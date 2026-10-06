"""Check the frozen AC operating point with the provisional SG/GFL split."""
import csv
import cmath
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "reports" / "experiment_E" / "tables"
def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

eq = read_csv(ROOT / "reports" / "experiment_A" / "matrices" / "bus33_baseline_equilibrium.csv")
voltage_parts = {}
for r in eq:
    if r["differential_or_algebraic"] != "algebraic" or "busbar" not in r["state_name"]:
        continue
    b = int(r["bus"])
    voltage_parts.setdefault(b, [None, None])[0 if "u_r" in r["state_name"] else 1] = float(r["equilibrium_value"])
V = {b: complex(*u) for b, u in voltage_parts.items()}
Iline = {b: 0j for b in V}
for r in read_csv(ROOT / "reports" / "experiment_D" / "inputs" / "branch.csv"):
    a, b = int(r["src_bus"]), int(r["dst_bus"])
    t = float(r["r_src"])
    z = complex(float(r["R"]), float(r["X"]))
    va, vb = t*V[a], V[b]
    i = (va-vb)/z
    Iline[a] += -t*(i+complex(float(r["G_src"]),float(r["B_src"]))*va)
    Iline[b] += i-complex(float(r["G_dst"]),float(r["B_dst"]))*vb

Yload = {int(r["bus"]): complex(float(r["initialized_admittance_real"]),
                                  float(r["initialized_admittance_imag"]))
         for r in read_csv(TABLES / "TABLE_E03_initialized_ZIP_loads.csv")}
design = {int(r["bus"]): r for r in read_csv(TABLES / "TABLE_E14_provisional_per_generator_design.csv")}
rows = []
for b in sorted(V):
    d = design.get(b)
    sg = complex(float(d["retained_SG_MW"]),float(d["retained_SG_Mvar"]))/100 if d else 0j
    gfl = complex(float(d["GFL_MW"]),float(d["GFL_Mvar"]))/100 if d else 0j
    igen = ((sg+gfl)/V[b]).conjugate()
    iload = Yload.get(b,0j)*V[b]
    residual = Iline[b]+iload+igen
    rows.append(dict(bus=b, voltage_real=V[b].real,voltage_imag=V[b].imag,
                     SG_MW=sg.real*100,GFL_MW=gfl.real*100,
                     load_MW=-(V[b]*iload.conjugate()).real*100,
                     KCL_real=residual.real,KCL_imag=residual.imag,
                     KCL_abs=abs(residual)))
with (TABLES / "TABLE_E15_provisional_electrical_balance.csv").open("w",newline="",encoding="utf-8") as f:
    w = csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
max_residual=max(r["KCL_abs"] for r in rows)
print("OPERATING_POINT_DEPENDENCE: FIXED")
print("MAX_NODAL_KCL_RESIDUAL_PU:",max_residual)
if max_residual>1e-8:
    raise SystemExit("electrical balance failed")
