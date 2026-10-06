import csv
import cmath
import math

ROOT = "reports/experiment_D/inputs/"
with open("reports/experiment_A/matrices/bus33_baseline_equilibrium.csv", encoding="utf-8") as f:
    eq = list(csv.DictReader(f))
V = {}
for r in eq:
    if r["differential_or_algebraic"] != "algebraic":
        continue
    b = int(r["bus"])
    if "busbar" not in r["state_name"]:
        continue
    V.setdefault(b, [None, None])
    k = 0 if "u_r" in r["state_name"] else 1
    V[b][k] = float(r["equilibrium_value"])
V = {b: complex(*v) for b, v in V.items()}
Iline = {b: 0j for b in V}
with open(ROOT + "branch.csv", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        a, b = int(r["src_bus"]), int(r["dst_bus"])
        t = float(r["r_src"])
        z = complex(float(r["R"]), float(r["X"]))
        va, vb = t * V[a], V[b]
        i = (va - vb) / z
        Iline[a] += -t * (i + complex(float(r["G_src"]), float(r["B_src"])) * va)
        Iline[b] += i - complex(float(r["G_dst"]), float(r["B_dst"])) * vb
with open(ROOT + "load.csv", encoding="utf-8") as f:
    L = {int(r["bus"]): complex(float(r["Pset"]), float(r["Qset"])) for r in csv.DictReader(f)}
with open(ROOT + "bus.csv", encoding="utf-8") as f:
    B = {int(r["bus"]): r for r in csv.DictReader(f)}
with open(ROOT + "machine.csv", encoding="utf-8") as f:
    M = {int(r["bus"]): r for r in csv.DictReader(f)}
X = {}
for r in eq:
    if r["differential_or_algebraic"] == "differential" and int(r["bus"]) in M:
        X.setdefault(int(r["bus"]), []).append((int(r["state_index"]), float(r["equilibrium_value"])))
def local_sg_current(b):
    m = M[b]
    x = [v for _, v in sorted(X[b])][-6:]
    psi2q, psi2d, ed, eq, w, delta = x
    xd1 = (float(m["X″_d"]) - float(m["X_ls"])) / (float(m["X′_d"]) - float(m["X_ls"]))
    xq1 = (float(m["X″_q"]) - float(m["X_ls"])) / (float(m["X′_q"]) - float(m["X_ls"]))
    cd = xd1 * eq + (1 - xd1) * psi2d
    cq = -xq1 * ed + (1 - xq1) * psi2q
    d, q = math.sin(delta), math.cos(delta)
    vd = d * V[b].real - q * V[b].imag
    vq = q * V[b].real + d * V[b].imag
    rs, xqpp, xdpp = [float(m[k]) for k in ("R_s", "X″_q", "X″_d")]
    a, c = -vd - w * cq, -vq + w * cd
    determinant = rs * rs + w * w * xqpp * xdpp
    id_ = (rs * a + w * xqpp * c) / determinant
    iq = (-w * xdpp * a + rs * c) / determinant
    factor = float(m["Sn"]) / 100
    return factor * complex(d * id_ + q * iq, -q * id_ + d * iq)
for b in range(30, 40):
    iload = L.get(b, 0j).conjugate() / V[b].conjugate() if False else (L.get(b, 0j) / V[b]).conjugate()
    igen = -Iline[b] - iload
    sgen = V[b] * igen.conjugate()
    print(b, "sg_MW_from_KCL", round(sgen.real * 100, 4),
          "sg_Mvar_from_KCL", round(sgen.imag * 100, 4),
          "sg_MW_from_state", round((V[b] * local_sg_current(b).conjugate()).real * 100, 4),
          "busrow_P_MW", round(float(B[b]["P"]) * 100, 4),
          "load_MW", round(-L.get(b, 0j).real * 100, 4))
    iload_actual = -Iline[b] - local_sg_current(b)
    print("  inferred_load_MW", round((V[b] * iload_actual.conjugate()).real * 100, 4))
