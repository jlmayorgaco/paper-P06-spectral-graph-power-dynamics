"""Compare direct PD and unchanged analytical models in physical state coordinates."""
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_M"
TABLES = OUT / "tables"
CASES = ("all_SG", "ExpG_candidate", "ExpK_nominal")


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pd_key(raw: str) -> tuple[int, str, str]:
    match = re.fullmatch(r"VIndex\((\d+), :(.*)\)", raw)
    if not match:
        raise ValueError(f"unknown PD symbol: {raw}")
    bus, name = int(match[1]), match[2]
    tail = name.split("₊")[-1]
    if name.startswith("busbar₊"):
        return bus, "network", "busbar_" + tail
    if "₊gov₊" in name:
        return bus, "SG", "gov_" + tail
    if "₊avr₊" in name:
        return bus, "SG", "avr_" + tail
    if "machine₊" in name:
        machine = {"ψ″_q": "psi2q", "ψ″_d": "psi2d", "E′_d": "Epd",
                   "E′_q": "Epq", "ω": "omega", "δ": "delta"}
        return bus, "SG", "machine_" + machine[tail]
    if "₊cc1₊" in name:
        return bus, "GFL", {"γ_q": "gamma_q", "γ_d": "gamma_d"}[tail]
    if "₊pll₊" in name:
        return bus, "GFL", {"θ": "theta", "Δω_rad_s": "delta_omega_rad_s",
                             "Δω_i_rad_s": "delta_omega_i_rad_s"}[tail]
    if "₊filter₊" in name or "₊v_dc_" in name:
        return bus, "GFL", tail
    raise ValueError(f"unmapped physical state: {raw}")


def an_key(row: dict) -> tuple[int, str, str]:
    return int(row["bus"]), row["kind"], row["state_name"]


def family(key: tuple[int, str, str]) -> str:
    _, kind, name = key
    if kind == "network":
        return "network"
    if kind == "SG":
        return "SG_governor" if name.startswith("gov_") else "SG_AVR" if name.startswith("avr_") else "SG_machine"
    if name in ("theta", "delta_omega_rad_s", "delta_omega_i_rad_s"):
        return "GFL_PLL"
    if name.startswith("gamma_"):
        return "GFL_current_control"
    if name.startswith("i_f_"):
        return "GFL_filter"
    return "GFL_DC"


def run_case(case: str) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict]]:
    path = OUT / "matrices" / case
    pd_states = read_csv(path / "PD_state_map.csv")
    an_states = read_csv(path / "AN_state_map.csv")
    pd_keys = [pd_key(row["state_name"]) for row in pd_states]
    an_keys = [an_key(row) for row in an_states]
    inventory = [dict(case=case, PD_states=len(pd_keys), AN_states=len(an_keys),
        PD_differential=sum(x["differential"] == "true" for x in pd_states),
        AN_differential=sum(x["differential"] == "true" for x in an_states),
        PD_algebraic=sum(x["differential"] == "false" for x in pd_states),
        AN_algebraic=sum(x["differential"] == "false" for x in an_states),
        matched_physical_states=len(set(pd_keys) & set(an_keys)),
        missing_in_analytic=len(set(pd_keys) - set(an_keys)),
        extra_in_analytic=len(set(an_keys) - set(pd_keys)))]
    pd_lookup = {key: i for i, key in enumerate(pd_keys)}
    alignment = []
    for j, key in enumerate(an_keys):
        i = pd_lookup.get(key)
        alignment.append(dict(case=case, AN_index=j+1, PD_index="" if i is None else i+1,
            bus=key[0],kind=key[1],physical_state=key[2],
            status="MISSING_IN_PD" if i is None else "MATCHED_RENAMED" if pd_states[i]["state_name"] != key[2] else "MATCHED_EXACT",
            PD_equilibrium="" if i is None else pd_states[i]["equilibrium_value"],
            AN_equilibrium=an_states[j]["equilibrium_value"],
            equilibrium_difference="" if i is None else float(pd_states[i]["equilibrium_value"])-float(an_states[j]["equilibrium_value"])))
    write_csv(TABLES / f"TABLE_M10_state_alignment_{case}.csv", alignment)
    if set(pd_keys) != set(an_keys):
        return inventory, [], [], [], []
    a_pd=np.loadtxt(path/"PD_A.csv",delimiter=",")
    m_pd=np.loadtxt(path/"PD_M.csv",delimiter=",")
    a_an=np.loadtxt(path/"AN_Ared.csv",delimiter=",")
    pd_d=np.flatnonzero(np.diag(m_pd)==1)
    pd_z=np.flatnonzero(np.diag(m_pd)==0)
    pd_red=a_pd[np.ix_(pd_d,pd_d)]-a_pd[np.ix_(pd_d,pd_z)]@solve(a_pd[np.ix_(pd_z,pd_z)],a_pd[np.ix_(pd_z,pd_d)])
    an_d=[j for j,r in enumerate(an_states) if r["differential"]=="true"]
    assert len(an_d)==a_an.shape[0]
    reduced_position={pd_keys[i]:k for k,i in enumerate(pd_d)}
    permutation=[reduced_position[an_keys[j]] for j in an_d]
    aligned=pd_red[np.ix_(permutation,permutation)]
    error=aligned-a_an
    maxidx=np.unravel_index(np.abs(error).argmax(),error.shape)
    rel=np.linalg.norm(error)/max(np.linalg.norm(a_an),1e-300)
    matrices=[dict(case=case,matrix="reduced_A",dimension=error.shape[0],
        absolute_frobenius=float(np.linalg.norm(error)),relative_frobenius=float(rel),
        max_absolute=float(np.max(np.abs(error))),max_row=int(maxidx[0])+1,
        max_column=int(maxidx[1])+1,PD_minus_AN=float(error[maxidx]),
        row_state=str(an_keys[an_d[maxidx[0]]]),column_state=str(an_keys[an_d[maxidx[1]]]),
        interpretation="Dynamical Schur reduction in physically matched coordinates")]
    full_order=[pd_lookup[key] for key in an_keys]
    for label,pd_matrix,an_matrix in (("descriptor_M",m_pd,np.loadtxt(path/"AN_M.csv",delimiter=",")),
                                      ("descriptor_A",a_pd,np.loadtxt(path/"AN_A.csv",delimiter=","))):
        delta=pd_matrix[np.ix_(full_order,full_order)]-an_matrix
        extremum=np.unravel_index(np.abs(delta).argmax(),delta.shape)
        matrices.append(dict(case=case,matrix=label,dimension=delta.shape[0],
            absolute_frobenius=float(np.linalg.norm(delta)),
            relative_frobenius=float(np.linalg.norm(delta)/max(np.linalg.norm(an_matrix),1e-300)),
            max_absolute=float(np.max(np.abs(delta))),max_row=int(extremum[0])+1,
            max_column=int(extremum[1])+1,PD_minus_AN=float(delta[extremum]),
            row_state=str(an_keys[extremum[0]]),column_state=str(an_keys[extremum[1]]),
            interpretation="Raw descriptor rows; algebraic equation normalization may differ"))
    # The official closed-network linearization is autonomous (M,A).  The
    # legacy analytical B/C/D are open bus-port maps, so a direct numerical
    # subtraction would compare different input/output definitions.
    for label in ("descriptor_B", "descriptor_C", "descriptor_D"):
        matrices.append(dict(case=case,matrix=label,dimension="NOT_DEFINED_FOR_AUTONOMOUS_PD_NETWORK",
            absolute_frobenius="",relative_frobenius="",max_absolute="",max_row="",
            max_column="",PD_minus_AN="",row_state="",column_state="",
            interpretation="NOT_COMPARABLE: official closed-network PD output is (M,A); legacy AN B/C/D are open bus-port channels"))
    groups=sorted(set(family(an_keys[j]) for j in an_d))
    blocks=[]
    for gr in groups:
        rows=[k for k,j in enumerate(an_d) if family(an_keys[j])==gr]
        for gc in groups:
            cols=[k for k,j in enumerate(an_d) if family(an_keys[j])==gc]
            sub=error[np.ix_(rows,cols)]
            ref=a_an[np.ix_(rows,cols)]
            blocks.append(dict(case=case,row_block=gr,column_block=gc,
                row_count=len(rows),column_count=len(cols),
                absolute_frobenius=float(np.linalg.norm(sub)),
                relative_frobenius=float(np.linalg.norm(sub)/max(np.linalg.norm(ref),1e-300)),
                max_absolute=float(np.max(np.abs(sub)))))
    pd_poles=np.loadtxt(path/"PD_poles.csv",delimiter=",",skiprows=1)
    an_poles=np.loadtxt(path/"AN_poles.csv",delimiter=",",skiprows=1)
    pd_l=pd_poles[:,0]+1j*pd_poles[:,1]
    an_l=an_poles[:,0]+1j*an_poles[:,1]
    # The closure is gauge quotiented once. A second pole can approach zero
    # physically, so remove exactly one smallest-magnitude PD gauge eigenvalue.
    gauge_index=int(np.argmin(np.abs(pd_l)))
    pd_l=np.delete(pd_l,gauge_index)
    if len(pd_l)!=len(an_l):
        return inventory, matrices, blocks, [], []
    rr,cc=linear_sum_assignment(np.abs(pd_l[:,None]-an_l[None,:]))
    poles=[dict(case=case,lambda_PD_real=float(pd_l[i].real),lambda_PD_imag=float(pd_l[i].imag),
        lambda_AN_real=float(an_l[j].real),lambda_AN_imag=float(an_l[j].imag),
        absolute_error=float(abs(pd_l[i]-an_l[j]))) for i,j in zip(rr,cc)]
    write_csv(TABLES/f"TABLE_M16_pole_matching_{case}.csv",poles)
    summary=[dict(case=case,PD_physical_poles=len(pd_l),AN_physical_poles=len(an_l),
        PD_alpha=float(np.max(pd_l.real)),AN_alpha=float(np.max(an_l.real)),
        max_pole_error=float(max(row["absolute_error"] for row in poles)),
        median_pole_error=float(np.median([row["absolute_error"] for row in poles])))]
    return inventory, matrices, blocks, summary, poles


def main() -> None:
    inv=[];mat=[];blocks=[];summ=[]
    for case in CASES:
        if not (OUT/"matrices"/case/"PD_A.csv").exists() or not (OUT/"matrices"/case/"AN_Ared.csv").exists():
            continue
        i,m,b,s,_=run_case(case)
        inv+=i;mat+=m;blocks+=b;summ+=s
    write_csv(TABLES/"TABLE_M09_state_inventory_comparison.csv",inv)
    write_csv(TABLES/"TABLE_M11_matrix_identity.csv",mat)
    write_csv(TABLES/"TABLE_M12_block_error_ledger.csv",blocks)
    write_csv(TABLES/"TABLE_M16_pole_matching_summary.csv",summ)
    for row in mat+summ:
        print(row)


if __name__=="__main__":
    main()
