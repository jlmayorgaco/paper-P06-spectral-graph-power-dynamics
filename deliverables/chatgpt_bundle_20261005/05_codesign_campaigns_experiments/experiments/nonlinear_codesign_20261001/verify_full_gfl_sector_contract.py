"""Independent rational proof for the full local GFL polytopic certificate.

All 32 vertices and model residual bounds are reconstructed from scalars. A
small input-radius derating covers finite-precision trim offsets of the actual
frozen model. No network closure or global design-optimality claim is made.
"""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import time
import tomllib
import numpy as np
from verify_pll_sector_contract import (rational as rat, zeros, transpose, multiply,
    exact_ldl_pivots, trace, dual_psd, atan_interval)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/"reports/nonlinear_codesign_20261001/full_gfl_sector_contract"


def reconstruct(spec, device, candidate):
    s = {k:rat(v) for k,v in spec.items()}
    d = {k:rat(v) for k,v in device.items() if isinstance(v, (int, float))}
    p = {k:rat(v) for k,v in device["parameters"].items()}
    kp, ki = rat(candidate["Kp"]), rat(candidate["Ki"])
    omega_scale = rat(2*np.pi*spec["frequency_max_hz"])
    scale = [s["angle_box"]]+[omega_scale]*3+[s["gamma_scale"]]*2+[
             s["id_box"], s["iq_box"], s["integral_current_scale"], s["vdc_box"]]
    pi_lower = F("3.14159265358979323846264338327950288")
    lo5, hi5 = atan_interval(F(1, 5), 40)
    lo239, hi239 = atan_interval(F(1, 239), 12)
    assert pi_lower <= 16*lo5-4*hi239
    assert omega_scale <= 2*pi_lower*s["frequency_max_hz"]
    bounds = [(1-s["angle_box"]**2/6,F(1)),(-s["angle_box"]/2,s["angle_box"]/2),
              (-s["id_box"],s["id_box"]),(-s["iq_box"],s["iq_box"]),
              (1/(p["Vdc"]+s["vdc_box"]),1/(p["Vdc"]-s["vdc_box"]))]
    assert 0 < s["angle_box"]**2 < 6 and 0 < s["vdc_box"] < p["Vdc"]
    cd, cq = [F(0)]*10, [F(0)]*10
    cd[5],cd[6],cd[8],cd[9] = p["cc_ki"],-p["cc_kp"],p["cc_kp"],p["cc_kp"]*p["dc_kp"]
    cq[4],cq[7] = p["cc_ki"],-p["cc_kp"]
    vv=[]
    for ss,cc,did,diq,beta in product(*bounds):
        A,B=zeros(10,10),zeros(10,4)
        ac,wb=p["omega_base"]/p["Xf"],p["omega_base"]*p["omega_frame"]
        A[0][1]=F(1);B[0][2]=-s["frame_rate_per_voltage"]
        A[1][0]=-kp*d["V0"]*ss/p["pll_tau"];A[1][1]=-1/p["pll_tau"];A[1][2]=1/p["pll_tau"]
        B[1][1]=kp/p["pll_tau"]
        A[2][0]=-ki*d["V0"]*ss;B[2][1]=ki
        A[3][1]=1/s["measurement_tau"];A[3][3]=-1/s["measurement_tau"]
        A[4][7]=F(-1)
        A[5][6]=F(-1);A[5][8]=F(1);A[5][9]=p["dc_kp"]
        A[6]=[ac*c for c in cd];A[6][6]-=ac*p["Rf"];A[6][7]+=wb
        A[6][0]-=ac*d["V0"]*cc;A[6][1]+=d["iq0"]+diq;B[6][0]=-ac
        A[7]=[ac*c for c in cq];A[7][7]-=ac*p["Rf"];A[7][6]-=wb
        A[7][0]+=ac*d["V0"]*ss;A[7][1]-=d["id0"]+did;B[7][1]=-ac
        A[8][9]=p["dc_ki"]
        hh=[d["id0"]*cd[j]+d["iq0"]*cq[j] for j in range(10)]
        hh[6]+=d["vid0"];hh[7]+=d["viq0"]
        A[9]=[-beta/p["Cdc"]*(hh[j]+did*cd[j]+diq*cq[j]) for j in range(10)]
        B[9][3]=beta/p["Cdc"]*s["dc_power_per_voltage"]
        A=[[s["time_scale"]*A[i][j]*scale[j]/scale[i] for j in range(10)] for i in range(10)]
        B=[[s["time_scale"]*B[i][j]/scale[i] for j in range(4)] for i in range(10)]
        vv.append((A,B))
    C=zeros(6,10)
    for k,i in enumerate([0,1,6,7,9]): C[k][i]=F(1)
    C[5][1]=s["frequency_max_hz"]/(s["measurement_tau"]*s["rocof_max_hz_s"])
    C[5][3]=-C[5][1]
    return s,d,p,scale,vv,C


def verify(spec, device, candidate):
    start=time.monotonic()
    s,d,p,scale,vv,C=reconstruct(spec,device,candidate)
    X=[[rat(a) for a in row] for row in candidate["X_scaled"]]
    eps=rat(candidate["epsilon_sdp"])
    g=s["matrix_guard"]; alpha=s["decay_rate"]*s["time_scale"]
    xpiv=exact_ldl_pivots([[X[i][j]-(g if i==j else 0) for j in range(10)] for i in range(10)])
    outmargins=[1-multiply(multiply([c],X),transpose([c]))[0][0] for c in C]
    assert all(a>0 for a in outmargins)
    assert trace(X)<=s["trace_bound"] and eps<=s["epsilon_bound"]
    minpiv=None
    guarded=True
    for A,B in vv:
        AX,XAT=multiply(A,X),multiply(X,transpose(A))
        neg=zeros(14,14)
        for i in range(10):
            for j in range(10):neg[i][j]=-(AX[i][j]+XAT[i][j]+2*alpha*X[i][j])
            for j in range(4):neg[i][10+j]=neg[10+j][i]=-eps*B[i][j]
        for j in range(4):neg[10+j][10+j]=2*alpha
        piv=exact_ldl_pivots(neg)
        minpiv=min(piv) if minpiv is None else min(minpiv,min(piv))
        try:
            exact_ldl_pivots([[neg[i][j]-(g if i==j else 0) for j in range(14)] for i in range(14)])
        except AssertionError:
            guarded=False

    # The ideal incremental model subtracts its nominal constants. Bound the
    # actual model's tiny trim offsets for every state in the declared box.
    x0=[rat(a) for a in device["x0"]]
    ed0=x0[7]-d["id0"];eq0=p["iset_q"]-d["iq0"]
    vd0=p["cc_kp"]*ed0+p["cc_ki"]*x0[1]
    vq0=p["cc_kp"]*eq0+p["cc_ki"]*x0[0]
    ac,wb=p["omega_base"]/p["Xf"],p["omega_base"]*p["omega_frame"]
    residual=[F(0)]*10
    residual[4]=abs(eq0);residual[5]=abs(ed0)
    residual[6]=abs(ac*(vd0-d["V0"]-p["Rf"]*d["id0"])+wb*d["iq0"])
    residual[7]=abs(ac*(vq0-p["Rf"]*d["iq0"])-wb*d["id0"])
    numerator=(abs(p["Pdc"]-vd0*d["id0"]-vq0*d["iq0"])
               +abs(vd0-d["vid0"])*s["id_box"]+abs(vq0-d["viq0"])*s["iq_box"])
    residual[9]=numerator/(p["Cdc"]*(p["Vdc"]-s["vdc_box"]))
    r2=sum((r*r/(sc*sc) for r,sc in zip(residual,scale)),F(0))
    beta=F(999,1000)
    # M <= ||S^-1 residual||/sqrt(g), so compare squared nonnegative bounds.
    inward=s["decay_rate"]*(1-beta*beta)
    assert r2/g < inward*inward
    eps_actual=beta*eps

    # Exact upper bound on the ideal (possibly guarded) fixed-gain SDP.
    ZX,sx=dual_psd(candidate["dual_X_lower"])
    lam=[max(F(0),rat(a)) for a in candidate["dual_outputs"]]
    lamT=max(F(0),rat(candidate["dual_trace"]))
    lamE=max(F(0),rat(candidate["dual_epsilon_upper"]))
    R=[[ZX[i][j]-sum((lam[k]*C[k][i]*C[k][j] for k in range(len(C))),F(0))
        -(lamT if i==j else 0) for j in range(10)] for i in range(10)]
    se=1-lamE
    const=-g*trace(ZX)+(1-g)*sum(lam)+lamT*s["trace_bound"]+lamE*s["epsilon_bound"]
    shifts=[sx]
    for (A,B),rawZ in zip(vv,candidate["dual_vertices"]):
        Z,shift=dual_psd(rawZ);shifts.append(shift)
        Z11=[row[:10] for row in Z[:10]]
        ZA,ATZ=multiply(Z11,A),multiply(transpose(A),Z11)
        for i in range(10):
            for j in range(10):R[i][j]-=ZA[i][j]+ATZ[i][j]+2*alpha*Z11[i][j]
        se-=2*sum((Z[i][10+j]*B[i][j] for i in range(10) for j in range(4)),F(0))
        const+=2*alpha*sum((Z[i][i] for i in range(10,14)),F(0))-g*trace(Z)
    rnorm=max(sum((abs(a) for a in row),F(0)) for row in R)
    upper=const+s["trace_bound"]*rnorm+s["epsilon_bound"]*max(F(0),se)
    if guarded and all(a>=g for a in outmargins): assert upper>=eps
    return dict(status="FULL_LOCAL_GFL_INVARIANCE_EXACTLY_VERIFIED_WITH_TRIM_RESIDUAL",
        bus=device["bus"],Kp=candidate["Kp"],Ki=candidate["Ki"],
        epsilon_centered_verified=float(eps),epsilon_actual_model_verified=float(eps_actual),
        epsilon_actual_model_exact=str(eps_actual),input_derating=float(beta),
        vertices_verified=len(vv),min_negative_LDL_pivot=float(minpiv),
        min_X_guarded_LDL_pivot=float(min(xpiv)),
        physical_trim_residual_bounds=[float(a) for a in residual],
        storage_trim_disturbance_bound=math.sqrt(float(r2/g)),
        inward_margin_reserved=float(inward),
        output_margins=[float(a) for a in outmargins],
        guarded_SDP_primal_verified=guarded and all(a>=g for a in outmargins),
        ideal_fixed_gain_SDP_upper=float(upper),
        ideal_fixed_gain_relative_gap=float((upper-eps)/upper),
        exact_dual_upper=str(upper),
        max_dual_psd_repair=float(max(shifts)),
        elapsed=time.monotonic()-start,
        scope="All local GFL states + instrument; joint input ball; network interconnection, SGs, hardware ratings and replacement capacity not certified")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--label",default="all_network_gains")
    args=parser.parse_args()
    source=OUT/f"synthesis_{args.label}.json"
    syn=json.loads(source.read_text(encoding="utf-8"))
    models_path=OUT/"models.toml"
    assert syn["models_sha256"]==hashlib.sha256(models_path.read_bytes()).hexdigest()
    models=tomllib.loads(models_path.read_text(encoding="utf-8"))
    devices={d["bus"]:d for d in models["devices"]}
    result=dict(status="VERIFICATION_IN_PROGRESS",certificates=[],
                synthesis_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for candidate in syn["candidates"]:
        try:
            assert "X_scaled" in candidate, "No primal SDP candidate; solver status: "+candidate["status"]
            r=verify(syn["specification"],devices[candidate["bus"]],candidate)
        except AssertionError as exc:
            r=dict(status="EXACT_VERIFICATION_FAILED",bus=candidate["bus"],error=str(exc))
        result["certificates"].append(r)
        print(json.dumps({k:r.get(k) for k in ("bus","status","epsilon_actual_model_verified","guarded_SDP_primal_verified","ideal_fixed_gain_relative_gap","elapsed","error")}),flush=True)
        (OUT/f"verification_{args.label}.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    result["status"]="ALL_LOCAL_CERTIFICATES_VERIFIED" if all(r["status"].startswith("FULL_LOCAL") for r in result["certificates"]) else "SOME_CERTIFICATES_FAILED"
    (OUT/f"verification_{args.label}.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
