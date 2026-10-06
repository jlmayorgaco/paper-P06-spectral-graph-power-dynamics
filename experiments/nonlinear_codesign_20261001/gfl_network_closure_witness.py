"""Exact-current-slice obstruction to composing independent local GFL envelopes.

Constructs states inside each local ellipsoid whose algebraically consistent
network voltages violate at least one local input assumption. This is a failure
of that proposed certificate composition, NOT a physical instability witness.
"""
from pathlib import Path
import csv
import hashlib
import json
import tomllib
import numpy as np
from scipy.optimize import minimize_scalar
from full_gfl_sector_contract import Specification, state_scale, STATE_NAMES

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports/nonlinear_codesign_20261001/full_gfl_sector_contract"
syn=json.loads((OUT/"synthesis_final.json").read_text(encoding="utf-8"))
ver=json.loads((OUT/"verification_final.json").read_text(encoding="utf-8"))
net=tomllib.loads((OUT/"network_closure.toml").read_text(encoding="utf-8"))
assert ver["status"]=="ALL_LOCAL_CERTIFICATES_VERIFIED"
assert ver["synthesis_sha256"]==hashlib.sha256((OUT/"synthesis_final.json").read_bytes()).hexdigest()
buses=net["buses"]
contracts={r["bus"]:r for r in syn["candidates"]}
verified={r["bus"]:r for r in ver["certificates"]}
s=Specification(**syn["specification"]);sc=state_scale(s)
Z=np.array(net["current_to_voltage"])
G=np.array(net["G"])
assert np.linalg.norm(G@Z+np.eye(len(G)),ord=np.inf)<1e-10
C=np.zeros((2,10));C[0,6]=sc[6];C[1,7]=sc[7]
rows=[];witnesses=[];matrix=np.zeros((len(buses),len(buses)))
conditional=[];transfers=[[None for j in buses] for i in buses]
for j,donor in enumerate(buses):
    X=np.array(contracts[donor]["X_scaled"])
    P=np.linalg.inv(X)
    X0=X-np.outer(X[:,0],X[0,:])/X[0,0]
    X0=(X0+X0.T)/2;X0[0,:]=0;X0[:,0]=0
    conditional.append(X0)
    theta=net["theta0"][j]
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    injected=net["rho"][j]*R@C
    for i,receiver in enumerate(buses):
        T=Z[2*i:2*i+2,2*j:2*j+2]@injected
        H=T@X0@T.T
        transfers[i][j]=T
        eig,Q=np.linalg.eigh(H)
        norm=np.sqrt(max(0,eig[-1]))
        yy=.999*(X0@T.T@Q[:,-1])/norm
        yy[0]=0
        storage=float(yy@P@yy)
        assert storage<1 and np.max(np.abs(yy[[0,1,6,7,9]]))<1
        eps=verified[receiver]["epsilon_actual_model_verified"]
        predicted=T@yy
        ratio=np.linalg.norm(predicted)/eps
        matrix[i,j]=norm/eps
        row=dict(donor_bus=donor,receiver_bus=receiver,local_storage=storage,
                 one_donor_exact_slice_gain=norm/eps,
                 witness_voltage_norm_pu=float(np.linalg.norm(predicted)),
                 receiver_input_radius_pu=eps,witness_input_usage=float(ratio))
        rows.append(row)
        witnesses.append({**row,**dict(zip(STATE_NAMES,(sc*yy).tolist())),
                          "predicted_delta_ur":float(predicted[0]),"predicted_delta_ui":float(predicted[1])})
for name,data in [("closure_matrix_entries.csv",rows),("closure_witnesses.csv",witnesses)]:
    with (OUT/name).open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
worst=max(rows,key=lambda r:r["witness_input_usage"])
joint_scale=.9/float(matrix.max())
joint_best=None
for i,receiver in enumerate(buses):
    HH=[transfers[i][j]@conditional[j]@transfers[i][j].T for j in range(len(buses))]
    def support(theta):
        direction=np.array([np.cos(theta),np.sin(theta)])
        return sum(np.sqrt(max(0,direction@H@direction)) for H in HH)
    angles=np.linspace(0,2*np.pi,2049)[:-1]
    kk=int(np.argmax([support(t) for t in angles]));tt=angles[kk];step=2*np.pi/2048
    opt=minimize_scalar(lambda t:-support(t),bounds=(tt-step,tt+step),method="bounded")
    direction=np.array([np.cos(opt.x),np.sin(opt.x)])
    yy=[]
    for j,donor in enumerate(buses):
        T=transfers[i][j];X0=conditional[j]
        y=joint_scale*X0@T.T@direction/np.sqrt(direction@HH[j]@direction)
        y[0]=0;yy.append(y)
    delta=sum(transfers[i][j]@yy[j] for j in range(len(buses)))
    ratio=np.linalg.norm(delta)/verified[receiver]["epsilon_actual_model_verified"]
    single_max=max(np.linalg.norm(transfers[ii][jj]@yy[jj])/verified[bb]["epsilon_actual_model_verified"]
                   for ii,bb in enumerate(buses) for jj in range(len(buses)))
    stores=[float(y@np.linalg.solve(np.array(contracts[bb]["X_scaled"]),y)) for bb,y in zip(buses,yy)]
    assert single_max<.901 and max(stores)<1
    if joint_best is None or ratio>joint_best["joint_input_usage"]:
        joint_best=dict(receiver_bus=receiver,joint_input_usage=float(ratio),
                        largest_single_input_usage=float(single_max),
                        common_state_scale=joint_scale,max_local_storage=max(stores),
                        predicted_delta_ur=float(delta[0]),predicted_delta_ui=float(delta[1]),
                        receiver_input_radius_pu=verified[receiver]["epsilon_actual_model_verified"])
        joint_states=yy
joint_rows=[dict(donor_bus=bb,receiver_bus=joint_best["receiver_bus"],
                 **dict(zip(STATE_NAMES,(sc*y).tolist()))) for bb,y in zip(buses,joint_states)]
with (OUT/"joint_closure_witness.csv").open("w",newline="",encoding="utf-8") as f:
    writer=csv.DictWriter(f,fieldnames=list(joint_rows[0]));writer.writeheader();writer.writerows(joint_rows)
result=dict(status="LOCAL_ENVELOPES_DO_NOT_COMPOSE_ON_THIS_NETWORK" if worst["witness_input_usage"]>1 else "SINGLE_DONOR_CHECK_INCONCLUSIVE",
            worst=worst,violating_donor_receiver_pairs=sum(r["witness_input_usage"]>1 for r in rows),
            total_pairs=len(rows),matrix=matrix.tolist(),buses=buses,joint_witness=joint_best,
            synthesis_sha256=hashlib.sha256((OUT/"synthesis_final.json").read_bytes()).hexdigest(),
            verification_sha256=hashlib.sha256((OUT/"verification_final.json").read_bytes()).hexdigest(),
            network_closure_sha256=hashlib.sha256((OUT/"network_closure.toml").read_bytes()).hexdigest(),
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            scope="Counterexample to closure of the independent local input contracts; no unstable trajectory or infeasible replacement capacity is asserted")
(OUT/"network_closure_witness.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps({k:result[k] for k in ("status","worst","violating_donor_receiver_pairs","total_pairs")},indent=2))
print(json.dumps(joint_best,indent=2))
