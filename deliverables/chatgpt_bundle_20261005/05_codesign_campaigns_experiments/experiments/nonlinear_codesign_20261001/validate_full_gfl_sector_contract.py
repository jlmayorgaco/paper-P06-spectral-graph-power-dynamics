"""Nonlinear algebra, convex-hull and trajectory audit, independent of the proof."""
from pathlib import Path
from itertools import product
import hashlib
import json
import tomllib
import numpy as np
from scipy.integrate import solve_ivp
from full_gfl_sector_contract import (Specification, state_scale, state_parameters,
    parameter_bounds, physical_matrices, exact_rhs, output_rows)

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports/nonlinear_codesign_20261001/full_gfl_sector_contract"
label="final"
syn=json.loads((OUT/f"synthesis_{label}.json").read_text(encoding="utf-8"))
verification=json.loads((OUT/f"verification_{label}.json").read_text(encoding="utf-8"))
models=tomllib.loads((OUT/"models.toml").read_text(encoding="utf-8"))
assert verification["status"]=="ALL_LOCAL_CERTIFICATES_VERIFIED"
assert verification["synthesis_sha256"]==hashlib.sha256((OUT/f"synthesis_{label}.json").read_bytes()).hexdigest()
devices={d["bus"]:d for d in models["devices"]}
certificates={d["bus"]:d for d in verification["certificates"]}
s=Specification(**syn["specification"])
scale=state_scale(s)
rng=np.random.default_rng(20261001)


def actual_rhs(d,kp,ki,z,epsilon,w):
    p=d["parameters"]
    de,om,xi,chi,gq,gd,di,dq,dvi,dv=z
    id_=d["id0"]+di;iq=d["iq0"]+dq
    vdc=p["Vdc"]+dv;vdi=d["x0"][7]+dvi
    ed=p["dc_kp"]*dv+vdi-id_;eq=p["iset_q"]-iq
    vd=p["cc_kp"]*ed+p["cc_ki"]*(d["x0"][1]+gd)
    vq=p["cc_kp"]*eq+p["cc_ki"]*(d["x0"][0]+gq)
    ud=d["V0"]*np.cos(de)+epsilon*w[0]
    uq=-d["V0"]*np.sin(de)+epsilon*w[1]
    a=p["omega_base"]/p["Xf"];frequency=p["omega_base"]*p["omega_frame"]+om
    return np.array([om-epsilon*s.frame_rate_per_voltage*w[2],
        (xi+kp*uq-om)/p["pll_tau"],ki*uq,(om-chi)/s.measurement_tau,
        eq,ed,a*(vd-ud-p["Rf"]*id_)+frequency*iq,
        a*(vq-uq-p["Rf"]*iq)-frequency*id_,p["dc_ki"]*dv,
        (p["Pdc"]+epsilon*s.dc_power_per_voltage*w[3]-vd*id_-vq*iq)/(p["Cdc"]*vdc)])


points=[];traces=[];histories={}
for candidate in syn["candidates"]:
    bus=candidate["bus"];d=devices[bus];cert=certificates[bus]
    kp,ki=candidate["Kp"],candidate["Ki"]
    eps=cert["epsilon_actual_model_verified"]
    X=np.array(candidate["X_scaled"]);P=np.linalg.inv(X);L=np.linalg.cholesky(X)
    bounds=parameter_bounds(s,d)
    raw_vertices=[physical_matrices(s,d,kp,ki,pp) for pp in product(*bounds)]
    max_rhs=max_hull=max_trim=0.;max_boundary=-np.inf
    for repeat in range(200):
        z=scale*rng.uniform(-.9,.9,10)
        pp=state_parameters(s,d,z)
        A,B=physical_matrices(s,d,kp,ki,pp)
        w=rng.normal(size=4);w/=max(1,np.linalg.norm(w))
        direct=exact_rhs(s,d,kp,ki,z,eps,w)
        max_rhs=max(max_rhs,float(np.max(np.abs(direct-(A@z+eps*B@w)))))
        coordinates=[(p-lo)/(hi-lo) for p,(lo,hi) in zip(pp,bounds)]
        assert min(coordinates)>=-1e-10 and max(coordinates)<=1+1e-10
        weights=[np.prod([u if k else 1-u for u,k in zip(coordinates,bits)])
                 for bits in product([0,1],repeat=5)]
        Ah=sum(weight*mat[0] for weight,mat in zip(weights,raw_vertices))
        Bh=sum(weight*mat[1] for weight,mat in zip(weights,raw_vertices))
        max_hull=max(max_hull,float(np.max(np.abs(A-Ah))),float(np.max(np.abs(B-Bh))))
        residual=actual_rhs(d,kp,ki,z,eps,w)-direct
        residual_bound=np.array(cert["physical_trim_residual_bounds"])
        # Direct floating-point subtraction has an arithmetic floor; the formal
        # residual enclosure itself was proved separately with exact fractions.
        max_trim=max(max_trim,float(np.max(np.abs(residual)-residual_bound)))
        direction=rng.normal(size=10);direction/=np.linalg.norm(direction)
        y=L@direction
        fa=actual_rhs(d,kp,ki,scale*y,eps,w)
        max_boundary=max(max_boundary,float(2*y@P@(fa/scale)))
    assert max_rhs<2e-8 and max_hull<2e-8 and max_trim<2e-8 and max_boundary<1e-6
    points.append(dict(bus=bus,samples=200,max_exact_factorization_error=max_rhs,
                       max_convex_hull_error=max_hull,max_trim_enclosure_excess_float=max_trim,
                       max_sampled_boundary_derivative=max_boundary))

    def forcing(t):
        w=np.array([np.sin(3*t),np.sin(5*t+t*t),np.cos(2*t),np.sin(7*t+.4)])
        return w/max(1,np.linalg.norm(w))
    for kind in (["origin","current_boundary"] if bus in (30,39) else ["origin"]):
        c=np.eye(10)[6]
        y0=np.zeros(10) if kind=="origin" else .999*(X@c)/np.sqrt(c@X@c)
        rhs=lambda t,y:actual_rhs(d,kp,ki,scale*y,eps,forcing(t))/scale
        sol=solve_ivp(rhs,(0,2),y0,method="Radau",rtol=2e-8,atol=1e-10,max_step=.01,dense_output=True)
        assert sol.success
        t=np.linspace(0,2,2001);yy=sol.sol(t);zz=scale[:,None]*yy
        storage=np.einsum("ij,ji->i",yy.T@P,yy)
        norm_outputs=output_rows(s)@yy
        vdc=d["parameters"]["Vdc"]+zz[9]
        current=np.hypot(d["id0"]+zz[6],d["iq0"]+zz[7])
        metric=dict(bus=bus,initial=kind,max_storage=float(storage.max()),
                    max_normalized_output=float(np.max(np.abs(norm_outputs))),
                    min_vdc=float(vdc.min()),max_vdc=float(vdc.max()),max_current=float(current.max()))
        assert metric["max_storage"]<=1+1e-6 and metric["max_normalized_output"]<=1+1e-6
        traces.append(metric)
        histories[f"bus{bus}_{kind}"]=np.column_stack([t,storage,current,vdc,zz[1]/(2*np.pi),(zz[1]-zz[3])/(2*np.pi*s.measurement_tau)])
    print(json.dumps(dict(bus=bus,status="NONLINEAR_CHECK_PASSED",factorization=max_rhs,hull=max_hull)),flush=True)
np.savez_compressed(OUT/"nonlinear_histories.npz",**histories)
result=dict(status="FULL_GFL_NONLINEAR_IMPLEMENTATION_CHECKS_PASSED",point_checks=points,trajectories=traces,
            synthesis_sha256=hashlib.sha256((OUT/f"synthesis_{label}.json").read_bytes()).hexdigest(),
            verification_sha256=hashlib.sha256((OUT/f"verification_{label}.json").read_bytes()).hexdigest(),
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            scope="Implementation checks; exact rational certificates, not these finite tests, establish the conditional input guarantee.")
(OUT/"nonlinear_validation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
