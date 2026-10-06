"""Inspect graph coupling in the computed storage; a scoped ablation only."""
import hashlib
import json
import numpy as np
import scipy.linalg as la
from threadpoolctl import threadpool_limits
from coupled_interval_model import Model,OUT

with threadpool_limits(limits=1):
    m=Model();s=np.load(OUT/"certificate_modal.npz");H=s["H"];A=np.array(m.raw["A"])
    P=H.T@H;blocks=[];bd=np.zeros_like(P)
    for bus in range(39):
        indices=[len(m.keep)+bus,len(m.keep)+39+bus]
        if bus>=29:
            old=np.r_[m.sg[bus-29],m.gf[bus-29]]
            indices += [int(np.where(m.keep==k)[0][0]) for k in old if k!=m.ig]
        blocks.append(indices);bd[np.ix_(indices,indices)]=P[np.ix_(indices,indices)]
    assert sum(map(len,blocks))==m.nx
    ablation=la.eigvalsh((A.T@bd+bd@A)/2,bd)
    # Structural graph only. The actual dynamics continue using the full AC Y.
    W=np.zeros((39,39))
    for i in range(39):
        for j in range(i+1,39):
            W[i,j]=W[j,i]=(la.norm(m.Y[2*i:2*i+2,2*j:2*j+2],"fro")+
                           la.norm(m.Y[2*j:2*j+2,2*i:2*i+2],"fro"))/(2*np.sqrt(2))
    L=np.diag(W.sum(axis=1))-W
    Lg=L[29:,29:]-L[29:,:29]@np.linalg.solve(L[:29,:29],L[:29,29:])
    lam,U=la.eigh((Lg+Lg.T)/2)
    report=dict(status="GRAPH_AND_STORAGE_ABLATION_DIAGNOSTIC",
        graph_nodes=39,graph_edges=int(np.count_nonzero(np.triu(W,1))),
        generator_structural_laplacian_eigenvalues=lam.tolist(),
        structural_laplacian_row_sum_residual=float(np.abs(Lg.sum(axis=1)).max()),
        dropped_cross_bus_storage_logarithmic_norm=float(ablation[-1]),
        nominal_full_metric_decay_reported=float(json.loads((OUT/"synthesis_modal.json").read_text())["nominal_decay"]),
        interpretation="Removing cross-bus blocks from this P loses its Lyapunov property if the ablation norm is positive; this does not exclude all possible nodal certificates. The structural Laplacian is an interpretation coordinate, not a replacement for AC Y or a diagonalization of nonlinear dynamics.",
        model_sha256=hashlib.sha256((OUT/"model.toml").read_bytes()).hexdigest(),
        metric_sha256=hashlib.sha256((OUT/"certificate_modal.npz").read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(open(__file__,"rb").read()).hexdigest())
    np.savez_compressed(OUT/"graph_structure.npz",W=W,L=L,L_generators=Lg,U_generators=U,eigenvalues=lam)
    (OUT/"graph_structure_audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
