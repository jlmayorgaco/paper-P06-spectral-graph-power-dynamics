"""Independent exported-Jacobian check; does not change the sealed prediction."""
import sys
sys.dont_write_bytecode=True
from predict import *

def main():
    lock=json.loads((OUT/"PREDICTION_LOCK.json").read_text())
    for name,expected in lock["trial_hashes"].items():
        assert sha(OUT/name)==expected,f"Prediction changed: {name}"
    base=Model(pack(tomllib.loads((OUT/"designs/baseline.toml").read_text())))
    predictions=pd.read_csv(OUT/"PREDICTED_ROOTS.csv")
    rows=[]
    target=lock["target_root"]
    for name in ["analytic","fixed"]:
        m=Model(pack(tomllib.loads((OUT/"designs"/(name+".toml")).read_text())))
        for mat in ["A","A0","B","C"]:
            setattr(m,mat,pd.read_csv(OUT/"independent"/name/(mat+".csv")).to_numpy())
        for row in predictions[predictions.design==name].itertuples():
            old=complex(row.baseline_real,row.baseline_imag)
            sb,vb,_,_=pair(base,old);qb=vb[base.ports.pll_angle_index.to_numpy(int)]
            sp,v,_,res=pair(m,complex(row.predicted_real,row.predicted_imag))
            qp=v[m.ports.pll_angle_index.to_numpy(int)]
            mac=abs(np.vdot(qb,qp))**2/(np.vdot(qb,qb).real*np.vdot(qp,qp).real)
            movement=sp-old;prediction=complex(row.predicted_real,row.predicted_imag)-old
            rows.append(dict(design=name,root=row.root,real=sp.real,imag=sp.imag,full_residual=res,
                tracked_pattern_MAC=mac,predicted_delta_real=prediction.real,predicted_delta_imag=prediction.imag,
                actual_delta_real=movement.real,actual_delta_imag=movement.imag,
                absolute_prediction_error=abs(movement-prediction),
                relative_prediction_error=abs(movement-prediction)/max(abs(movement),1e-12),
                target=(row.root==target),target_pole_error=abs(movement) if row.root==target else None))
    frame("TABLE_09_INDEPENDENT_ROOTS.csv",rows)
    row=next(r for r in rows if r["design"]=="analytic" and r["target"])
    assert row["target_pole_error"]<=1e-5 and row["tracked_pattern_MAC"]>=.999999
    assert row["full_residual"]<=1e-9
    print(json.dumps({"independent_target_error":row["target_pole_error"],
                     "target_MAC":row["tracked_pattern_MAC"],
                     "minimum_tracked_MAC":min(r["tracked_pattern_MAC"] for r in rows)}))
if __name__=="__main__":main()
