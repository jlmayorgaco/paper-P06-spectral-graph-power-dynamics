"""Interval checks for explicitly gauge-restored binary64 export, not physics uncertainty."""
from moments import *
sys.path.insert(0,str(ROOT/'experiments/regional_paper_closure_20261003/vendor'))
from flint import arb,arb_mat,ctx
ctx.prec=192
def mat(a):return arb_mat([[arb(float(x)) for x in row] for row in np.asarray(a)])
def main():
    m=Moments();A=mat(m.A);B=mat(m.B);C=mat(m.C);D=mat(m.D)
    Ai=A.inv();X=Ai*B;G0=D+C*X;G1=-(C*Ai*X)
    one=arb_mat([[1] for _ in range(10)]);projector=arb_mat(10,10)
    for i in range(10):
        for j in range(10):projector[i,j]=(1 if i==j else 0)-arb(1)/10
    K0=-(G0*projector);K1=-G1
    minor=arb_mat([[K0[i,j] for j in range(9)] for i in range(9)])
    determinant=minor.det();assert not determinant.contains(0)
    L=arb_mat([[K0[j,i] for j in range(10)] for i in range(9)]+[[1]*10])
    ell=L.inv()*arb_mat([[0] for _ in range(9)]+[[1]])
    d=(ell.transpose()*K1*one)[0,0];assert not d.contains(0)
    ws=[ell[i,0]/d for i in range(10)]
    finite=all(Ai[i,j].is_finite() for i in range(174) for j in range(174))
    assert finite
    rows=[dict(bus=i+30,interval=w.str(35),sign='positive' if w>0 else 'negative' if w<0 else 'unresolved',
               midpoint=float(w.mid()),floating_difference=float(w.mid())-m.weights[i]) for i,w in enumerate(ws)]
    pd.DataFrame(rows).to_csv(OUT/'TABLE_16_INTERVAL_WEIGHTS.csv',index=False)
    result={'scope':'exact binary64 matrices exported by Julia; G0 restored by right multiplication I-11T/10',
        'precision_bits':ctx.prec,'hidden_inverse_enclosed':finite,'gauge_nine_minor_nonzero':True,
        'gauge_minor_determinant':determinant.str(35),'d_nonzero':True,'d_interval':d.str(35),
        'all_weight_signs_enclosed':all(x['sign']!='unresolved' for x in rows),
        'max_float_weight_error':max(abs(x['floating_difference']) for x in rows),
        'physical_uncertainty_certificate':False,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dump('MOMENT_ASSUMPTION_CERTIFICATE.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
