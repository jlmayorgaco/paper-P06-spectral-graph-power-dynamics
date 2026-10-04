from close_design import *
from certify_roots import certify
PAIR_DISPATCH=float(PORTS.P0.iloc[0]+PORTS.P0.iloc[7])

def evaluate(h):
    if h==0:return ANCHOR.copy(),TARGET,0.
    pi,zi,_=local(0,2*h,h);pj,zj,_=local(7,2*h,h)
    p=pi+pj-ANCHOR;z=TARGET
    for f in np.linspace(.05,1,20):z=Model(ANCHOR+f*(p-ANCHOR)).refine(z)[0]
    return p,z,(z-zi-zj+TARGET).real

def main():
    import hashlib
    save('ADDENDUM_03_LOCK.json',{'sha256':hashlib.sha256((OUT/'ADDENDUM_03_POLICY_PATH.md').read_bytes()).hexdigest()})
    deriv=[]
    for dh in [1e-4,5e-5]:
        _,zp,ip=evaluate(dh);_,zm,im=evaluate(-dh);curv=(ip+im)/(2*dh*dh)
        deriv.append(dict(h=dh,quadratic_coefficient=curv,linear=(ip-im)/(2*dh),
            predicted_h=np.sqrt((-.05-TARGET.real)/curv) if curv>0 else None))
    save('POLICY_CURVATURE.json',deriv)
    rows=[]
    for h in np.linspace(0,.012,13):
        try:
            p,z,ii=evaluate(float(h));rows.append(dict(h=h,added_MW=PAIR_DISPATCH*h,alpha=z.real,imag=z.imag,interaction=ii,valid=True,error=''))
        except Exception as ex:rows.append(dict(h=h,valid=False,error=str(ex)))
    pd.DataFrame(rows).to_csv(OUT/'TABLE_13_POLICY_PATH.csv',index=False)
    for a,b in zip(rows,rows[1:]):
        if a['valid'] and b['valid'] and a['alpha']<-.05<b['alpha']:break
    else:save('POLICY_CROSSING.json',{'status':'NO_RESOLVED_CROSSING'});return
    lo,hi=a['h'],b['h'];history=[]
    while hi-lo>1e-8:
        mid=(lo+hi)/2;p,z,ii=evaluate(mid);history.append(dict(h=mid,alpha=z.real))
        if z.real<-.05:lo=mid
        else:hi=mid
    # Outward endpoints separated enough for root-box resolution1e-7.
    lo-=1e-6;hi+=1e-6
    plo,zlo,_=evaluate(lo);phi,zhi,_=evaluate(hi)
    certify('policy_below',plo,zlo);certify('policy_above',phi,zhi)
    save('POLICY_CROSSING.json',{'status':'NUMERICAL_FIRST_CROSSING_WITH_CERTIFIED_ENDPOINT_SIGNS',
        'h_lower':lo,'h_upper':hi,'added_MW_lower':PAIR_DISPATCH*lo,'added_MW_upper':PAIR_DISPATCH*hi,
        'alpha_lower':zlo.real,'alpha_upper':zhi.real,'curvature_prediction_h':deriv[-1]['predicted_h'],
        'gain_policy':'Kp_i=Kp_anchor*(1+2h); singleton alpha fixed via Ki_i; rho_i+=h at buses30,37',
        'scope':'not a global or local co-design maximum; interval-wide monotonicity not certified','history':history})
    print((OUT/'POLICY_CROSSING.json').read_text())
if __name__=='__main__':main()
