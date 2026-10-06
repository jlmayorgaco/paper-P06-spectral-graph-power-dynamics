"""Exploratory actual-pole sensitivity attribution; see THEORY section 6."""
from analyze_ports import *
import copy

def port_derivative(m,s,t,Dp):
    D=m.delta(s,t);n,h=m.n,m.h
    cc=D[np.ix_(h,h)];cn=D[np.ix_(h,n)];nc=D[np.ix_(n,h)]
    X=la.solve(cc,cn)
    return m.M[:,None]*(Dp[np.ix_(n,n)]-Dp[np.ix_(n,h)]@X-
        nc@la.solve(cc,Dp[np.ix_(h,n)]-Dp[np.ix_(h,h)]@X))

def run():
    table=pd.read_csv(OUT/'TABLE_03_TRACKED_MODES.csv');rows=[]
    for key in ['N','T']:
        m=Model(key);t=.04
        row=table[(table.design==key)&(table.tau_ms==40)].sort_values('real').iloc[-1]
        s=complex(row.real,row.imag);v,w,res=m.mode(s,t)
        Z=m.impedance(s,t);U,sv,Vh=la.svd(Z);ell=U[:,-1];r=v[m.n];r/=la.norm(r)
        Ds=m.I+t*np.exp(-s*t)*(m.B@m.C)
        Zs=port_derivative(m,s,t,Ds);den=np.vdot(ell,Zs@r)
        for p in ['tau']+[f'{gain}_{j}' for gain in ['logKp','logKi'] for j in range(10)]:
            if p=='tau':
                Dp=s*np.exp(-s*t)*(m.B@m.C);h=1e-7
                sp=m.refine(s,t+h)[0];sm=m.refine(s,t-h)[0]
            else:
                gain,j=p.split('_');j=int(j)
                # Only the Kp row (PLL frequency) or Ki row (integrator) changes.
                idx=int(m.ports.pll_frequency_index.iloc[j])+(gain=='logKi')
                dB=np.zeros_like(m.B);dB[idx,j]=m.B[idx,j]
                Dp=-np.exp(-s*t)*(dB@m.C);h=1e-5
                mp,mm=copy.copy(m),copy.copy(m)
                mp.B=m.B+(np.exp(h)-1)*dB;mm.B=m.B+(np.exp(-h)-1)*dB
                sp=mp.refine(s,t)[0];sm=mm.refine(s,t)[0]
            full=-np.vdot(w,Dp@v)/np.vdot(w,Ds@v)
            Zp=port_derivative(m,s,t,Dp)
            port=-np.vdot(ell,Zp@r)/den
            node=-np.vdot(ell,np.diag(np.diag(Zp))@r)/den
            cross=port-node;fd=(sp-sm)/(2*h)
            rows.append(dict(design=key,parameter=p,real_s=s.real,imag_s=s.imag,
                full_real=full.real,full_imag=full.imag,port_real=port.real,port_imag=port.imag,
                nodal_real=node.real,collective_real=cross.real,fd_real=fd.real,fd_imag=fd.imag,
                port_full_relative=abs(port-full)/max(abs(full),1e-10),
                fd_absolute=abs(full-fd),fd_relative=abs(full-fd)/max(abs(full),1e-10),
                collective_fraction=abs(cross.real)/max(abs(node.real)+abs(cross.real),1e-20),
                status='EXPLORATORY_EXACT_SENSITIVITY_ATTRIBUTION'))
    pd.DataFrame(rows).to_csv(OUT/'TABLE_08_ACTION_ATTRIBUTION.csv',index=False)
    print(pd.DataFrame(rows).groupby('design')[['port_full_relative','fd_relative']].max().to_string())

if __name__=='__main__': run()
