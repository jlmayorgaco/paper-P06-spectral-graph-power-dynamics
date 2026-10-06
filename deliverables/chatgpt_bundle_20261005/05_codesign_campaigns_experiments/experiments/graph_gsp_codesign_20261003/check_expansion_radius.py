from graph_design import *
d=tomllib.loads((OUT/'baseline.toml').read_text())
m=Model(np.r_[d['rho'],np.log(d['Kp']),np.log(d['Ki'])])
L=read('L');lam,U=la.eigh(L)
rows=[]
for tau in [0,.02,.04]:
    for f in [.1,.5,1,2,4,5,8,12]:
        Z=U.T@m.impedance(2j*np.pi*f,tau)@U
        K=np.diag(1/np.diag(Z))@(Z-np.diag(np.diag(Z)))
        radius=max(abs(la.eigvals(K)))
        rows.append(dict(tau_ms=1000*tau,frequency_hz=f,neumann_spectral_radius=radius,
            infinite_series_converges=bool(radius<1),norm_condition=la.norm(K,2)<1))
df=pd.DataFrame(rows)
df.to_csv(OUT/'TABLE_05B_NEUMANN_CONVERGENCE.csv',index=False)
print(df.to_string(index=False))
