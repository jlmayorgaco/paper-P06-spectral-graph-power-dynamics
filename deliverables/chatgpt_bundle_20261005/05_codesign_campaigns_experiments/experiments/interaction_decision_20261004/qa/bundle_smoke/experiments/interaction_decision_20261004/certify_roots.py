"""Continuous Arb root counts, with balanced open-loop resolvent Taylor balls.

Only the computational enclosure uses a convergent resolvent series. The physical
delay remains exp(-s*tau). The source CSV/gain floats are exact binary64 inputs.
"""
from model import *
import sys,time
sys.path.insert(0,str(ROOT/'experiments/regional_paper_closure_20261003/vendor'))
from flint import acb,arb,acb_mat,ctx
ctx.prec=160

def conv(a):return acb_mat([[acb(float(x)) for x in row] for row in np.asarray(a)])
def norminf(A):return max(sum((abs(A[i,j]) for j in range(A.ncols())),arb(0)).upper() for i in range(A.nrows()))
def bounds(z):return {'real':z.real.str(30),'imag':z.imag.str(30)}

class BallModel:
    def __init__(self,p):
        F={k:conv(v) for k,v in M.items()};self.p=p
        G=acb_mat(20,20);N=acb_mat(20,204)
        for i in range(20):
            w=1-arb(float(p[i//2]))
            for j in range(20):G[i,j]=F['Y'][i,j]+w*F['Ds'][i,j]
            for j in range(204):N[i,j]=w*F['Cs'][i,j]+(1-w)*F['Cf'][i,j]
        vx=-(G.inv()*N);A=F['Adev']+F['Bv']*vx;C=F['Etheta']+F['Hv']*vx
        B=acb_mat(204,10)
        for i in range(204):
            for j in range(10):B[i,j]=F['Bp'][i,j]*arb(float(p[10+j]))+F['Bi'][i,j]*arb(float(p[20+j]))
        # Binary powers make this similarity exact and avoid physical unit norms.
        _,(sc,_)=la.matrix_balance(Model(p).A,permute=False,separate=True)
        self.A=acb_mat([[A[i,j]*arb(float(sc[j]))/arb(float(sc[i])) for j in range(204)] for i in range(204)])
        self.B=acb_mat([[B[i,j]/arb(float(sc[i])) for j in range(10)] for i in range(204)])
        self.C=acb_mat([[C[i,j]*arb(float(sc[j])) for j in range(204)] for i in range(10)])
        self.tau=[arb(float(v)) for v in TAU]
        self.G_regular=True
    def prepare(self,center,radius,order=5):
        self.center=acb(float(center.real),float(center.imag));self.radius=arb(float(radius));self.order=order
        DD=-self.A
        for i in range(204):DD[i,i]+=self.center
        Q=DD.inv();self.qnorm=arb(norminf(Q));self.rb=[Q*self.B]
        for _ in range(order+1):self.rb.append(Q*self.rb[-1])
        self.rnorm=arb(norminf(self.rb[0]))
        assert arb(2).sqrt()*self.radius*self.qnorm<1,'open-loop chart not enclosed'
        # Pre-multiplication by C is exact ball arithmetic; remainder bound also
        # receives the C infinity norm. This reduces every panel to 10x10.
        self.cr=[self.C*r for r in self.rb];self.cnorm=arb(norminf(self.C))
    def matrix(self,s):
        d=s-self.center;q=abs(d)*self.qnorm
        assert q<1
        R=acb_mat(10,10);Rp=acb_mat(10,10);power=acb(1)
        for k in range(self.order+1):
            R+=self.cr[k]*power
            Rp-=self.cr[k+1]*(power*(k+1))
            power*=-d
        er=self.cnorm*self.rnorm*q**(self.order+1)/(1-q)
        ed=self.cnorm*self.rnorm*self.qnorm*q**(self.order+1)*((self.order+2)/(1-q)+q/(1-q)**2)
        K=acb_mat(10,10);Ks=acb_mat(10,10)
        for i in range(10):
            for j in range(10):
                rr=R[i,j]+acb(arb(0,er.upper()),arb(0,er.upper()))
                rp=Rp[i,j]+acb(arb(0,ed.upper()),arb(0,ed.upper()))
                e=(-s*self.tau[j]).exp()
                K[i,j]=(1 if i==j else 0)-rr*e
                Ks[i,j]=(-rp+self.tau[j]*rr)*e
        return K,Ks
    def count(self,n):
        r=self.radius;center=self.center;total=acb(0);panels=[]
        corners=[(-1,-1),(1,-1),(1,1),(-1,1),(-1,-1)]
        for a,b in zip(corners,corners[1:]):
            for k in range(n):
                ar=arb(a[0])+arb(b[0]-a[0])*k/n;ai=arb(a[1])+arb(b[1]-a[1])*k/n
                br=arb(a[0])+arb(b[0]-a[0])*(k+1)/n;bi=arb(a[1])+arb(b[1]-a[1])*(k+1)/n
                ds=r*acb(br-ar,bi-ai)
                s=center+r*acb((ar+br)/2+arb(0,(abs(br-ar)/2).upper()),(ai+bi)/2+arb(0,(abs(bi-ai)/2).upper()))
                K,Ks=self.matrix(s)
                midpoint=acb_mat([[K[i,j].mid() for j in range(10)] for i in range(10)])
                R=midpoint.inv()
                W=(R*K).inv()*(R*Ks)
                val=sum((W[i,i] for i in range(10)),acb(0))*ds
                total+=val;panels.append(bounds(val))
        count=total/(2*arb.pi()*acb(0,1))
        proved=bool(count.contains(acb(1)) and count.real>arb('.5') and count.real<arb('1.5'))
        return count,proved,panels

def certify(name,p,seed):
    z,res=Model(p).refine(seed);start=time.perf_counter();bm=BallModel(p)
    radius=1e-7;bm.prepare(z,radius)
    for n in [8,16,32,64,128]:
        try:count,ok,panels=bm.count(n)
        except (ValueError,ZeroDivisionError):continue
        if ok:break
    else:raise RuntimeError('root count not isolated')
    result={'name':name,'status':'PROVED_EXPORTED_MODEL_ROOT_ENCLOSURE','center':[z.real,z.imag],
        'radius':radius,'real_lower':(arb(float(z.real))-arb(float(radius))).str(30),
        'real_upper':(arb(float(z.real))+arb(float(radius))).str(30),
        'count':bounds(count),'count_one_proved':ok,'panels':4*n,
        'open_loop_resolvent_norm_upper':bm.qnorm.str(20),
        'source_semantics':'CSV and design parsed as exact binary64; delays exact binary64 0.04, no physical uncertainty',
        'seconds':time.perf_counter()-start,'p':p.tolist()}
    save('CERT_'+name+'.json',result);save('CERT_PANELS_'+name+'.json',panels)
    print(name,'CERTIFIED',count,'seconds',result['seconds'],flush=True)
    return result

def main():
    from compensated import ANCHOR,TARGET,combined
    selected=json.loads((OUT/'damping_candidate_6.json').read_text());p=np.array(selected['p'])
    singles=[]
    for i in [0,7]:
        pp=ANCHOR.copy();ix=[i,10+i,20+i];pp[ix]=p[ix];singles.append(pp)
    targets=[('anchor',ANCHOR,TARGET),('single30',singles[0],-.06570892701739738+30.93216223759933j),
        ('single37',singles[1],-.06570892701739738+30.55842598246854j),
        ('joint',p,complex(*selected['root'])),('complex_pair',combined((0,7),.01),TARGET)]
    for name,pp,z in targets:
        if not (OUT/('CERT_'+name+'.json')).exists():certify(name,pp,z)
if __name__=='__main__':main()
