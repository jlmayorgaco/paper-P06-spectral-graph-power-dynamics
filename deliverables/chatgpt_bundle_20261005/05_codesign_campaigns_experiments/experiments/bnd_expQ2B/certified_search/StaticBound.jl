module StaticBound
using LinearAlgebra
import ..StaticFrequency as S
function bernstein_transform!(b)
    for axis in 1:10
        stride=3^(axis-1);block=3stride
        for start in 0:block:length(b)-1, j in 1:stride
            i=start+j
            b[i+stride]=2b[i+stride]-(b[i]+b[i+2stride])/2
        end
    end
    b
end
function coefficients(m,weights)
    n=3^10;D=zeros(n);N=zeros(n);cost=zeros(n)
    for k in 0:n-1
        e=[Float64((k÷3^(i-1))%3)/2 for i in 1:10]
        val=S.determinants(m,e);D[k+1]=val.D;N[k+1]=val.N;cost[k+1]=dot(weights,e)
    end
    scale=max(maximum(abs.(D)),maximum(abs.(N)))
    (;D=bernstein_transform!(D./scale),N=bernstein_transform!(N./scale),cost,scale)
end
# Every nonnegative multiplier gives a valid dual lower bound for exact
# Bernstein coefficients. Numerical coefficients are explicitly provisional.
function dual_bound(c,P;cycles=15)
    mu=zeros(size(P,2));best=minimum(c)
    for _ in 1:cycles, j in eachindex(mu)
        base=c-P*mu+P[:,j]*mu[j]
        at(t)=begin
            vals=base-t*P[:,j];k=argmin(vals)
            (vals[k],-P[k,j])
        end
        lo=0.;hi=max(1.,mu[j]);value,slope=at(0.)
        if slope<=0;mu[j]=0.;continue;end
        for k in 1:60
            at(hi)[2]<=0 && break
            hi*=2
        end
        for k in 1:65
            mid=(lo+hi)/2;value,slope=at(mid)
            if slope>0;lo=mid;else;hi=mid;end
        end
        mu[j]=(lo+hi)/2
        best=max(best,minimum(c-P*mu))
    end
    (;bound=best,multipliers=mu)
end
function root_bounds(b)
    P=hcat(.5b.D+b.N,.5b.D-b.N)
    pos=dual_bound(b.cost,P);neg=dual_bound(b.cost,-P)
    (;positive=pos,negative=neg,lower=min(pos.bound,neg.bound))
end
end
