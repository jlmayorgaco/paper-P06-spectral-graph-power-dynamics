module LocalOracle
using LinearAlgebra, TOML, SHA, CSV, DataFrames
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
const N=PDExactDesignN
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
const CTX=N.design_context(ROOT)
const BETA=1.6991206999182038e-6
const EVALUATIONS=Ref(0)
const ROBUST_FACTORIZATIONS=Ref(0)
const OUT=get(ENV,"EXPQ2B_CERT_OUT",joinpath(ROOT,"reports","experiment_Q2B","CERTIFIED_SEARCH"))
mkpath(OUT)
BLAS.set_num_threads(1)

struct Architecture
    support::Vector{Int}
    Q::Matrix{Float64}
    angle::Matrix{Float64}
    input::Vector{Float64}
    c::Vector{Float64}
end
function encode(s,d)
    e=1 .-Float64.(d["rho"])
    vcat(e[s.-29],(d["Kp"].-CTX.kpmin)./(CTX.kpmax.-CTX.kpmin),
        (d["Ki"].-CTX.kimin)./(CTX.kimax.-CTX.kimin))
end
function decode(a,x)
    n=length(a.support);rho=ones(10);rho[a.support.-29].=1 .-x[1:n]
    kp=CTX.kpmin.+x[n+1:n+10].*(CTX.kpmax.-CTX.kpmin)
    ki=CTX.kimin.+x[n+11:n+20].*(CTX.kimax.-CTX.kimin)
    (;rho,kp,ki)
end
function architecture(s,d)
    s=sort(Int.(s));m=N.descriptor(CTX,d["rho"],d["Kp"],d["Ki"])
    g=N.gauge_vector(m);Q=nullspace(reshape(g/norm(g),1,:));angle=zeros(10,78)
    for bus in 30:39
        v=CTX.net.voltage[bus]
        angle[bus-29,2bus-1]=-imag(v)/abs2(v)
        angle[bus-29,2bus]=real(v)/abs2(v)
    end
    input=LinearSecurity.load_input_vector(CTX,m,16)
    Architecture(s,Q,angle,input,vcat(CTX.power[s.-29]./1000,zeros(20)))
end

function matrices(a,x;derivatives=false)
    p=decode(a,x);m=N.descriptor(CTX,p.rho,p.kp,p.ki);Q=a.Q
    gy=lu(m.Gy);Y=gy\m.C;L=gy\a.input;bf=-m.B*L
    A=Q'*m.Ared*Q;B=Q'*bf;C=-a.angle*Y*m.Ared*Q/(2pi)
    D=-a.angle*Y*bf/(2pi);jump=-a.angle*L
    base=(;A,B,C,D,jump,m,p)
    derivatives || return base
    ds=NamedTuple[];ns=length(a.support)
    for j in eachindex(x)
        da=zeros(size(m.A));db=zeros(size(m.B));dc=zeros(size(m.C));dd=zeros(size(m.D))
        if j<=ns
            bus=a.support[j];yi=2bus-1:2bus
            for (b,r) in zip(m.blocks,eachrow(m.state_map))
                b.bus==bus || continue
                xi=Int(r.first):Int(r.last);sgn=b.kind=="SG" ? 1. : -1.
                dc[yi,xi].+=sgn.*b.J.C;dd[yi,yi].+=sgn.*b.J.D
            end
        else
            gain=j-ns;ii=mod1(gain,10);bus=ii+29
            k=findfirst(b->b.bus==bus&&b.kind=="GFL",m.blocks)
            r=m.state_map[k,:];xi=Int(r.first):Int(r.last);yi=2bus-1:2bus
            op=N.trim_gfl(CTX,bus);theta=op.x[3]
            row=gain<=10 ? 4 : 5
            scale=gain<=10 ? (CTX.kpmax[ii]-CTX.kpmin[ii])/op.pars.pll_tau : CTX.kimax[ii]-CTX.kimin[ii]
            da[xi[row],xi[3]]=scale*(-cos(theta)*op.u[1]-sin(theta)*op.u[2])
            db[xi[row],yi[1]]=-scale*sin(theta);db[xi[row],yi[2]]=scale*cos(theta)
        end
        dy=gy\(dc-dd*Y);dl=-(gy\(dd*L))
        dar=da-db*Y-m.B*dy;dbf=-db*L-m.B*dl
        push!(ds,(A=Q'*dar*Q,B=Q'*dbf,
            C=-a.angle*(dy*m.Ared+Y*dar)*Q/(2pi),
            D=-a.angle*(dy*bf+Y*dbf)/(2pi),jump=-a.angle*dl))
    end
    merge(base,(;ds))
end

function phi2(z,t)
    t<=0 && return zero(z)
    u=z*t
    abs(u)<1e-3 ? t^2*(.5+u/6+u^2/24+u^3/120+u^4/720+u^5/5040) : (expm1(u)-u)/z^2
end
function dphi2(z,t)
    t<=0 && return zero(z)
    u=z*t
    abs(u)<1e-3 ? t^3*(1/6+u/12+u^2/40+u^3/180+u^4/1008) : ((u-2)*exp(u)+u+2)/z^3
end
phi1(z,t)=t<=0 ? zero(z) : expm1(z*t)/z
function coeff(lam,t,kind;derivative=false)
    if isinf(t)
        return derivative ? 1 ./lam.^2 : -1 ./lam
    end
    T=.5;fn=derivative ? dphi2 : phi2
    kind==:F ? (fn.(lam,t).-fn.(lam,t-T))./T :
        (fn.(lam,t).-2fn.(lam,t-T).+fn.(lam,t-2T))./T^2
end
function coefficients(t,kind)
    T=.5
    isinf(t) && return (1.,0.)
    kind==:F ? (min(t,T)/T, Float64(t<T)/(2pi*T)) :
        ((t-2max(t-T,0.)+max(t-2T,0.))/T^2,
            (Float64(t<T)-Float64(T<=t<2T))/(2pi*T^2))
end
function coeff_t(lam,t,kind)
    T=.5
    kind==:F ? (phi1.(lam,t).-phi1.(lam,t-T))./T :
        (phi1.(lam,t).-2phi1.(lam,t-T).+phi1.(lam,t-2T))./T^2
end
function output(v,t,kind,ch)
    if kind==:F && t>=.5
        return v.dc[ch]+100real(sum(v.R[ch,:].*tailcoeff(v.lam,t,kind)))
    elseif kind==:R && t>=1.
        return 100real(sum(v.R[ch,:].*tailcoeff(v.lam,t,kind)))
    end
    d,j=coefficients(t,kind)
    100real(sum(v.R[ch,:].*coeff(v.lam,t,kind))+v.D[ch]*d+v.jump[ch]*j)
end
function output_dt(v,t,kind,ch)
    if (kind==:F && t>=.5) || (kind==:R && t>=1.)
        return 100real(sum(v.R[ch,:].*v.lam.*tailcoeff(v.lam,t,kind)))
    end
    T=.5;d=kind==:F ? Float64(t<T)/T : (Float64(t<T)-Float64(T<=t<2T))/T^2
    100real(sum(v.R[ch,:].*coeff_t(v.lam,t,kind))+v.D[ch]*d)
end
function tailcoeff(lam,t,kind;derivative=false)
    T=.5;em=expm1.(lam*T)
    if kind==:F
        f=exp.(lam*(t-T)).*em./(T*lam.^2)
        derivative ? f.*(t-T .-2 ./lam).+exp.(lam*t)./lam.^2 : f
    else
        f=exp.(lam*(t-2T)).*em.^2 ./ (T^2*lam.^2)
        derivative ? f.*(t-2T .-2 ./lam).+2exp.(lam*(t-T)).*em./(T*lam.^2) : f
    end
end
function bisect(f,a,b)
    fa=f(a)
    for _ in 1:45
        c=(a+b)/2;fc=f(c)
        (b-a)<1e-9 && return c
        if signbit(fc)==signbit(fa);a=c;fa=fc;else;b=c;end
    end
    (a+b)/2
end
function extrema(v,kind;horizon=60.,fine=false)
    grid=sort(unique(vcat(0.,10. .^range(-6,log10(.02),length=32),
        collect(.02:(fine ? .01 : .05):horizon),.5-1e-10,.5,.5+1e-10,1-1e-10,1.,1+1e-10,horizon)))
    H=reduce(hcat,[coeff(v.lam,t,kind) for t in grid]);HT=reduce(hcat,[coeff_t(v.lam,t,kind) for t in grid])
    cd=[coefficients(t,kind)[1] for t in grid];cj=[coefficients(t,kind)[2] for t in grid]
    ddt=kind==:F ? Float64.(grid.<.5)./.5 : (Float64.(grid.<.5).-Float64.((grid.>=.5).&(grid.<1)))./.25
    vals=100real.(v.R*H+v.D*cd'+v.jump*cj');dvals=100real.(v.R*HT+v.D*ddt')
    late=findall(t->t>=(kind==:F ? .5 : 1.),grid)
    tail=reduce(hcat,[tailcoeff(v.lam,grid[k],kind) for k in late])
    vals[:,late]=100real.(v.R*tail)
    kind==:F && (vals[:,late].+=v.dc)
    dvals[:,late]=100real.(v.R*(tail.*v.lam))
    result=NamedTuple[]
    for ch in 1:10
        # Keep the finite extrema and the asymptotic plateau as separate
        # inequalities. A last sample at the arbitrary horizon is not an
        # extremum and must not become an almost-duplicate active surface.
        endpoints=findall(t->t==0. || abs(t-.5)<2e-10 || abs(t-1.)<2e-10,grid)
        idx=endpoints[argmax(abs.(vals[ch,endpoints]))]
        best=abs(vals[ch,idx]);bt=grid[idx];sign=signbit(vals[ch,idx]) ? -1. : 1.
        for k in 1:length(grid)-1
            grid[k]<.5<=grid[k+1] && continue
            grid[k]<1<=grid[k+1] && continue
            signbit(dvals[ch,k])==signbit(dvals[ch,k+1]) && continue
            max(abs(vals[ch,k]),abs(vals[ch,k+1]))<best*.85 && continue
            t=bisect(t->output_dt(v,t,kind,ch),grid[k],grid[k+1]);f=output(v,t,kind,ch)
            if abs(f)>best;best=abs(f);bt=t;sign=signbit(f) ? -1. : 1.;end
        end
        push!(result,(;value=best,time=bt,sign,channel=ch,kind))
    end
    result
end

function robust_minimum(A,ws)
    function at(w)
        ROBUST_FACTORIZATIONS[]+=1
        F=svd(inv(Matrix(im*w*I-A-.05I)))
        (;omega=w,beta=1/F.S[1],F,domega=real(im*dot(F.V[:,1],F.U[:,1])))
    end
    points=at.(ws);initial=minimum(p.beta for p in points)
    for j in eachindex(ws)
        w=ws[j];w<=1e-8 && continue
        points[j].beta>max(5initial,10BETA) && continue
        width=max(.1,.03w);left=at(max(0.,w-width));right=at(w+width)
        for _ in 1:5
            left.domega<0 && right.domega>0 && break
            width*=2;left=at(max(0.,w-width));right=at(w+width)
        end
        if left.domega<0 && right.domega>0
            for _ in 1:40
                mid=at((left.omega+right.omega)/2)
                if mid.domega<0;left=mid;else;right=mid;end
                right.omega-left.omega<1e-9 && break
            end
            push!(points,at((left.omega+right.omega)/2))
        end
    end
    points[argmin([p.beta for p in points])]
end

function structural_rows(v)
    # In a connected synchronized equilibrium, all bus phase rates are the
    # same. Remove only this known DC redundancy; retain the full vector in
    # every feasibility evaluation. A transient peak at infinity is the same
    # DC inequality and is not counted twice for LICQ.
    uniform=maximum(v.dc)-minimum(v.dc)<1e-8
    rows=uniform ? [1,2,2+argmax(abs.(v.dc))] : collect(1:12)
    append!(rows,13:22)
    append!(rows,23:32)
    rows
end

function evaluate(a,x;derivatives=false,fine=false)
    EVALUATIONS[]+=1
    b=matrices(a,x;derivatives);E=eigen(b.A);lam=E.values;V=E.vectors;Vi=inv(V);Z=Vi*b.B
    CV=b.C*V;R=CV.*reshape(Z,1,:);alpha=maximum(real.(lam))
    ws=unique(vcat(0.,abs.(imag.(lam[sortperm(real.(lam);rev=true)[1:min(8,length(lam))]]))))
    # Inverse iteration/SVD avoids loss of relative accuracy in the tiny last
    # singular value of a matrix whose largest entry is O(10^6).
    robust=robust_minimum(b.A,ws);beta=robust.beta
    v=merge(b,(;lam,V,Vi,Z,CV,R,alpha,beta,omega=robust.omega,Fs=robust.F,
        beta_frequency_stationarity=abs(robust.domega)))
    dc=100real.(-b.C*(b.A\b.B)+b.D)
    v=merge(v,(;dc))
    fp=extrema(v,:F;fine);rp=extrema(v,:R;fine)
    # Every local output has its own inequality; structural ties stay explicit.
    g=vcat((alpha+.05)/.05,(BETA-beta)/BETA,abs.(dc)./.5 .-1,
        [p.value/.5-1 for p in fp],[p.value/.5-1 for p in rp])
    v=merge(v,(;g,fp,rp,dc,J=1000dot(a.c,x),
        Fpeak=max(maximum(abs.(dc)),maximum(p.value for p in fp)),Rpeak=maximum(p.value for p in rp)))
    derivatives || return v
    n=length(x);Jac=zeros(length(g),n);j=argmax(real.(lam));r=V[:,j];l=Vi[j,:]
    sf=robust.F;u=sf.V[:,1];w=sf.U[:,1];ab=b.A\b.B
    dcJ=zeros(10,n)
    for col in 1:n
        z=b.ds[col]
        Jac[1,col]=real(sum(l.*(z.A*r)))/.05
        Jac[2,col]=real(dot(u,z.A*w))/BETA
        dc_z=100real.(-z.C*ab+b.C*(b.A\(z.A*ab-z.B))+z.D)
        dcJ[:,col]=dc_z
        Jac[3:12,col].=sign.(dc).*dc_z./.5
    end
    # Frechet divided differences of the matrix functions at each active peak.
    cache=Dict{Tuple{Symbol,Float64},Matrix{ComplexF64}}()
    for (row,p) in enumerate(vcat(fp,rp))
        if isinf(p.time)
            Jac[12+row,:].=Jac[2+p.channel,:];continue
        end
        late=(p.kind==:F && p.time>=.5)||(p.kind==:R && p.time>=1.)
        fn=late ? tailcoeff : coeff
        fs=fn(lam,p.time,p.kind);dfs=fn(lam,p.time,p.kind;derivative=true)
        divided=get!(cache,(p.kind,p.time)) do
            [abs(lam[i]-lam[j])<1e-7*max(1.,abs(lam[i])) ? (dfs[i]+dfs[j])/2 :
                (fs[i]-fs[j])/(lam[i]-lam[j]) for i in eachindex(lam),j in eachindex(lam)]
        end
        # Contract the Frechet operator before differentiating each parameter.
        L=(reshape(CV[p.channel,:],:,1).*divided).*reshape(Z,1,:)
        weight=transpose(Vi)*L*transpose(V)
        d,jump=coefficients(p.time,p.kind)
        late && (d=0.;jump=0.)
        for col in 1:n
            z=b.ds[col]
            val=sum((z.C[p.channel,:]'*V)[:].*fs.*Z)+sum(CV[p.channel,:].*fs.*(Vi*z.B))+
                sum(weight.*z.A)+z.D[p.channel]*d+z.jump[p.channel]*jump
            Jac[12+row,col]=p.sign*(100real(val)+(late && p.kind==:F ? dcJ[p.channel,col] : 0.))/.5
        end
    end
    merge(v,(;Jac))
end

end
