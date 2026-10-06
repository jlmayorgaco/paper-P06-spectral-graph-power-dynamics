# Reuse the explicit QP and NNLS algebra (no external nonlinear optimizer).
legacy=read(joinpath(L.ROOT,"experiments/bnd_expQ2B/local_certificate/solve.jl"),String)
include_string(Main,"function nnls"*split(split(legacy,"function nnls";limit=2)[2],"function main()";limit=2)[1],"reused_active_set_algebra")
function audit(a,x,v;tol=2e-7)
    rows=Vector{Float64}[];vals=Float64[];names=String[];ids=Int[]
    for k in eachindex(v.g)
        v.g[k]>=-tol || continue
        push!(rows,v.Jac[k,:]);push!(vals,v.g[k]);push!(names,v.names[k]);push!(ids,k)
    end
    for j in eachindex(x)
        if x[j]<=tol
            r=zeros(length(x));r[j]=-1;push!(rows,r);push!(vals,-x[j]);push!(names,"lower$j");push!(ids,-j)
        end
        if x[j]>=1-tol
            r=zeros(length(x));r[j]=1;push!(rows,r);push!(vals,x[j]-1);push!(names,"upper$j");push!(ids,-j-length(x))
        end
    end
    A=isempty(rows) ? zeros(0,length(x)) : reduce(vcat,transpose.(rows))
    scales=max.([norm(r) for r in rows],1e-12);As=A./scales
    mu=isempty(rows) ? Float64[] : nnls(As',-a.c;tol=1e-12)./scales
    (;stationarity=norm(a.c+A'*mu,Inf),primal=max(maximum(v.g),maximum(-x),maximum(x.-1),0.),
        complementarity=maximum(abs.(mu.*vals);init=0.),mu,A,As,names,ids,values=vals,
        rank=rank(As;rtol=1e-9),active_count=length(rows))
end
function save_point(name,a,x,v,au;status="EXPLORATORY")
    p=L.decode(a,x)
    data=Dict("support"=>a.support,"rho"=>p.rho,"Kp"=>p.kp,"Ki"=>p.ki,"J_MW"=>v.J,
        "alpha"=>v.alpha,"beta_observed"=>v.beta,"Fpeak"=>v.Fpeak,"Rpeak"=>v.Rpeak,
        "stationarity"=>au.stationarity,"primal"=>au.primal,"complementarity"=>au.complementarity,
        "local_certified"=>false,"status"=>status)
    open(joinpath(L.OUT,name),"w") do io;TOML.print(io,data);end
end
function retract(a,xt,vref,qmu;iterations=4)
    xx=copy(xt);active=findall(qmu.>1e-10)
    isempty(active) && return xx
    free=findall((xx.>1e-10).&(xx.<1-1e-10))
    for _ in 1:iterations
        vv=L.multi_evaluate(a,xx;derivatives=true)
        residual=matched_values(vv,vref)[active]
        maximum(abs.(residual))<2e-9 && return xx
        J=L.matched_jacobian(vv,vref)[active,free]
        sc=max.(sqrt.(sum(abs2,J;dims=2))[:],1e-12)
        correction=pinv(J./sc;rtol=1e-11)*(-residual./sc)
        fraction=1.
        for (k,j) in enumerate(free)
            correction[k]>0 && (fraction=min(fraction,(1-xx[j])/correction[k]))
            correction[k]<0 && (fraction=min(fraction,-xx[j]/correction[k]))
        end
        xx[free].+=min(1.,.99fraction).*correction
    end
    xx
end
function solve_multi(a,x0;maxiter=250,label="MULTIPEAK_SOC")
    x=copy(x0);ns=length(a.support);n=length(x);lo=vcat(fill(1e-9,ns),zeros(20));hi=vcat(fill(1-1e-9,ns),ones(20))
    v=L.multi_evaluate(a,x;derivatives=true)
    H=Matrix(Diagonal(vcat(fill(.2,ns),fill(.001,20))));trust=.01;penalty=1000.
    history=NamedTuple[];started=time();status="MAX_ITERATIONS"
    for it in 1:maxiter
        au=audit(a,x,v)
        push!(history,(;iteration=it,retained_MW=v.J,primal=au.primal,stationarity=au.stationarity,
            complementarity=au.complementarity,alpha=v.alpha,beta=v.beta,F=v.Fpeak,R=v.Rpeak,
            trust,elapsed_s=time()-started,active=join(au.names,";")))
        CSV.write(joinpath(L.OUT,"$(label)_HISTORY.csv"),DataFrame(history))
        save_point("$(label)_CURRENT.toml",a,x,v,au)
        println("MULTI it=$it J=$(v.J) primal=$(au.primal) KKT=$(au.stationarity) nactive=$(au.active_count) trust=$trust");flush(stdout)
        if au.primal<2e-8 && au.stationarity<1e-6 && au.complementarity<1e-8
            status="FIRST_ORDER_CANDIDATE";break
        end
        if any(x[1:ns].<2e-9)||any(x[1:ns].>1-2e-9)
            status="ARCHITECTURE_EVENT";break
        end
        selected=findall(v.g.>-.15);G0=v.Jac[selected,:];g0=v.g[selected]
        scales=max.(sqrt.(sum(abs2,G0;dims=2))[:],1.)
        G=G0./scales;gs=g0./scales;m=length(gs)
        Iq=Matrix{Float64}(I,n,n)
        AH=vcat(hcat(G,-ones(m)),hcat(Iq,zeros(n)),hcat(-Iq,zeros(n)),hcat(zeros(1,n),[-1.;;]))
        bb=vcat(-gs,min.(hi-x,trust),min.(x-lo,trust),0.)
        HQ=zeros(n+1,n+1);HQ[1:n,1:n]=H;HQ[end,end]=1.
        sol=qp(HQ,vcat(a.c,penalty),AH,bb,vcat(zeros(n),max(0.,maximum(gs))+1e-10))
        if sol.status!="OK"
            println("QP $(sol.status)");trust*=.5
            trust<1e-8 && (status=sol.status;break)
            continue
        end
        step=sol.z[1:n];qmu=zeros(length(v.g))
        for (k,id) in enumerate(sol.W)
            id<=m && (qmu[selected[id]]=max(0.,sol.mu[k])/scales[id])
        end
        # A fixed set of row scales is used for an entire line search. All
        # temporal surfaces are matched by continuation, with new violated
        # surfaces included in the exact feasibility test.
        phi0=dot(a.c,x)+penalty*sum(max.(gs,0.))
        pred=dot(a.c,step)+penalty*(sum(max.(gs+G*step,0.))-sum(max.(gs,0.)))
        accepted=false
        for power in 0:18
            t=2. ^(-power);xt=clamp.(x+t*step,lo,hi)
            if power<=3
                xr=try retract(a,xt,v,qmu) catch;xt end
                norm(xr-xt,Inf)<=max(.5norm(t*step,Inf),1e-6) && (xt=xr)
            end
            vt=try L.multi_evaluate(a,xt) catch;continue end
            gt=matched_values(vt,v)
            phit=dot(a.c,xt)+penalty*sum(max.(gt[selected]./scales,0.))
            if phit<=phi0+1e-4*t*min(pred,-1e-12) && maximum(vt.g)<=max(.02,2au.primal)
                vn=L.multi_evaluate(a,xt;derivatives=true);s=xt-x
                y=(L.matched_jacobian(vn,v)-v.Jac)'*qmu
                hs=H*s;ss=dot(s,hs);sy=dot(s,y)
                if ss>1e-19
                    theta=sy>=.2ss ? 1. : .8ss/(ss-sy);r=theta*y+(1-theta)*hs
                    H=H-hs*hs'/ss+r*r'/dot(s,r);H=(H+H')/2
                    mineig=eigmin(Symmetric(H));mineig<1e-9 && (H+=(1e-9-mineig)*I)
                end
                x=xt;v=vn;accepted=true
                power==0 && (trust=min(.15,1.3trust))
                power>4 && (trust=max(1e-7,.5trust))
                break
            end
        end
        if !accepted
            trust*=.5
            if trust<1e-8;status="MERIT_STAGNATION";break;end
        end
    end
    au=audit(a,x,v);save_point("$(label)_FINAL.toml",a,x,v,au;status)
    open(joinpath(L.OUT,"$(label)_COUNTERS.toml"),"w") do io
        TOML.print(io,Dict("model_evaluations"=>L.EVALUATIONS[],"robust_factorizations"=>L.ROBUST_FACTORIZATIONS[],"elapsed_s"=>time()-started))
    end
    println("MULTI_FINAL $status J=$(v.J) KKT=$(au.stationarity) primal=$(au.primal)");flush(stdout)
    (;a,x,v,audit=au,status,elapsed_s=time()-started)
end
function matched_values(vnew,vold)
    gs=copy(vold.g);gs[1:3]=vnew.g[1:3]
    for (i,p) in enumerate(vold.surfaces)
        choices=findall(q->q.channel==p.channel&&q.kind==p.kind,vnew.surfaces)
        if isempty(choices)
            gs[i+3]=abs(L.output(vnew,p.time,p.kind,p.channel))/.5-1
        else
            k=choices[argmin([abs(vnew.surfaces[k].time-p.time) for k in choices])]
            gs[i+3]=vnew.g[k+3]
        end
    end
    gs
end
