include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, CSV, DataFrames, TOML, SHA
const L=LocalOracle

function nnls(A,b;tol=1e-10)
    n=size(A,2);x=zeros(n);P=Int[]
    for _ in 1:10n+100
        w=A'*(b-A*x);inactive=setdiff(1:n,P)
        isempty(inactive) && return x
        j=inactive[argmax(w[inactive])];w[j]<=tol && return x;push!(P,j)
        for _ in 1:n+2
            z=zeros(n);z[P]=pinv(A[:,P];rtol=1e-11)*b
            all(z[P].>0) && (x=z;break)
            bad=[j for j in P if z[j]<=0];alpha=minimum(x[j]/(x[j]-z[j]+eps()) for j in bad)
            x+=alpha*(z-x);P=[j for j in P if x[j]>tol]
        end
    end
    x
end
function audit(a,x,v;tol=1e-6)
    rows=Vector{Float64}[];values=Float64[];names=String[]
    for k in L.structural_rows(v)
        if v.g[k]>=-tol
            push!(rows,collect(v.Jac[k,:]));push!(values,v.g[k]);push!(names,"g$k")
        end
    end
    for j in eachindex(x)
        if x[j]<=tol
            r=zeros(length(x));r[j]=-1;push!(rows,r);push!(values,-x[j]);push!(names,"lower$j")
        end
        if x[j]>=1-tol
            r=zeros(length(x));r[j]=1;push!(rows,r);push!(values,x[j]-1);push!(names,"upper$j")
        end
    end
    A=isempty(rows) ? zeros(0,length(x)) : reduce(vcat,transpose.(rows))
    scales=[max(norm(r),1e-12) for r in rows];As=A./reshape(scales,:,1)
    mu=isempty(rows) ? Float64[] : nnls(As',-a.c)./scales
    (;stationarity=norm(a.c+A'*mu,Inf),complementarity=isempty(mu) ? 0. : maximum(abs.(mu.*values)),
        primal=max(maximum(v.g),maximum(-x),maximum(x.-1),0.),mu,A,names,values,
        rank=rank(As;rtol=1e-8),active_count=size(A,1))
end

function qp(H,c,A,b,z)
    W=Int[];n=length(c)
    for it in 1:500
        if isempty(W)
            p=-(H\(H*z+c));mu=Float64[]
        else
            AW=A[W,:];Z=nullspace(AW;rtol=1e-12)
            p=-Z*((Z'*H*Z)\(Z'*(H*z+c)))
            mu=pinv(AW';rtol=1e-12)*(-(H*(z+p)+c))
        end
        if norm(p,Inf)<2e-9
            if isempty(W)||minimum(mu)>=-1e-8
                return (;z,mu,W,status="OK",iterations=it)
            end
            deleteat!(W,argmin(mu));continue
        end
        alpha=1.;block=0
        for i in eachindex(b)
            i in W && continue
            den=dot(A[i,:],p)
            den<=1e-10 && continue
            t=(b[i]-dot(A[i,:],z))/den
            if t<alpha;alpha=max(0.,t);block=i;end
        end
        z+=alpha*p
        if block!=0
            if isempty(W)||rank(A[vcat(W,block),:];rtol=1e-10)>length(W)
                push!(W,block)
            else
                # Duplicate inactive surfaces can block only within roundoff.
                return (;z,mu,W,status="DEPENDENT_BLOCK",iterations=it)
            end
        end
    end
    (;z,mu=Float64[],W,status="QP_LIMIT",iterations=500)
end

function main()
    gate=CSV.read(joinpath(L.OUT,"TABLE_gradient_gate.csv"),DataFrame)
    all(gate.pass) || error("Derivative validation did not pass; optimizer blocked")
    path=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(L.ROOT,"reports","codesign_validation_20261001","candidate_repaired.toml")
    d=TOML.parsefile(path);a=L.architecture(d["support"],d);x=L.encode(a.support,d)
    ns=length(a.support);n=length(x);lo=vcat(fill(1e-8,ns),zeros(20));hi=vcat(fill(1-1e-8,ns),ones(20))
    maxiter=length(ARGS)>1 ? parse(Int,ARGS[2]) : 100
    H=Matrix(Diagonal(vcat(fill(.1,ns),fill(1e-3,20))));trust=.025;penalty=100.;hist=NamedTuple[]
    v=L.evaluate(a,x;derivatives=true);bestx=copy(x);best=v;status="MAX_ITERATIONS";start=time()
    for it in 1:maxiter
        au=audit(a,x,v)
        println("LOCAL it=",it," J=",v.J," maxg=",maximum(v.g)," KKT=",au.stationarity," F=",v.Fpeak," beta=",v.beta," trust=",trust);flush(stdout)
        push!(hist,(;iteration=it,J_MW=v.J,maxg=maximum(v.g),stationarity=au.stationarity,
            complementarity=au.complementarity,alpha=v.alpha,beta=v.beta,
            Fpeak=v.Fpeak,Rpeak=v.Rpeak,
            active=join(au.names,";"),trust,elapsed_s=time()-start))
        CSV.write(joinpath(L.OUT,"TABLE_local_sqp_history.csv"),DataFrame(hist))
        p=L.decode(a,x)
        open(joinpath(L.OUT,"CURRENT_LOCAL_CANDIDATE.toml"),"w") do io
            TOML.print(io,Dict("support"=>a.support,"rho"=>p.rho,"Kp"=>p.kp,"Ki"=>p.ki,
                "J_MW"=>v.J,"stationarity"=>au.stationarity,"iteration"=>it,"alpha"=>v.alpha,
                "beta_sampled"=>v.beta,"Fpeak"=>v.Fpeak,
                "Rpeak"=>v.Rpeak,"finite_peak_times"=>[p.time for p in v.fp],
                "dc_output_spread"=>maximum(v.dc)-minimum(v.dc),
                "eigenvector_condition"=>cond(v.V),"local_certified"=>false))
        end
        if au.stationarity<1e-6 && au.primal<1e-7 && au.complementarity<1e-7
            status="KKT_NUMERICAL_CANDIDATE";break
        end
        if any(x[1:ns].<=2e-8) || any(x[1:ns].>=1-2e-8)
            status="ARCHITECTURE_BOUNDARY_REBUILD_REQUIRED";break
        end
        rows=[j for j in L.structural_rows(v) if v.g[j]>-.3]
        # Remove duplicate local affine inequalities for the QP only.
        selected=Int[]
        for j in rows
            any(k->norm(v.Jac[j,:]-v.Jac[k,:])<1e-9 && abs(v.g[j]-v.g[k])<1e-9,selected) || push!(selected,j)
        end
        G=v.Jac[selected,:];gs=v.g[selected];scales=max.(sqrt.(sum(abs2,G;dims=2))[:],1.)
        G=G./scales;gs=gs./scales
        m=length(gs);Iq=Matrix{Float64}(I,n,n)
        AH=vcat(hcat(G,-ones(m)),hcat(Iq,zeros(n)),hcat(-Iq,zeros(n)),hcat(zeros(1,n),[-1.;;]))
        bb=vcat(-gs,min.(hi-x,trust),min.(x-lo,trust),0.)
        HQ=zeros(n+1,n+1);HQ[1:n,1:n]=H;HQ[end,end]=1.
        sol=qp(HQ,vcat(a.c,penalty),AH,bb,vcat(zeros(n),max(0.,maximum(gs))+1e-8))
        if sol.status!="OK"
            trust*=.5
            println("QP ",sol.status);flush(stdout)
            trust<1e-7 && (status=sol.status;break)
            continue
        end
        step=sol.z[1:n]
        phi0=dot(a.c,x)+penalty*sum(max.(gs,0.))
        predicted=dot(a.c,step)+penalty*(sum(max.(gs+G*step,0.))-sum(max.(gs,0.)))
        committed=false
        for power in 0:16
            alpha=2. ^(-power);xt=clamp.(x+alpha*step,lo,hi)
            vt=try L.evaluate(a,xt) catch;continue end
            # Merit includes all constraints, using the same row normalization.
            phit=dot(a.c,xt)+penalty*sum(max.(vt.g[selected]./scales,0.))
            inactive=setdiff(L.structural_rows(vt),selected)
            if phit<=phi0+1e-4*alpha*min(predicted,-1e-10) && maximum(vt.g[inactive];init=-Inf)<=1e-7
                vn=L.evaluate(a,xt;derivatives=true);s=xt-x
                qpmu=zeros(length(v.g))
                for (k,idx) in enumerate(sol.W)
                    idx<=m && (qpmu[selected[idx]]=sol.mu[k]/scales[idx])
                end
                y=(vn.Jac-v.Jac)'*qpmu
                hs=H*s;sHs=dot(s,hs);sy=dot(s,y)
                if sHs>1e-16
                    theta=sy>=.2sHs ? 1. : .8sHs/(sHs-sy)
                    r=theta*y+(1-theta)*hs
                    H=H-hs*hs'/sHs+r*r'/dot(s,r);H=(H+H')/2
                    emin=eigmin(Symmetric(H));emin<1e-8 && (H+=(1e-8-emin)*I)
                end
                x=xt;v=vn;committed=true
                if maximum(v.g)<1e-7 && (maximum(best.g)>1e-7 || v.J<best.J);bestx=copy(x);best=v;end
                power==0 && (trust=min(.2,trust*1.25))
                power>=4 && (trust=max(1e-6,trust*.5))
                break
            end
        end
        if !committed
            println("MERIT_REJECT predicted=",predicted," step_inf=",norm(step,Inf));flush(stdout)
            trust*=.5
            trust<1e-7 && (status="MERIT_STAGNATION";break)
        end
    end
    au=audit(a,x,v)
    pf=L.decode(a,x)
    open(joinpath(L.OUT,"FINAL_EXPLORATORY_CANDIDATE.toml"),"w") do io
        TOML.print(io,Dict("support"=>a.support,"rho"=>pf.rho,"Kp"=>pf.kp,"Ki"=>pf.ki,
            "J_MW"=>v.J,"alpha"=>v.alpha,"beta_observed"=>v.beta,
            "Fpeak"=>v.Fpeak,"Rpeak"=>v.Rpeak,"stationarity"=>au.stationarity,
            "primal"=>au.primal,"local_certified"=>false,"frozen_for_validation"=>false))
    end
    pb=L.decode(a,bestx)
    has_feasible=maximum(best.g)<1e-7
    open(joinpath(L.OUT,"BEST_TRACKED_POINT.toml"),"w") do io
        TOML.print(io,Dict("support"=>a.support,"rho"=>pb.rho,"Kp"=>pb.kp,"Ki"=>pb.ki,
            "J_MW"=>best.J,"alpha"=>best.alpha,"beta_sampled"=>best.beta,
            "Fpeak"=>best.Fpeak,"Rpeak"=>best.Rpeak,
            "max_normalized_constraint"=>maximum(best.g),"numerically_feasible"=>has_feasible,
            "local_certified"=>false))
    end
    open(joinpath(L.OUT,"LOCAL_ATTEMPT_STATUS.toml"),"w") do io
        TOML.print(io,Dict("status"=>status,"stationarity"=>au.stationarity,"primal"=>au.primal,
            "complementarity"=>au.complementarity,"LICQ_rank"=>au.rank,"active_rows"=>au.names,
            "active_count"=>au.active_count,"multipliers"=>au.mu,"elapsed_s"=>time()-start,
            "retained_MW"=>v.J,"best_feasible_MW"=>(has_feasible ? best.J : "NONE"),
            "has_numerically_feasible_incumbent"=>has_feasible,"local_certified"=>false,
            "global_certified"=>false,"event_bus"=>16,"disturbance_MW"=>100.,"window_s"=>.5,
            "model_evaluations"=>L.EVALUATIONS[],"robust_frequency_factorizations"=>L.ROBUST_FACTORIZATIONS[],
            "scope"=>"fixed-support analytical problem; stationary frequency refinement pending full-band audit"))
    end
    println("FINAL ",status," J=",v.J," KKT=",au.stationarity," primal=",au.primal);flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
