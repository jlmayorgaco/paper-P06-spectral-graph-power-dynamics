module JointKKT

using LinearAlgebra, DataFrames, Statistics
using ..ExpP
const N=ExpP.PDExactDesignN
const SIGMA=0.05

export solve_fixed_support, candidate_kkt_audit, graph_feature_audit

function _point(ctx,buses,epsv,kp,ki)
    rho=ones(10)
    for (j,bus) in enumerate(buses)
        rho[bus-29]=1-epsv[j]
    end
    N.spectrum(ctx,rho,kp,ki)
end

function _active_jacobian(ctx,buses,epsv,kp,ki,sp,eps_scale;tol_active=1e-6)
    active=findall(real.(sp.lambda).>=sp.alpha-tol_active)
    G=zeros(Float64,length(active),1+20)
    diagnostics=NamedTuple[]
    rho=ones(10)
    for (j,bus) in enumerate(buses);rho[bus-29]=1-epsv[j];end
    for (r,idx) in enumerate(active)
        d=N.simple_mode_sensitivities(ctx,rho,kp,ki;mode=idx)
        mismatch=minimum(abs.(sp.lambda.-d.lambda))
        mismatch<1e-7 || error("active pole branch/order mismatch $mismatch")
        G[r,1]=-real(d.rho[buses[1]-29])*eps_scale
        for i in 1:10
            G[r,1+i]=real(d.Kp[i])*N.K0P
            G[r,11+i]=real(d.Ki[i])*N.K0I
        end
        push!(diagnostics,(;index=idx,lambda=d.lambda,condition=d.condition,
            pairing_mismatch=mismatch))
    end
    active,G,DataFrame(diagnostics)
end

"""Primal active-set QP for a convex, positive-definite regularized SQP model."""
function _solve_qp_active(H,c,A,b;tol=1e-10,maxiter=500)
    n=length(c);d=zeros(Float64,n)
    minimum(b-A*d)>=-1e-9 || error("QP start is not feasible")
    slack=b-A*d
    W=findall(abs.(slack).<=1e-13)
    history=NamedTuple[]
    for it in 1:maxiter
        grad=H*d+c
        if isempty(W)
            p=-(H\grad);λ=Float64[]
        else
            Aw=A[W,:]
            K=[H transpose(Aw);Aw zeros(length(W),length(W))]
            sol=K\vcat(-grad,zeros(length(W)))
            p=sol[1:n];λ=sol[n+1:end]
        end
        if norm(p,Inf)<=tol
            if isempty(λ) || minimum(λ)>=-1e-9
                stat=norm(H*d+c+(isempty(W) ? zeros(n) : transpose(A[W,:])*λ),Inf)
                push!(history,(;iteration=it,step=norm(d,Inf),working_set=length(W),
                    minimum_multiplier=isempty(λ) ? NaN : minimum(λ),stationarity=stat,
                    outcome="QP_OPTIMAL"))
                return (;d,working_set=W,multipliers=λ,history=DataFrame(history),
                    stationarity=stat,status="QP_OPTIMAL")
            end
            k=argmin(λ)
            deleteat!(W,k)
            push!(history,(;iteration=it,step=norm(d,Inf),working_set=length(W),
                minimum_multiplier=minimum(λ),stationarity=NaN,outcome="DROP_NEGATIVE_MULTIPLIER"))
            continue
        end
        α=1.0;blocker=0
        inW=Set(W)
        for k in eachindex(b)
            k in inW && continue
            den=dot(view(A,k,:),p)
            if den>tol
                ak=(b[k]-dot(view(A,k,:),d))/den
                if ak<α
                    α=max(0.0,ak);blocker=k
                end
            end
        end
        d .+= α.*p
        blocker!=0 && !(blocker in W) && push!(W,blocker)
        push!(history,(;iteration=it,step=norm(d,Inf),working_set=length(W),
            minimum_multiplier=isempty(λ) ? NaN : minimum(λ),stationarity=NaN,
            outcome=blocker==0 ? "FULL_STEP" : "ADD_BLOCKER"))
    end
    (;d,working_set=W,multipliers=Float64[],history=DataFrame(history),
      stationarity=NaN,status="QP_ITERATION_LIMIT")
end

function _all_bounds(ctx,buses,eps_scale)
    n=21
    lo=fill(-Inf,n);hi=fill(Inf,n)
    for j in eachindex(buses)
        lo[j]=1e-8/eps_scale;hi[j]=1/eps_scale
    end
    for i in 1:10
        lo[1+i]=ctx.kpmin[i]/N.K0P;hi[1+i]=ctx.kpmax[i]/N.K0P
        lo[11+i]=ctx.kimin[i]/N.K0I;hi[11+i]=ctx.kimax[i]/N.K0I
    end
    lo,hi
end

function _decode(y,buses,eps_scale)
    epsv=y[1:length(buses)].*eps_scale
    kp=y[1 .+ (1:10)].*N.K0P
    ki=y[11 .+ (1:10)].*N.K0I
    epsv,kp,ki
end

function _encode(epsv,kp,ki,eps_scale)
    vcat(epsv./eps_scale,kp./N.K0P,ki./N.K0I)
end

"""Joint ε/Kp/Ki active-set trust-region continuation on one fixed support."""
function solve_fixed_support(ctx,buses,eps0,kp0,ki0;guard=1e-9,maxiter=12,
                             tol_active=1e-6,trust_radius=0.05)
    length(buses)==length(eps0) || throw(DimensionMismatch("support/epsilon mismatch"))
    eps_scale=minimum(eps0)
    y=_encode(eps0,kp0,ki0,eps_scale)
    lo,hi=_all_bounds(ctx,buses,eps_scale)
    all(lo.-1e-10 .<= y .<= hi.+1e-10) || error("initial point outside fixed-support bounds")
    weight=sum(ctx.power[bus-29]*eps0[j] for (j,bus) in enumerate(buses))
    c=zeros(length(y))
    for (j,bus) in enumerate(buses);c[j]=ctx.power[bus-29]*eps_scale;end
    H=Matrix{Float64}(I,length(y),length(y))
    history=NamedTuple[];qp_hist=DataFrame[]
    accepted=0;termination="MAX_ITERATIONS"
    for it in 1:maxiter
        epsv,kp,ki=_decode(y,buses,eps_scale)
        sp=_point(ctx,buses,epsv,kp,ki)
        active,G,branch=_active_jacobian(ctx,buses,epsv,kp,ki,sp,eps_scale;tol_active)
        # Every complete-spectrum pole is checked after every corrector step;
        # the local QP keeps all modes currently within tol_active active.
        g=real.(sp.lambda[active]).+SIGMA
        g .+= guard
        n=length(y);Id=Matrix{Float64}(I,n,n)
        A=vcat(G,Id,-Id,Id,-Id)
        b=vcat(-g,hi-y,y-lo,fill(trust_radius,n),fill(trust_radius,n))
        qp=_solve_qp_active(H,c,A,b)
        append!(qp_hist,[qp.history])
        normstep=norm(qp.d,Inf)
        push!(history,(;iteration=it,cost=weight,alpha=sp.alpha,
            guarded_constraint=maximum(real.(sp.lambda))+SIGMA+guard,
            active_modes=join(active,";"),active_mode_count=length(active),
            step_inf=normstep,qp_status=qp.status,full_spectrum_pass=sp.alpha<=-SIGMA-guard+1e-10,
            accepted=false,line_search_scale=0.0,stationarity_qp=qp.stationarity))
        if normstep<1e-9
            termination="QP_STATIONARY";break
        end
        if qp.status!="QP_OPTIMAL"
            termination=qp.status;break
        end
        committed=false
        for scale in (1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625)
            yt=y+scale*qp.d
            all(lo.-1e-11 .<= yt .<= hi.+1e-11) || continue
            et,kpt,kit=_decode(yt,buses,eps_scale)
            spt=_point(ctx,buses,et,kpt,kit)
            cst=sum(ctx.power[bus-29]*et[j] for (j,bus) in enumerate(buses))
            feas=spt.alpha<=-N.SIGMA-guard+1e-10
            if feas && cst<weight-1e-13
                y=yt;weight=cst;accepted+=1;committed=true
                history[end]=merge(history[end],(;accepted=true,line_search_scale=scale,
                    cost=weight,alpha=spt.alpha,full_spectrum_pass=true))
                break
            end
        end
        if !committed
            termination="LINE_SEARCH_NO_FEASIBLE_DESCENT";break
        end
    end
    epsv,kp,ki=_decode(y,buses,eps_scale)
    sp=_point(ctx,buses,epsv,kp,ki)
    active,G,branch=_active_jacobian(ctx,buses,epsv,kp,ki,sp,eps_scale;tol_active)
    (;buses,epsilon=epsv,rho=1 .- epsv,Kp=kp,Ki=ki,
      retained_SG_MW=sum(ctx.power[b-29]*epsv[j] for (j,b) in enumerate(buses)),
      converted_GFL_MW=sum(ctx.power)-sum(ctx.power[b-29]*epsv[j] for (j,b) in enumerate(buses)),
      alpha=sp.alpha,lambda=sp.lambda,active,G,branch,history=DataFrame(history),
      qp_history=isempty(qp_hist) ? DataFrame() : vcat(qp_hist...),accepted_steps=accepted,
      termination,guard,full_spectrum_pass=sp.alpha<=-SIGMA-guard+1e-10)
end

"""Numerical KKT, LICQ and critical-cone/SOSC audit at a fixed-support point."""
function candidate_kkt_audit(ctx,bus,epsilon,kp,ki;guard=1e-9,tol_active=1e-6)
    rho=ones(10);rho[bus-29]=1-epsilon
    sp=N.spectrum(ctx,rho,kp,ki)
    active,G,branch=_active_jacobian(ctx,[bus],[epsilon],kp,ki,sp,epsilon;tol_active)
    # ExpN's candidate is expected to have one simple rightmost finite pole.
    length(active)==1 || return (;status="ACTIVE_MODE_SWITCH",active,G,branch,alpha=sp.alpha)
    grow=vec(G[1,:]);g=sp.alpha+SIGMA+guard
    grad=zeros(21);grad[1]=ctx.power[bus-29]*epsilon
    μ=-grad[1]/grow[1]
    lag=grad+μ*grow
    lower=[ctx.kpmin./N.K0P;ctx.kimin./N.K0I]
    upper=[ctx.kpmax./N.K0P;ctx.kimax./N.K0I]
    y=[1.0;kp./N.K0P;ki./N.K0I]
    lower_mult=Float64[];upper_mult=Float64[];stationarity=Float64[]
    for j in 1:10
        k=1+j
        if abs(y[k]-lower[j])<1e-8
            push!(lower_mult,lag[k]);push!(upper_mult,0.0);push!(stationarity,abs(lag[k]-max(lag[k],0.0)))
        elseif abs(y[k]-upper[j])<1e-8
            push!(lower_mult,0.0);push!(upper_mult,-lag[k]);push!(stationarity,abs(lag[k]-min(lag[k],0.0)))
        else
            push!(lower_mult,0.0);push!(upper_mult,0.0);push!(stationarity,abs(lag[k]))
        end
    end
    for j in 1:10
        k=11+j
        if abs(y[k]-lower[10+j])<1e-8
            push!(lower_mult,lag[k]);push!(upper_mult,0.0);push!(stationarity,abs(lag[k]-max(lag[k],0.0)))
        elseif abs(y[k]-upper[10+j])<1e-8
            push!(lower_mult,0.0);push!(upper_mult,-lag[k]);push!(stationarity,abs(lag[k]-min(lag[k],0.0)))
        else
            push!(lower_mult,0.0);push!(upper_mult,0.0);push!(stationarity,abs(lag[k]))
        end
    end
    Jbounds=zeros(20,21)
    for j in 1:20;Jbounds[j,j+1]=1;end
    licq=rank(vcat(transpose(grow),Jbounds))==21
    strict=all((lower_mult .> 1e-8) .| (upper_mult .> 1e-8))
    # Every gain is active at one bound. Strict box multipliers and nonzero
    # spectral ε derivative leave only dε=0 in the critical cone.
    comp=abs(μ*g)
    (;status=μ>=0 && licq && strict ? "LOCAL_KKT_CERTIFIED" : "KKT_FAIL",
      alpha=sp.alpha,active_count=length(active),active_lambda=sp.lambda[active],
      guard_constraint=g,primal_residual=max(0.0,g),multiplier_spectral=μ,
      lower_multipliers=lower_mult,upper_multipliers=upper_mult,
      stationarity_residual=max(abs(grad[1]+μ*grow[1]),maximum(stationarity)),
      complementarity_residual=comp,LICQ=licq,active_jacobian_rank=rank(G),
      SOSC=strict ? "VACUOUS_CRITICAL_CONE_STRICT_GAIN_BOUNDS" : "NOT_CERTIFIED",
      pole_condition=only(branch.condition),branch,full_spectrum_pass=sp.alpha<=-SIGMA-guard+1e-10)
end

"""Offline graph feature basis; no remote states/signals enter any PLL."""
function graph_feature_audit(ctx,kp,ki,rho)
    ys=ctx.net.y_static;n=39
    G=zeros(n,n);B=zeros(n,n)
    for i in 1:n,j in 1:n
        G[i,j]=ys[2i-1,2j-1]
        B[i,j]=ys[2i,2j-1]
    end
    Goff=copy(G);Boff=copy(B)
    for i in 1:n;Goff[i,i]=0;Boff[i,i]=0;end
    Wg=abs.(Goff);Wb=abs.(Boff)
    one=ones(n)
    raw=hcat(one,Wg*one,Wb*one,Wg*(Wb*one),Wb*(Wg*one))
    Phi=raw[30:39,:]
    for j in 2:size(Phi,2)
        col=Phi[:,j];centered=col.-mean(col)
        norm(centered)>1e-12 ? (Phi[:,j].=centered./norm(centered)) : (Phi[:,j].=0)
    end
    sv=svdvals(Phi);r=count(sv.>maximum(sv)*1e-10)
    targetp=kp./N.K0P.-1;targeti=ki./N.K0I.-1
    θp=all(targetp .== first(targetp)) ? vcat(first(targetp),zeros(size(Phi,2)-1)) : Phi\targetp
    θi=all(targeti .== first(targeti)) ? vcat(first(targeti),zeros(size(Phi,2)-1)) : Phi\targeti
    fitp=Phi*θp;fiti=Phi*θi
    residual=max(norm(targetp-fitp),norm(targeti-fiti))
    kpfit=(1 .+ fitp).*N.K0P;kifit=(1 .+ fiti).*N.K0I
    sp=N.spectrum(ctx,rho,kpfit,kifit)
    (;feature_rank=r,feature_count=size(Phi,2),parameters_per_channel=r,
      total_parameters=2r,free_parameters=20,singular_values=sv,
      Phi,theta_Kp=θp,theta_Ki=θi,fit_relative_residual=residual,
      replacement_MW=sum(ctx.power.*(1 .- rho)),alpha=sp.alpha,
      feasible=sp.alpha<=-SIGMA+1e-10)
end

end
