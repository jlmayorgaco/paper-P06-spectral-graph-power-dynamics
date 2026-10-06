module RobustAudit

using LinearAlgebra, DataFrames
using ..ExpP
const N=ExpP.PDExactDesignN

export shifted_quotient, sigma_min_at, certify_radius, robust_boundary,
       robust_kkt_audit, hinf_upper_certificate, observed_radius

function _rho(bus,epsilon)
    rho=ones(10);rho[bus-29]=1-epsilon;rho
end

function shifted_quotient(ctx,bus,epsilon,kp,ki;sigma=0.05)
    m=N.descriptor(ctx,_rho(bus,epsilon),kp,ki)
    g=N.gauge_vector(m)
    norm(g)>0 || error("empty gauge direction")
    Q=nullspace(reshape(g/norm(g),1,:))
    Aq=transpose(Q)*m.Ared*Q
    As=Matrix{Float64}(Aq)+sigma*I
    (;As,Q,model=m,Aq,epsilon,bus)
end

function sigma_min_at(As,omega)
    M=ComplexF64.(im*omega*I(size(As,1))-As)
    minimum(svdvals(M))
end

"""Deterministic observed minimum; this is a lower-bound observation only."""
function observed_radius(As;points=260,refine_steps=55)
    λ=eigvals(As)
    wmax=max(1.0,4maximum(abs.(imag.(λ)),init=0.0),20maximum(abs.(λ),init=0.0))
    grid=unique(vcat(0.0,10.0.^range(-7,log10(wmax),length=160),
        collect(range(0.0,wmax,length=points))))
    values=[sigma_min_at(As,w) for w in grid]
    best=argmin(values);wbest=grid[best];vbest=values[best];refinements=0
    phi=(sqrt(5.0)-1)/2
    for j in 2:length(grid)-1
        values[j]<=values[j-1] && values[j]<=values[j+1] || continue
        a=grid[j-1];b=grid[j+1];c=b-phi*(b-a);d=a+phi*(b-a)
        fc=sigma_min_at(As,c);fd=sigma_min_at(As,d)
        for _ in 1:refine_steps
            if fc<fd
                b=d;d=c;fd=fc;c=b-phi*(b-a);fc=sigma_min_at(As,c)
            else
                a=c;c=d;fc=fd;d=a+phi*(b-a);fd=sigma_min_at(As,d)
            end
        end
        w=(a+b)/2;v=sigma_min_at(As,w);refinements+=1
        if v<vbest;wbest=w;vbest=v;end
    end
    (;status="OBSERVED_NOT_CERTIFIED",beta_upper_observed=vbest,
      gamma_lower_observed=1/vbest,omega_observed=wbest,
      pointwise_min_at_zero=values[1],frequency_samples=length(grid),refinements,
      frequencies=grid,sigma_min=values)
end

"""Rigorous Lipschitz interval lower bound for minω σmin(jωI−Aσ).

For every real ω, the smallest singular value is 1-Lipschitz in ω. The
unsearched tail is bounded with σmin(jωI−Aσ) ≥ |ω|−||Aσ||F. A reported
frequency-grid minimum is retained separately as an observed upper bound on
the exact radius; it is never labelled a norm certificate.
"""
function certify_radius(As,beta;max_nodes=50000,max_depth=64)
    n=size(As,1);normF=norm(As)
    W=normF+beta+1.0
    fcache=Dict{Float64,Float64}()
    function f(w)
        get!(fcache,Float64(w)) do
            sigma_min_at(As,w)
        end
    end
    f0=f(0.0);fW=f(W)
    f0<beta && return (;status="FAIL_POINT_WITNESS",certified=false,
        beta_lower_cert=NaN,beta_upper_observed=f0,gamma_lower_observed=1/f0,
        omega_witness=0.0,nodes=length(fcache),max_depth=0,tail_lower=W-normF)
    fW<beta && return (;status="FAIL_POINT_WITNESS",certified=false,
        beta_lower_cert=NaN,beta_upper_observed=fW,gamma_lower_observed=1/fW,
        omega_witness=W,nodes=length(fcache),max_depth=0,tail_lower=W-normF)
    stack=[(0.0,W,f0,fW,0)]
    lower_bounds=Float64[];observed=min(f0,fW);maxseen=0
    while !isempty(stack)
        a,b,fa,fb,depth=pop!(stack);maxseen=max(maxseen,depth)
        mid=(a+b)/2;fm=f(mid);observed=min(observed,fm)
        if fm<beta
            return (;status="FAIL_POINT_WITNESS",certified=false,
                beta_lower_cert=NaN,beta_upper_observed=observed,
                gamma_lower_observed=1/observed,omega_witness=mid,
                nodes=length(fcache),max_depth=maxseen,tail_lower=W-normF)
        end
        lb=fm-(b-a)/2
        if lb>=beta
            push!(lower_bounds,lb)
        else
            if depth>=max_depth || length(fcache)>=max_nodes
                return (;status="INCOMPLETE_INTERVAL_BOUND",certified=false,
                    beta_lower_cert=isempty(lower_bounds) ? -Inf : minimum(lower_bounds),
                    beta_upper_observed=observed,gamma_lower_observed=1/observed,
                    omega_witness=NaN,nodes=length(fcache),max_depth=maxseen,
                    tail_lower=W-normF)
            end
            push!(stack,(a,mid,fa,fm,depth+1))
            push!(stack,(mid,b,fm,fb,depth+1))
        end
    end
    tail=W-normF
    certlower=min(tail,minimum(lower_bounds))
    (;status=certlower>=beta ? "CERTIFIED" : "INCOMPLETE_INTERVAL_BOUND",
      certified=certlower>=beta,beta_lower_cert=certlower,
      beta_upper_observed=observed,gamma_lower_observed=1/observed,
      omega_witness=NaN,nodes=length(fcache),max_depth=maxseen,tail_lower=tail)
end

"""Bounded-real Riccati upper certificate for ||(sI−Aσ)⁻¹||∞.

The stabilizing Hamiltonian subspace yields P solving
Aσ'P+PAσ+I+P²/γ²=0. Positivity of P, stability of Aσ+P/γ²,
and a small independently checked CARE residual certify the requested γ upper
bound. This is separate from any frequency-grid observation.
"""
function hinf_upper_certificate(As,gamma;residual_tol=5e-7)
    n=size(As,1);gamma>0 || error("gamma must be positive")
    Iₙ=Matrix{Float64}(I,n,n)
    # A symplectic similarity balances the γ^-2 and unit off-diagonal blocks.
    H=[As Iₙ./gamma; -Iₙ./gamma -transpose(As)]
    F=schur(H);stablemask=real.(F.values).<0
    nstable=count(stablemask)
    nstable==n || return (;status="NO_STABILIZING_HAMILTONIAN_SUBSPACE",
        certified=false,gamma,stable_eigenvalue_count=nstable,
        care_relative_residual=Inf,min_P_eigenvalue=NaN,
        closed_loop_abscissa=NaN,hamiltonian_abscissa=maximum(real.(F.values)))
    Fs=ordschur(F,stablemask)
    U1=Fs.Z[1:n,1:n];U2=Fs.Z[n+1:2n,1:n]
    Praw=gamma.*(U2/U1)
    P=real.((Praw+adjoint(Praw))/2)
    # Refine the invariant-subspace result by Newton–Kleinman Sylvester steps.
    for _ in 1:8
        Acl=As+P/gamma^2
        Fp=transpose(As)*P+P*As+Iₙ+(P*P)./gamma^2
        Pnew=sylvester(transpose(Acl),Acl,(P*P)./gamma^2-Iₙ)
        Pnew=real.((Pnew+transpose(Pnew))/2)
        Fnew=transpose(As)*Pnew+Pnew*As+Iₙ+(Pnew*Pnew)./gamma^2
        norm(Fnew)<norm(Fp) || break
        P=Pnew
        norm(Fnew)<=1e-13*max(norm(P),1.0) && break
    end
    care_of(Px)=transpose(As)*Px+Px*As+Iₙ+(Px*Px)./gamma^2
    scales=sort!(unique(vcat(1.0,1 .+ 10.0.^(-8:-1),1 .- 10.0.^(-8:-1),
        collect(range(0.9,1.1,length=41)))))
    candidates=[(scale=s,care=care_of(s.*P)) for s in scales]
    chosen=candidates[argmin([maximum(eigvals(Symmetric((x.care+transpose(x.care))/2))) for x in candidates])]
    P=chosen.scale.*P;care=chosen.care
    rel=norm(care)/max(norm(transpose(As)*P)+norm(P*As)+norm(Iₙ)+norm(P*P)/gamma^2,eps())
    careeig=eigvals(Symmetric((care+transpose(care))/2))
    mineig=minimum(eigvals(Symmetric(P)))
    clalpha=maximum(real.(eigvals(As+P/gamma^2)))
    lmi_max=maximum(careeig)
    valid=mineig>0 && clalpha<0 && lmi_max < -residual_tol
    (;status=valid ? "BOUNDED_REAL_CERTIFIED" : "BOUNDED_REAL_NUMERICAL_FAILURE",
      certified=valid,gamma,gamma_upper=gamma,beta_lower=1/gamma,
      stable_eigenvalue_count=nstable,care_relative_residual=rel,
      P_scale=chosen.scale,
      min_P_eigenvalue=mineig,care_max_eigenvalue=lmi_max,
      care_min_eigenvalue=minimum(careeig),closed_loop_abscissa=clalpha,
      hamiltonian_abscissa=maximum(real.(F.values)))
end

function _beta0(ctx,bus,epsilon,kp,ki)
    sh=shifted_quotient(ctx,bus,epsilon,kp,ki)
    sigma_min_at(sh.As,0.0)
end

"""Safeguarded scalar active-frequency continuation for a real critical mode."""
function robust_boundary(ctx,bus,epsilon0,kp,ki,beta_req;reserve=1e-3,
                         max_bisections=70)
    target=beta_req*(1+reserve)
    b0=_beta0(ctx,bus,epsilon0,kp,ki)
    b0<target || error("starting point unexpectedly satisfies robust target")
    lo=epsilon0;hi=epsilon0+max(1e-8,epsilon0*1e-4)
    bh=_beta0(ctx,bus,hi,kp,ki);expansions=0
    while bh<target
        hi=epsilon0+2*(hi-epsilon0);expansions+=1
        hi<1 || error("robust root is outside epsilon ∈ (0,1]")
        bh=_beta0(ctx,bus,hi,kp,ki)
        expansions<40 || error("could not bracket robust root")
    end
    for _ in 1:max_bisections
        mid=(lo+hi)/2;bm=_beta0(ctx,bus,mid,kp,ki)
        bm<target ? (lo=mid) : (hi=mid)
        hi-lo<=max(1e-13,epsilon0*1e-10) && break
    end
    epsilon=(lo+hi)/2
    (;epsilon,target_beta=target,beta0=_beta0(ctx,bus,epsilon,kp,ki),
      bracket=(lo,hi),initial_beta=b0,high_beta=bh,expansions,
      bisections=max_bisections,reserve)
end

"""One-sided gain derivatives and local KKT test at a real-frequency active radius."""
function robust_kkt_audit(ctx,bus,epsilon,kp,ki,beta_target;step_normalized=1e-3)
    epsscale=epsilon;β=_beta0(ctx,bus,epsilon,kp,ki)
    de=max(1e-9,epsilon*1e-5)
    βeps=(_beta0(ctx,bus,epsilon+de,kp,ki)-
          _beta0(ctx,bus,epsilon-de,kp,ki))/(2de)*epsscale
    βg=zeros(20);locations=String[];stepabs=Float64[]
    for i in 1:10
        for (kind,vec,base,lo,hi) in (("Kp",kp,N.K0P,ctx.kpmin[i],ctx.kpmax[i]),
                                     ("Ki",ki,N.K0I,ctx.kimin[i],ctx.kimax[i]))
            y=vec[i]/base;h=step_normalized
            plus=y+h<=hi/base+1e-12
            minus=y-h>=lo/base-1e-12
            vp=copy(vec);vm=copy(vec)
            if plus && minus
                vp[i]+=h*base;vm[i]-=h*base
                der=(_beta0(ctx,bus,epsilon,kp===vec ? vp : kp,ki===vec ? vp : ki)-
                     _beta0(ctx,bus,epsilon,kp===vec ? vm : kp,ki===vec ? vm : ki))/(2h)
                push!(locations,"INTERIOR")
            elseif plus
                vp[i]+=h*base
                der=(_beta0(ctx,bus,epsilon,kp===vec ? vp : kp,ki===vec ? vp : ki)-β)/h
                push!(locations,"LOWER")
            elseif minus
                vm[i]-=h*base
                der=(β-_beta0(ctx,bus,epsilon,kp===vec ? vm : kp,ki===vec ? vm : ki))/h
                push!(locations,"UPPER")
            else
                der=NaN;push!(locations,"NO_FEASIBLE_FD")
            end
            βg[kind=="Kp" ? i : 10+i]=der
            push!(stepabs,h)
        end
    end
    power=ctx.power[bus-29]
    μ=power*epsscale/βeps
    grad_lag=-μ.*βg
    lower=zeros(20);upper=zeros(20);residual=zeros(20)
    for i in 1:20
        if locations[i]=="LOWER"
            lower[i]=grad_lag[i];residual[i]=max(0.0,-grad_lag[i])
        elseif locations[i]=="UPPER"
            upper[i]=-grad_lag[i];residual[i]=max(0.0,grad_lag[i])
        else
            residual[i]=abs(grad_lag[i])
        end
    end
    active_row=zeros(1,21);active_row[1]=-βeps;active_row[2:end].=-βg
    Jbounds=zeros(20,21);for j in 1:20;Jbounds[j,j+1]=1;end
    licq=rank(vcat(active_row,Jbounds))==21
    g=beta_target-β
    stationarity=max(abs(power*epsscale-μ*βeps),maximum(residual))
    strict=all((lower .> 1e-10).|(upper .>1e-10))
    (;beta=β,beta_target,g,robust_multiplier=μ,beta_gradient_epsilon=βeps,
      beta_gradient_gains=βg,gain_bound_side=locations,
      lower_multipliers=lower,upper_multipliers=upper,
      primal_residual=max(0.0,g),stationarity_residual=stationarity,
      complementarity_residual=abs(μ*g),LICQ=licq,active_jacobian_rank=rank(active_row),
      SOSC=strict ? "VACUOUS_CRITICAL_CONE_STRICT_GAIN_BOUNDS" : "NOT_CERTIFIED",
      status=μ>0 && licq && strict && stationarity<1e-5 ? "LOCAL_KKT_CERTIFIED" : "KKT_NOT_CERTIFIED",
      singular_gap=NaN,finite_difference_normalized_step=step_normalized)
end

end
