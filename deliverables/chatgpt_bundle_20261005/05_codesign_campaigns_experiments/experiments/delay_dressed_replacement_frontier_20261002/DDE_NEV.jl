module DDENEV
using LinearAlgebra
include("DelayCharacteristic.jl")
const D=DelayCharacteristic

"Bordered Newton correction for a simple nonlinear eigenvalue root."
function newton_root(L,s0,tau,v0;tol=2e-10,maxiter=30,max_step=2.0)
    s=ComplexF64(s0);v=ComplexF64.(v0);v./=norm(v);anchor=copy(v)
    n=length(v);lastres=Inf
    for k in 1:maxiter
        M=D.delta_matrix(L,s,tau);Ms=D.delta_s(L,s,tau)
        r=M*v;lastres=norm(r)/max(1,norm(M)*norm(v))
        lastres<=tol && return (;s,v,residual=lastres,iterations=k-1,converged=true)
        K=Matrix{ComplexF64}(undef,n+1,n+1)
        K[1:n,1:n].=M;K[1:n,n+1].=Ms*v;K[n+1,1:n].=conj.(anchor);K[n+1,n+1]=0
        step=try K\vcat(-r,0.0+0im) catch; return (;s,v,residual=lastres,iterations=k,converged=false) end
        ds=step[end];fac=abs(ds)>max_step ? max_step/abs(ds) : 1.0
        s+=fac*ds;v.+=fac.*step[1:n]
        scale=dot(anchor,v)
        abs(scale)>1e-10 || return (;s,v,residual=lastres,iterations=k,converged=false)
        v./=scale
    end
    M=D.delta_matrix(L,s,tau);lastres=norm(M*v)/max(1,norm(M)*norm(v))
    (;s,v,residual=lastres,iterations=maxiter,converged=lastres<=tol)
end

"Continue all zero-delay ODE roots in a declared real-part strip."
function track_roots(L,tau;real_cut=-1.0,continuation_steps=8,max_roots=typemax(Int),maxiter=30)
    e=eigen(ComplexF64.(L.A));eligible=findall(real.(e.values).>=real_cut)
    ordered=eligible[sortperm(real.(e.values[eligible]),rev=true)]
    ids=ordered[1:min(length(ordered),max_roots)]
    out=NamedTuple[];path=collect(range(0.0,1.0,length=continuation_steps+1))[2:end]
    for idx in ids
        s=e.values[idx];v=e.vectors[:,idx];ok=true;res=Inf;iters=0
        for q in path
            r=newton_root(L,s,q.*tau,v;maxiter)
            if !r.converged;ok=false;res=r.residual;iters+=r.iterations;break;end
            s=r.s;v=r.v;res=r.residual;iters+=r.iterations
        end
        push!(out,(;root_id="ode_$idx",initial=e.values[idx],s,v,converged=ok,residual=res,iterations=iters))
    end
    out
end

function _contour(gamma,radius,step)
    corners=ComplexF64[gamma-im*radius,radius-im*radius,radius+im*radius,gamma+im*radius,gamma-im*radius]
    z=ComplexF64[]
    for k in 1:4
        a=corners[k];b=corners[k+1];n=max(2,ceil(Int,abs(b-a)/step))
        append!(z,[a+(b-a)*(j/n) for j in 0:n-1])
    end
    z
end

"Numerical argument-principle root count in Re(s)>gamma, with a norm-derived finite contour."
function count_rightmost_roots(L,tau;gamma=-1.0,initial_step=20.0,max_step_halvings=5)
    bound=opnorm(L.A0,2)+sum(opnorm(L.Ai[j],2)*exp(-gamma*tau[j]) for j in eachindex(tau))
    radius=bound+1.0
    F=schur(ComplexF64.(L.A0));T=F.T;Z=F.Z;mu=diag(T)
    left=L.C'*Z;right=Z'*L.B;step=initial_step;maxjump=Inf;winding=NaN;points=0;min_logabs=Inf
    for refinement in 0:max_step_halvings
        contour=_contour(gamma,radius,step);phi=Vector{Float64}(undef,length(contour));min_logabs=Inf
        for (j,s) in enumerate(contour)
            Q=UpperTriangular(s*I-T);X=Q\right
            small=I-Diagonal(exp.(-s.*tau))*left*X
            la,sgn=logabsdet(small);min_logabs=min(min_logabs,la)
            iszero(sgn) && return (;count=missing,gamma,radius,npoints=length(contour),max_phase_increment=Inf,
                min_logabs_small_det=la,refinements=refinement,status="CONTOUR_SINGULAR")
            phi[j]=sum(angle.(s .- mu))+angle(sgn)
        end
        inc=[angle(exp(im*(phi[mod1(j+1,length(phi))]-phi[j]))) for j in eachindex(phi)]
        maxjump=maximum(abs,inc);winding=round(Int,sum(inc)/(2pi));points=length(contour)
        maxjump<=pi/3 && break
        step/=2
    end
    status=maxjump<=pi/3 ? "NUMERICAL_CONTOUR_COUNT" : "CONTOUR_RESOLUTION_UNRESOLVED"
    (;count=winding,gamma,radius,npoints=points,max_phase_increment=maxjump,min_logabs_small_det=min_logabs,
      refinements=round(Int,log2(initial_step/step)),status)
end

end
