module SingleSupportKKT

using CSV, DataFrames, LinearAlgebra
include(joinpath(@__DIR__,"..","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN

export solve_single_support_branches

const SIGMA=0.05
const STRICT_GUARD=1e-9

function spectrum_at(ctx,bus,eps,kp,ki)
    rho=ones(10);rho[bus-29]=1-eps
    PDExactDesignN.spectrum(ctx,rho,kp,ki)
end

"""Bracket and correct every sampled infeasible→feasible interval.

The architecture is fixed while eps>0. Each correction uses the complete
physical spectrum and deterministic bisection; it never crosses eps=0.
"""
function solve_single_support_branches(ctx,bus,kp,ki,grid)
    samples=[(;eps,alpha=spectrum_at(ctx,bus,eps,kp,ki).alpha) for eps in grid]
    roots=NamedTuple[]
    target=-SIGMA-STRICT_GUARD
    for j in 2:length(samples)
        a=samples[j-1];b=samples[j]
        a.eps==0 && b.eps==0 && continue
        a.alpha>target && b.alpha<=target || continue
        lo=a.eps;hi=b.eps
        for iteration in 1:60
            mid=(lo+hi)/2
            α=spectrum_at(ctx,bus,mid,kp,ki).alpha
            if α<=target
                hi=mid
            else
                lo=mid
            end
            hi-lo<1e-14 && break
        end
        eps=hi
        rho=ones(10);rho[bus-29]=1-eps
        sp=PDExactDesignN.spectrum(ctx,rho,kp,ki)
        d=PDExactDesignN.simple_mode_sensitivities(ctx,rho,kp,ki)
        g_eps=-real(d.rho[bus-29])
        mu= -ctx.power[bus-29]/g_eps
        gp=real.(d.Kp);gi=real.(d.Ki)
        stationarity=abs(ctx.power[bus-29]+mu*g_eps)
        # At the lower Kp box, grad L must be nonnegative. At the upper
        # Ki box it must be nonpositive. Interior gain coordinates need zero.
        gain_res=Float64[];strict=Bool[]
        for i in 1:10
            for (val,der,lower,upper) in ((kp[i],mu*gp[i],ctx.kpmin[i],ctx.kpmax[i]),
                                          (ki[i],mu*gi[i],ctx.kimin[i],ctx.kimax[i]))
                if abs(val-lower)<1e-9
                    push!(gain_res,max(0.0,-der));push!(strict,der>1e-10)
                elseif abs(val-upper)<1e-9
                    push!(gain_res,max(0.0,der));push!(strict,der< -1e-10)
                else
                    push!(gain_res,abs(der));push!(strict,false)
                end
            end
        end
        fullres=max(stationarity,maximum(gain_res))
        active=findall(real.(sp.lambda).>=sp.alpha-1e-6)
        kkt=mu>=0 && fullres<1e-7 && length(active)==1
        # Strict box multipliers and the one nonzero spectral gradient leave
        # no nonzero critical direction. This is the empty-cone SOSC case.
        sosc=kkt && all(strict)
        push!(roots,(;bus,epsilon=eps,rho=1-eps,
            retained_SG_MW=ctx.power[bus-29]*eps,
            GFL_MW=sum(ctx.power)-ctx.power[bus-29]*eps,
            alpha=sp.alpha,active_real=real(d.lambda),
            active_imag=imag(d.lambda),active_count=length(active),
            d_alpha_d_epsilon=g_eps,mu,stationarity_residual=fullres,
            complementarity_residual=abs(mu*(sp.alpha+SIGMA)),
            primal_residual=max(0.0,sp.alpha+SIGMA),
            LICQ=abs(g_eps)>1e-8,SOSC=sosc,
            KKT=kkt && sosc,pole_condition=d.condition,
            lower_bracket=a.eps,upper_bracket=b.eps))
    end
    samples,roots
end

end
