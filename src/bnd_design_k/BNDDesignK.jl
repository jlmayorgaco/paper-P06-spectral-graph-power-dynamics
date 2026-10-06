module BNDDesignK

using LinearAlgebra, CSV, DataFrames, SHA, TOML, Statistics

# The design layer contains no PowerDynamics dependency.
include(joinpath(@__DIR__, "..", "bnd_design_g", "BNDDesignG.jl"))
using .BNDDesignG

const SIGMA_REQUIRED = 0.05
const ACTIVE_TOL = 1e-6
const SUPPORT_EPS = 1e-9

export DesignContext, design_context, support_mask, support_buses, support_epsilon,
    evaluate, pole_family, fixed_support_gradient, gain_box_support,
    minimum_gain_correction, derivative_check, exact_beta_star

struct DesignContext
    root::String
    net::Any
    buses::Vector{Int}
    power::Vector{Float64}
    kp0::Vector{Float64}
    ki0::Vector{Float64}
    kpmin::Vector{Float64}
    kpmax::Vector{Float64}
    kimin::Vector{Float64}
    kimax::Vector{Float64}
end

function design_context(root::AbstractString)
    net = BNDDesignG.CollectiveModel.frozen_network(root)
    discovered = BNDDesignG.PhysicalData.discover_generators(root)
    rows = discovered[discovered.replaceable .== true, :]
    buses = Int.(rows.bus)
    buses == sort(collect(keys(net.sg))) || error("discovered SGs disagree with network")
    power = Float64.(rows.SG_dispatch_initial_MW)
    nominal = BNDDesignG.PhysicalData.nominal_pll_gains()
    n = length(buses)
    DesignContext(String(root), net, buses, power,
        fill(nominal.Kp,n), fill(nominal.Ki,n),
        fill(0.25nominal.Kp,n), fill(4nominal.Kp,n),
        fill(0.25nominal.Ki,n), fill(4nominal.Ki,n))
end

support_mask(epsilon::AbstractVector) = sum((epsilon[i] > 0 ? 1 : 0) << (i-1) for i in eachindex(epsilon))
support_buses(ctx::DesignContext, mask::Integer) = [ctx.buses[i] for i in eachindex(ctx.buses) if (mask >> (i-1)) & 1 == 1]
function support_epsilon(ctx::DesignContext, mask::Integer, value::Real)
    [((mask >> (i-1)) & 1 == 1) ? Float64(value) : 0.0 for i in eachindex(ctx.buses)]
end

"""Exact physical finite spectrum after projection of the known angular gauge.

The state architecture is rebuilt at every endpoint. Derivatives must never be
evaluated across eps=0 (SG removal) or eps=1 (GFL removal).
"""
function evaluate(ctx::DesignContext, epsilon, kp=ctx.kp0, ki=ctx.ki0; vectors=false)
    n=length(ctx.buses)
    length(epsilon)==length(kp)==length(ki)==n || throw(DimensionMismatch("ten triples required"))
    all(0 .<= epsilon .<= 1) || throw(ArgumentError("epsilon outside [0,1]"))
    all(ctx.kpmin .<= kp .<= ctx.kpmax) || throw(ArgumentError("Kp outside frozen box"))
    all(ctx.kimin .<= ki .<= ctx.kimax) || throw(ArgumentError("Ki outside frozen box"))
    model=BNDDesignG.CollectiveModel.mixed_jacobian(ctx.net,1 .- epsilon,kp,ki)
    gauge=BNDDesignG.gauge_vector(ctx.net,model)
    norm(gauge)>0 || error("no angular gauge vector")
    residual=norm(model.Ared*gauge)/max(norm(model.Ared)*norm(gauge),eps(Float64))
    residual<=1e-8 || error("gauge residual $residual")
    Q=nullspace(reshape(gauge/norm(gauge),1,:))
    Aq=transpose(Q)*model.Ared*Q
    if vectors
        F=eigen(Aq)
        λ=F.values
        R=F.vectors
        L=inv(R)'
    else
        λ=eigvals(Aq)
        R=nothing; L=nothing
    end
    α=maximum(real.(λ))
    active=findall(z->real(z)>=α-ACTIVE_TOL,λ)
    return (;epsilon=Float64.(epsilon),kp=Float64.(kp),ki=Float64.(ki),
        model,gauge,Q,Aq,lambda=λ,right=R,left=L,alpha=α,active,
        support=support_mask(epsilon),gauge_residual=residual,
        retained_mw=dot(ctx.power,epsilon),gfl_mw=sum(ctx.power)-dot(ctx.power,epsilon))
end

function pole_family(ev, j)
    r=ev.Q*ev.right[:,j]
    energy=abs2.(r); total=sum(energy)
    sg=0.0; gfl=0.0
    for row in eachrow(ev.model.state_map)
        e=sum(energy[Int(row.first):Int(row.last)])
        if row.kind=="SG"; sg+=e else gfl+=e end
    end
    λ=ev.lambda[j]
    family=abs(λ)<1e-2 ? "COLLECTIVE_ZERO_BRANCH" :
        sg/total>0.75 ? "SG_INTERNAL_BRANCH" :
        gfl/total>0.75 ? "GFL_INTERNAL_BRANCH" :
        max(sg,gfl)/total<0.75 ? "NETWORK_INTERMODAL_BRANCH" : "OTHER"
    return (;family,sg_energy=sg/total,gfl_energy=gfl/total)
end

"""Ordinary reduced-matrix derivative on one fixed SG/GFL architecture."""
function pole_gradient(ctx::DesignContext, ev, j::Integer)
    all(0 .< ev.epsilon .< 1) || error("fixed-architecture derivative requires 0<eps<1 for each bus; use active coordinates")
    m=ev.model; GyinvC=m.Gy\m.C; n=length(ctx.buses)
    r=ev.right[:,j]; l=ev.left[:,j]; den=dot(l,r)
    ge=zeros(ComplexF64,n); gp=zeros(ComplexF64,n); gi=zeros(ComplexF64,n)
    for (i,b) in enumerate(ctx.buses)
        yi=(2b-1):(2b)
        sgrow=findfirst(row->row.bus==b && row.kind=="SG",eachrow(m.state_map))
        gfrow=findfirst(row->row.bus==b && row.kind=="GFL",eachrow(m.state_map))
        sgrow===nothing && error("SG block absent at active epsilon")
        gfrow===nothing && error("GFL block absent at active epsilon")
        sr=m.state_map[sgrow,:]; fr=m.state_map[gfrow,:]
        sx=Int(sr.first):Int(sr.last); fx=Int(fr.first):Int(fr.last)
        Js=ctx.net.sg[b].J
        op=ctx.net.gfl[b].op
        Jf=BNDDesignG.CollectiveModel.AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=ev.kp[i],ki=ev.ki[i])
        dC=zeros(size(m.C)); dD=zeros(size(m.D))
        dC[yi,sx].=Js.C; dC[yi,fx].=-Jf.C
        dD[yi,yi].=Js.D-Jf.D
        # d(A-B G^-1 C) = B G^-1 dG G^-1 C - B G^-1 dC.
        dAe=m.B*(m.Gy\(dD*GyinvC-dC))
        ge[i]=dot(l,transpose(ev.Q)*dAe*ev.Q*r)/den
        gd=BNDDesignG.CollectiveModel.AnalyticGFLPLL.gain_derivatives(op.x,op.u,op.parameters)
        for (dst,dA,dB) in ((gp,gd.A_kp,gd.B_kp),(gi,gd.A_ki,gd.B_ki))
            dAr=zeros(size(m.A)); dBr=zeros(size(m.B))
            dAr[fx,fx].=dA; dBr[fx,yi].=dB
            dAz=dAr-dBr*GyinvC
            dst[i]=dot(l,transpose(ev.Q)*dAz*ev.Q*r)/den
        end
    end
    (;epsilon=ge,kp=gp,ki=gi,denominator=den)
end

"""Derivative with respect to coordinates that remain interior in the support."""
function fixed_support_gradient(ctx::DesignContext, ev, j::Integer)
    m=ev.model; GyinvC=m.Gy\m.C; n=length(ctx.buses)
    r=ev.right[:,j]; l=ev.left[:,j]; den=dot(l,r)
    ge=fill(ComplexF64(NaN),n); gp=zeros(ComplexF64,n); gi=zeros(ComplexF64,n)
    for (i,b) in enumerate(ctx.buses)
        yi=(2b-1):(2b)
        sgrow=findfirst(row->row.bus==b && row.kind=="SG",eachrow(m.state_map))
        gfrow=findfirst(row->row.bus==b && row.kind=="GFL",eachrow(m.state_map))
        if sgrow!==nothing && gfrow!==nothing
            sr=m.state_map[sgrow,:]; fr=m.state_map[gfrow,:]
            sx=Int(sr.first):Int(sr.last); fx=Int(fr.first):Int(fr.last)
            Js=ctx.net.sg[b].J
            op=ctx.net.gfl[b].op
            Jf=BNDDesignG.CollectiveModel.AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=ev.kp[i],ki=ev.ki[i])
            dC=zeros(size(m.C)); dD=zeros(size(m.D))
            dC[yi,sx].=Js.C; dC[yi,fx].=-Jf.C
            dD[yi,yi].=Js.D-Jf.D
            dAe=m.B*(m.Gy\(dD*GyinvC-dC))
            ge[i]=dot(l,transpose(ev.Q)*dAe*ev.Q*r)/den
        end
        if gfrow!==nothing
            fr=m.state_map[gfrow,:]; fx=Int(fr.first):Int(fr.last)
            op=ctx.net.gfl[b].op
            gd=BNDDesignG.CollectiveModel.AnalyticGFLPLL.gain_derivatives(op.x,op.u,op.parameters)
            for (dst,dA,dB) in ((gp,gd.A_kp,gd.B_kp),(gi,gd.A_ki,gd.B_ki))
                dAr=zeros(size(m.A)); dBr=zeros(size(m.B))
                dAr[fx,fx].=dA; dBr[fx,yi].=dB
                dst[i]=dot(l,transpose(ev.Q)*(dAr-dBr*GyinvC)*ev.Q*r)/den
            end
        end
    end
    (;epsilon=ge,kp=gp,ki=gi,denominator=den)
end

"Exact support function for the frozen rectangular gain movement box."
function gain_box_support(v,x,lo,hi)
    length(v)==length(x)==length(lo)==length(hi) || throw(DimensionMismatch())
    delta=[v[i]>0 ? hi[i]-x[i] : v[i]<0 ? lo[i]-x[i] : 0.0 for i in eachindex(v)]
    (;value=dot(v,delta),delta)
end

"Minimum weighted gain movement for one active modal inequality, with box clipping."
function minimum_gain_correction(g,residual,x,lo,hi,R)
    # Solve min 1/2 d'R d subject to g'd <= residual; an explicit active set
    # clips any violated coordinates and repeats the one-row pseudoinverse.
    d=zeros(length(x)); free=trues(length(x))
    for _ in 1:length(x)+1
        remaining=residual-dot(g,d)
        remaining>=0 && return (;delta=d,feasible=true)
        ids=findall(free)
        isempty(ids) && break
        invR=1.0 ./ R[ids]
        denom=sum(abs2.(g[ids]).*invR)
        denom>0 || break
        trial=remaining .* (invR.*g[ids])./denom
        violations=[k for k in eachindex(ids) if x[ids[k]]+trial[k]<lo[ids[k]] || x[ids[k]]+trial[k]>hi[ids[k]]]
        if isempty(violations)
            d[ids].=trial
            return (;delta=d,feasible=true)
        end
        k=first(violations); i=ids[k]
        d[i]=clamp(x[i]+trial[k],lo[i],hi[i])-x[i]
        free[i]=false
    end
    (;delta=d,feasible=dot(g,d)<=residual+1e-10)
end

function derivative_check(ctx,ev,j,group,i;h=1e-5)
    g=fixed_support_gradient(ctx,ev,j)
    ep=copy(ev.epsilon); em=copy(ev.epsilon)
    kp=copy(ev.kp); km=copy(ev.kp); ip=copy(ev.ki); im=copy(ev.ki)
    if group==:epsilon
        ep[i]+=h; em[i]-=h; analytical=g.epsilon[i]
    elseif group==:kp
        kp[i]+=h; km[i]-=h; analytical=g.kp[i]
    elseif group==:ki
        ip[i]+=h; im[i]-=h; analytical=g.ki[i]
    else error("unknown group") end
    lp=evaluate(ctx,ep,kp,ip).lambda
    lm=evaluate(ctx,em,km,im).lambda
    target=ev.lambda[j]
    fd=(lp[argmin(abs.(lp.-target))]-lm[argmin(abs.(lm.-target))])/(2h)
    (;analytical,finite_difference=fd,absolute_error=abs(analytical-fd))
end

function exact_beta_star(ev;scale=1.0)
    BNDDesignG.resolvent_peak(ev.Aq,SIGMA_REQUIRED,scale)
end

end
