module CollectiveModel

using CSV, DataFrames, LinearAlgebra
include(joinpath(@__DIR__, "..", "bnd_design", "AnalyticSG.jl"))
include(joinpath(@__DIR__, "..", "bnd_design", "AnalyticGFLPLL.jl"))

export FrozenNetwork, frozen_network, mixed_jacobian, spectrum, port_admittance

struct FrozenNetwork
    root::String
    buses::Vector{Int}
    voltage::Vector{ComplexF64}
    y_static::Matrix{Float64}
    sg::Dict{Int,Any}
    gfl::Dict{Int,Any}
    load_audit::DataFrame
    no_load_kcl_inf::Float64
end

"Convert a complex nodal admittance to real and imaginary coordinates."
function realify(Y::AbstractMatrix{<:Complex})
    n = size(Y, 1)
    R = zeros(Float64, 2n, 2n)
    for i in 1:n, j in 1:n
        z = Y[i,j]
        R[2i-1,2j-1] = real(z)
        R[2i-1,2j] = -imag(z)
        R[2i,2j-1] = imag(z)
        R[2i,2j] = real(z)
    end
    R
end

"""Reconstruct the fixed branch/ZIP network from frozen data and states.

All archived loads are constant impedance. Their initialized effective
admittances are recovered from nodal current balance because the archived CSV
setpoints need not equal the initialized injections at generator-plus-load buses.
"""
function frozen_network(root::AbstractString)
    buses_df = CSV.read(joinpath(root,"reports","experiment_D","inputs","bus.csv"),DataFrame)
    branch_df = CSV.read(joinpath(root,"reports","experiment_D","inputs","branch.csv"),DataFrame)
    load_df = CSV.read(joinpath(root,"reports","experiment_D","inputs","load.csv"),DataFrame)
    eq = CSV.read(joinpath(root,"reports","experiment_A","matrices","bus33_baseline_equilibrium.csv"),DataFrame)
    buses = sort(Int.(buses_df.bus)); n=length(buses)
    buses == collect(1:n) || error("expected contiguous archived bus numbering")
    voltage = Vector{ComplexF64}(undef,n)
    for b in buses
        rows = filter(r -> Int(r.bus)==b && r.differential_or_algebraic=="algebraic" &&
            occursin("busbar",String(r.state_name)), eq)
        length(rows.state_index)==2 || error("missing frozen busbar voltage at bus $b")
        ur=only(Float64.(filter(r->occursin("u_r",String(r.state_name)),rows).equilibrium_value))
        ui=only(Float64.(filter(r->occursin("u_i",String(r.state_name)),rows).equilibrium_value))
        voltage[b]=complex(ur,ui)
    end
    Yline=zeros(ComplexF64,n,n)
    for r in eachrow(branch_df)
        a,b=Int(r.src_bus),Int(r.dst_bus)
        t=Float64(r.r_src); z=complex(Float64(r.R),Float64(r.X))
        y=inv(z); ysrc=complex(Float64(r.G_src),Float64(r.B_src))
        ydst=complex(Float64(r.G_dst),Float64(r.B_dst))
        Yline[a,a]-=t^2*(y+ysrc)
        Yline[a,b]+=t*y
        Yline[b,a]+=t*y
        Yline[b,b]-=y+ydst
    end
    sg=Dict{Int,Any}(); gfl=Dict{Int,Any}()
    Igen=zeros(ComplexF64,n)
    for r in eachrow(buses_df)
        Bool(r.has_gen) || continue
        b=Int(r.bus)
        op=AnalyticSG.frozen_bus_operating_point(root,b)
        J=AnalyticSG.jacobians(op.x,op.u,op.parameters)
        d=AnalyticSG.steady_state_diagnostics(op.x,op.u,op.parameters)
        d.rhs_norm<1e-8 || error("SG state residual at bus $b")
        sg[b]=(op=op,J=J,diagnostic=d)
        Igen[b]=complex(d.port_current...)
        gop=AnalyticGFLPLL.operating_point_for_injection(op.u,d.network_power.p,d.network_power.q)
        gd=AnalyticGFLPLL.steady_state_diagnostics(gop.x,gop.u,gop.parameters)
        norm(gd.rhs,Inf)<1e-8 || error("GFL state residual at bus $b")
        gfl[b]=(op=gop,)
    end
    Iline=Yline*voltage
    Yload=zeros(ComplexF64,n)
    audit=NamedTuple[]
    for r in eachrow(load_df)
        b=Int(r.bus)
        (r.KpZ==1 && r.KqZ==1 && r.KpI==0 && r.KqI==0 && r.KpC==0 && r.KqC==0) ||
            error("non-impedance ZIP load at bus $b needs explicit linearization")
        iload=-Iline[b]-Igen[b]
        Yload[b]=iload/voltage[b]
        s=voltage[b]*conj(iload)
        push!(audit,(bus=b,initialized_load_MW=-100real(s),initialized_load_Mvar=-100imag(s),
            CSV_setpoint_load_MW=-100Float64(r.Pset),CSV_setpoint_load_Mvar=-100Float64(r.Qset),
            initialized_admittance_real=real(Yload[b]),initialized_admittance_imag=imag(Yload[b])))
    end
    no_load=maximum(abs(Iline[b]+Igen[b]) for b in buses if !any(load_df.bus .== b))
    Ystatic=Yline+Diagonal(Yload)
    FrozenNetwork(String(root),buses,voltage,realify(Ystatic),sg,gfl,DataFrame(audit),no_load)
end

"""Assemble exact physical endpoint removal and independent-gain DAE Schur reduction."""
function mixed_jacobian(net::FrozenNetwork,rho::AbstractVector,Kp::AbstractVector,Ki::AbstractVector)
    gens=sort(collect(keys(net.sg))); ng=length(gens)
    length(rho)==length(Kp)==length(Ki)==ng || throw(DimensionMismatch("one independent triple per SG"))
    blocks=NamedTuple[]
    for (j,b) in enumerate(gens)
        0<=rho[j]<=1 || throw(ArgumentError("rho outside [0,1] at bus $b"))
        if rho[j]<1
            J=net.sg[b].J
            push!(blocks,(bus=b,kind="SG",share=1-rho[j],J=J))
        end
        if rho[j]>0
            op=net.gfl[b].op
            J=AnalyticGFLPLL.jacobians(op.x,op.u,op.parameters;kp=Kp[j],ki=Ki[j])
            push!(blocks,(bus=b,kind="GFL",share=rho[j],J=J))
        end
    end
    nx=sum(size(t.J.A,1) for t in blocks); ny=size(net.y_static,1)
    A=zeros(nx,nx); B=zeros(nx,ny); C=zeros(ny,nx); D=zeros(ny,ny)
    state_map=NamedTuple[]; off=0
    for t in blocks
        J=t.J; d=size(J.A,1); xi=(off+1):(off+d); yi=(2t.bus-1):(2t.bus)
        A[xi,xi].=J.A; B[xi,yi].=J.B
        C[yi,xi].+=t.share.*J.C
        D[yi,yi].+=t.share.*J.D
        push!(state_map,(bus=t.bus,kind=t.kind,first=first(xi),last=last(xi)))
        off+=d
    end
    Gy=net.y_static+D
    Ared=A-B*(Gy\C)
    return (Ared=Ared,A=A,B=B,C=C,D=D,Gy=Gy,state_map=DataFrame(state_map),
        n_dynamic=nx,n_algebraic=ny,condition_Gy=cond(Gy))
end

function spectrum(net::FrozenNetwork,rho,Kp,Ki)
    model=mixed_jacobian(net,rho,Kp,Ki)
    lambda=eigvals(model.Ared)
    return (;model,lambda,spectral_abscissa=maximum(real.(lambda)),
        least_damped=lambda[argmax(real.(lambda))])
end

port_admittance(J,s)=J.D+J.C*((s*I(size(J.A,1))-J.A)\J.B)

end
