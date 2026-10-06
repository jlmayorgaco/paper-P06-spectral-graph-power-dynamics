module AnalyticGFLPLL

using LinearAlgebra
using ForwardDiff
using CSV
using DataFrames

export GFLParameters, frozen_bus33_parameters, frozen_bus33_operating_point,
       state_names, output_names, rhs, output, jacobians, gain_derivatives,
       affine_update, pll_error, steady_state_diagnostics, rebase_to_rating,
       rated_output, rated_port_transfer, operating_point_for_injection

Base.@kwdef struct GFLParameters{T<:Real}
    Rf::T
    Xf::T
    omega_base::T
    omega_frame::T
    pll_tau::T
    pll_kp::T
    pll_ki::T
    cc_kp::T
    cc_ki::T
    cc_f::T
    cc_fcoupl::T
    dc_capacitance::T
    dc_voltage_ref::T
    dc_kp::T
    dc_ki::T
    iset_q::T
    p_dc::T
end

const STATE_NAMES = ["gamma_q","gamma_d","theta","delta_omega_rad_s",
    "delta_omega_i_rad_s","i_f_i","i_f_r","v_dc_i","v_dc_state"]
const OUTPUT_NAMES = ["terminal_i_r","terminal_i_i"]
state_names() = copy(STATE_NAMES)
output_names() = copy(OUTPUT_NAMES)

"PLL phase detector from the installed SimpleGFLDC equation, independent of gains."
pll_error(theta, u_r, u_i) = -sin(theta)*u_r + cos(theta)*u_i

function rotated_to_dq(xr, xi, theta)
    return cos(theta)*xr + sin(theta)*xi, -sin(theta)*xr + cos(theta)*xi
end
function rotated_to_ri(xd, xq, theta)
    return cos(theta)*xd - sin(theta)*xq, sin(theta)*xd + cos(theta)*xq
end

"""
    rhs(x, u, p; kp=p.pll_kp, ki=p.pll_ki)

Independent nine-state SimpleGFLDC differential model. State and equation
ordering follow the installed component: CC1 integrators, PLL states, L-filter
currents, then DC-link PI/capacitor states. The terminal voltage input and
terminal current output are in the common real/imaginary per-unit convention.
"""
function rhs(x::AbstractVector, u::AbstractVector, p::GFLParameters;
             kp=p.pll_kp, ki=p.pll_ki)
    length(x)==9 || throw(DimensionMismatch("SimpleGFLDC has nine differential states"))
    length(u)==2 || throw(DimensionMismatch("terminal voltage input must be [u_r,u_i]"))
    gamma_q,gamma_d,theta,domega,domega_i,i_f_i,i_f_r,vdc_i,vdc = x
    u_r,u_i = u

    i_d,i_q = rotated_to_dq(i_f_r,i_f_i,theta)
    v_d,v_q = rotated_to_dq(u_r,u_i,theta)
    iref_d = (p.dc_voltage_ref-vdc)*p.dc_kp+vdc_i
    iref_q = p.iset_q
    err_d,err_q = iref_d-i_d,iref_q-i_q
    # These are the installed CC1 equations, including its declared feedforward
    # and cross-coupling settings (both are zero in the frozen nominal model).
    v_i_d = -p.cc_fcoupl*p.Xf*i_q + p.cc_kp*err_d + p.cc_ki*gamma_d + p.cc_f*v_d
    v_i_q =  p.cc_fcoupl*p.Xf*i_d + p.cc_kp*err_q + p.cc_ki*gamma_q + p.cc_f*v_q
    v_i_r,v_i_i = rotated_to_ri(v_i_d,v_i_q,theta)

    di_f_r = (p.omega_base/p.Xf) * (v_i_r-u_r-p.Rf*i_f_r+p.omega_frame*p.Xf*i_f_i)
    di_f_i = (p.omega_base/p.Xf) * (v_i_i-u_i-p.Rf*i_f_i-p.omega_frame*p.Xf*i_f_r)
    p_ac = v_i_d*i_d + v_i_q*i_q
    e = pll_error(theta,u_r,u_i)

    return [iref_q-i_q,
            iref_d-i_d,
            domega,
            (domega_i+kp*e-domega)/p.pll_tau,
            ki*e,
            di_f_i,
            di_f_r,
            (p.dc_voltage_ref-vdc)*p.dc_ki,
            (p_ac-p.p_dc)/(p.dc_capacitance*vdc)]
end

output(x::AbstractVector, u::AbstractVector, p::GFLParameters) = [x[7],x[6]]

"""Convert a system-base GFL realization to its own MVA base exactly.

The component has no `Sn` parameter. This conversion changes the current,
filter, current-loop, and DC-link quantities so that a device rating `Sn`
shares the same voltage base while its terminal current is later referred back
to the system base. PLL states and PLL gains are base-independent.
"""
function rebase_to_rating(op,rating_mva,system_base_mva)
    α=Float64(rating_mva)/Float64(system_base_mva)
    α>0 || throw(ArgumentError("device rating must be positive"))
    p=op.parameters
    pown=GFLParameters(Rf=p.Rf*α,Xf=p.Xf*α,omega_base=p.omega_base,
        omega_frame=p.omega_frame,pll_tau=p.pll_tau,pll_kp=p.pll_kp,pll_ki=p.pll_ki,
        cc_kp=p.cc_kp*α,cc_ki=p.cc_ki*α,cc_f=p.cc_f,cc_fcoupl=p.cc_fcoupl,
        dc_capacitance=p.dc_capacitance/α,dc_voltage_ref=p.dc_voltage_ref,
        dc_kp=p.dc_kp/α,dc_ki=p.dc_ki/α,iset_q=p.iset_q/α,p_dc=p.p_dc/α)
    xown=Float64.(op.x)
    xown[[1,2,6,7,8]]./=α
    return (x=xown,u=copy(op.u),parameters=pown,rating_ratio=α,
        source_operating_point=op)
end

"Terminal current converted from the device MVA base to the common system base."
rated_output(x,u,p,rating_ratio)=rating_ratio.*output(x,u,p)

"Port transfer in system-base current units for a rated device realization."
rated_port_transfer(J,s,rating_ratio)=rating_ratio.*(J.D+J.C*((s*I(size(J.A,1))-J.A)\J.B))

function jacobians(x,u,p::GFLParameters; kp=p.pll_kp,ki=p.pll_ki)
    fz=z->rhs(z,u,p;kp,ki)
    fu=v->rhs(x,v,p;kp,ki)
    hz=z->output(z,u,p)
    hu=v->output(x,v,p)
    return (A=ForwardDiff.jacobian(fz,x), B=ForwardDiff.jacobian(fu,u),
            C=ForwardDiff.jacobian(hz,x), D=ForwardDiff.jacobian(hu,u))
end

"Exact parameter derivatives: only the PLL frequency and integrator rows move."
function gain_derivatives(x,u,p::GFLParameters)
    theta=x[3]; u_r,u_i=u
    e_x=zeros(promote_type(eltype(x),eltype(u)),9)
    e_x[3]=-cos(theta)*u_r-sin(theta)*u_i
    e_u=[-sin(theta),cos(theta)]
    A_kp=zeros(eltype(e_x),9,9); A_ki=similar(A_kp); fill!(A_ki,0)
    B_kp=zeros(eltype(e_x),9,2); B_ki=similar(B_kp); fill!(B_ki,0)
    A_kp[4,:].=e_x./p.pll_tau
    A_ki[5,:].=e_x
    B_kp[4,:].=e_u./p.pll_tau
    B_ki[5,:].=e_u
    C_kp=zeros(eltype(e_x),2,9); C_ki=similar(C_kp); fill!(C_ki,0)
    D_kp=zeros(eltype(e_x),2,2); D_ki=similar(D_kp); fill!(D_ki,0)
    return (A_kp=A_kp,A_ki=A_ki,B_kp=B_kp,B_ki=B_ki,
            C_kp=C_kp,C_ki=C_ki,D_kp=D_kp,D_ki=D_ki,
            error_gradient_x=e_x,error_gradient_u=e_u)
end

"Construct only the documented stock SimpleGFLDC nominal constants."
function frozen_bus33_parameters(; iset_q, p_dc)
    fbase=60.0; fpll=5.0; ftau=300.0; fi=600.0
    xf=0.03; vdc=2.5; cdc=1.25
    return GFLParameters(Rf=0.01,Xf=xf,omega_base=2pi*fbase,omega_frame=1.0,
        pll_tau=1/(ftau*2pi),pll_kp=fpll*2pi,pll_ki=(fpll*2pi)^2/4,
        cc_kp=(xf/(2pi*fbase))*(fi*2pi),
        cc_ki=(xf/(2pi*fbase))*(fi*2pi)^2/4,
        cc_f=0.0,cc_fcoupl=0.0,dc_capacitance=cdc,dc_voltage_ref=vdc,
        dc_kp=vdc*cdc*(5*2pi),dc_ki=vdc*cdc*(5*2pi)^2/4,
        iset_q=Float64(iset_q),p_dc=Float64(p_dc))
end

"""Construct the exact frozen SimpleGFLDC equilibrium for a specified AC injection.

`P` and `Q` are positive network injections in system-base pu. This derives the
PLL angle, dq current references, controller integrator states, and DC power
from the device equations without an initialization solver.
"""
function operating_point_for_injection(u::AbstractVector,P::Real,Q::Real)
    length(u)==2 || throw(DimensionMismatch("voltage must be [u_r,u_i]"))
    V=complex(u[1],u[2]); S=complex(P,Q)
    Iinj=conj(S/V)
    theta=angle(V)
    id,iq=rotated_to_dq(real(Iinj),imag(Iinj),theta)
    p0=frozen_bus33_parameters(iset_q=iq,p_dc=0.0)
    vd,vq=rotated_to_dq(u[1],u[2],theta)
    # Steady L-filter equation in the PLL frame, with the stock SimpleGFLDC
    # feed-forward and cross-coupling settings (both zero in A-C).
    vid=vd+p0.Rf*id-p0.omega_frame*p0.Xf*iq
    viq=vq+p0.Rf*iq+p0.omega_frame*p0.Xf*id
    gamma_d=vid/p0.cc_ki
    gamma_q=viq/p0.cc_ki
    x=zeros(Float64,9)
    x[1]=gamma_q; x[2]=gamma_d; x[3]=theta
    x[6]=imag(Iinj); x[7]=real(Iinj); x[8]=id; x[9]=p0.dc_voltage_ref
    pac=vid*id+viq*iq
    p=frozen_bus33_parameters(iset_q=iq,p_dc=pac)
    return (x=x,u=Float64.(u),parameters=p,injection=(p=Float64(P),q=Float64(Q)),
        current=Float64[real(Iinj),imag(Iinj)])
end

function steady_state_diagnostics(x,u,p::GFLParameters)
    e=pll_error(x[3],u[1],u[2])
    return (pll_error=e, rhs=rhs(x,u,p), output=output(x,u,p),
            relative_voltage_angle_error=abs(x[3]-atan(u[2],u[1])))
end

"""
    frozen_bus33_operating_point(root)

Read the explicitly frozen ExpA bus-33 operating point. This is a static
operating-point input; the design equations and Jacobians below are assembled
from the component equations above and do not call PowerDynamics.
"""
function frozen_bus33_operating_point(root)
    statefile=joinpath(root,"reports","experiment_A","tables","TABLE_A03_state_partition.csv")
    eqfile=joinpath(root,"reports","experiment_A","matrices","bus33_mixed_equilibrium.csv")
    states=CSV.read(statefile,DataFrame)
    eqrows=CSV.read(eqfile,DataFrame)
    vals=Dict(Int(r.state_index)=>Float64(r.equilibrium_value) for r in eachrow(eqrows))
    names=["gfl₊cc1₊γ_q","gfl₊cc1₊γ_d","gfl₊pll₊θ","gfl₊pll₊Δω_rad_s",
           "gfl₊pll₊Δω_i_rad_s","gfl₊filter₊i_f_i","gfl₊filter₊i_f_r",
           "gfl₊v_dc_i","gfl₊v_dc_state"]
    function value_for(suffix;bus="33",kind="differential")
        # State names are rendered as `VIndex(bus, :component₊variable)`.
        idx=findfirst(r->string(r.bus)==bus && string(r.differential_or_algebraic)==kind &&
            endswith(String(r.state_name),suffix*")"),eachrow(states))
        idx===nothing && error("ExpA state map is missing $suffix at bus $bus")
        row=states[idx,:]
        return vals[Int(row.state_index)]
    end
    x=[value_for(n) for n in names]
    ur=value_for("busbar₊u_r";kind="algebraic")
    ui=value_for("busbar₊u_i";kind="algebraic")
    # At the frozen point the d-axis is aligned with the bus voltage and the q
    # current reference equals the measured q-axis filter current. P_dc is
    # solved from the same steady-state DC-link balance, not by a gain sweep.
    x[3]=atan(ui,ur); x[4]=0.0; x[5]=0.0
    theta=x[3]
    i_d,i_q=rotated_to_dq(x[7],x[6],theta)
    provisional=frozen_bus33_parameters(iset_q=i_q,p_dc=0.0)
    iref_d=(provisional.dc_voltage_ref-x[9])*provisional.dc_kp+x[8]
    err_d,err_q=iref_d-i_d,provisional.iset_q-i_q
    v_i_d=provisional.cc_kp*err_d+provisional.cc_ki*x[2]
    v_i_q=provisional.cc_kp*err_q+provisional.cc_ki*x[1]
    p_dc=v_i_d*i_d+v_i_q*i_q
    p=frozen_bus33_parameters(iset_q=i_q,p_dc=p_dc)
    return (x=Float64.(x),u=[ur,ui],parameters=p,state_names=names,
            frozen_source=eqfile,full_state_source=statefile)
end

end
