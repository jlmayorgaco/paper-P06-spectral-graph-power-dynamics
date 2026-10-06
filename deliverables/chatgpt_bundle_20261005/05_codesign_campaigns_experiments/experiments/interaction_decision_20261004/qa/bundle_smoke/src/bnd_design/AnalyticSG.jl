module AnalyticSG

using LinearAlgebra
using ForwardDiff
using CSV
using DataFrames

export SGParameters, rating_fractions, scale_current_on_system_base,
       frozen_parameters, frozen_bus_operating_point, state_names,
       uncontrolled_state_names, machine_algebraic, rhs, rhs_uncontrolled,
       output, jacobians, port_transfer,
       steady_state_diagnostics

"Physical rating shares for a co-located aggregate; internal equations are untouched."
function rating_fractions(rho::Real)
    0<=rho<=1 || throw(ArgumentError("rho must be in [0,1]"))
    return (sg=1-Float64(rho),gfl=Float64(rho))
end

"Convert a per-unit device current to the common system base by its rating share."
scale_current_on_system_base(i_device_pu,rating_mva,system_base_mva) =
    (rating_mva/system_base_mva).*i_device_pu

Base.@kwdef struct SGParameters{T<:Real}
    rs::T
    xd::T
    xq::T
    xdp::T
    xqp::T
    xdpp::T
    xqpp::T
    xls::T
    tdp::T
    tdpp::T
    tqp::T
    tqpp::T
    inertia::T
    damping::T
    rating_mva::T
    system_base_mva::T
    omega_base::T
    omega_frame::T
    avr_ka::T
    avr_ke::T
    avr_kf::T
    avr_ta::T
    avr_tf::T
    avr_te::T
    avr_tr::T
    avr_vr_min::T
    avr_vr_max::T
    avr_e1::T
    avr_se1::T
    avr_e2::T
    avr_se2::T
    avr_vref::T
    gov_vmin::T
    gov_vmax::T
    gov_r::T
    gov_t1::T
    gov_t2::T
    gov_t3::T
    gov_dt::T
    gov_omega_ref::T
    gov_p_ref::T
    controlled::Bool
    vf_set::T
    tau_m_set::T
end

const STATE_NAMES = ["gov_xg1","gov_xg2","avr_v_fb","avr_vfout","avr_vr","avr_vm",
    "machine_psi2q","machine_psi2d","machine_Epd","machine_Epq","machine_omega","machine_delta"]
state_names()=copy(STATE_NAMES)
uncontrolled_state_names()=["machine_psi2q","machine_psi2d","machine_Epd","machine_Epq","machine_omega","machine_delta"]

function csv_row(path, key, value)
    table=CSV.read(path,DataFrame)
    idx=findfirst(==(value),table[!,Symbol(key)])
    idx===nothing && error("missing $key=$value in $path")
    return table[idx,:]
end

function optional_csv_row(path,key,value)
    table=CSV.read(path,DataFrame)
    idx=findfirst(==(value),table[!,Symbol(key)])
    return idx===nothing ? nothing : table[idx,:]
end

"Load bus-specific, copied IEEE-39 parameter rows. This does not import PowerDynamics."
function frozen_parameters(root,bus::Integer;system_base_mva=100.0,omega_base=2pi*60,omega_frame=1.0)
    data=joinpath(root,"reports","experiment_D","inputs")
    m=csv_row(joinpath(data,"machine.csv"),"bus",bus)
    a=optional_csv_row(joinpath(data,"avr.csv"),"bus",bus)
    g=optional_csv_row(joinpath(data,"gov.csv"),"bus",bus)
    controlled=(a!==nothing&&g!==nothing)
    isfile(joinpath(root,"reports","experiment_A","matrices","bus33_baseline_equilibrium.csv")) ||
        error("missing frozen ExpA all-SG baseline state file")
    T=Float64
    SGParameters{T}(
        rs=Float64(m[Symbol("R_s")]), xd=Float64(m[Symbol("X_d")]), xq=Float64(m[Symbol("X_q")]),
        xdp=Float64(m[Symbol("X′_d")]), xqp=Float64(m[Symbol("X′_q")]),
        xdpp=Float64(m[Symbol("X″_d")]), xqpp=Float64(m[Symbol("X″_q")]), xls=Float64(m[Symbol("X_ls")]),
        tdp=Float64(m[Symbol("T′_d0")]), tdpp=Float64(m[Symbol("T″_d0")]),
        tqp=Float64(m[Symbol("T′_q0")]), tqpp=Float64(m[Symbol("T″_q0")]),
        inertia=Float64(m[:H]), damping=Float64(m[:D]), rating_mva=Float64(m[:Sn]),
        system_base_mva=Float64(system_base_mva),omega_base=Float64(omega_base),omega_frame=Float64(omega_frame),
        avr_ka=controlled ? Float64(a[:Ka]) : 1.0,avr_ke=controlled ? Float64(a[:Ke]) : 0.0,
        avr_kf=controlled ? Float64(a[:Kf]) : 0.0,avr_ta=controlled ? Float64(a[:Ta]) : 1.0,
        avr_tf=controlled ? Float64(a[:Tf]) : 1.0,avr_te=controlled ? Float64(a[:Te]) : 1.0,
        avr_tr=controlled ? Float64(a[:Tr]) : 1.0,
        avr_vr_min=controlled ? Float64(a[:vr_min]) : -Inf,avr_vr_max=controlled ? Float64(a[:vr_max]) : Inf,
        avr_e1=controlled ? Float64(a[:E1]) : 1.0,avr_se1=controlled ? Float64(a[:Se1]) : 0.0,
        avr_e2=controlled ? Float64(a[:E2]) : 2.0,avr_se2=controlled ? Float64(a[:Se2]) : 0.0,
        avr_vref=0.0,gov_vmin=controlled ? Float64(g[:V_min]) : -Inf,
        gov_vmax=controlled ? Float64(g[:V_max]) : Inf,gov_r=controlled ? Float64(g[:R]) : 1.0,
        gov_t1=controlled ? Float64(g[:T1]) : 1.0,gov_t2=controlled ? Float64(g[:T2]) : 0.0,
        gov_t3=controlled ? Float64(g[:T3]) : 1.0,gov_dt=controlled ? Float64(g[:DT]) : 0.0,
        gov_omega_ref=controlled ? Float64(g[:ω_ref]) : 1.0,gov_p_ref=0.0,
        controlled=controlled,vf_set=0.0,tau_m_set=0.0)
end

function frozen_bus_operating_point(root,bus::Integer,p::SGParameters=frozen_parameters(root,bus))
    # ExpA stores the full-network all-SG equilibrium in this bus-33-named file.
    eqfile=joinpath(root,"reports","experiment_A","matrices","bus33_baseline_equilibrium.csv")
    eqrows=CSV.read(eqfile,DataFrame)
    vals=Dict(Int(r.state_index)=>Float64(r.equilibrium_value) for r in eachrow(eqrows))
    suffixes=p.controlled ? ["ctrld_gen₊gov₊xg1","ctrld_gen₊gov₊xg2","ctrld_gen₊avr₊v_fb",
        "ctrld_gen₊avr₊vfout","ctrld_gen₊avr₊vr","ctrld_gen₊avr₊vm",
        "ctrld_gen₊machine₊ψ″_q","ctrld_gen₊machine₊ψ″_d","ctrld_gen₊machine₊E′_d",
        "ctrld_gen₊machine₊E′_q","ctrld_gen₊machine₊ω","ctrld_gen₊machine₊δ"] :
        ["machine₊ψ″_q","machine₊ψ″_d","machine₊E′_d","machine₊E′_q","machine₊ω","machine₊δ"]
    function find_value(suffix;kind="differential")
        idx=findfirst(r->Int(r.bus)==bus && String(r.differential_or_algebraic)==kind &&
            endswith(String(r.state_name),suffix*")"),eachrow(eqrows))
        idx===nothing && error("ExpA state map is missing bus $bus state $suffix")
        return vals[Int(eqrows[idx,:state_index])]
    end
    x=Float64[find_value(s) for s in suffixes]
    ur=find_value("busbar₊u_r";kind="algebraic"); ui=find_value("busbar₊u_i";kind="algebraic")
    u=[ur,ui]
    # References are recovered from the frozen equilibrium equations rather than
    # treated as free fitting parameters.
    if p.controlled
        vref=x[6]+x[3]+x[5]/p.avr_ka
        # TGOV1 uses ref_sig=(p_ref-DeltaOmega)/R; DeltaOmega=omega-omega_ref.
        pref=p.gov_r*x[1]+(x[11]-p.gov_omega_ref)
        p2=SGParameters{Float64}(; (name=>(name==:avr_vref ? vref : name==:gov_p_ref ? pref : getfield(p,name))
            for name in fieldnames(typeof(p)))...)
    else
        m=machine_algebraic(x,u,p)
        gd1=(p.xdpp-p.xls)/(p.xdp-p.xls); gd2=(p.xdp-p.xdpp)/(p.xdp-p.xls)^2
        id=m.id; psi2d=x[2]; ed=x[3]; eq=x[4]
        vf=eq+(p.xd-p.xdp)*(id-gd2*psi2d-(1-gd1)*id+gd2*eq)
        taue=m.psid*m.iq-m.psiq*m.id
        taum=taue+p.damping*(x[5]-1)
        p2=SGParameters{Float64}(; (name=>(name==:vf_set ? vf : name==:tau_m_set ? taum : getfield(p,name))
            for name in fieldnames(typeof(p)))...)
    end
    return (x=x,u=[ur,ui],parameters=p2,bus=bus,
        frozen_state_source=eqfile,frozen_state_map=eqfile,
        parameter_sources=[joinpath(root,"reports","experiment_D","inputs",f) for f in ("machine.csv","avr.csv","gov.csv")])
end

function machine_algebraic(xm,u,p::SGParameters)
    psi2q,psi2d,ed,eq,omega,delta=xm
    xd1=(p.xdpp-p.xls)/(p.xdp-p.xls); xq1=(p.xqpp-p.xls)/(p.xqp-p.xls)
    cd=xd1*eq+(1-xd1)*psi2d
    cq=-xq1*ed+(1-xq1)*psi2q
    sn,cs=sin(delta),cos(delta)
    # PowerDynamics defines terminal voltage = T_to_glob(delta)*[Vd,Vq];
    # therefore the input port uses the inverse transform T_to_loc.
    vd=sn*u[1]-cs*u[2]
    vq=cs*u[1]+sn*u[2]
    mat=[p.rs -omega*p.xqpp; omega*p.xdpp p.rs]
    id,iq=mat\[-vd-omega*cq,-vq+omega*cd]
    psid=-p.xdpp*id+cd
    psiq=-p.xqpp*iq+cq
    iglobal=(p.rating_mva/p.system_base_mva).*[sn*id+cs*iq,-cs*id+sn*iq]
    return (id=id,iq=iq,vd=vd,vq=vq,psid=psid,psiq=psiq,
        i_terminal=iglobal,p_machine=vd*id+vq*iq,q_machine=vq*id-vd*iq)
end

function quad_se(u,se1,se2,e1,e2)
    a=sqrt(se1*e1/(se2*e2)); A=e2-(e1-e2)/(a-1)
    u<=A && return zero(u)
    B=se2*e2*(a-1)^2/(e1-e2)^2
    return B*(u-A)^2/u
end

"Independent Sauer–Pai + Type-I AVR + TGOV1 differential equations."
function rhs(x,u,p::SGParameters)
    length(x)==12 || throw(DimensionMismatch("controlled Sauer-Pai bus has twelve differential states"))
    xg1,xg2,vfb,vfout,vr,vm,psi2q,psi2d,ed,eq,omega,delta=x
    ma=machine_algebraic(view(x,7:12),u,p)
    id,iq,vd,vq,psid,psiq=ma.id,ma.iq,ma.vd,ma.vq,ma.psid,ma.psiq
    gd1=(p.xdpp-p.xls)/(p.xdp-p.xls); gd2=(p.xdp-p.xdpp)/(p.xdp-p.xls)^2
    gq1=(p.xqpp-p.xls)/(p.xqp-p.xls); gq2=(p.xqp-p.xqpp)/(p.xqp-p.xls)^2
    vmag=sqrt(vd^2+vq^2)
    vfceil=abs(vfout)*quad_se(abs(vfout),p.avr_se1,p.avr_se2,p.avr_e1,p.avr_e2)
    vfoutdot=(vr-vfceil-p.avr_ke*vfout)/p.avr_te
    amp=p.avr_ka*(p.avr_vref-vm-vfb)
    vrdot=(amp-vr)/p.avr_ta
    vrdot=ifelse(((vr>p.avr_vr_max)&&(amp>vr))||((vr<p.avr_vr_min)&&(amp<vr)),zero(vrdot),vrdot)
    vfbdot=(p.avr_kf*vfoutdot-vfb)/p.avr_tf
    vmdot=(vmag-vm)/p.avr_tr
    ref=(p.gov_p_ref-(omega-p.gov_omega_ref))/p.gov_r
    xg1dot=(ref-xg1)/p.gov_t1
    xg1dot=ifelse(((xg1>p.gov_vmax)&&(ref>xg1))||((xg1<p.gov_vmin)&&(ref<xg1)),zero(xg1dot),xg1dot)
    xg2dot=(xg1+p.gov_t2*xg1dot-xg2)/p.gov_t3
    taum=(xg2-p.gov_dt*(omega-p.gov_omega_ref))/omega
    taue=psid*iq-psiq*id
    omegadot=(taum-taue-p.damping*(omega-1))/(2p.inertia)
    deltadot=p.omega_base*(omega-p.omega_frame)
    eqdot=(-eq-(p.xd-p.xdp)*(id-gd2*psi2d-(1-gd1)*id+gd2*eq)+vfout)/p.tdp
    eddot=(-ed+(p.xq-p.xqp)*(iq-gq2*psi2q-(1-gq1)*iq-gq2*ed))/p.tqp
    psi2ddot=(-psi2d+eq-(p.xdp-p.xls)*id)/p.tdpp
    psi2qdot=(-psi2q-ed-(p.xqp-p.xls)*iq)/p.tqpp
    return [xg1dot,xg2dot,vfbdot,vfoutdot,vrdot,vmdot,psi2qdot,psi2ddot,eddot,eqdot,omegadot,deltadot]
end

"Independent six-state Sauer–Pai model with fixed field voltage and torque inputs."
function rhs_uncontrolled(x,u,p::SGParameters)
    length(x)==6 || throw(DimensionMismatch("uncontrolled Sauer-Pai machine has six differential states"))
    psi2q,psi2d,ed,eq,omega,delta=x
    ma=machine_algebraic(x,u,p)
    id,iq=ma.id,ma.iq
    gd1=(p.xdpp-p.xls)/(p.xdp-p.xls); gd2=(p.xdp-p.xdpp)/(p.xdp-p.xls)^2
    gq1=(p.xqpp-p.xls)/(p.xqp-p.xls); gq2=(p.xqp-p.xqpp)/(p.xqp-p.xls)^2
    eqdot=(-eq-(p.xd-p.xdp)*(id-gd2*psi2d-(1-gd1)*id+gd2*eq)+p.vf_set)/p.tdp
    eddot=(-ed+(p.xq-p.xqp)*(iq-gq2*psi2q-(1-gq1)*iq-gq2*ed))/p.tqp
    psi2ddot=(-psi2d+eq-(p.xdp-p.xls)*id)/p.tdpp
    psi2qdot=(-psi2q-ed-(p.xqp-p.xls)*iq)/p.tqpp
    taue=ma.psid*iq-ma.psiq*id
    omegadot=(p.tau_m_set-taue-p.damping*(omega-1))/(2p.inertia)
    deltadot=p.omega_base*(omega-p.omega_frame)
    return [psi2qdot,psi2ddot,eddot,eqdot,omegadot,deltadot]
end

output(x,u,p::SGParameters)=machine_algebraic(length(x)==12 ? view(x,7:12) : x,u,p).i_terminal
function jacobians(x,u,p::SGParameters)
    f=z->(length(z)==12 ? rhs(z,u,p) : rhs_uncontrolled(z,u,p))
    fu=v->(length(x)==12 ? rhs(x,v,p) : rhs_uncontrolled(x,v,p))
    fz=z->f(z); hz=z->output(z,u,p); hu=v->output(x,v,p)
    return (A=ForwardDiff.jacobian(fz,x),B=ForwardDiff.jacobian(fu,u),
        C=ForwardDiff.jacobian(hz,x),D=ForwardDiff.jacobian(hu,u))
end
port_transfer(A,B,C,D,s)=C*((s*I-A)\B)+D

function steady_state_diagnostics(x,u,p::SGParameters)
    m=machine_algebraic(length(x)==12 ? view(x,7:12) : x,u,p)
    f=length(x)==12 ? rhs(x,u,p) : rhs_uncontrolled(x,u,p)
    return (rhs=f,rhs_norm=norm(f,Inf),port_current=output(x,u,p),
        machine_power=(p=m.p_machine,q=m.q_machine),
        network_power=(p=real((u[1]+im*u[2])*conj(output(x,u,p)[1]+im*output(x,u,p)[2])),
            q=imag((u[1]+im*u[2])*conj(output(x,u,p)[1]+im*output(x,u,p)[2]))),
        stator_residual=norm([p.rs*m.id+(length(x)==12 ? x[11] : x[5])*m.psiq+m.vd,
            p.rs*m.iq-(length(x)==12 ? x[11] : x[5])*m.psid+m.vq],Inf))
end

end
