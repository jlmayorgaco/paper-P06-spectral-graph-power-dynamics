module PDReferenceN

using CSV, DataFrames, Graphs, LinearAlgebra, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__, "..", "pd39", "PD39.jl"))
using .PD39

export FrozenBaseline, frozen_baseline, build_architecture, trim_state,
    direct_power_audit, residual_audit

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BASE_MVA = 100.0
const NOMINAL_KP = 2pi * 5
const NOMINAL_KI = (2pi * 5)^2 / 4

struct FrozenBaseline
    nw
    state
    rows::Dict{Int,NamedTuple}
end

function frozen_baseline()
    nw = PD39.baseline_network()
    state = PD39.initialize_equilibrium(nw; sparse=false, check=:error).state
    table = CSV.read(joinpath(ROOT,"reports","experiment_N",
        "TABLE_N01_original_operating_point.csv"), DataFrame)
    rows = Dict{Int,NamedTuple}()
    for row in eachrow(table)
        bus=Int(row.bus)
        r=(;V_real=Float64(row.V_real_pu),V_imag=Float64(row.V_imag_pu),
            P=Float64(row.P_gen_MW),Q=Float64(row.Q_gen_Mvar),
            Sn=Float64(row.Sn_original_MVA),H=Float64(row.H_seconds))
        prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
        measured_P=Float64(state[VIndex(bus,Symbol(prefix*"Sn"))]) *
            Float64(state[VIndex(bus,Symbol(prefix*"P"))])
        abs(measured_P-r.P)<1e-8 || error("N01 P mismatch at bus $bus")
        rows[bus]=r
    end
    FrozenBaseline(nw,state,rows)
end

function _gfl(;rho,kp,ki)
    xf=0.03;vdc=2.5;cdc=1.25;ftau=300.0;fi=600.0
    PD39.PD39Model.WeightedSimpleGFLDC(
        Xf=xf,Rf=0.01,PLL_Kp=kp,PLL_Ki=ki,
        PLL_τ_lpf=1/(ftau*2pi),
        CC1_KP=(xf/(2pi*60))*(fi*2pi),
        CC1_KI=(xf/(2pi*60))*(fi*2pi)^2/4,
        CC1_F=0,CC1_Fcoupl=0,C_dc=cdc,V_dc=vdc,
        kp_v_dc=vdc*cdc*(5*2pi),
        ki_v_dc=vdc*cdc*(5*2pi)^2/4,
        port_scale=rho,name=:gfl)
end

"""Build the compiled PD reference architecture, retaining ZIP loads separately."""
function build_architecture(base::FrozenBaseline,rho::AbstractVector,
                            kp::AbstractVector=fill(NOMINAL_KP,10),
                            ki::AbstractVector=fill(NOMINAL_KI,10))
    length(rho)==length(kp)==length(ki)==10 || throw(DimensionMismatch("ten generator triples required"))
    all(0 .<= rho .<= 1) || throw(ArgumentError("rho outside [0,1]"))
    mdl=PD39.PD39Model
    data=mdl.ieee39_data()
    vertices,edges=mdl.copy_network_components(base.nw)
    for bus in 30:39
        i=bus-29; r=Float64(rho[i]); r==0 && continue
        busrow=data.bus[findfirst(==(bus),data.bus.bus),:]
        parts=Any[]
        if r<1
            sg=deepcopy(bus==39 ? mdl.OfficialIEEE39.uncontrolled_machine :
                mdl.OfficialIEEE39.controlled_machine)
            push!(parts,sg)
        end
        mdl.set_Sbase!(100.0);mdl.set_fbase!(60.0)
        push!(parts,_gfl(;rho=r,kp=Float64(kp[i]),ki=Float64(ki[i])))
        Bool(busrow.has_load) && push!(parts,deepcopy(mdl.OfficialIEEE39.load))
        pf=mdl.get_pfmodel(vertices[bus])
        v=mdl.compile_bus(mdl.MTKBus(parts...);
            name=Symbol("bus$(bus)"),vidx=bus,pf)
        mdl.set_Sbase!();mdl.set_fbase!()
        mdl.set_default!(v,Symbol("busbar₊Vbase"),Float64(busrow.base_kv))
        if Bool(busrow.has_load)
            mdl.OfficialIEEE39.apply_csv_params!(v,data.load,bus)
            vset=Float64(base.state[VIndex(bus,:ZIPLoad₊Vset)])
            mdl.set_default!(v,r"ZIPLoad₊Vset$",vset)
        end
        if r<1
            mdl.OfficialIEEE39.apply_csv_params!(v,data.machine,bus)
            Bool(busrow.has_avr) && mdl.OfficialIEEE39.apply_csv_params!(v,data.avr,bus)
            Bool(busrow.has_gov) && mdl.OfficialIEEE39.apply_csv_params!(v,data.gov,bus)
            mdl.set_default!(v,r"Sn$",(1-r)*base.rows[bus].Sn)
        end
        vertices[bus]=v
    end
    nw=mdl.Network(vertices,edges)
    mdl.set_jac_prototype!(nw)
    return nw
end

"""Closed-form trim from installed SG, PLL, CC1, filter and DC equations."""
function trim_state(nw,base::FrozenBaseline,rho::AbstractVector,
                    kp::AbstractVector=fill(NOMINAL_KP,10),
                    ki::AbstractVector=fill(NOMINAL_KI,10))
    length(rho)==length(kp)==length(ki)==10 || throw(DimensionMismatch("ten generator triples required"))
    s=NWState(nw)
    mixedvars=Dict(string(x)=>x for x in NetworkDynamics.SII.variable_symbols(nw))
    mixedpars=Dict(string(x)=>x for x in NetworkDynamics.SII.parameter_symbols(nw))
    for (sym,val) in zip(NetworkDynamics.SII.variable_symbols(base.nw),uflat(base.state))
        key=string(sym)
        haskey(mixedvars,key) && (s[mixedvars[key]]=val)
    end
    for (sym,val) in zip(NetworkDynamics.SII.parameter_symbols(base.nw),pflat(base.state))
        key=string(sym)
        haskey(mixedpars,key) && (s[mixedpars[key]]=val)
    end
    for bus in 30:39
        i=bus-29;r=Float64(rho[i]);r==0 && continue
        op=base.rows[bus]
        ur=op.V_real;ui=op.V_imag;v2=ur^2+ui^2;vm=sqrt(v2)
        p=op.P/BASE_MVA;q=op.Q/BASE_MVA
        ir=(p*ur+q*ui)/v2;ii=(p*ui-q*ur)/v2
        id=p/vm;iq=-q/vm;theta=atan(ui,ur)
        rf=0.01;xf=0.03
        vid=vm+rf*id-xf*iq;viq=rf*iq+xf*id
        cki=(xf/(2pi*60))*(600*2pi)^2/4
        for (name,value) in (
            ("gfl₊filter₊i_f_r",ir),("gfl₊filter₊i_f_i",ii),
            ("gfl₊pll₊θ",theta),("gfl₊pll₊Δω_rad_s",0.0),
            ("gfl₊pll₊Δω_i_rad_s",0.0),
            ("gfl₊cc1₊γ_d",vid/cki),("gfl₊cc1₊γ_q",viq/cki),
            ("gfl₊v_dc_state",2.5),("gfl₊v_dc_i",id))
            s[VIndex(bus,Symbol(name))]=value
        end
        s[VIndex(bus,:gfl₊iset_q)]=iq
        s[VIndex(bus,:gfl₊P_dc)]=vid*id+viq*iq
        s[VIndex(bus,:gfl₊filter₊PortScale)]=r
        s[VIndex(bus,:gfl₊PLL_Kp)]=Float64(kp[i])
        s[VIndex(bus,:gfl₊PLL_Ki)]=Float64(ki[i])
        if r<1
            prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
            s[VIndex(bus,Symbol(prefix*"Sn"))]=(1-r)*op.Sn
        end
    end
    _refine_sg_roundoff!(nw,s,rho)
    s
end

"""Remove amplified floating-point residuals in surviving SG states.

At very small retained ratings, the network current is a sum of SG and GFL
contributions. The resulting Float64 cancellation can leave O(1e-10) SG
residuals even though the frozen P/Q sharing is accurate to O(1e-12) pu.
This solves the three local machine/AVR equilibrium rows using their own
physical states. It is invoked only beyond the stated 1e-10 residual gate;
voltage, power setpoints, and design gains are held fixed.
"""
function _refine_sg_roundoff!(nw,s,rho)
    u=uflat(s);p=pflat(s)
    du=similar(u);nw(du,u,p,0.0)
    maximum(abs,du)<1e-10 && return s
    syms=string.(NetworkDynamics.SII.variable_symbols(nw))
    for bus in 30:38
        0<rho[bus-29]<1 || continue
        names=("ctrld_gen₊avr₊vm","ctrld_gen₊machine₊ψ″_q",
               "ctrld_gen₊machine₊ψ″_d")
        idx=[findfirst(==(string(VIndex(bus,Symbol(n)))),syms) for n in names]
        any(isnothing,idx) && continue
        j=Int.(idx)
        maximum(abs,du[j])<1e-10 && continue
        J=zeros(3,3);h=1e-7
        for (col,k) in enumerate(j)
            up=copy(u);um=copy(u);up[k]+=h;um[k]-=h
            fp=similar(du);fm=similar(du)
            nw(fp,up,p,0.0);nw(fm,um,p,0.0)
            J[:,col]=(fp[j]-fm[j])/(2h)
        end
        δ= -J\du[j]
        maximum(abs,δ)<1e-8 || error("SG roundoff correction too large at bus $bus")
        u[j].+=δ
        for k in j
            s[NetworkDynamics.SII.variable_symbols(nw)[k]]=u[k]
        end
        nw(du,u,p,0.0)
    end
    maximum(abs,du)<1e-10 || error("trim residual remains above gate")
    s
end

function residual_audit(nw,s)
    du=similar(uflat(s));nw(du,uflat(s),pflat(s),0.0)
    (;maximum=maximum(abs,du),norm=norm(du),finite=all(isfinite,du))
end

function direct_power_audit(s,base::FrozenBaseline,rho)
    rows=NamedTuple[]
    for bus in 30:39
        i=bus-29;r=Float64(rho[i]);op=base.rows[bus]
        ur=Float64(s[VIndex(bus,:busbar₊u_r)])
        ui=Float64(s[VIndex(bus,:busbar₊u_i)])
        sgp=0.0;sgq=0.0;gfp=0.0;gfq=0.0
        if r<1
            prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
            sn=Float64(s[VIndex(bus,Symbol(prefix*"Sn"))])
            sgp=sn*Float64(s[VIndex(bus,Symbol(prefix*"P"))])
            sgq=sn*Float64(s[VIndex(bus,Symbol(prefix*"Q"))])
        end
        if r>0
            ir=Float64(s[VIndex(bus,:gfl₊filter₊i_f_r)])
            ii=Float64(s[VIndex(bus,:gfl₊filter₊i_f_i)])
            gfp=BASE_MVA*r*(ur*ir+ui*ii)
            gfq=BASE_MVA*r*(ui*ir-ur*ii)
        end
        pl=bus in (31,39) ? BASE_MVA*Float64(s[VIndex(bus,:ZIPLoad₊P)]) : 0.0
        ql=bus in (31,39) ? BASE_MVA*Float64(s[VIndex(bus,:ZIPLoad₊Q)]) : 0.0
        pn=Float64(s[VIndex(bus,:busbar₊P_MW)])
        qn=Float64(s[VIndex(bus,:busbar₊Q_MVAr)])
        tau=if r<1
            prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
            Float64(s[VIndex(bus,Symbol(prefix*"τ_m"))])
        else
            NaN
        end
        bound_checks=Bool[]
        if r<1
            prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
            for tail in ("τ_m","τ_e","vf")
                push!(bound_checks,Float64(s[VIndex(bus,Symbol(prefix*tail))])>=-1e-12)
            end
            if bus!=39
                vfout=Float64(s[VIndex(bus,:ctrld_gen₊avr₊vfout)])
                vr=Float64(s[VIndex(bus,:ctrld_gen₊avr₊vr)])
                vrmin=Float64(s[VIndex(bus,:ctrld_gen₊avr₊vr_min)])
                vrmax=Float64(s[VIndex(bus,:ctrld_gen₊avr₊vr_max)])
                xg1=Float64(s[VIndex(bus,:ctrld_gen₊gov₊xg1)])
                vmin=Float64(s[VIndex(bus,:ctrld_gen₊gov₊V_min)])
                vmax=Float64(s[VIndex(bus,:ctrld_gen₊gov₊V_max)])
                append!(bound_checks,(vfout>=-1e-12,vrmin-1e-12<=vr<=vrmax+1e-12,
                    vmin-1e-12<=xg1<=vmax+1e-12))
            end
        end
        r>0 && push!(bound_checks,Float64(s[VIndex(bus,:gfl₊v_dc_state)])>0)
        push!(rows,(;bus,rho=r,SG_P_MW=sgp,SG_Q_Mvar=sgq,
            GFL_P_MW=gfp,GFL_Q_Mvar=gfq,ZIP_P_MW=pl,ZIP_Q_Mvar=ql,
            net_P_MW=pn,net_Q_Mvar=qn,
            target_SG_P_MW=(1-r)*op.P,target_SG_Q_Mvar=(1-r)*op.Q,
            target_GFL_P_MW=r*op.P,target_GFL_Q_Mvar=r*op.Q,
            max_P_error_pu=max(abs(sgp-(1-r)*op.P),abs(gfp-r*op.P))/BASE_MVA,
            max_Q_error_pu=max(abs(sgq-(1-r)*op.Q),abs(gfq-r*op.Q))/BASE_MVA,
            KCL_P_error_MW=sgp+gfp+pl-pn,
            KCL_Q_error_Mvar=sgq+gfq+ql-qn,
            SG_tau_m_pu=tau,
            bounds_status=all(bound_checks) ? "PASS" : "TRIM_INFEASIBLE"))
    end
    DataFrame(rows)
end

end
