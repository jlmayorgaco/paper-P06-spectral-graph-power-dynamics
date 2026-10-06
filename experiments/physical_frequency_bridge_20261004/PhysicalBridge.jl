# Generated only in this new experiment; frozen sources remain unchanged.
module PhysicalBridge
using Graphs, LinearAlgebra, NetworkDynamics, PowerDynamics, PowerDynamics.Library
using ModelingToolkitBase
using ModelingToolkitBase: @component, @variables, @parameters, t_nounits, D_nounits, System
const P = Main.PDReferenceN
const ScaledPortLFilter = P.PD39.PD39Model.ScaledPortLFilter
const NOMINAL_KP = P.NOMINAL_KP
const NOMINAL_KI = P.NOMINAL_KI
@component function PhysicalGFLDC(; name, port_scale=1.0, defaults...)
    pars = @parameters begin
        Rf=0.01
        Xf=0.03
        PLL_Kp=2pi*5
        PLL_Ki=(2pi*5)^2/4
        PLL_τ_lpf=1/(2pi*300)
        CC1_KP=0.36
        CC1_KI=135.0
        CC1_F=0.0
        CC1_Fcoupl=0.0
        C_dc=1.25
        V_dc=2.5
        kp_v_dc
        ki_v_dc
        iset_q, [guess=0.0]
        P_dc, [guess=0.0]
    end
    vars = @variables begin
        v_dc_state(t_nounits), [guess=2.5]
        v_dc_i(t_nounits), [guess=0.0]
    end
    @named filter=ScaledPortLFilter(;Rf,Xf,PortScale=port_scale)
    @named pll=ComposableInverter.PLL_LPF(;Kp=PLL_Kp,Ki=PLL_Ki,τ_lpf=PLL_τ_lpf)
    @named cc1=ComposableInverter.CC1(;Xf,F=CC1_F,Fcoupl=CC1_Fcoupl,KP=CC1_KP,KI=CC1_KI)
    systems=@named begin
        terminal=Terminal()
    end
    push!(systems,filter); push!(systems,pll); push!(systems,cc1)
    iset_d_dc=(v_dc_state-V_dc)*kp_v_dc+v_dc_i
    p_ac=cc1.V_I_d*cc1.i_f_d+cc1.V_I_q*cc1.i_f_q
    eqs=Equation[
        connect(terminal,filter.terminal),
        pll.u_r ~ terminal.u_r,
        pll.u_i ~ terminal.u_i,
        cc1.i_f_ref_d ~ iset_d_dc,
        cc1.i_f_ref_q ~ iset_q,
        C_dc*D_nounits(v_dc_state) ~ (P_dc-p_ac)/v_dc_state,
        D_nounits(v_dc_i) ~ (v_dc_state-V_dc)*ki_v_dc,
    ]
    append!(eqs,[cc1.i_f_d,cc1.i_f_q] .~ ComposableInverter._ri_to_dq(filter.i_f_r,filter.i_f_i,pll.θ))
    append!(eqs,[cc1.V_C_d,cc1.V_C_q] .~ ComposableInverter._ri_to_dq(terminal.u_r,terminal.u_i,pll.θ))
    append!(eqs,[filter.V_I_r,filter.V_I_i] .~ ComposableInverter._dq_to_ri(cc1.V_I_d,cc1.V_I_q,pll.θ))
    sys=System(eqs,t_nounits,vars,pars;name,systems)
    return set_mtk_defaults(sys,defaults)
end

function _gfl(;rho,kp,ki)
    xf=0.03;vdc=2.5;cdc=1.25;ftau=300.0;fi=600.0
    PhysicalGFLDC(
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
function build_architecture(base::P.FrozenBaseline,rho::AbstractVector,
                            kp::AbstractVector=fill(NOMINAL_KP,10),
                            ki::AbstractVector=fill(NOMINAL_KI,10))
    length(rho)==length(kp)==length(ki)==10 || throw(DimensionMismatch("ten generator triples required"))
    all(0 .<= rho .<= 1) || throw(ArgumentError("rho outside [0,1]"))
    mdl=P.PD39.PD39Model
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


end
