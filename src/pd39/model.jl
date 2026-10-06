module PD39Model

using CSV
using DataFrames
using Graphs
using ModelingToolkitBase
using NetworkDynamics
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase: @component, @variables, @parameters, t_nounits, D_nounits, System

export POWERDYNAMICS_IEEE39_EXAMPLE, CANDIDATE_SG_BUSES, SLACK_BUS,
    BASE_MVA, BASE_FREQ, GFL_MODEL_ID, ieee39_data, copy_network_components,
    baseline_network, simple_gfldc_template, replace_bus, replace_buses,
    weighted_replacement_network, candidate_table

const POWERDYNAMICS_IEEE39_EXAMPLE = joinpath(
    pkgdir(PowerDynamics),
    "docs",
    "examples",
    "ieee39_part1.jl",
)

"""
Private namespace containing the maintained PowerDynamics IEEE-39 example.

The example is included verbatim so that the baseline is the library's own
qualified reference model. New PD39 code composes copies of its compiled
components instead of modifying that source file.
"""
module OfficialIEEE39
    using PowerDynamics
    using PowerDynamics.Library
    using ModelingToolkitBase
    using NetworkDynamics
    using DataFrames
    using CSV

    include(joinpath(pkgdir(PowerDynamics), "docs", "examples", "ieee39_part1.jl"))
end

const CANDIDATE_SG_BUSES = [30, 32, 33, 34, 35, 36, 37, 38]
const SLACK_BUS = 31
const BASE_MVA = 100.0
const BASE_FREQ = 60.0
const GFL_MODEL_ID = "PowerDynamics.Library.ComposableInverter.SimpleGFLDC"

"LFilter whose external current port has a fixed share; state equations stay stock."
@component function ScaledPortLFilter(; name, port_scale=1.0, defaults...)
    pars=@parameters begin
        ωbase, [bound_to=:systembase₊ωbase]
        ωframe, [bound_to=:systembase₊ωframe]
        Rf=0.01
        Xf=0.03
        PortScale=1.0
    end
    vars=@variables begin
        i_f_r(t_nounits), [guess=0.0]
        i_f_i(t_nounits), [guess=0.0]
        i_f_mag(t_nounits)
        V_I_r(t_nounits)
        V_I_i(t_nounits)
    end
    @named terminal=Terminal()
    eqs=Equation[
        (Xf/ωbase)*D_nounits(i_f_r) ~ V_I_r-terminal.u_r-Rf*i_f_r+ωframe*Xf*i_f_i,
        (Xf/ωbase)*D_nounits(i_f_i) ~ V_I_i-terminal.u_i-Rf*i_f_i-ωframe*Xf*i_f_r,
        terminal.i_r ~ PortScale*i_f_r,
        terminal.i_i ~ PortScale*i_f_i,
        i_f_mag ~ sqrt(i_f_r^2+i_f_i^2),
    ]
    return set_mtk_defaults(System(eqs,t_nounits,vars,pars;name,systems=[terminal]),defaults)
end

"PowerDynamics SimpleGFLDC with its terminal current multiplied by a fixed share."
@component function WeightedSimpleGFLDC(; name, port_scale=1.0, defaults...)
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
    iset_d_dc=(V_dc-v_dc_state)*kp_v_dc+v_dc_i
    p_ac=cc1.V_I_d*cc1.i_f_d+cc1.V_I_q*cc1.i_f_q
    eqs=Equation[
        connect(terminal,filter.terminal),
        pll.u_r ~ terminal.u_r,
        pll.u_i ~ terminal.u_i,
        cc1.i_f_ref_d ~ iset_d_dc,
        cc1.i_f_ref_q ~ iset_q,
        C_dc*D_nounits(v_dc_state) ~ (p_ac-P_dc)/v_dc_state,
        D_nounits(v_dc_i) ~ (V_dc-v_dc_state)*ki_v_dc,
    ]
    append!(eqs,[cc1.i_f_d,cc1.i_f_q] .~ ComposableInverter._ri_to_dq(filter.i_f_r,filter.i_f_i,pll.θ))
    append!(eqs,[cc1.V_C_d,cc1.V_C_q] .~ ComposableInverter._ri_to_dq(terminal.u_r,terminal.u_i,pll.θ))
    append!(eqs,[filter.V_I_r,filter.V_I_i] .~ ComposableInverter._dq_to_ri(cc1.V_I_d,cc1.V_I_q,pll.θ))
    sys=System(eqs,t_nounits,vars,pars;name,systems)
    return set_mtk_defaults(sys,defaults)
end

"Return copies of the official CSV tables used by the IEEE-39 example."
function ieee39_data()
    return (
        bus = copy(OfficialIEEE39.bus_df),
        branch = copy(OfficialIEEE39.branch_df),
        load = copy(OfficialIEEE39.load_df),
        machine = copy(OfficialIEEE39.machine_df),
        avr = copy(OfficialIEEE39.avr_df),
        gov = copy(OfficialIEEE39.gov_df),
    )
end

"Copy the compiled network components without sharing mutable initialization metadata."
function copy_network_components(nw)
    vertices = [copy(nw[VIndex(i)]) for i in 1:Graphs.nv(nw)]
    edges = [copy(nw[EIndex(i)]) for i in 1:Graphs.ne(nw)]
    return vertices, edges
end

"Construct an isolated copy of the official PowerDynamics IEEE-39 network."
function baseline_network()
    vertices, edges = copy_network_components(OfficialIEEE39.nw)
    return Network(vertices, edges)
end

"Construct the preregistered stock GFL model used for replacement buses."
function simple_gfldc_template(; pll_scale = 1.0, filter_scale = 1.0, current_control_scale = 1.0)
    set_Sbase!(BASE_MVA)
    set_fbase!(BASE_FREQ)

    V_dc = 2.5
    C_dc = 1.25
    f_v_dc = 5.0
    Xf = 0.03 * filter_scale
    Rf = 0.01
    f_pll = 5.0 * pll_scale
    f_tau_pll = 300.0
    f_i_dq = 600.0 * current_control_scale

    @named gfl = ComposableInverter.SimpleGFLDC(
        Xf = Xf,
        Rf = Rf,
        PLL_Kp = f_pll * 2 * pi,
        PLL_Ki = (f_pll * 2 * pi)^2 / 4,
        PLL_τ_lpf = 1 / (f_tau_pll * 2 * pi),
        CC1_KP = (Xf / (2 * pi * BASE_FREQ)) * (f_i_dq * 2 * pi),
        CC1_KI = (Xf / (2 * pi * BASE_FREQ)) * (f_i_dq * 2 * pi)^2 / 4,
        CC1_F = 0,
        CC1_Fcoupl = 0,
        C_dc = C_dc,
        V_dc = V_dc,
        kp_v_dc = V_dc * C_dc * (f_v_dc * 2 * pi),
        ki_v_dc = V_dc * C_dc * (f_v_dc * 2 * pi)^2 / 4,
    )
    template = compile_bus(MTKBus(gfl); name = :SimpleGFLDC)

    # The component captured the bases at construction time; leave the global
    # setters in their default state for callers that build other networks.
    set_Sbase!()
    set_fbase!()
    return template
end

"Replace one compiled IEEE-39 bus while preserving its original PF model."
function replace_bus(nw, bus_idx::Integer; template = simple_gfldc_template())
    1 <= bus_idx <= Graphs.nv(nw) || throw(ArgumentError("bus_idx must be in 1:nv(nw)"))
    vertices, edges = copy_network_components(nw)
    original_pfmodel = get_pfmodel(vertices[bus_idx])
    vertices[bus_idx] = compile_bus(
        template;
        pf = original_pfmodel,
        vidx = bus_idx,
        name = Symbol("bus$(bus_idx)"),
    )
    replaced = Network(vertices, edges)
    set_jac_prototype!(replaced)
    return replaced
end

"Replace several candidate buses in one network reconstruction."
function replace_buses(nw, bus_idxs; template = simple_gfldc_template())
    vertices, edges = copy_network_components(nw)
    for bus_idx in bus_idxs
        1 <= bus_idx <= Graphs.nv(nw) || throw(ArgumentError("bus_idx must be in 1:nv(nw)"))
        original_pfmodel = get_pfmodel(vertices[bus_idx])
        vertices[bus_idx] = compile_bus(
            template;
            pf = original_pfmodel,
            vidx = bus_idx,
            name = Symbol("bus$(bus_idx)"),
        )
    end
    replaced = Network(vertices, edges)
    set_jac_prototype!(replaced)
    return replaced
end

"""Build a PowerDynamics cross-fade of the original SG and stock GFL port.

Both devices share the same bus voltage. Their current injections are weighted
by `(1-rho)` and `rho`, respectively, while their differential equations remain
unscaled. At the endpoints the zero-share device is removed structurally to
avoid counting canceled dormant-state poles.
"""
function weighted_replacement_network(nw, bus_idx::Integer, rho::Real, beta::Real=1.0)
    0.0 <= rho <= 1.0 || throw(ArgumentError("rho must be in [0,1]"))
    0.9 <= beta <= 1.1 || throw(ArgumentError("PLL scale must be in [0.9,1.1]"))
    rho == 0.0 && return nw
    rho == 1.0 && return replace_bus(nw, bus_idx;
        template=simple_gfldc_template(pll_scale=beta))

    data=ieee39_data()
    busrow=data.bus[findfirst(==(bus_idx),data.bus.bus),:]
    has_avr=Bool(busrow.has_avr); has_gov=Bool(busrow.has_gov)
    (has_avr&&has_gov) || throw(ArgumentError("weighted replacement currently supports controlled SG buses"))

    # Put both stock devices on the same bus. The machine rating scales its
    # terminal current while leaving its ODE untouched; the copied GFL model
    # applies rho only to the external terminal-current equations.
    sg=deepcopy(OfficialIEEE39.controlled_machine)
    set_Sbase!(BASE_MVA); set_fbase!(BASE_FREQ)
    fpll=5.0*beta; ftau=300.0; fi=600.0; xf=0.03; vdc=2.5; cdc=1.25
    gfl=WeightedSimpleGFLDC(
        Xf=xf,Rf=0.01,
        PLL_Kp=fpll*2pi,PLL_Ki=(fpll*2pi)^2/4,
        PLL_τ_lpf=1/(ftau*2pi),
        CC1_KP=(xf/(2pi*BASE_FREQ))*(fi*2pi),
        CC1_KI=(xf/(2pi*BASE_FREQ))*(fi*2pi)^2/4,
        CC1_F=0,CC1_Fcoupl=0,C_dc=cdc,V_dc=vdc,
        kp_v_dc=vdc*cdc*(5*2pi),ki_v_dc=vdc*cdc*(5*2pi)^2/4,
        port_scale=rho,name=:gfl)

    vertices,edges=copy_network_components(nw)
    pf=get_pfmodel(vertices[bus_idx])
    replacement=compile_bus(MTKBus(sg,gfl);name=Symbol("bus$(bus_idx)"),vidx=bus_idx,pf)
    # The network component stores bases at compile time, so restore defaults
    # only after compiling this bus at the case's 100 MVA / 60 Hz bases.
    set_Sbase!(); set_fbase!()
    set_default!(replacement,Symbol("busbar₊Vbase"),Float64(busrow.base_kv))
    OfficialIEEE39.apply_csv_params!(replacement,data.machine,bus_idx)
    OfficialIEEE39.apply_csv_params!(replacement,data.avr,bus_idx)
    OfficialIEEE39.apply_csv_params!(replacement,data.gov,bus_idx)
    mrow=data.machine[findfirst(==(bus_idx),data.machine.bus),:]
    set_default!(replacement,r"Sn$",(1-rho)*Float64(mrow.Sn))
    vertices[bus_idx]=replacement
    replaced=Network(vertices,edges)
    set_jac_prototype!(replaced)
    return replaced
end

"Clone a compiled interior mixture and change only its frozen share and PLL defaults."
function set_weighted_replacement_parameters(nw,bus_idx::Integer,rho::Real,beta::Real)
    0.0 < rho < 1.0 || throw(ArgumentError("parameter updates require an interior rho"))
    0.9 <= beta <= 1.1 || throw(ArgumentError("PLL scale must be in [0.9,1.1]"))
    vertices,edges=copy_network_components(nw)
    v=vertices[bus_idx]
    data=ieee39_data()
    mrow=data.machine[findfirst(==(bus_idx),data.machine.bus),:]
    kp=2pi*5.0*beta
    set_default!(v,r"Sn$",(1-rho)*Float64(mrow.Sn))
    set_default!(v,r"PortScale$",rho)
    set_default!(v,r"PLL_Kp$",kp)
    set_default!(v,r"PLL_Ki$",kp^2/4)
    out=Network(vertices,edges)
    set_jac_prototype!(out)
    return out
end

"Extract the static candidate metadata used in the preregistered campaign."
function candidate_table()
    data = ieee39_data()
    rows = NamedTuple[]
    for bus_idx in CANDIDATE_SG_BUSES
        row_idx = findfirst(data.bus.bus .== bus_idx)
        row_idx === nothing && error("missing bus $bus_idx in official IEEE-39 data")
        row = data.bus[row_idx, :]
        machine_idx = findfirst(data.machine.bus .== bus_idx)
        machine_idx === nothing && error("missing machine $bus_idx in official IEEE-39 data")
        machine = data.machine[machine_idx, :]
        push!(rows, (
            bus = bus_idx,
            category = String(row.category),
            bus_type = String(row.bus_type),
            dispatch_p_pu = Float64(row.P),
            dispatch_p_mw = 100.0 * Float64(row.P),
            voltage_setpoint_pu = Float64(row.V),
            machine_rating_mva = Float64(machine.Sn),
        ))
    end
    return DataFrame(rows)
end

end
