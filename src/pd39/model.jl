module PD39Model

using CSV
using DataFrames
using Graphs
using ModelingToolkitBase
using NetworkDynamics
using PowerDynamics
using PowerDynamics.Library

export POWERDYNAMICS_IEEE39_EXAMPLE, CANDIDATE_SG_BUSES, SLACK_BUS,
    BASE_MVA, BASE_FREQ, GFL_MODEL_ID, ieee39_data, copy_network_components,
    baseline_network, simple_gfldc_template, replace_bus, replace_buses,
    candidate_table

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
