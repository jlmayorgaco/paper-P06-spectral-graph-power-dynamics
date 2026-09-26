using CSV
using DataFrames
using Graphs
using NetworkDynamics
using PowerDynamics

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
const OUT = joinpath(ROOT, "results", "PD39_CONFIRMATORY_OUTAGE_DIAGNOSTICS.csv")
const OUTAGE_VOLTAGE_MIN_PU = 0.90
const OUTAGE_VOLTAGE_MAX_PU = 1.10
mkpath(dirname(OUT))

function voltage_range(s)
    v = Float64[]
    for bus in 1:39
        try
            ur = Float64(s[PowerDynamics.VIndex(bus, :busbar₊u_r)])
            ui = Float64(s[PowerDynamics.VIndex(bus, :busbar₊u_i)])
            push!(v, hypot(ur, ui))
        catch
        end
    end
    return isempty(v) ? (NaN, NaN) : (minimum(v), maximum(v))
end

base = baseline_network()
orig = CSV.read(joinpath(ROOT, "results", "pd39", "link_outage", "link_outage_results.csv"), DataFrame)
data = ieee39_data()
rows = NamedTuple[]
for eidx in 1:Graphs.ne(base)
    vertices, edges = copy_network_components(base)
    deleteat!(edges, eidx)
    nw = Network(vertices, edges)
    src, dst = Int(data.branch[eidx, :src_bus]), Int(data.branch[eidx, :dst_bus])
    connected = try is_connected(SimpleGraph(nw.im.g)) catch; false end
    stage, type, message = "graph", "", ""
    pffinite, statefinite, fixed, alpha = false, false, false, NaN
    vmin, vmax = NaN, NaN
    if !connected
        type, message = "graph_disconnected", "undirected network graph disconnected"
    else
        local pf
        try
            pf = solve_powerflow(nw; verbose = false, sparse = false)
            pffinite = all(isfinite, uflat(pf))
        catch err
            stage = "powerflow"
            type = occursin("converg", lowercase(sprint(showerror, err))) ? "powerflow_divergence" : "powerflow_solver_failure"
            message = sprint(showerror, err)
        end
        if isempty(type) && !pffinite
            stage, type, message = "powerflow", "powerflow_nonfinite", "PF state is nonfinite"
        end
        if isempty(type)
            local s
            try
                s = initialize_from_pf!(nw; pfs = pf, verbose = false, sparsepf = false)
                statefinite = all(isfinite, uflat(s))
            catch err
                stage, type, message = "dynamic_initialization", "dynamic_initialization_failure", sprint(showerror, err)
            end
            if isempty(type) && !statefinite
                stage, type, message = "dynamic_initialization", "dynamic_initialization_nonfinite", "dynamic state is nonfinite"
            end
            if isempty(type)
                fixed = isfixpoint(s; tol = 1e-8)
                if !fixed
                    stage, type, message = "fixed_point", "fixed_point_failure", "fixed-point residual exceeds tolerance"
                end
                vmin, vmax = voltage_range(s)
                if isempty(type) && (vmin < OUTAGE_VOLTAGE_MIN_PU || vmax > OUTAGE_VOLTAGE_MAX_PU)
                    stage, type, message = "voltage", "voltage_infeasibility", "voltage outside frozen confirmatory interval"
                end
            end
            if isempty(type)
                try
                    alpha = stability_audit(s).max_real
                    stage = "spectrum"
                    type = alpha >= 0 ? "true_instability" : "qualified_stable"
                catch err
                    stage, type, message = "spectrum", "spectrum_failure", sprint(showerror, err)
                end
            end
        end
    end
    old = only(filter(r -> r.link == eidx, orig))
    push!(rows, (link = eidx, src_bus = src, dst_bus = dst,
        discovery_status = old.equilibrium_status, diagnostic_stage = stage,
        diagnostic_type = type, diagnostic_message = message,
        graph_connected = connected, powerflow_finite = pffinite,
        dynamic_state_finite = statefinite, fixed_point = fixed,
        voltage_min_pu = vmin, voltage_max_pu = vmax, alpha = alpha,
        discovery_dynamic_margin = old.dynamic_margin,
        discovery_static_weakness = old.static_weakness))
    println(eidx, " ", type)
end
CSV.write(OUT, DataFrame(rows))
println("wrote ", OUT, " rows=", length(rows))
