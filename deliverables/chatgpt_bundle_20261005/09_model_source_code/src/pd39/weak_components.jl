using DataFrames
using Graphs
using NetworkDynamics

export static_weak_nodes, static_weak_links, dynamic_node_scores, dynamic_weak_links

"""
Static node proxy fixed before dynamic results: inverse incident electrical
strength, with branch strength defined as 1/|X| in the official CSV model.
"""
function static_weak_nodes(data = ieee39_data())
    buses = Int.(data.bus.bus)
    strength = Dict(bus => 0.0 for bus in buses)
    for row in eachrow(data.branch)
        branch_strength = 1 / max(abs(Float64(row.X)), eps(Float64))
        strength[Int(row.src_bus)] += branch_strength
        strength[Int(row.dst_bus)] += branch_strength
    end
    return DataFrame(
        bus = buses,
        electrical_strength = [strength[bus] for bus in buses],
        static_weakness = [1 / max(strength[bus], eps(Float64)) for bus in buses],
    )
end

"Static link proxy fixed before dynamic results: larger series reactance is weaker."
function static_weak_links(data = ieee39_data())
    return DataFrame(
        link = collect(1:nrow(data.branch)),
        src_bus = Int.(data.branch.src_bus),
        dst_bus = Int.(data.branch.dst_bus),
        series_reactance = abs.(Float64.(data.branch.X)),
        static_weakness = abs.(Float64.(data.branch.X)),
    )
end

"Compute dynamic single-replacement margin loss from a campaign table."
function dynamic_node_scores(results; baseline_margin)
    required = (:bus, :dynamic_margin)
    all(name -> name in propertynames(results), required) ||
        throw(ArgumentError("results must contain bus and dynamic_margin columns"))
    out = copy(results)
    out[!, :dynamic_margin_loss] = baseline_margin .- out.dynamic_margin
    return sort(out, :dynamic_margin_loss, rev = true)
end

"""
Evaluate declared one-at-a-time branch-removal diagnostics.

This is intentionally separate from the node replacement campaign: each
branch is removed, the PF is re-solved, and failures remain explicit.
"""
function dynamic_weak_links(base_nw; base_margin = nothing)
    rows = NamedTuple[]
    base_margin_value = if isnothing(base_margin)
        base_eq = initialize_equilibrium(base_nw; sparse = false)
        stability_audit(base_eq.state).dynamic_margin
    else
        Float64(base_margin)
    end
    for link_idx in 1:Graphs.ne(base_nw)
        vertices, edges = copy_network_components(base_nw)
        deleteat!(edges, link_idx)
        row = (
            link = link_idx,
            equilibrium_status = "equilibrium_failed",
            dynamic_margin = NaN,
            max_real = NaN,
            stable = false,
            margin_loss = NaN,
        )
        try
            outage_nw = Network(vertices, edges)
            eq = initialize_equilibrium(outage_nw; sparse = false)
            if eq.powerflow_finite && eq.state_finite && eq.fixed_point
                stab = stability_audit(eq.state)
                row = (
                    link = link_idx,
                    equilibrium_status = "ok",
                    dynamic_margin = stab.dynamic_margin,
                    max_real = stab.max_real,
                    stable = stab.stable,
                    margin_loss = base_margin_value - stab.dynamic_margin,
                )
            end
        catch
            # Keep the declared failure label; the campaign script should add
            # exception text if detailed error provenance is required.
        end
        push!(rows, row)
    end
    return DataFrame(rows)
end
