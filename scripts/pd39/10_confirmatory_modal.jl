using CSV
using DataFrames
using Dates
using Printf

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT = joinpath(@__DIR__, "..", "..", "results", "PD39_7OF8_TO_8OF8_MODAL_ANALYSIS.csv")
mkpath(dirname(OUT))

v8 = collect(CANDIDATE_SG_BUSES)
predecessors = Dict(bus => filter(!=(bus), v8) for bus in v8)
rows = NamedTuple[]
mode_shapes = Dict{Tuple{String,String},Vector{ComplexF64}}()

for scenario in uncertainty_scenarios()
    controller_delta = (scenario.pll_scale - 1, scenario.filter_scale - 1,
                        scenario.current_control_scale - 1)
    cases = [("V8", v8)]
    append!(cases, [("missing_$(bus)", predecessors[bus]) for bus in v8])
    for (portfolio_id, portfolio) in cases
        print("modal ", portfolio_id, " ", scenario.id, " ... ")
        nw = try
            build_confirmatory_network(portfolio; controller_delta = controller_delta,
                                       bounds = :discovery)
        catch err
            println("build failed")
            push!(rows, (portfolio = portfolio_string(portfolio), portfolio_id = portfolio_id,
                missing_bus = portfolio_id == "V8" ? 0 : parse(Int, split(portfolio_id, "_")[2]),
                scenario = scenario.id, status = "failed", error_type = "build_failure",
                error_message = sprint(showerror, err), alpha = NaN, margin = NaN,
                stable = false, robust = false, critical_eigenvalue = "",
                critical_frequency_hz = NaN, damping_ratio = NaN,
                critical_mode_family = "", critical_state_labels = "",
                critical_participation = "", equilibrium_residual = NaN,
                jacobian_condition = NaN, g_z_condition = NaN,
                smallest_singular_value = NaN, g_z_smallest_singular_value = NaN,
                eigenvector_condition = NaN, mac_to_v8 = NaN))
            continue
        end
        c = run_confirmatory_case(nw; label = portfolio_id * "/" * scenario.id)
        d = diagnostic_report(c; case_id = portfolio_id * "/" * scenario.id)
        mac = NaN
        if c.modal !== nothing
            try
                mode_shapes[(portfolio_id, scenario.id)] = common_mode_shape(c.state)
            catch err
                d = merge(d, (error_type = "common_observable_projection_failure",
                              error_message = sprint(showerror, err)))
            end
        end
        if portfolio_id != "V8" && haskey(mode_shapes, ("V8", scenario.id)) && haskey(mode_shapes, (portfolio_id, scenario.id))
            mac = modal_assurance(mode_shapes[("V8", scenario.id)], mode_shapes[(portfolio_id, scenario.id)])
        end
        println(c.status, " alpha=", c.modal === nothing ? "NaN" : c.modal.critical_eigenvalue)
        push!(rows, merge(d, (portfolio = portfolio_string(portfolio), portfolio_id = portfolio_id,
                              missing_bus = portfolio_id == "V8" ? 0 : parse(Int, split(portfolio_id, "_")[2]),
                              scenario = scenario.id,
                              robust = c.modal !== nothing && real(c.modal.critical_eigenvalue) <= -0.05,
                              mac_to_v8 = mac)))
    end
end

out = DataFrame(rows)
CSV.write(OUT, out)

# Add the discovery exact comparison in a separate compact table.
disc_path = joinpath(@__DIR__, "..", "..", "results", "pd39", "portfolio_campaign", "portfolio_scenario_results.csv")
disc = CSV.read(disc_path, DataFrame)
cmp = NamedTuple[]
for r in eachrow(out)
    d = filter(x -> x.portfolio == r.portfolio && x.scenario == r.scenario, disc)
    isempty(d) && continue
    dr = only(d)
    push!(cmp, (portfolio = r.portfolio, portfolio_id = r.portfolio_id,
        missing_bus = r.missing_bus, scenario = r.scenario,
        alpha_confirmatory = r.alpha, alpha_discovery = dr.max_real,
        delta_alpha_confirmatory_minus_discovery = r.alpha - dr.max_real,
        margin_confirmatory = r.margin, margin_discovery = dr.dynamic_margin,
        residual = r.equilibrium_residual, mac_to_v8 = r.mac_to_v8))
end
CSV.write(joinpath(dirname(OUT), "PD39_MODAL_DISCOVERY_COMPARISON.csv"), DataFrame(cmp))
println("wrote ", OUT, " rows=", nrow(out))
