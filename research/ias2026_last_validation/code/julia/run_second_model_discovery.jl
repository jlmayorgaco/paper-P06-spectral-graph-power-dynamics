"""Second-model policy discovery for official PowerDynamics SimpleGFLDC.

The frozen operating-point policy is retained.  Only the PLL and current-loop
bandwidth multipliers are varied over a preregistered physically interpretable
grid, and all 16 replacement subsets are evaluated at each point.
"""

using CSV
using DataFrames
using LinearAlgebra
using NetworkDynamics
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase

include(joinpath(@__DIR__, "run_p2_pd39_portfolios.jl"))

const DISC_RAW = joinpath(CAMPAIGN, "raw", "second_model_search")
const DISC_REPORT = joinpath(CAMPAIGN, "reports", "SECOND_MODEL_POLICY_DISCOVERY.md")
mkpath(DISC_RAW)

const POLICY_SCALES = (0.125, 0.5, 1.0, 2.0, 8.0)
const DISCOVERY_PORTFOLIOS = (Set([30, 33, 35]), Set([30, 33, 35, 37]))

function discovery_case(replaced::Set{Int}, axis::String, scale::Float64)
    kwargs = axis == "pll" ? (simple_pll_scale=scale, simple_current_scale=1.0) : (simple_pll_scale=1.0, simple_current_scale=scale)
    pf_net = build_portfolio(replaced, Dict{Int,Float64}(); model=:simplegfldc, kwargs...)
    pf_state = solve_powerflow(pf_net; pfnw=powerflow_model(pf_net), verbose=false)
    interface = interface_values(pf_state)
    vrefs = Dict{Int,Float64}()
    for bus in replaced
        grid_index = bus + count(x -> x <= bus, replaced)
        ur = interface[VIndex(grid_index, :busbar₊u_r)]
        ui = interface[VIndex(grid_index, :busbar₊u_i)]
        vrefs[bus] = hypot(ur, ui)
    end
    net = build_portfolio(replaced, vrefs; model=:simplegfldc, kwargs...)
    state = initialize_from_pf(net; verbose=false, subverbose=false,
                               check=:none, tol=INIT_TOL, nwtol=NETWORK_TOL)
    residual = state_residual(net, state)
    lin = reduce_dae(linearize_network(state))
    values = filter(isfinite, collect(eigen(lin.A).values))
    gauge = abs.(values) .< 1e-8
    transverse = values[.!gauge]
    alpha = isempty(transverse) ? NaN : maximum(real, transverse)
    critical = isempty(transverse) ? values[argmax(real.(values))] : transverse[argmax(real.(transverse))]
    (powerflow_status="PASS", initialization_status="PASS", state_count=length(uflat(state)),
     residual=residual, alpha_transverse=alpha, critical_frequency_hz=abs(imag(critical))/(2π),
     stable=isfinite(alpha) && alpha < 0.0, gauge_modes_removed=count(gauge),
     spectrum_count=length(values), error="")
end

function run_discovery()
    rows = NamedTuple[]
    for axis in ("pll", "current")
        for scale in POLICY_SCALES
            for replaced in DISCOVERY_PORTFOLIOS
                portfolio = portfolio_key(replaced)
                result = try
                    discovery_case(replaced, axis, scale)
                catch err
                    (powerflow_status="FAIL", initialization_status="FAIL", state_count=0,
                     residual=Inf, alpha_transverse=NaN, critical_frequency_hz=NaN,
                     stable=false, gauge_modes_removed=0, spectrum_count=0,
                     error=sprint(showerror, err))
                end
                push!(rows, merge((axis=axis, scale=scale, portfolio=portfolio,
                                  replaced_count=length(replaced), model="PowerDynamics.ComposableInverter.SimpleGFLDC",
                                  matched_dispatch_PQ_system_base=true, matched_rating=false), result))
            end
            println("SECOND_MODEL_DISCOVERY axis=", axis, " scale=", scale,
                    " pass=", count(r -> r.axis == axis && r.scale == scale && r.initialization_status == "PASS", rows),
                    " unstable=", count(r -> r.axis == axis && r.scale == scale && r.initialization_status == "PASS" && !r.stable, rows))
        end
    end
    path = joinpath(DISC_RAW, "simplegfldc_policy_discovery.csv")
    CSV.write(path, DataFrame(rows))
    groups = Dict{Tuple{String,Float64},Vector{NamedTuple}}()
    for row in rows
        key = (row.axis, row.scale)
        push!(get!(groups, key, NamedTuple[]), row)
    end
    mixed = Tuple{String,Float64}[]
    open(DISC_REPORT, "w") do io
        println(io, "# Second-model policy discovery")
        println(io)
        println(io, "status: `FRESH_SECOND_MODEL_POLICY_DISCOVERY`")
        println(io, "model: `PowerDynamics.ComposableInverter.SimpleGFLDC`")
        println(io, "policy coordinates: PLL bandwidth multiplier and current-loop bandwidth multiplier")
        println(io, "grid: ", join(POLICY_SCALES, ", "))
        println(io, "representative portfolios: 30+33+35 (proper stable reference), 30+33+35+37 (flagship blocker candidate)")
        println(io)
        println(io, "| coordinate | scale | initialized | stable | unstable | mixed |")
        println(io, "|---|---:|---:|---:|---:|---|")
        for axis in ("pll", "current")
            for scale in POLICY_SCALES
                group = groups[(axis, scale)]
                pass = filter(r -> r.initialization_status == "PASS", group)
                stable = count(r -> r.stable, pass)
                unstable = length(pass) - stable
                is_mixed = stable > 0 && unstable > 0
                is_mixed && push!(mixed, (axis, scale))
                println(io, "| ", axis, " | ", scale, " | ", length(pass), " | ", stable, " | ", unstable, " | ", is_mixed, " |")
            end
        end
        println(io)
        if isempty(mixed)
            println(io, "No mixed stable/unstable policy region was found on this frozen second-model grid. This negative search result is recorded rather than converted into a mixed holdout.")
        else
            println(io, "Mixed policy groups: ", mixed)
        end
        println(io)
        println(io, "All classifications use gauge-aware transverse spectra, a common matched-dispatch system-base operating point, and the corrected rating scope: matched_rating=false and rating_equivalence_established=false.")
    end
    open(joinpath(DISC_RAW, "discovery_summary.json"), "w") do io
        println(io, "{\"status\":\"FRESH_SECOND_MODEL_POLICY_DISCOVERY\",\"mixed_policy_groups\":\"", join(string.(mixed), ";"), "\",\"scales\":\"", join(POLICY_SCALES, ","), "\"}")
    end
    isempty(mixed) ? 0 : 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(run_discovery())
end
