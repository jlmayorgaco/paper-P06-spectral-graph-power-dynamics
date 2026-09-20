"""Execute one phase of the preregistered alternative-model P6 holdout."""

using CSV
using DataFrames
using LinearAlgebra

include(joinpath(@__DIR__, "run_p2_pd39_portfolios.jl"))

function argvalue(flag, default)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("$flag requires a value")
    return ARGS[idx + 1]
end

function run_phase(phase)
    input_path = joinpath(CAMPAIGN, "prereg", "p6_genuine_holdout_inputs.csv")
    inputs = CSV.read(input_path, DataFrame)
    selected = filter(r -> String(r.phase) == phase, eachrow(inputs))
    output_path = joinpath(CAMPAIGN, "raw", "p6", "p6_genuine_$(phase)_labels.csv")
    rows = NamedTuple[]
    for row in selected
        replaced = Set(parse.(Int, split(String(row.portfolio), "+")))
        gain = Float64(row.gfl_gain)
        key = String(row.portfolio)
        status = "PASS"
        residual = Inf
        alpha = NaN
        frequency = NaN
        message = ""
        try
            pf_net = build_portfolio(replaced, Dict{Int,Float64}(); model=:gfl11, gfl_gain=gain)
            pf_state = solve_powerflow(pf_net; pfnw=powerflow_model(pf_net), verbose=false)
            interface = interface_values(pf_state)
            vrefs = Dict{Int,Float64}()
            for bus in replaced
                grid_index = bus + count(x -> x <= bus, replaced)
                vrefs[bus] = hypot(interface[VIndex(grid_index, :busbar₊u_r)], interface[VIndex(grid_index, :busbar₊u_i)])
            end
            net = build_portfolio(replaced, vrefs; model=:gfl11, gfl_gain=gain)
            state = initialize_from_pf(net; verbose=false, subverbose=false, check=:none, tol=INIT_TOL, nwtol=NETWORK_TOL)
            residual = state_residual(net, state)
            vals = filter(isfinite, collect(jacobian_eigenvals(state)))
            transverse = vals[.! (abs.(vals) .< 1e-8)]
            critical = isempty(transverse) ? vals[argmax(real.(vals))] : transverse[argmax(real.(transverse))]
            alpha = real(critical)
            frequency = abs(imag(critical))/(2π)
        catch err
            status = "STOPPED_BY_GATE"
            message = sprint(showerror, err)
        end
        push!(rows, (case_id=String(row.case_id), phase=phase, portfolio=key, gfl_gain=gain,
                     stable=status == "PASS" && alpha < 0, alpha_transverse=alpha,
                     critical_frequency_hz=frequency, residual=residual, status=status,
                     message=message))
        println("P6_GENUINE ", phase, " ", row.case_id, " status=", status, " alpha=", alpha)
    end
    CSV.write(output_path, DataFrame(rows))
end

phase = argvalue("--phase", "discovery")
phase in ("discovery", "holdout") || error("phase must be discovery or holdout")
run_phase(phase)
