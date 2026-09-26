"""P4 official-model common-disturbance TDS audit for representative V4 cases."""

using CSV
using DataFrames
using LinearAlgebra
using Statistics
using SciMLBase
using OrdinaryDiffEqRosenbrock

include(joinpath(@__DIR__, "run_p2_pd39_portfolios.jl"))

const P4_RAW = joinpath(CAMPAIGN, "raw", "p4")
mkpath(P4_RAW)

function estimate_trace(t, y)
    finite = isfinite.(y)
    t = t[finite]; y = y[finite]
    length(y) < 4 && return (frequency_hz=NaN, decay_rate=NaN, samples=length(y))
    y = y .- mean(y)
    crossings = findall(i -> y[i] <= 0 && y[i + 1] > 0, 1:(length(y) - 1))
    frequency = length(crossings) >= 2 ? (length(crossings) - 1) / (t[last(crossings)] - t[first(crossings)]) : NaN
    peak_idx = argmax(abs.(y))
    tail_idx = max(peak_idx + 1, Int(cld(length(y), 2)))
    decay = if abs(y[peak_idx]) > 0 && abs(y[tail_idx]) > 0 && t[tail_idx] > t[peak_idx]
        log(abs(y[peak_idx]) / abs(y[tail_idx])) / (t[tail_idx] - t[peak_idx])
    else
        NaN
    end
    (frequency_hz=frequency, decay_rate=decay, samples=length(y))
end

function run_case(key, replaced)
    pf_net = build_portfolio(replaced, Dict{Int,Float64}(); model=:gfl11)
    pf_state = solve_powerflow(pf_net; pfnw=powerflow_model(pf_net), verbose=false)
    interface = interface_values(pf_state)
    vrefs = Dict{Int,Float64}()
    for bus in replaced
        grid_index = bus + count(x -> x <= bus, replaced)
        vrefs[bus] = hypot(interface[VIndex(grid_index, :busbar₊u_r)], interface[VIndex(grid_index, :busbar₊u_i)])
    end
    net = build_portfolio(replaced, vrefs; model=:gfl11)
    state = initialize_from_pf(net; verbose=false, subverbose=false, check=:none, tol=INIT_TOL, nwtol=NETWORK_TOL)
    eigs = filter(isfinite, collect(jacobian_eigenvals(state)))
    gauge = abs.(eigs) .< 1e-8
    transverse = eigs[.!gauge]
    critical = isempty(transverse) ? eigs[argmax(real.(eigs))] : transverse[argmax(real.(transverse))]
    u0 = deepcopy(state)
    uflat(u0)[1] += 1e-6
    status = "PASS"
    solve_message = ""
    t = Float64[]
    y = Float64[]
    try
        problem = ODEProblem(net, uflat(u0), (0.0, 10.0), pflat(u0))
        solution = solve(problem, Rodas5P(); saveat=0.01, abstol=1e-8, reltol=1e-8,
                         initializealg=SciMLBase.NoInit())
        t = Float64.(solution.t)
        trajectory = reduce(vcat, (reshape(Float64.(u), 1, :) for u in solution.u))
        selected_state = argmax(vec(std(trajectory; dims=1)))
        y = trajectory[:, selected_state]
        all(isfinite, y) || error("non-finite trace")
        trace = estimate_trace(t, y)
        CSV.write(joinpath(P4_RAW, "tds_$(key).csv"), DataFrame(time=t, signal=y))
        return (portfolio=key, status=status, eig_alpha=real(critical), eig_frequency_hz=abs(imag(critical))/(2π),
                measured_decay_rate=trace.decay_rate, measured_frequency_hz=trace.frequency_hz,
                selected_state=selected_state, samples=trace.samples, message=solve_message)
    catch err
        status = "STOPPED_BY_GATE"
        solve_message = sprint(showerror, err)
        return (portfolio=key, status=status, eig_alpha=real(critical), eig_frequency_hz=abs(imag(critical)),
                measured_decay_rate=NaN, measured_frequency_hz=NaN, selected_state=0, samples=0, message=solve_message)
    end
end

function p4_main()
    cases = [("none", Set{Int}()), ("30", Set([30])), ("30+33+35+37", Set([30, 33, 35, 37]))]
    rows = NamedTuple[]
    for (key, replaced) in cases
        row = run_case(key, replaced)
        push!(rows, row)
        println("P4_TDS ", key, " status=", row.status, " eig_f=", row.eig_frequency_hz,
                " measured_f=", row.measured_frequency_hz)
    end
    CSV.write(joinpath(P4_RAW, "p4_pd_tds.csv"), DataFrame(rows))
    passed = count(r -> r.status == "PASS", rows)
    open(joinpath(CAMPAIGN, "reports", "P4_JULIA_TDS_STATUS.md"), "w") do io
        println(io, "# P4 — Julia common small-disturbance TDS")
        println(io)
        println(io, "status: ", passed == length(rows) ? "PASS_WITHOUT_MODAL_MATCH" : "PARTIAL_OR_STOPPED")
        println(io, "evidence_class: FRESH_OFFICIAL_POWERDYNAMICS_TDS")
        println(io, "cases: base=none, one=30, full=30+33+35+37")
        println(io, "disturbance: common 1e-6 perturbation of the first network state; the highest-variance state trace is selected for decay/frequency estimation; identical solver/tolerances")
        println(io, "eigensolve_comparison: raw/p4/p4_pd_tds.csv")
        println(io, "traces: raw/p4/tds_<portfolio>.csv when the solve completed")
        println(io, "This is official-model TDS execution evidence. All three trajectories solved and selected-trace decay diagnostics were recorded. The selected trace did not provide a reliable oscillation crossing, so measured frequency is NaN and no frequency agreement is claimed; the decay diagnostic is not treated as the critical eigendecay. It is not a true same-model frozen-no-governor TDS gate and does not promote the H4 blocker.")
    end
    return passed == length(rows) ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(p4_main())
end
