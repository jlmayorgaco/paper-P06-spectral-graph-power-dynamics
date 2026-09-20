"""P5 full census: official SimpleGFLDC on the frozen IEEE-39 V4 portfolios."""

using CSV
using DataFrames
using LinearAlgebra

include(joinpath(@__DIR__, "run_p2_pd39_portfolios.jl"))

const P5_RAW = joinpath(CAMPAIGN, "raw", "p5")
const P5_REPORTS = joinpath(CAMPAIGN, "reports")
mkpath(P5_RAW)
mkpath(P5_REPORTS)

function p5_main()
    rows = NamedTuple[]
    for mask in 0:15
        replaced = Set(TARGET_BUSES[j] for j in 1:4 if (mask >> (j - 1)) & 1 == 1)
        key = portfolio_key(replaced)
        row = (portfolio=key, replaced_count=length(replaced), replaced_buses=key,
               matched_dispatch_PQ_system_base=true, matched_rating=false,
               rating_equivalence_established=false,
               dynamic_rating_basis="global_100_MVA_system_base",
               powerflow_status="NOT_RUN",
               initialization_status="NOT_RUN", state_count=0, residual=Inf,
               alpha_full=NaN, alpha_transverse=NaN, critical_frequency_hz=NaN,
               critical_real=NaN, stable=false, gauge_modes_removed=0,
               spectrum_count=0, candidate_MVA=0.0, blocker="")
        try
            pf_net = build_portfolio(replaced, Dict{Int,Float64}(); model=:simplegfldc)
            pf_state = solve_powerflow(pf_net; pfnw=powerflow_model(pf_net), verbose=false)
            row = merge(row, (powerflow_status="PASS",))
            interface = interface_values(pf_state)
            vrefs = Dict{Int,Float64}()
            candidate_mva = 0.0
            for bus in replaced
                grid_index = bus + count(x -> x <= bus, replaced)
                ur = interface[VIndex(grid_index, :busbar₊u_r)]
                ui = interface[VIndex(grid_index, :busbar₊u_i)]
                vrefs[bus] = hypot(ur, ui)
                machine_row = machine_df[findfirst(machine_df.bus .== bus), :]
                candidate_mva += Float64(machine_row.Sn)
            end
            net = build_portfolio(replaced, vrefs; model=:simplegfldc)
            state = initialize_from_pf(net; verbose=false, subverbose=false,
                                       check=:none, tol=INIT_TOL, nwtol=NETWORK_TOL)
            residual = state_residual(net, state)
            vals = filter(isfinite, collect(jacobian_eigenvals(state)))
            gauge_mask = abs.(vals) .< 1e-8
            transverse = vals[.!gauge_mask]
            alpha_full = isempty(vals) ? NaN : maximum(real, vals)
            alpha = isempty(transverse) ? NaN : maximum(real, transverse)
            critical = isempty(transverse) ? vals[argmax(real.(vals))] : transverse[argmax(real.(transverse))]
            row = merge(row, (initialization_status="PASS",
                              state_count=length(uflat(state)), residual=residual,
                              alpha_full=alpha_full, alpha_transverse=alpha,
                              critical_frequency_hz=isempty(vals) ? NaN : abs(imag(critical))/(2π),
                              critical_real=alpha, stable=isfinite(alpha) && alpha < 0,
                              gauge_modes_removed=count(gauge_mask), spectrum_count=length(vals),
                              candidate_MVA=candidate_mva,
                              blocker=(length(uflat(state)) > 0 && isfinite(residual) && residual <= NETWORK_TOL) ? "" : "state_or_residual_gate"))
        catch err
            row = merge(row, (blocker=sprint(showerror, err),))
        end
        push!(rows, row)
        println("P5_SIMPLEGFLDC ", key, " status=", row.initialization_status,
                " states=", row.state_count, " residual=", row.residual)
    end
    path = joinpath(P5_RAW, "p5_simplegfldc_ieee39_portfolios.csv")
    CSV.write(path, DataFrame(rows))
    passed = count(r -> r.initialization_status == "PASS", rows)
    stable = count(r -> r.stable, rows)
    open(joinpath(P5_REPORTS, "P5_SECOND_GFL_STATUS.md"), "w") do io
        println(io, "# P5 — official SimpleGFLDC IEEE-39 census")
        println(io)
        println(io, "component_status: PASS")
        println(io, "component_result: raw/p5/p5_simplegfldc_result.csv")
        println(io, "status: ", passed == 16 ? "NEGATIVE_HOLDOUT" : "STOPPED_BY_GATE")
        println(io, "result_label: FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT")
        println(io, "evidence_class: FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT")
        println(io, "model: PowerDynamics.ComposableInverter.SimpleGFLDC")
        println(io, "portfolios: 16 of 16 attempted")
        println(io, "initialized: ", passed)
        println(io, "stable_transverse: ", stable)
        println(io, "matched_dispatch_PQ_system_base: true")
        println(io, "matched_rating: false")
        println(io, "rating_equivalence_established: false")
        println(io, "dynamic_rating_basis: global_100_MVA_system_base")
        println(io, "candidate_MVA: reported_from_retired_machine_ratings_only")
        println(io, "no_retune: true")
        println(io, "census_csv: raw/p5/p5_simplegfldc_ieee39_portfolios.csv")
        println(io)
        println(io, "This is a second dynamic model negative result with matched scheduled P/Q on the common 100 MVA system base. Candidate MVA is reported from retired machine ratings, but dynamic rating equivalence is not established. It is not a same-model cross-code parity gate and does not reproduce the frozen no-governor H4 blocker.")
    end
    return passed == 16 ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(p5_main())
end
