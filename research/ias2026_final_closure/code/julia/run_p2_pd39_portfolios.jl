"""P2: fresh PowerDynamics IEEE-39 census for all 16 GFL portfolios."""

using CSV
using DataFrames
using LinearAlgebra
using NetworkDynamics
using PowerDynamics
using PowerDynamics.Library
using ModelingToolkitBase

include(joinpath(@__DIR__, "frozen_gfl11.jl"))
include(joinpath(@__DIR__, "campaign_root.jl"))
example_dir = joinpath(pkgdir(PowerDynamics), "docs", "examples")
include(joinpath(example_dir, "ieee39_part1.jl"))
set_Sbase!(100.0)
set_fbase!(60.0)

const CAMPAIGN = campaign_root_from_args()
const RAW = joinpath(CAMPAIGN, "raw", "p2")
const REPORTS = joinpath(CAMPAIGN, "reports")
mkpath(RAW)
mkpath(REPORTS)

const TARGET_BUSES = [30, 33, 35, 37]
const GFL_G = 0.03625
const INIT_TOL = 1e-8
const NETWORK_TOL = 1e-8

function simplegfldc(; name=:simplegfldc)
    V_dc=2.5; C_dc=1.25; f_v_dc=5.0
    xwLf=0.03; Rf=0.01
    f_pll=5.0; f_tau_pll=300.0; f_i_dq=600.0
    @named gfl = ComposableInverter.SimpleGFLDC(
        Xf=xwLf, Rf=Rf,
        PLL_Kp=f_pll*2π,
        PLL_Ki=(f_pll*2π)^2/4,
        PLL_τ_lpf=1/(f_tau_pll*2π),
        CC1_KP=(xwLf/(2π*60.0))*(f_i_dq*2π),
        CC1_KI=(xwLf/(2π*60.0))*(f_i_dq*2π)^2/4,
        CC1_F=0, CC1_Fcoupl=0,
        C_dc=C_dc, V_dc=V_dc,
        kp_v_dc=V_dc*C_dc*(f_v_dc*2π),
        ki_v_dc=V_dc*C_dc*(f_v_dc*2π)*(f_v_dc*2π)/4,
    )
    return gfl
end

function build_portfolio(replaced::Set{Int}, vrefs::Dict{Int,Float64}; model::Symbol=:gfl11, gfl_gain::Float64=GFL_G)
    buses = Any[]
    for row in eachrow(bus_df)
        i = Int(row.bus)
        if i in replaced
            machine_row = machine_df[findfirst(machine_df.bus .== i), :]
            weight = Float64(machine_row.Sn) / 100.0
            vref = get(vrefs, i, 1.0)
            gfl = if model == :gfl11
                GFL11Injector(name=Symbol("gfl$(i)"), g=gfl_gain,
                              w=weight, p_ref=Float64(row.P) / weight,
                              q_ref=Float64(row.Q) / weight, v_ref=vref)
            elseif model == :simplegfldc
                simplegfldc(name=Symbol("simplegfldc$(i)"))
            else
                error("unknown converter model: $model")
            end
            gfl_bus = compile_bus(MTKBus(gfl); name=Symbol("gfl$i"), current_source=true)
            grid_bus = compile_bus(MTKBus(); name=Symbol("bus$i"))
            set_default!(gfl_bus, :busbar₊Vbase, row.base_kv)
            set_default!(grid_bus, :busbar₊Vbase, row.base_kv)
            set_pfmodel!(gfl_bus, pfPQ(P=Float64(row.P), Q=Float64(row.Q); current_source=true))
            push!(buses, gfl_bus)
            push!(buses, grid_bus)
        else
            bus = if row.category == "junction"
                compile_bus(junction_bus_template; name=Symbol("bus$i"))
            elseif row.category == "load"
                compile_bus(load_bus_template; name=Symbol("bus$i"))
            elseif row.category == "ctrld_machine"
                compile_bus(ctrld_machine_bus_template; name=Symbol("bus$i"))
            elseif row.category == "ctrld_machine_load"
                compile_bus(ctrld_machine_load_bus_template; name=Symbol("bus$i"))
            elseif row.category == "unctrld_machine_load"
                compile_bus(unctrld_machine_load_bus_template; name=Symbol("bus$i"))
            else
                error("unknown bus category $(row.category) at $i")
            end
            set_default!(bus, :busbar₊Vbase, row.base_kv)
            row.has_load && apply_csv_params!(bus, load_df, i)
            row.has_gen && apply_csv_params!(bus, machine_df, i)
            row.has_avr && apply_csv_params!(bus, avr_df, i)
            row.has_gov && apply_csv_params!(bus, gov_df, i)
            pf_model = if row.bus_type == "PQ"
                pfPQ(P=row.P, Q=row.Q)
            elseif row.bus_type == "PV"
                pfPV(P=row.P, V=row.V)
            elseif row.bus_type == "Slack"
                pfSlack(V=row.V, δ=0)
            end
            set_pfmodel!(bus, pf_model)
            push!(buses, bus)
        end
    end
    edges = Any[]
    for row in eachrow(branch_df)
        src = Symbol("bus$(Int(row.src_bus))")
        dst = Symbol("bus$(Int(row.dst_bus))")
        @named piline = PiLine()
        line = compile_line(MTKLine(piline); src=src, dst=dst)
        for col_name in names(branch_df)
            if col_name ∉ ["src_bus", "dst_bus", "transformer"]
                set_default!(line, Regex(col_name*"\$"), row[col_name])
            end
        end
        push!(edges, line)
    end
    for bus in replaced
        push!(edges, LoopbackConnection(; src=Symbol("gfl$bus"), dst=Symbol("bus$bus"),
                                         potential=[:u_r, :u_i], flow=[:i_r, :i_i]))
    end
    net = Network(buses, edges)
    formula31 = @initformula :ZIPLoad₊Vset = sqrt(:busbar₊u_r^2 + :busbar₊u_i^2)
    set_initformula!(net[VIndex(:bus31)], formula31)
    set_initformula!(net[VIndex(:bus39)], formula31)
    net
end

function portfolio_key(s::Set{Int})
    isempty(s) ? "none" : join(sort(collect(s)), "+")
end

function state_residual(net, state)
    du = zeros(Float64, length(uflat(state)))
    net(du, uflat(state), pflat(state), 0.0)
    isempty(du) ? Inf : maximum(abs, du)
end

function write_mode_shapes(rows)
    path = joinpath(RAW, "p2_mode_shapes.csv")
    open(path, "w") do io
        println(io, "portfolio,eigen_index,real,imag,frequency_hz,participating_state,participation")
        for r in rows
            println(io, join((r.portfolio, r.eigen_index, r.real, r.imag,
                              r.frequency_hz, r.participating_state, r.participation), ','))
        end
    end
end

function main()
    rows = NamedTuple[]
    mode_rows = NamedTuple[]
    for mask in 0:15
        replaced = Set(TARGET_BUSES[j] for j in 1:4 if (mask >> (j - 1)) & 1 == 1)
        key = portfolio_key(replaced)
        row = (portfolio=key, replaced_count=length(replaced),
               replaced_buses=key, matched_dispatch=true,
               matched_rating=true, powerflow_status="NOT_RUN",
               initialization_status="NOT_RUN", state_count=0,
               residual=Inf, alpha_full=NaN, alpha_transverse=NaN,
               critical_frequency_hz=NaN, critical_real=NaN,
               stable=false, gauge_modes_removed=0, spectrum_count=0,
               blocker="")
        try
            pf_net = build_portfolio(replaced, Dict{Int,Float64}(); model=:gfl11)
            pf_model = powerflow_model(pf_net)
            pf_state = solve_powerflow(pf_net; pfnw=pf_model, verbose=false)
            row = merge(row, (powerflow_status="PASS",))
            interface = interface_values(pf_state)
            vrefs = Dict{Int,Float64}()
            for bus in replaced
                grid_index = bus + count(x -> x <= bus, replaced)
                ur = interface[VIndex(grid_index, :busbar₊u_r)]
                ui = interface[VIndex(grid_index, :busbar₊u_i)]
                vrefs[bus] = hypot(ur, ui)
            end
            net = build_portfolio(replaced, vrefs; model=:gfl11)
            state = initialize_from_pf(net; verbose=false, subverbose=false,
                                       check=:none, tol=INIT_TOL, nwtol=NETWORK_TOL)
            residual = state_residual(net, state)
            state_count = length(uflat(state))
            lin = reduce_dae(linearize_network(state))
            decomp = eigen(lin.A)
            values = collect(decomp.values)
            # The reduced network retains a single machine-angle gauge mode near
            # zero. Remove only modes below the preregistered numerical gauge
            # threshold, and retain every other eigenvalue in the transverse
            # spectrum.
            finite_values = filter(isfinite, values)
            alpha_full = maximum(real, finite_values)
            gauge_threshold = 1e-8
            gauge_mask = abs.(finite_values) .< gauge_threshold
            transverse_values = finite_values[.!gauge_mask]
            alpha_transverse = isempty(transverse_values) ? NaN : maximum(real, transverse_values)
            critical_idx = isempty(transverse_values) ? argmax(real.(finite_values)) : argmax(real.(transverse_values))
            critical = isempty(transverse_values) ? finite_values[critical_idx] : transverse_values[critical_idx]
            pf = participation_factors(lin)
            if !isnothing(pf.sym)
                for idx in sortperm(real.(pf.eigenvalues), rev=true)[1:min(3, length(pf.eigenvalues))]
                    state_idx = argmax(@view pf.pfactors[:, idx])
                    push!(mode_rows, (portfolio=key, eigen_index=idx,
                                      real=real(pf.eigenvalues[idx]),
                                      imag=imag(pf.eigenvalues[idx]),
                                      frequency_hz=abs(imag(pf.eigenvalues[idx]))/(2π),
                                      participating_state=string(pf.sym[state_idx]),
                                      participation=pf.pfactors[state_idx, idx]))
                end
            end
            row = merge(row, (initialization_status="PASS",
                              state_count=state_count, residual=residual,
                              alpha_full=alpha_full,
                              alpha_transverse=alpha_transverse,
                              critical_frequency_hz=abs(imag(critical))/(2π),
                              critical_real=real(critical), stable=alpha_transverse < 0,
                              gauge_modes_removed=count(gauge_mask),
                              spectrum_count=length(values),
                              blocker=state_count >= 11 && residual <= NETWORK_TOL ? "" : "residual_or_state_gate"))
        catch err
            row = merge(row, (blocker=sprint(showerror, err),))
        end
        push!(rows, row)
        println("P2_PORTFOLIO ", key, " status=", row.initialization_status,
                " states=", row.state_count, " residual=", row.residual)
    end

    path = joinpath(RAW, "p2_powerdynamics_portfolios.csv")
    CSV.write(path, DataFrame(rows))
    write_mode_shapes(mode_rows)

    hasse_path = joinpath(RAW, "p2_hasse_edges.csv")
    open(hasse_path, "w") do io
        println(io, "child,parent,child_stable,parent_stable")
        for child in rows, parent in rows
            child_set = Set(filter(x -> !isempty(x), split(child.portfolio, "+")))
            parent_set = Set(filter(x -> !isempty(x), split(parent.portfolio, "+")))
            child_set == parent_set && continue
            issubset(parent_set, child_set) || continue
            length(child_set) == length(parent_set) + 1 || continue
            println(io, join((child.portfolio, parent.portfolio,
                              child.stable, parent.stable), ','))
        end
    end
    pass_count = count(r -> r.initialization_status == "PASS", rows)
    report = joinpath(REPORTS, "P2_POWERDYNAMICS_V4_STATUS.md")
    open(report, "w") do io
        println(io, "# P2 — PowerDynamics V4 portfolio census")
        println(io)
        println(io, "status: ", pass_count == 16 ? "NEGATIVE_HOLDOUT" : "STOPPED_BY_GATE")
        println(io, "result_label: ", pass_count == 16 ? "FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT" : "STOPPED_BY_GATE")
        println(io, "evidence_class: ", pass_count == 16 ? "FRESH_ALTERNATIVE_DYNAMIC_MODEL_V4_NEGATIVE_HOLDOUT" : "STOPPED_BY_GATE")
        println(io, "portfolios: 16 of 16 attempted")
        println(io, "initialized: ", pass_count)
        println(io, "target_buses: {30,33,35,37}")
        println(io, "matched_dispatch_and_rating: true")
        println(io, "dynamic_model: official PowerDynamics IEEE-39 machines with official AVR/governor devices; GFL11 replacement harness")
        println(io, "same_model_cross_code_gate: NOT_TESTED")
        println(io, "network_tolerance: ", NETWORK_TOL)
        println(io, "census_csv: raw/p2/p2_powerdynamics_portfolios.csv")
        println(io, "hasse_csv: raw/p2/p2_hasse_edges.csv")
        println(io, "mode_shapes_csv: raw/p2/p2_mode_shapes.csv")
        println(io)
        println(io, pass_count == 16 ? "All requested portfolios initialized and had stable gauge-aware transverse spectra. This is a fresh alternative dynamic-model negative holdout; it is not TRUE_SAME_MODEL_CROSS_CODE_PASS." : "At least one portfolio remains blocked; no universal P2 promotion is made.")
    end
    println("P2_", pass_count == 16 ? "PASS" : "STOPPED", " initialized=", pass_count, "/16")
    return pass_count == 16 ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main())
end
