module PhysicalData

using CSV
using DataFrames
using SHA
using TOML
include(joinpath(@__DIR__, "..", "bnd_design", "AnalyticSG.jl"))

export discover_generators, freeze_design_domain, nominal_pll_gains

const SYSTEM_BASE_MVA = 100.0

nominal_pll_gains() = (Kp = 5.0 * 2pi, Ki = (5.0 * 2pi)^2 / 4)
file_sha256(path) = bytes2hex(sha256(read(path)))

"""Discover the SG components from frozen IEEE-39 tables, including local loads.

`bus.P,Q` are net injections. For buses with a colocated ZIP load, its
initialized injection can differ substantially from the CSV setpoint. The
actual SG dispatch comes from its frozen states and independent machine port.
"""
function discover_generators(root::AbstractString)
    input = joinpath(root, "reports", "experiment_D", "inputs")
    buses = CSV.read(joinpath(input, "bus.csv"), DataFrame)
    machines = CSV.read(joinpath(input, "machine.csv"), DataFrame)
    loads = CSV.read(joinpath(input, "load.csv"), DataFrame)
    allunique(Int.(machines.bus)) || error("machine.csv contains duplicate generator buses")
    rows = NamedTuple[]
    for machine in eachrow(sort(machines, :bus))
        bus = Int(machine.bus)
        busidx = findfirst(==(bus), Int.(buses.bus))
        busidx === nothing && error("machine bus $bus is absent from bus.csv")
        b = buses[busidx, :]
        colocated = loads[loads.bus .== bus, :]
        load_p = sum(Float64.(colocated.Pset))
        load_q = sum(Float64.(colocated.Qset))
        net_p = Float64(b.P)
        net_q = Float64(b.Q)
        op = AnalyticSG.frozen_bus_operating_point(root, bus)
        measured = AnalyticSG.steady_state_diagnostics(op.x, op.u, op.parameters)
        measured.rhs_norm < 1e-8 || error("frozen SG equations fail at bus $bus")
        gen_p = measured.network_power.p
        gen_q = measured.network_power.q
        actual_load_p = gen_p - net_p
        actual_load_q = gen_q - net_q
        controlled = Bool(b.has_avr) && Bool(b.has_gov)
        uncontrolled = !Bool(b.has_avr) && !Bool(b.has_gov)
        supported = controlled || uncontrolled
        rating = Float64(machine.Sn)
        replaceable = Bool(b.has_gen) && supported && rating > 0 && gen_p > 0
        exclusion = replaceable ? "" :
            (!Bool(b.has_gen) ? "no generator flag" :
             !supported ? "unsupported partial AVR/governor composition" :
             rating <= 0 ? "nonpositive machine rating" : "nonpositive SG dispatch")
        push!(rows, (
            bus=bus,
            component=String(b.category),
            machine_type=controlled ? "SauerPai_AVRTypeI_TGOV1" : "SauerPai_uncontrolled",
            bus_type=String(b.bus_type),
            SG_dispatch_initial_MW=SYSTEM_BASE_MVA * gen_p,
            SG_reactive_initial_Mvar=SYSTEM_BASE_MVA * gen_q,
            machine_rating_MVA=rating,
            colocated_load_MW=SYSTEM_BASE_MVA * actual_load_p,
            colocated_load_Mvar=SYSTEM_BASE_MVA * actual_load_q,
            nominal_load_setpoint_MW=-SYSTEM_BASE_MVA * load_p,
            nominal_load_setpoint_Mvar=-SYSTEM_BASE_MVA * load_q,
            load_initialization_shift_MW=SYSTEM_BASE_MVA * (actual_load_p + load_p),
            net_bus_injection_MW=SYSTEM_BASE_MVA * net_p,
            net_bus_injection_Mvar=SYSTEM_BASE_MVA * net_q,
            replaceable=replaceable,
            exclusion_reason=exclusion,
        ))
    end
    return DataFrame(rows)
end

"""Freeze independent Kp/Ki rectangles before any ExpE spectral computation."""
function freeze_design_domain(root::AbstractString)
    report = joinpath(root, "reports", "experiment_E")
    tables = joinpath(report, "tables")
    configdir = joinpath(root, "experiments", "bnd_expE", "configs")
    mkpath(tables)
    mkpath(configdir)
    discovered = discover_generators(root)
    CSV.write(joinpath(tables, "TABLE_E01_replaceable_generators.csv"), discovered)
    candidates = discovered[discovered.replaceable .== true, :]
    nrow(candidates) > 0 || error("no replaceable SGs discovered")
    nominal = nominal_pll_gains()
    domain_table = DataFrame(
        bus=Int.(candidates.bus),
        rho_min=fill(0.0, nrow(candidates)),
        rho_max=fill(1.0, nrow(candidates)),
        Kp_nom=fill(nominal.Kp, nrow(candidates)),
        Kp_min=fill(0.25 * nominal.Kp, nrow(candidates)),
        Kp_max=fill(4.0 * nominal.Kp, nrow(candidates)),
        Ki_nom=fill(nominal.Ki, nrow(candidates)),
        Ki_min=fill(0.25 * nominal.Ki, nrow(candidates)),
        Ki_max=fill(4.0 * nominal.Ki, nrow(candidates)),
        gains_independent=fill(true, nrow(candidates)),
    )
    CSV.write(joinpath(tables, "TABLE_E02_design_domain.csv"), domain_table)

    path = joinpath(configdir, "DESIGN_DOMAIN_FROZEN.toml")
    sidecar = path * ".sha256"
    if isfile(path) || isfile(sidecar)
        isfile(path) && isfile(sidecar) || error("incomplete ExpE domain freeze")
        strip(read(sidecar, String)) == file_sha256(path) || error("frozen ExpE domain hash mismatch")
        old = TOML.parsefile(path)
        Int.(old["replaceable_buses"]) == Int.(candidates.bus) ||
            error("frozen generator set differs from current frozen input tables")
        old["Kp_nom"] == nominal.Kp && old["Ki_nom"] == nominal.Ki ||
            error("frozen nominal PLL gains differ from analytic model")
        input = joinpath(root,"reports","experiment_D","inputs")
        for f in ("bus.csv","machine.csv","load.csv")
            old[replace(f,".csv"=>"_csv_sha256")] == file_sha256(joinpath(input,f)) ||
                error("frozen input $f changed after the ExpE design domain was frozen")
        end
        return (path=path, sha256=file_sha256(path), generators=discovered,
            domain=domain_table, status="EXISTING_FROZEN")
    end

    input = joinpath(root, "reports", "experiment_D", "inputs")
    open(path, "w") do io
        println(io, "experiment = \"BND_EXP_E\"")
        println(io, "domain_status = \"FROZEN_BEFORE_EXP_E_STABILITY\"")
        println(io, "system_base_MVA = ", repr(SYSTEM_BASE_MVA))
        println(io, "sigma_req_per_s = 0.05")
        println(io, "gauge_tolerance = 1.0e-8")
        println(io, "replaceable_buses = [", join(Int.(candidates.bus), ", "), "]")
        println(io, "Kp_nom = ", repr(nominal.Kp))
        println(io, "Ki_nom = ", repr(nominal.Ki))
        println(io, "Kp_factor_min = 0.25")
        println(io, "Kp_factor_max = 4.0")
        println(io, "Ki_factor_min = 0.25")
        println(io, "Ki_factor_max = 4.0")
        println(io, "independent_PLL_gains = true")
        println(io, "primary_objective = \"sum(actual_SG_dispatch_MW[i] * rho[i])\"")
        println(io, "epsilon_metric = 1.0e-6")
        println(io, "eta_PLL_metric = 1.0")
        println(io, "gamma_graph_metric = 1.0")
        println(io, "closure_residual_tolerance = 1.0e-9")
        println(io, "KKT_residual_tolerance = 1.0e-7")
        println(io, "active_constraint_tolerance = 1.0e-9")
        println(io, "step_norm_tolerance = 1.0e-8")
        println(io, "objective_relative_tolerance = 1.0e-8")
        println(io, "bus_csv_sha256 = \"", file_sha256(joinpath(input, "bus.csv")), "\"")
        println(io, "machine_csv_sha256 = \"", file_sha256(joinpath(input, "machine.csv")), "\"")
        println(io, "load_csv_sha256 = \"", file_sha256(joinpath(input, "load.csv")), "\"")
        for row in eachrow(domain_table)
            println(io)
            println(io, "[[generator]]")
            println(io, "bus = ", row.bus)
            println(io, "rho_min = ", repr(row.rho_min))
            println(io, "rho_max = ", repr(row.rho_max))
            println(io, "Kp_min = ", repr(row.Kp_min))
            println(io, "Kp_max = ", repr(row.Kp_max))
            println(io, "Ki_min = ", repr(row.Ki_min))
            println(io, "Ki_max = ", repr(row.Ki_max))
        end
    end
    write(sidecar, file_sha256(path) * "\n")
    return (path=path, sha256=file_sha256(path), generators=discovered,
        domain=domain_table, status="NEWLY_FROZEN")
end

end
