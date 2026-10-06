using CSV, DataFrames, LinearAlgebra, NetworkDynamics, PowerDynamics, SHA

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_N")
include(joinpath(ROOT, "src", "pd39", "PD39.jl"))
using .PD39

function write_once(path, table)
    isfile(path) && error("immutable ExpN original operating point already exists: $path")
    CSV.write(path, DataFrame(table))
    write(path * ".sha256", bytes2hex(sha256(read(path))) * "\n")
end

function main()
    nw = PD39.baseline_network()
    eq = PD39.initialize_equilibrium(nw; sparse=false, check=:error)
    sys = linearize_network(eq.state)
    n = size(sys.A, 1)
    mass = sys.M isa UniformScaling ? ones(n) : diag(sys.M)
    names = string.(sys.sym)
    values = Float64.(uflat(eq.state))
    psyms = string.(NetworkDynamics.SII.parameter_symbols(nw))
    pvals = Float64.(pflat(eq.state))
    length(psyms) == length(pvals) || error("parameter map length mismatch")
    parameters = Dict(zip(psyms, pvals))
    param(bus, suffix) = get(parameters, "VIndex($bus, :$suffix)", NaN)
    op = NamedTuple[]
    states = NamedTuple[]
    refs = NamedTuple[]
    for bus in 30:39
        prefix = bus == 39 ? "machine₊" : "ctrld_gen₊machine₊"
        vr = Float64(eq.state[VIndex(bus, :busbar₊u_r)])
        vi = Float64(eq.state[VIndex(bus, :busbar₊u_i)])
        sn = Float64(eq.state[VIndex(bus, Symbol(prefix * "Sn"))])
        h = Float64(eq.state[VIndex(bus, Symbol(prefix * "H"))])
        pg = sn * Float64(eq.state[VIndex(bus, Symbol(prefix * "P"))])
        qg = sn * Float64(eq.state[VIndex(bus, Symbol(prefix * "Q"))])
        pl = bus in (31, 39) ? 100 * Float64(eq.state[VIndex(bus, :ZIPLoad₊P)]) : 0.0
        ql = bus in (31, 39) ? 100 * Float64(eq.state[VIndex(bus, :ZIPLoad₊Q)]) : 0.0
        pn = Float64(eq.state[VIndex(bus, :busbar₊P_MW)])
        qn = Float64(eq.state[VIndex(bus, :busbar₊Q_MVAr)])
        abs(pg + pl - pn) < 1e-8 || error("P split invalid at bus $bus")
        abs(qg + ql - qn) < 1e-8 || error("Q split invalid at bus $bus")
        push!(op, (;bus, V_real_pu=vr, V_imag_pu=vi, V_abs_pu=hypot(vr, vi),
            P_gen_MW=pg, Q_gen_Mvar=qg, S_gen_real_MW=pg, S_gen_imag_Mvar=qg,
            Sn_original_MVA=sn, H_seconds=h, kinetic_HSn_MVA_s=h*sn,
            P_load_MW=pl, Q_load_Mvar=ql, P_net_MW=pn, Q_net_Mvar=qn,
            avr_vref_pu=param(bus, "ctrld_gen₊avr₊vref"),
            governor_p_ref_pu=param(bus, "ctrld_gen₊gov₊p_ref"),
            governor_omega_ref_pu=param(bus, "ctrld_gen₊gov₊ω_ref"),
            machine_torque_set_pu=param(bus, "machine₊τ_m_set"),
            zip_Pset_pu=param(bus, "ZIPLoad₊Pset"),
            zip_Qset_pu=param(bus, "ZIPLoad₊Qset"),
            zip_Vset_pu=param(bus, "ZIPLoad₊Vset"),
            source="independent_PowerDynamics_initialized_component"))
        for i in eachindex(names)
            occursin("VIndex($bus, :", names[i]) || continue
            push!(states, (;bus, state_index=i, state_name=names[i],
                differential=mass[i] == 1.0, equilibrium_value=values[i]))
        end
        for (symbol, value) in parameters
            startswith(symbol, "VIndex($bus, :") || continue
            (occursin("₊vref)", symbol) || occursin("₊p_ref)", symbol) ||
             occursin("₊ω_ref)", symbol) || occursin("₊τ_m_set)", symbol) ||
             occursin("ZIPLoad₊", symbol)) || continue
            push!(refs, (;bus, parameter=symbol, initialized_value=value))
        end
    end
    total = sum(row.P_gen_MW for row in op)
    abs(total - 5402.761089978776) < 1e-8 ||
        @warn("New component-level SG total differs from previous", total)
    mkpath(OUT)
    write_once(joinpath(OUT, "TABLE_N01_original_operating_point.csv"), op)
    write_once(joinpath(OUT, "TABLE_N01_dynamic_states.csv"), states)
    write_once(joinpath(OUT, "TABLE_N01_references.csv"), refs)
    println("N01_FROZEN total_original_SG_MW=", total,
        " state_rows=", length(states), " reference_rows=", length(refs))
end

main()
