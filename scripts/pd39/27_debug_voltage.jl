using CSV
using DataFrames
using PowerDynamics
include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const ROOT = joinpath(@__DIR__, "..", "..")
c = CSV.read(joinpath(ROOT, "results", "PD39_HOLDOUT_CONDITIONS.csv"), DataFrame)[1, :]
bd = zeros(Float64, 46)
bd[1] = 1 / 1.091552734375 - 1
nw = build_confirmatory_network(collect(CANDIDATE_SG_BUSES);
    controller_delta = (c.delta_pll, c.delta_xf, c.delta_cc),
    load_delta = (c.delta_load_p, c.delta_load_q), ibr_delta = c.delta_ibr_p,
    branch_delta = bd, bounds = :repair)
eq = initialize_equilibrium(nw; sparse = false)
println("state_type=", typeof(eq.state), " finite=", eq.state_finite, " fixed=", eq.fixed_point)
rr = run_margin_case(nw)
println("after_margin status=", rr.status, " alpha=", rr.alpha)
eq2 = initialize_equilibrium(nw; sparse = false)
println("after_margin_eq finite=", eq2.state_finite, " fixed=", eq2.fixed_point)
vals = Float64[]
for bus in 1:39
    try
        push!(vals, hypot(Float64(eq2.state[PowerDynamics.VIndex(bus, :busbar₊u_r)]),
                          Float64(eq2.state[PowerDynamics.VIndex(bus, :busbar₊u_i)])))
    catch err
        println("voltage error bus=", bus, " ", sprint(showerror, err))
    end
end
println("voltage_count=", length(vals), " min=", isempty(vals) ? NaN : minimum(vals), " max=", isempty(vals) ? NaN : maximum(vals))
for bus in 1:3
    for sym in (:busbar₊u_r, :busbar₊u_i)
        idx = PowerDynamics.VIndex(bus, sym)
        try
            println(bus, " ", sym, " idx=", idx, " value=", eq.state[idx])
        catch err
            println(bus, " ", sym, " idx=", idx, " error=", sprint(showerror, err))
        end
    end
end
