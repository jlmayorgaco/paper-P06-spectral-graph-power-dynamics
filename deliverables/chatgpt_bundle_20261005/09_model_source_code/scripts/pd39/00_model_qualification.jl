using CSV
using DataFrames
using Dates
using Graphs
using NetworkDynamics
using PowerDynamics
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT_DIR = joinpath(@__DIR__, "..", "..", "results", "pd39", "model_qualification")
mkpath(OUT_DIR)

nw = baseline_network()
data = ieee39_data()
pfs = solve_powerflow(nw; verbose = false, sparse = false)
s0 = initialize_from_pf!(nw; pfs = pfs, verbose = false, sparsepf = false)

gfl_template = simple_gfldc_template()
gfl_nw = replace_bus(nw, 32; template = gfl_template)
gfl_pfs = solve_powerflow(gfl_nw; verbose = false, sparse = false)
gfl_s0 = initialize_from_pf!(gfl_nw; pfs = gfl_pfs, verbose = false, sparsepf = false)

qualification = Dict(
    "timestamp_utc" => string(now(UTC)),
    "julia_version" => string(VERSION),
    "powerdynamics_source" => pkgdir(PowerDynamics),
    "powerdynamics_version" => "5.0.0",
    "powerdynamics_model" => GFL_MODEL_ID,
    "network_buses" => Graphs.nv(nw),
    "network_branches" => Graphs.ne(nw),
    "official_bus_rows" => nrow(data.bus),
    "official_branch_rows" => nrow(data.branch),
    "candidate_sg_buses" => CANDIDATE_SG_BUSES,
    "slack_bus" => SLACK_BUS,
    "baseline_state_dim" => length(s0),
    "gfl_replacement_state_dim" => length(gfl_s0),
    "baseline_bus39_u_r" => s0[VIndex(39, :busbar₊u_r)],
    "baseline_bus39_u_i" => s0[VIndex(39, :busbar₊u_i)],
    "gfl_bus32_u_r" => gfl_s0[VIndex(32, :busbar₊u_r)],
    "baseline_powerflow_finite" => all(isfinite, uflat(pfs)),
    "baseline_dynamic_initialization_finite" => all(isfinite, uflat(s0)),
    "gfl_powerflow_finite" => all(isfinite, uflat(gfl_pfs)),
    "gfl_dynamic_initialization_finite" => all(isfinite, uflat(gfl_s0)),
)

CSV.write(joinpath(OUT_DIR, "candidate_table.csv"), candidate_table())
open(joinpath(OUT_DIR, "qualification.toml"), "w") do io
    TOML.print(io, qualification)
end

println("PD39 model qualification complete")
println("baseline_state_dim=$(length(s0))")
println("gfl_replacement_state_dim=$(length(gfl_s0))")
println("candidate_buses=$(join(CANDIDATE_SG_BUSES, ','))")
