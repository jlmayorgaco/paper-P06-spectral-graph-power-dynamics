using CSV
using DataFrames
using Dates
using NetworkDynamics
using PowerDynamics
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))
using .PD39

const OUT_DIR = joinpath(@__DIR__, "..", "..", "results", "pd39", "baseline")
mkpath(OUT_DIR)

nw = baseline_network()
eq = initialize_equilibrium(nw; sparse = false)
eq.powerflow_finite || error("baseline PF is not finite")
eq.state_finite || error("baseline dynamic state is not finite")
eq.fixed_point || error("baseline fixed-point check failed")
stab = stability_audit(eq.state)

top = dominant_modes(stab; n = min(20, length(stab.nontrivial_eigenvalues)))
CSV.write(
    joinpath(OUT_DIR, "dominant_modes.csv"),
    DataFrame(mode = 1:length(top), real = real.(top), imag = imag.(top), magnitude = abs.(top)),
)

summary = Dict(
    "timestamp_utc" => string(now(UTC)),
    "julia_version" => string(VERSION),
    "powerdynamics_source" => pkgdir(PowerDynamics),
    "state_dim" => length(eq.state),
    "raw_mode_count" => length(stab.eigenvalues),
    "gauge_mode_count" => length(stab.gauge_eigenvalues),
    "dynamic_margin" => stab.dynamic_margin,
    "max_real_nontrivial" => stab.max_real,
    "stable" => stab.stable,
)
open(joinpath(OUT_DIR, "baseline.toml"), "w") do io
    TOML.print(io, summary)
end

println("PD39 baseline complete: margin=$(stab.dynamic_margin), stable=$(stab.stable)")
