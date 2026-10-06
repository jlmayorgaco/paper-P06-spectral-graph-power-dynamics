# Phase M: the frozen five-event nonlinear campaign (DelayedEvents.CASES, same adaptive method-of-steps solver and acceptance guards) for a design TOML.
# usage: julia --project=. events.jl CASE.toml LABEL   -> raw/events/events_LABEL.csv   (nothing is written into historical directories)
using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const DE = DelayedEvents; const R = DE.R
BLAS.set_num_threads(1)
d = TOML.parsefile(ARGS[1]); label = ARGS[2]
rho, kp, ki = Float64.(d["rho"]), Float64.(d["Kp"]), Float64.(d["Ki"]); tau = Float64(get(d, "tau", 0.04))
ctx = R.N.design_context(R.ROOT); rows = NamedTuple[]
dir = joinpath(@__DIR__, "..", "raw", "events"); mkpath(dir)
for (bus, delta) in [(16, 100.0)]
    println("EVENT_START ", label, " ", bus, " ", delta); flush(stdout)
    m = R.model(ctx, rho, kp, ki; bus, delta, dc_convention=:physical_supply)
    try
        met, r = DE.adaptive(m; tau)
        trajdir = joinpath(dir, "traj"); mkpath(trajdir)
        R.metrics(m, r; dt=.01, window=.5, monitor_buses=collect(1:39), savepath=joinpath(trajdir, "traj_$(label)_bus$(bus)_$(Int(delta)).csv"))
        push!(rows, (; label, tau_ms=1000tau, bus, delta, F=met.Fpeak_Hz, R=met.Rpeak_Hz_s, Vmin=met.Vmin, Vmax=met.Vmax, slack=met.limiter_fraction, runtime=met.runtime_s,
            pass=met.Fpeak_Hz <= .5 && met.Rpeak_Hz_s <= .5 && met.Vmin >= .9 && met.Vmax <= 1.1 && met.limiter_fraction >= .002, error=""))
    catch err
        push!(rows, (; label, tau_ms=1000tau, bus, delta, F=NaN, R=NaN, Vmin=NaN, Vmax=NaN, slack=NaN, runtime=NaN, pass=false, error=sprint(showerror, err)[1:min(end, 100)]))
    end
    CSV.write(joinpath(dir, "events_$(label).csv"), DataFrame(rows)); println("EVENT_DONE ", last(rows)); flush(stdout)
end
