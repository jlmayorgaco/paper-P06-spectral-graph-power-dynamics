# T15: save frequency trajectories of one event for a retuning rule (method of steps, exact delay).
# usage: julia --project=. t15_trajectories.jl LABEL RULE NODE_RE NODE_IM H_MS RHO HORIZON EVENT_IDX
#   RULE = exact | lowfreq | none ; EVENT_IDX 1..5 (2 = bus 16 +100 MW)
using LinearAlgebra, CSV, DataFrames, TOML
const EXP = joinpath(@__DIR__, "..", "..", "experiments")
include(joinpath(EXP, "graph_gsp_codesign_20261003", "DelayedEvents.jl"))
const R = DelayedEvents.R
BLAS.set_num_threads(1)

function transport(kp, ki, lam, h)
    a, w = real(lam), imag(lam); e = exp(a * h)
    (e * (kp * cos(w * h) + (ki + a * kp) / w * sin(w * h)),
     e * (ki * cos(w * h) - (abs2(lam) * kp + a * ki) / w * sin(w * h)))
end

function main()
    label, rule = ARGS[1], ARGS[2]
    lam = complex(parse(Float64, ARGS[3]), parse(Float64, ARGS[4]))
    h = parse(Float64, ARGS[5]) / 1000; rhou = parse(Float64, ARGS[6]); horizon = parse(Float64, ARGS[7]); idx = parse(Int, ARGS[8])
    ctx = R.N.design_context(R.ROOT)
    base = TOML.parsefile(joinpath(EXP, "graph_gsp_codesign_20261003", "baseline.toml"))
    kp0, ki0 = Float64.(base["Kp"]), Float64.(base["Ki"]); rho = fill(rhou, 10); tau = 0.040 + h
    kp, ki = copy(kp0), copy(ki0)
    if rule == "exact"
        g = [transport(kp0[i], ki0[i], lam, h) for i in 1:10]; kp, ki = first.(g), last.(g)
    elseif rule == "lowfreq"
        kp = kp0 .+ h .* ki0
    end
    bus, delta = DelayedEvents.CASES[idx]
    println(label, " tau=", tau, " event=", bus, " ", delta); flush(stdout)
    m = R.model(ctx, rho, kp, ki; bus, delta, dc_convention=:physical_supply)
    try
        met, r = DelayedEvents.adaptive(m; tau=tau, horizon=horizon, dtmax=.01, tol=1e-9)
        R.metrics(m, r; dt=.01, window=.5, monitor_buses=collect(1:39), savepath=joinpath(@__DIR__, "T15_" * label * "_trajectory.csv"))
        println("DONE ", label, " F=", met.Fpeak_Hz, " R=", met.Rpeak_Hz_s); flush(stdout)
    catch err
        println("ABORTED ", label, " ", first(sprint(showerror, err), 200)); flush(stdout)
    end
end
main()
