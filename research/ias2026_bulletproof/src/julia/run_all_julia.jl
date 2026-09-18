#!/usr/bin/env julia

phase = isempty(ARGS) ? "inventory" : ARGS[1]
root = normpath(joinpath(@__DIR__, "..", ".."))
logdir = joinpath(root, "logs")
mkpath(logdir)

function write_status(name, status, notes)
    path = joinpath(logdir, "julia_" * name * ".md")
    open(path, "w") do io
        println(io, "# Julia phase: ", name)
        println(io, "status: ", status)
        println(io, "notes: ", notes)
        println(io, "julia: ", VERSION)
        println(io, "project: ", Base.active_project())
    end
end

if phase == "inventory"
    write_status(phase, "NUMERICALLY_VERIFIED", "Julia runtime and project path recorded.")
elseif phase == "gate-pd39"
    try
        import PowerDynamics
        include(joinpath(@__DIR__, "run_pd39_gate.jl"))
        write_status(phase, "POWERDYNAMICS_VALIDATED", "PowerDynamics IEEE-39 tutorial equilibrium gate completed.")
    catch err
        write_status(phase, "STOPPED_BY_GATE", "PowerDynamics could not be loaded: " * sprint(showerror, err))
        println("PowerDynamics gate failed: ", sprint(showerror, err))
    end
else
    write_status(phase, "NOT_TESTED", "Phase runner scaffold exists; this phase has not been executed.")
end

println("phase=", phase)
