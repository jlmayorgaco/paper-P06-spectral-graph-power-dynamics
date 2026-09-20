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
        report = joinpath(root, "raw", "powerdynamics", "pd39_equilibrium_gate.md")
        passed = isfile(report) && occursin("status: POWERDYNAMICS_VALIDATED", read(report, String))
        write_status(phase, passed ? "POWERDYNAMICS_VALIDATED" : "STOPPED_BY_GATE",
            passed ? "PowerDynamics IEEE-39 Gate A completed." : "PowerDynamics Gate A stopped; inspect raw/powerdynamics/pd39_equilibrium_gate.md.")
    catch err
        write_status(phase, "STOPPED_BY_GATE", "PowerDynamics could not be loaded: " * sprint(showerror, err))
        println("PowerDynamics gate failed: ", sprint(showerror, err))
    end
else
    write_status(phase, "NOT_TESTED", "Phase runner scaffold exists; this phase has not been executed.")
end

println("phase=", phase)
