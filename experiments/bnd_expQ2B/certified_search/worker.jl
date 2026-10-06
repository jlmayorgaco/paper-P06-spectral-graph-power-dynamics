include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA
const L=LocalOracle
println("WORKER_READY");flush(stdout)
for line in eachline(stdin)
    strip(line)=="EXIT" && break
    try
        include_string(Main,line,"interactive_audit")
    catch err
        showerror(stdout,err,catch_backtrace());println()
    end
    println("COMMAND_DONE");flush(stdout)
end
