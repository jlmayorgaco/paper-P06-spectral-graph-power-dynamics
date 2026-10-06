include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA
const L=LocalOracle
Base.include(L,joinpath(@__DIR__,"MultiPeak.jl"))
include("SQP.jl");include("NewtonKKT.jl")
d=TOML.parsefile(abspath(ARGS[1]));a=L.architecture(d["support"],d);x=L.encode(a.support,d)
result=newton_kkt(a,x;maxiter=20)
println("NEWTON_DONE");flush(stdout)
for line in eachline(stdin)
    strip(line)=="EXIT" && break
    try;include_string(Main,line,"interactive_newton");catch err;showerror(stdout,err,catch_backtrace());println();end
    println("COMMAND_DONE");flush(stdout)
end
