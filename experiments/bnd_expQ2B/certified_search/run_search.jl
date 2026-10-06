include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA
const L=LocalOracle
Base.include(L,joinpath(@__DIR__,"MultiPeak.jl"))
include("SQP.jl")
include("StaticFrequency.jl")
d=TOML.parsefile(abspath(ARGS[1]));a=L.architecture(d["support"],d);x=L.encode(a.support,d)
sm=StaticFrequency.build()
println("STATIC zero=",StaticFrequency.frequency(sm,zeros(10))," point=",StaticFrequency.frequency(sm,1 .-d["rho"]));flush(stdout)
result=solve_multi(a,x;maxiter=length(ARGS)>1 ? parse(Int,ARGS[2]) : 180,label="MULTIPEAK_SOC")
println("SEARCH_DONE");flush(stdout)
for line in eachline(stdin)
    strip(line)=="EXIT" && break
    try;include_string(Main,line,"interactive_search");catch err;showerror(stdout,err,catch_backtrace());println();end
    println("COMMAND_DONE");flush(stdout)
end
